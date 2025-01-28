suppressPackageStartupMessages({
  library(tidyverse)
  library(foreach)
  library(iterators)
  library(parallel)
  
  PRINTdir <- "/home/dudekmf/local/src/PRINT/"
  source(paste0(PRINTdir, "code/utils.R"))
  source(paste0(PRINTdir, "code/getCounts.R"))
  source(paste0(PRINTdir, "code/getBias.R"))
  source(paste0(PRINTdir, "code/getFootprints.R"))
  source(paste0(PRINTdir, "code/getSubstructures.R"))
  source(paste0(PRINTdir, "code/visualization.R"))
  source(paste0(PRINTdir, "code/getGroupData.R"))
  source(paste0(PRINTdir, "code/footprintTracking.R"))
  source(paste0(PRINTdir, "code/getTFBS.R"))
})

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]

Sys.setenv(PYTHONPATH = "") # Loading modules populates this variable with other python installations
options(ignore.interactive = FALSE) # Prevent progress bar from being printed when run as a script

fpQTLs <- read.delim("../../regression/FP_methods/%s/regression_results/fpQTLs_covariates_FDR5.txt" %>% sprintf(FP_METHOD))

# Separate allelic_fragments file into multiple files based on variant and allele
allelic_fragments <- read.delim("../allelic_fragments/%s_allelic_fragments.txt" %>% sprintf(FP_METHOD))

allelic_fragments <- allelic_fragments %>%
  inner_join(fpQTLs %>% select(variant_id, variant_chrom), by = join_by(rsID == variant_id))

list_of_dfs <- allelic_fragments %>%
  group_split(rsID, allele) 

for (df in list_of_dfs) {
  rsID <- df$rsID[1]
  allele <- df$allele[1]
  
  rsID_allele <- paste(rsID, allele, sep = "_") %>% gsub(";", "_", ., fixed = TRUE)
  print(rsID_allele)
  
  outfile <- paste0("allelic_fragments/", rsID_allele, ".tsv")
  
  if (file.exists(outfile)) {
    next
  }
  
  # Put sample in the barcode slot so we can calculate per-sample FP scores
  df %>%
    select(variant_chrom, insertion_5prime, insertion_3prime, sample, count) %>%
    write.table(outfile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

# Test
# rsID_allele <- "rs61968021_alt"
# i <- which(fpQTLs$variant_id == rsID_allele)

h5Path <- paste0(PRINTdir, "data/TFBSPrediction/TFBS_model.h5")
h5Model <- loadTFBSModel(h5Path)

result <- data.frame()
for (i in 1:nrow(fpQTLs)) {
  for (allele in c("ref", "alt")) {
    
    rsID <- fpQTLs$variant_id[i]
    rsID_allele <- paste(rsID, allele, sep = "_") %>% gsub(";", "_", ., fixed = TRUE)
    cat("------------------------------------------------------------------------------------------------------\n")
    cat("i =", i, "- rsID_allele = ", rsID_allele, "\n")
    cat("------------------------------------------------------------------------------------------------------\n")
    
    # Initialize project
    project <- footprintingProject(projectName = rsID_allele,
                                   refGenome = "hg38")
    projectMainDir <- "PRINT_temp/"
    projectDataDir <- paste0(projectMainDir, "data/", rsID_allele, "/")
    dataDir(project) <- projectDataDir
    mainDir(project) <- projectMainDir
    
    # Set some project variables
    pathToFrags <- paste0("allelic_fragments/", rsID_allele, ".tsv")
    
    if (!file.exists(pathToFrags)) {
      cat("DUDEK skipping because frag file doesn't exist...\n")
      next
    }
    
    frags_df <- read.table(pathToFrags)
    het_samples <- frags_df %>% pull(V4) %>% unique()
    n <- length(het_samples)
    
    if (n < 2) {
      cat("DUDEK skipping because we have less than 2 het samples...\n")
      next
    }
    
    barcodeGroups <- data.frame(
      barcode = het_samples,
      group = 1:n
    )
    cat("n =", n, "\n")
    barcodeGrouping(project) <- barcodeGroups
    groups(project) <- mixedsort(unique(barcodeGroups$group))
    groupCellType(project) <- "liver"
    
    # Select region (100 bp around variant)
    w <- 100
    chrom <- fpQTLs$variant_chrom[i]
    start <- fpQTLs$variant_pos[i]-w
    end <- fpQTLs$variant_pos[i]+w
    
    # Create a second "dummy" region because PRINT complains if we only have one region
    regions <- GRanges(seqnames = c(chrom, chrom),
                       ranges = IRanges(start = c(start, start+1), 
                                        end = c(end, end+1)))
    regionRanges(project) <- regions
    
    
    tryCatch(
      {
        chunk_filename <- paste0(dataDir(project), "chunkedTFBSResults/chunk_1.rds")
        
        if (!file.exists(chunk_filename)) { # If we haven't previously calculated
          # Calculate bias
          project <- getPrecomputedBias(project, nCores = 1)
          
          # Get count tensors
          project <- getCountTensor(project, 
                                    pathToFrags, 
                                    barcodeGroups, 
                                    returnCombined = F,
                                    chunkSize = 2,
                                    nCores = 1)
          
          # Load dispersion model
          for(kernelSize in 2:100){
            dispModel(project, as.character(kernelSize)) <-
              readRDS(paste0(PRINTdir, "data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
          }
          
          # Load TFBS model
          TFBindingModel(project) <- h5Model
          
          # Get TF habitation scores
          project <- getTFBS(project, 
                             innerChunkSize = min(n, 2),
                             chunkSize = 2,
                             nCores = 1)
        } else {
          cat("Using previously calculated results...\n")
        }
        
        # Extract FP score
        TFBS_results <- readRDS(chunk_filename)
        # cat("TFBS_results=\n----------------------------------\n")
        # print(TFBS_results)
        # cat("----------------------------------\n")
        
        # Sometimes TFBS_results has a different number of dimensions? idk
        if (is.list(TFBS_results[[1]][[1]])) {
          FP_scores <- TFBS_results[[1]][[1]]$TFBSScores[1,]
        } else {
          FP_scores <- TFBS_results[[1]]$TFBSScores[1,]
        }
        
        
      },
      error = function(cond) {
        message("DUDEK ERROR:")
        message(conditionMessage(cond))
        FP_scores <- rep(NA,n)
      },
      finally = {
        cat("length(FP_scores) =", length(FP_scores), "\n")
        new_rows <- data.frame(
          variant_id = rep(rsID, n),
          variant_chrom = rep(chrom, n),
          variant_pos = rep(fpQTLs$variant_pos[i], n),
          sample = het_samples,
          allele = rep(allele, n),
          FP_score = FP_scores,
          row.names = NULL
        )
        result <- rbind(result, new_rows)
        
      }
    )
  }
}

result %>%
  write.table("allelic_FP_scores_per_sample/%s_allelic_FP_scores_per_sample.txt" %>% sprintf(FP_METHOD), 
              row.names = FALSE, quote = FALSE, sep = "\t")

