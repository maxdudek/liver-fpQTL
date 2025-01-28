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

cat("Separating allelic fragments file")
for (df in list_of_dfs) {
  rsID <- df$rsID[1]
  allele <- df$allele[1]
  
  barcode <- paste(rsID, allele, sep = "_") %>% gsub(";", "_", ., fixed = TRUE)
  print(barcode)
  
  outfile <- paste0("allelic_fragments/", barcode, ".tsv")
  
  if (file.exists(outfile)) {
    next
  }
  
  df %>%
    mutate(barcode = barcode) %>%
    select(variant_chrom, insertion_5prime, insertion_3prime, barcode, count) %>%
    write.table(outfile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

# Get FP scores
FP_scores <- data.frame()
for (i in 1:nrow(fpQTLs)) {
  for (allele in c("ref", "alt")) {
    
    rsID <- fpQTLs$variant_id[i]
    barcode <- paste(rsID, allele, sep = "_") %>% gsub(";", "_", ., fixed = TRUE)
    cat("Barcode = ", barcode, "\n")
    
    # Initialize project
    project <- footprintingProject(projectName = barcode,
                                   refGenome = "hg38")
    projectMainDir <- "PRINT_temp/"
    projectDataDir <- paste0(projectMainDir, "data/", barcode, "/")
    dataDir(project) <- projectDataDir
    mainDir(project) <- projectMainDir
    
    # Set some project variables
    pathToFrags <- paste0("allelic_fragments/", barcode, ".tsv")
    
    if (!file.exists(pathToFrags)) {
      cat("DUDEK skipping because frag file doesn't exist...\n")
      next
    }
    
    barcodeGroups <- data.frame(
      barcode = barcode,
      group = 1L
    )
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
    h5Path <- paste0(PRINTdir, "data/TFBSPrediction/TFBS_model.h5")
    TFBindingModel(project) <- loadTFBSModel(h5Path)
    
    tryCatch(
        {
          # Get TF habitation scores
  	    	project <- getTFBS(project, 
  				           innerChunkSize = 2,
  				           chunkSize = 2,
  				           nCores = 1)
    	    
    			# Extract FP score
    			TFBS_results <- readRDS(paste0(dataDir(project), "chunkedTFBSResults/chunk_1.rds"))
    			FP_score <- TFBS_results[[1]]$TFBSScores[1,1]
        },
          error = function(cond) {
    		    message("DUDEK ERROR:")
    		    message(conditionMessage(cond))
    		    FP_score <- NA
          },
          finally = {
    		    new_row <- data.frame(
    				variant_id = rsID,
    				variant_chrom = chrom,
    				variant_pos = fpQTLs$variant_pos[i],
    				allele = allele,
    				FP_score = FP_score,
    				row.names = NULL
			)
	    
	    FP_scores <- rbind(FP_scores, new_row)
        
        }
    )
  }
}

FP_scores %>%
  write.table("allelic_FP_scores/%s_allelic_FP_scores.txt" %>% sprintf(FP_METHOD), 
              row.names = FALSE, quote = FALSE, sep = "\t")

