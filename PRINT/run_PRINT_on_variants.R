suppressPackageStartupMessages({
    library(tidyverse)
    library(parallel)
    library(foreach)
    # library(reticulate)
    # use_condaenv("PRINT", required = TRUE)

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

Sys.setenv(PYTHONPATH = "") # Loading modules populates this variable with other python installations
options(ignore.interactive = FALSE) # Prevent progress bar from being printed when run as a script

args = commandArgs(trailingOnly=TRUE)
SAMPLE_ID <- args[1]

CORES <- 2
CHUNK_SIZE <- 1e5

# install.packages("tensorflow", repos='http://cran.us.r-project.org')
# library(tensorflow)
# install_tensorflow()
# library(reticulate)
# library(keras)
# install_keras()
# reticulate::py_install(
#     "numpy", pip = TRUE,
#     pip_options = c("--user", "--force-reinstall", "--upgrade"))
# numpy <- import("numpy")
# py_discover_config()
# reticulate::py_config() 

cat("DUDEK: Initializing project 3782-BW-", SAMPLE_ID, "...\n")
projectName <- paste0("3782-BW-", SAMPLE_ID)
project <- footprintingProject(projectName = projectName,
                               refGenome = "hg38")
projectMainDir <- "./"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Set some project variables
pathToFrags <- paste0("frags/", projectName, ".tsv")

barcodeGroups <- data.frame(
  barcode = projectName,
  group = 1L
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))
groupCellType(project) <- "liver"

out_filename <- paste0("results_by_sample/", projectName, "_variant_results.rds")
if (file.exists(out_filename)) {
    cat("WARNING DUDEK: The file ", out_filename, " already exists, so skipping.")
    quit()
}

cat("DUDEK: Loading regions of interest...\n")
variant_info_ocr <- readRDS("../../../raw_data/brandon_liver_ATAC/vcf/variant_info_maf5_ocr.Rds")
w <- 100

regions <- GRanges(seqnames = variant_info_ocr$variant_chrom,
                   ranges = IRanges(start = variant_info_ocr$variant_pos-w, 
                                    end = variant_info_ocr$variant_pos+w))
regionRanges(project) <- regions

cat("DUDEK: Loading region bias...\n")
regionBias(project) <- readRDS("variant_predBias.rds")

cat("DUDEK: Getting count tensors...\n")
project <- getCountTensor(project, 
                          pathToFrags, 
                          barcodeGroups, 
                          returnCombined = FALSE,
                          chunkSize = CHUNK_SIZE,
                          nCores = CORES)

cat("DUDEK: Loading dispersion and TFBS models...\n")
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <-
  readRDS(paste0(PRINTdir, "data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

h5Path <- paste0(PRINTdir, "data/TFBSPrediction/TFBS_model.h5")
TFBindingModel(project) <- loadTFBSModel(h5Path)

cat("DUDEK: Getting TFBS habitation scores...\n")
project <- getTFBS(project, 
                   innerChunkSize = 100,
                   chunkSize = CHUNK_SIZE,
                   nCores = CORES)

warnings()

cat("DUDEK: Extracting TFBS and insertion counts...\n")
TFBSDir <- paste0(dataDir(project), "chunkedTFBSResults/")
chunkDir <- paste0(dataDir(project), "chunkedCountTensor/")
TFBSChunkFiles <- gtools::mixedsort(list.files(TFBSDir))
nChunks <- length(TFBSChunkFiles)

scores <- c()
nInsertions <- c()
for (i in 1:nChunks) {
  cat("\tChunk ", i, "\n")
  TFBSChunkData <- readRDS(paste0(TFBSDir, "chunk_", i, ".rds"))
  CountChunkData <- readRDS(paste0(chunkDir, "chunk_", i, ".rds"))
  
  # Add scores
  scores <- c(scores, sapply(TFBSChunkData, function(x) {x$TFBSScores}))
  nInsertions <- c(nInsertions, sapply(CountChunkData, function(x) {sum(x$count)}))
}

cat("DUDEK: Writing results...\n")
data.frame(
  variant_id = variant_info_ocr$variant_id,
  TFBS = scores,
  local_insertion_count = as.integer(nInsertions)
) %>% 
    saveRDS(out_filename)

cat("DUDEK: Done\n")
