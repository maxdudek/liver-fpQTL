library(tidyverse)
ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

FP_METHOD <- "PRINT_no_gaussian"

DIR <- paste0("../../regression/FP_methods/", FP_METHOD)
fpscore_cov_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))
variant_info <- readRDS(paste0(DIR, "/data/variant_info_regression.Rds"))

zoom1 <- c(4.53e7, 4.66e7)
zoom2 <- c(4.62e7, 4.63e7)

brandon_peaks <- read.table(paste0(ROOT_DIR, "raw_data/brandon_liver_ATAC/peak_info/genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak_chr.bed"))

CHUNK_SIZE <- 1e5

variants_in_region <- fpscore_cov_regression %>%
  filter(variant_chrom == "chr17") %>%
  filter(variant_pos > zoom2[1] & variant_pos < zoom2[2]) %>%
  mutate(i = match(variant_id, variant_info$variant_id))

# Get list of samples to load
genotype_matrix <- readRDS(paste0(DIR, "/data/genotype_matrix_fpQTLs_with_covariates.Rds"))
samples <- colnames(genotype_matrix)


# Extract cutsites
cutsites <- data.frame()
for (s in samples) {
  print(s)
  countTensor <- readRDS(paste0("../../PRINT/data/", s, "/chunkedCountTensor/chunk_29.rds"))
  
  for (j in 1:nrow(variants_in_region)) {
    # Get variant info
    rs <- variants_in_region$variant_id[j]
    i <- variants_in_region$i[j]
    variant_pos <- variants_in_region$variant_pos[j]
    
    # Get position in chunk
    chunk_i <- i %% CHUNK_SIZE
    if (chunk_i == 0) {chunk_i <- CHUNK_SIZE} # R uses 1-indexing
    region <- countTensor[[chunk_i]]
    
    # If there are no local cutsites, add a dummy row so the sample is represented
    if (nrow(region) == 0) {
      region <- data.frame(region = i, position = 101, count = 0)
    }
    
    if (i != region$region[1]) {
      cat("ERROR: variant number i =", i, "for", rs, "does not match region",
          region$region[1], "in chunk", chunk, "\n")
      break
    }
    
    sample_cutsites <- data.frame(
      sample = s,
      position = variant_pos + region$position - 101,
      num_insertions = region$count,
      row.names = NULL
    )
    
    cutsites <- rbind(cutsites, sample_cutsites)
  }
}

cat("Removing duplicate rows...\n")
cutsites <- cutsites[!duplicated(cutsites),]

cat("Saving cutsites...\n")
cutsites %>% saveRDS("cutsites/chr17_peak_cutsites.Rds")
