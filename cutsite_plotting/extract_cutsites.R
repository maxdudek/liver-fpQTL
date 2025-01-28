suppressPackageStartupMessages({
  library(tidyverse)
})

# Usage:
# Rscript extract_cutsites.R [FP_METHOD]
# or
# Rscript extract_cutsites.R [FP_METHOD] rsXXXXX1 rsXXXXX2 ...
#
# If a list of rsIDs is provided, then cutsites for ONLY those 
# variants will be extracted
# Otherwise, cutsites for the top fpQTLs with q < 0.05 will be extracted

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]


DIR <- paste0("../regression/FP_methods/", FP_METHOD)
cat("FP_METHOD = ", FP_METHOD, "\n")

cat("Loading data...\n")
fpQTLs <- read.delim(paste0(DIR, "/regression_results/fpQTLs_covariates_FDR5.txt"))
variant_info <- variant_info <- readRDS("../regression/variant_info.Rds")
fpscore_matrix <- readRDS(paste0(DIR, "/data/fpscore_matrix_normalized_fpQTLs_with_covariates.Rds"))
genotype_matrix <- readRDS(paste0(DIR, "/data/genotype_matrix_fpQTLs_with_covariates.Rds"))
variant_predBias <- readRDS("../PRINT/variant_predBias.rds")
REGION_SIZE <- ncol(variant_predBias)
W <- REGION_SIZE %/% 2

CHUNK_SIZE <- 1e5
# NUM_VARIANTS <- min(nrow(fpQTLs), 25)
NUM_VARIANTS <- nrow(fpQTLs)

if (length(args) > 1) {
  variants <- args[2:length(args)]
} else {
  variants <- fpQTLs$variant_id[1:NUM_VARIANTS]
}

cutsite_dfs <- list()
for (rs in variants) {
    cat("Extracting cutsite data for ", rs, "...\n", sep="")

    # Get variant number
    i <- match(rs, variant_info$variant_id)

    # Get chunk number
    chunk <- ((i-1) %/% CHUNK_SIZE) + 1

    # Get position in chunk
    chunk_i <- i %% CHUNK_SIZE
    if (chunk_i == 0) {chunk_i <- CHUNK_SIZE} # R uses 1-indexing

    # Get samples in regression
    variant_genotypes <- genotype_matrix[rs, ] %>% na.omit()
    samples <- names(variant_genotypes)
    variant_cutsites <- data.frame()
    
    # Get bias for variant region
    region_bias <- variant_predBias[i, ]

    for (s in samples) {
        # cat("\tGetting chunked count tensor for sample", s, "\n")
        g <- variant_genotypes[s]
        fpscore <- fpscore_matrix[rs, s]

        countTensor <- readRDS(paste0("../PRINT/data/", s, "/chunkedCountTensor/chunk_", chunk, ".rds"))

        region <- countTensor[[chunk_i]]

        # If there are no local cutsites, add a dummy row so the sample is represented
        if (nrow(region) == 0) {
            region <- data.frame(region = i, position = 101, count = 0)
        }

        if (i != region$region[1]) {
            cat("ERROR: variant number i =", i, "for", rs, "does not match region", region$region[1], "in chunk", chunk, "\n")
            quit()
        }

        raw_insertions <- rep(0, REGION_SIZE)
        raw_insertions[region$position] <- region$count
        
        expected_insertions <- region_bias/sum(region_bias) * sum(raw_insertions)
        corrected_insertions <- raw_insertions - expected_insertions
        
        sample_cutsites <- data.frame(
          sample = s,
          genotype = g,
          position = -W:W,
          num_insertions = raw_insertions,
          bias = region_bias,
          expected_insertions = expected_insertions,
          corrected_insertions = corrected_insertions,
          row.names = NULL
        )

        variant_cutsites <- rbind(variant_cutsites, sample_cutsites)
    }

    # cat("\tSaving cutsites...\n")
    # variant_cutsites %>% saveRDS(paste0("cutsites/", FP_METHOD, "/", rs, ".Rds"))
    cutsite_dfs[[rs]] <- variant_cutsites

}

cat("Writing cutsites...\n")
cutsite_dfs %>% saveRDS(paste0("cutsites/", FP_METHOD, "/fpQTL_cutsites.Rds"))
