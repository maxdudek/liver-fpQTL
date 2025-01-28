suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(qvalue)
})



args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]
DIR <- paste0("FP_methods/", FP_METHOD)
cat("FP_METHOD = ", FP_METHOD, "\n")

# NCORES <- detectCores()
# NCORES <- 16
NCORES <- 1
if (length(args) > 1) {
  NCORES <- as.integer(args[2])
}
cat("N cores = %d\n" %>% sprintf(NCORES))

run_regression_on_variant <- function(i, covariates = FALSE, save_covariate_residuals = FALSE) {
  if (i %% 1000 == 0) {
    print(sprintf("Ran regression on %d variants", i))
  }
  
  x <- genotype_matrix[i,]
  y <- fpscore_matrix[i,]

  if (covariates) {
    if(save_covariate_residuals) {
      # Regress out covariates
      covariate_model <- lm(y ~ regression_covariates$sex + regression_covariates$batch + regression_covariates$PC1 +
                      regression_covariates$PC2 + regression_covariates$PC3)
      covariate_residuals <- residuals(covariate_model)

      # Store residuals 
      fpscore_covariate_residual_matrix[i,] <- covariate_residuals
    } 

    model <- lm(y ~ x + regression_covariates$sex + regression_covariates$batch + regression_covariates$PC1 +
                regression_covariates$PC2 + regression_covariates$PC3)          
  } else {
    model <- lm(y ~ x)
  }
  
  coef <- summary(model)$coefficients
  beta <- coef[2,1]
  pval <- coef[2,4]
  r_squared <- summary(model)$r.squared
  
  c(beta = beta, r_squared = r_squared, pval = pval)
}


get_regression_result <- function(covariates = FALSE, save_covariate_residuals = FALSE, ncores = 1L) {

  if (ncores > 1) {
    # Run regressions in parallel
    cl <- makeCluster(ncores, outfile="job_out/cluster_out.txt")
    f = file(); sink(file=f) # Silence output
    clusterEvalQ(cl, { suppressPackageStartupMessages({library(tidyverse)}) })
    sink(); close(f)
    clusterExport(cl=cl, varlist=c("genotype_matrix", "fpscore_matrix", "regression_covariates", "fpscore_covariate_residual_matrix"))
    regression_result <- parSapply(cl, 1:nrow(fpscore_matrix), run_regression_on_variant, 
                                   covariates=covariates, save_covariate_residuals=save_covariate_residuals)
    stopCluster(cl)
  } else {
    regression_result <- sapply(1:nrow(fpscore_matrix), run_regression_on_variant,
                                covariates=covariates, save_covariate_residuals=save_covariate_residuals)
  }
  
  # Turn results into data frame
  data.frame(
    variant_id = rownames(fpscore_matrix),
    beta = regression_result["beta", ],
    r_squared = regression_result["r_squared", ],
    pval = regression_result["pval", ]
  ) %>%
    cbind(variant_info %>% select(-variant_id)) ->
    fpscore_regression

  cat("Verifying that variant_ids line up in regression_results:\n")
  print(all(fpscore_regression$variant_id == variant_info$variant_id))
  
  # Calculate Storey q-value
  qval.object <- qvalue(fpscore_regression$pval)
  summary(qval.object)
  fpscore_regression$ST_qval <- qval.object$qvalues
  fpscore_regression$ST_lfdr <- qval.object$lfdr

  fpscore_regression
}


#----------------------------------------------------------------
# Regression WITHOUT missing variants removed (No PEER factors) |
#----------------------------------------------------------------
cat("Loading data without missing variants removed...\n")
fpscore_matrix <- readRDS(paste0(DIR, "/data/fpscore_matrix_normalized.Rds"))
genotype_matrix <- readRDS(paste0(DIR, "/data/genotype_matrix_regression.Rds"))
variant_info <- readRDS(paste0(DIR, "/data/variant_info_regression.Rds"))

# variant_info %>%
#   pull(variant_id) %>%
#   write(paste0(DIR, "/data/variant_info_rsIDs.txt"))

# Create an empty matrix to store the residuals from regressing out covariates
fpscore_covariate_residual_matrix <- matrix(, nrow = nrow(fpscore_matrix), ncol = ncol(fpscore_matrix))

regression_covariates <- read.table("../../../raw_data/brandon_liver_ATAC/sample_covariates/regression_covariates_170.txt", header = TRUE)

regression_covariates %>%
  filter(sample_id %in% colnames(fpscore_matrix)) %>%
  arrange(match(sample_id, colnames(fpscore_matrix))) ->
  regression_covariates

cat("Verifying sample order in covariates...\n")
all(colnames(fpscore_matrix) == regression_covariates$sample_id)

cat("Verifying variant order in data tables...\n")
all(rownames(fpscore_matrix) == variant_info$variant_id)
all(rownames(genotype_matrix) == variant_info$variant_id)

cat("Verifying sample order in data tables...\n")
all(colnames(fpscore_matrix) == colnames(genotype_matrix))

cat("Running regressions without covariates...\n")
fpscore_regression <- get_regression_result(covariates=FALSE, ncores=NCORES)

cat("Writing result without covariates...\n")
fpscore_regression %>% saveRDS(paste0(DIR, "/regression_results/fp_score_genotype_regression.Rds"))

# Get number of significant fpQTLs at this level
fpscore_regression %>%
  filter(ST_qval < 0.05) ->
  fpscore_regression_sig

fpscore_regression_sig %>%
  nrow() ->
  num_sig_fpQTLs
cat("num_sig_fpQTLs at FDR 0.05 = %d\n" %>% sprintf(num_sig_fpQTLs))

fpscore_regression_sig %>%
  arrange(pval) %>%
  head()

cat("Running regressions WITH covariates...\n")
fpscore_regression <- get_regression_result(covariates=TRUE, save_covariate_residuals=FALSE, ncores=NCORES)

cat("Writing result WITH covariates...\n")
fpscore_regression %>% saveRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))

cat("Writing fpscore residuals...\n")
fpscore_covariate_residual_matrix %>% saveRDS(paste0(DIR, "/data/fpscore_covariate_residual_matrix.Rds"))

# Get number of significant fpQTLs at this level
fpscore_regression %>%
  filter(ST_qval < 0.05) ->
  fpscore_regression_sig
  
fpscore_regression_sig %>%
  nrow() ->
  num_sig_fpQTLs
cat("num_sig_fpQTLs at FDR 0.05 = %d\n" %>% sprintf(num_sig_fpQTLs))

fpscore_regression_sig %>%
  arrange(pval) %>%
  head()

cat("Regression Done\n")
quit()

#------------------------------------------------------------
# Regression with missing variants removed (PEER residuals) |
#------------------------------------------------------------
# genotype_matrix <- readRDS("data/genotype_matrix_noMissing.Rds")

# peers <- 0:20
# num_fdr10_fpqtls <- rep(NA, length(peers))
# names(num_fdr10_fpqtls) <- peers

# for (peer in peers) {
#   cat("Loading data for PEER k = %d...\n" %>% sprintf(peer))
#   fpscore_matrix <- readRDS("data/fpscore_matrix_noMissing_normalized_PEER%d.Rds" %>% sprintf(peer))
#   rownames(fpscore_matrix) <- rownames(genotype_matrix)
#   colnames(fpscore_matrix) <- colnames(genotype_matrix)

#   cat("Running regressions for PEER k = %d...\n" %>% sprintf(peer))
#   fpscore_regression <- get_regression_result()
  
#   # Get number of significant fpQTLs at this level
#   fpscore_regression %>%
#     filter(ST_qval < 0.1) %>%
#     nrow() ->
#     num_sig_fpQTLs
#   num_fdr10_fpqtls[peer+1] <- num_sig_fpQTLs
  
#   cat("Writing results for PEER k = %d...\n" %>% sprintf(peer))
#   fpscore_regression %>%
#     write.table(file = "regression_results/fp_score_genotype_regression_PEER%d_chr21_test.txt.gz" %>% sprintf(peer) %>% gzfile(),
#                 quote = F,
#                 row.names = F,
#                 sep = "\t")
# }

# cat("Outputting results for all peers...\n")
# num_fdr10_fpqtls

# data.frame(
#   PEER = peers,
#   num_fdr10_fpqtls = num_fdr10_fpqtls
# ) %>%
#   ggplot(aes(x = PEER, y = num_fdr10_fpqtls)) +
#   geom_line() +
#   xlab("Number of PEER factors") +
#   ylab("Number of fpQTLs at FDR < 10%")

# ggsave("figures/num_fdr10_fpqtls_by_peers.png")



