suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(qvalue)
  library(betareg)
})

ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

args = commandArgs(trailingOnly=TRUE)
NCORES <- 1
if (length(args) > 0) {
  NCORES <- as.integer(args[1])
}
cat("N cores = %d\n" %>% sprintf(NCORES))

run_regression_on_variant <- function(i, covariates=TRUE, link="logit", phi_model=NULL) {
  if (i %% 10000 == 0) {
    print(sprintf("Ran regression on %d variants", i))
  }
  
  x <- genotype_matrix[i,]
  y <- fpscore_matrix[i,]
  
  model.formula <- "y ~ x"
  if (covariates) { model.formula <- paste0(model.formula, " + regression_covariates$sex + regression_covariates$batch + 
                                            regression_covariates$PC1 + regression_covariates$PC2 + regression_covariates$PC3") }
  
  if (!is.null(phi_model)) { model.formula <- paste0(model.formula, " | ", phi_model) }
  
  model <- betareg(as.formula(model.formula), link=link)
  
  
  coef <- summary(model)$coefficients$mean
  beta <- coef[2,1]
  pval <- coef[2,4]
  r_squared <- summary(model)$pseudo.r.squared
  
  c(beta = beta, r_squared = r_squared, pval = pval)
}


get_regression_result <- function(ncores = 1L, ...) {

  if (ncores > 1) {
    # Run regressions in parallel
    cl <- makeCluster(ncores, outfile="job_out/cluster_out_%.0f.txt" %>% sprintf(as.numeric(Sys.time())))
    f = file(); sink(file=f) # Silence output
    clusterEvalQ(cl, { suppressPackageStartupMessages({library(tidyverse); library(betareg)}) })
    sink(); close(f)
    clusterExport(cl=cl, varlist=c("genotype_matrix", "fpscore_matrix", "regression_covariates"))
    regression_result <- parSapply(cl, 1:nrow(fpscore_matrix), run_regression_on_variant, ...)
    stopCluster(cl)
  } else {
    regression_result <- sapply(1:nrow(fpscore_matrix), run_regression_on_variant, ...)
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
# Regression                                                    |
#----------------------------------------------------------------
DATA_DIR <- "FP_methods/PRINT_no_gaussian/data/"

cat("Loading data without missing variants removed...\n")
fpscore_matrix <- readRDS(paste0(DATA_DIR, "fpscore_matrix_normalized.Rds"))
genotype_matrix <- readRDS(paste0(DATA_DIR, "genotype_matrix_regression.Rds"))
variant_info <- readRDS(paste0(DATA_DIR, "variant_info_regression.Rds"))

# local_insertions <- readRDS("../PRINT/consolidated_results/PRINT_local_insertion_matrix.Rds")

regression_covariates <- read.delim(paste0(ROOT_DIR, "raw_data/brandon_liver_ATAC/sample_covariates/regression_covariates_170.txt"))

regression_covariates %>%
  filter(sample_id %in% colnames(fpscore_matrix)) %>%
  arrange(match(sample_id, colnames(fpscore_matrix))) ->
  regression_covariates

cat("Verifying sample order in covariates...\n")
all(colnames(fpscore_matrix) == regression_covariates$sample_id)

cat("Verifying variant order in data tables...\n")
all(rownames(fpscore_matrix) == variant_info$variant_id)
all(rownames(genotype_matrix) == variant_info$variant_id)
# all(rownames(local_insertions) == variant_info$variant_id)

cat("Verifying sample order in data tables...\n")
all(colnames(fpscore_matrix) == colnames(genotype_matrix))
# all(colnames(local_insertions) == colnames(genotype_matrix))

PHI_MODEL_COV="x + regression_covariates$sex + regression_covariates$batch +
                                            regression_covariates$PC1 + regression_covariates$PC2 + regression_covariates$PC3"

for (link in c("logit", "probit", "cloglog")) {
  for (covariates in c(TRUE, FALSE)) {
    for (phi_model in list(NULL)) {
      cat("Running regressions, covariates=%s, link=%s, phi_model=%s...\n" %>% sprintf(covariates, link, format(phi_model)))
      
      FP_method <- "PRINT_beta_%s" %>% sprintf(link)
      if (!is.null(phi_model)) {FP_method<-paste0(FP_method, "_phi")}
      DIR <- paste0("FP_methods/", FP_method)
      dir.create(file.path(DIR, "regression_results"), showWarnings = FALSE, recursive = TRUE)
      
      
      outfile <- ifelse(covariates, "fp_score_covariates_genotype_regression.Rds", "fp_score_genotype_regression.Rds")
      outfile <- paste0(DIR, "/regression_results/", outfile)
      if (!is.null(phi_model)) { phi_model <- ifelse(covariates, PHI_MODEL_COV, "x") }
      
      if (file.exists(outfile)) {
        cat("WARNING: skipping regression because outfile exists!\n")
        next
      }
      
      fpscore_regression <- get_regression_result(ncores=NCORES, covariates=covariates, link=link, phi_model=phi_model)
      
      cat("Writing results...\n")
      fpscore_regression %>% saveRDS(outfile)
      
      # Get number of significant fpQTLs at this level
      fpscore_regression_sig <- fpscore_regression %>%
        filter(ST_qval < 0.05)
      
      num_sig_fpQTLs <- fpscore_regression_sig %>% nrow()
      
      cat("num_sig_fpQTLs at FDR 0.05 = %d\n" %>% sprintf(num_sig_fpQTLs))
    }
  }
}

cat("Regression Done\n")
quit()


