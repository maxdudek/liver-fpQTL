suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(Rcpp)
  library(RcppArmadillo)
})

ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]
DIR <- paste0("FP_methods/", FP_METHOD, "/")
cat("FP_METHOD = ", FP_METHOD, "\n")

# temp
1:10 %>% saveRDS(paste0(DIR, "regression_results/residuals_with_covariates.Rds"))

NCORES <- 1
if (length(args) > 1) {
  NCORES <- as.integer(args[2])
}
cat("N cores = %d\n" %>% sprintf(NCORES))

cat("Loading data...\n")
fpscore_matrix <- readRDS(paste0(DIR, "data/fpscore_matrix_normalized.Rds"))
genotype_matrix <- readRDS(paste0(DIR, "data/genotype_matrix_regression.Rds"))
variant_info <- readRDS(paste0(DIR, "data/variant_info_regression.Rds"))

regression_covariates <- read.delim(paste0(ROOT_DIR, "raw_data/brandon_liver_ATAC/sample_covariates/regression_covariates_170.txt"))

regression_covariates_mat <- regression_covariates %>%
  filter(sample_id %in% colnames(fpscore_matrix)) %>%
  arrange(match(sample_id, colnames(fpscore_matrix))) %>%
  column_to_rownames("sample_id") %>%
  select(sex, batch, PC1, PC2, PC3) %>%
  as.matrix()

cat("Verifying sample order in covariates...\n")
all(colnames(fpscore_matrix) == rownames(regression_covariates_mat))

cat("Verifying variant order in data tables...\n")
all(rownames(fpscore_matrix) == variant_info$variant_id)
all(rownames(genotype_matrix) == variant_info$variant_id)

cat("Verifying sample order in data tables...\n")
all(colnames(fpscore_matrix) == colnames(genotype_matrix))

run_regression_on_variant_fast <- function(i) {
  y <- fpscore_matrix[i,]
  X <- cbind(I = 1,genotype_matrix[i,],regression_covariates_mat)

  fLmSEXP(X, y)$res
}

cat("Starting cluster...\n")
cl <- makeCluster(NCORES, outfile="job_out/cluster_out.txt")
f = file(); sink(file=f) # Silence output

clusterEvalQ(cl, { 
  suppressPackageStartupMessages({
    library(tidyverse)
    library(Rcpp)
    library(RcppArmadillo)
  }) 
  
  # Modify code to extract residuals
  src <- '
  Rcpp::List fLmSEXP(SEXP Xs, SEXP ys) {
  Rcpp::NumericMatrix Xr(Xs);
  Rcpp::NumericVector yr(ys);
  int n = Xr.nrow(), k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false);
  arma::colvec y(yr.begin(), yr.size(), false);
  // fit model y ~ X, extract residuals
  arma::colvec coef = arma::solve(X, y);
  arma::colvec res  = y - X*coef;
  // return the results
  return Rcpp::List::create(Rcpp::Named("coefficients")=coef,Rcpp::Named("res")=res);
  }
  '
  cppFunction(code=src, depends="RcppArmadillo")
})

sink(); close(f)
clusterExport(cl=cl, varlist=c("genotype_matrix", "fpscore_matrix", "regression_covariates_mat"))
cat("Getting residuals...\n")
residual_matrix <- parSapply(cl, 1:nrow(fpscore_matrix), run_regression_on_variant_fast) %>% t()
stopCluster(cl)

rownames(residual_matrix) <- variant_info$variant_id
colnames(residual_matrix) <- colnames(fpscore_matrix)

dim(residual_matrix)
residual_matrix[1:10, 1:10]

cat("Writing fpscore residuals...\n")
residual_matrix %>% saveRDS(paste0(DIR, "regression_results/residuals_with_covariates.Rds"))
