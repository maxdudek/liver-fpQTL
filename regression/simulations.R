
suppressPackageStartupMessages({
  library(tidyverse)
  library(betareg)
  library(parallel)
})

ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"
set.seed(0)

args = commandArgs(trailingOnly=TRUE)
NCORES <- 1
if (length(args) > 0) {
  NCORES <- as.integer(args[1])
}
cat("N cores = %d\n" %>% sprintf(NCORES))


DATA_DIR <- "FP_methods/PRINT_no_gaussian/data/"

cat("Reading in data...\n")
fpscore_matrix <- readRDS(paste0(DATA_DIR, "fpscore_matrix_normalized.Rds"))
genotype_matrix <- readRDS(paste0(DATA_DIR, "genotype_matrix_regression.Rds"))
variant_info <- readRDS(paste0(DATA_DIR, "variant_info_regression.Rds"))

regression_covariates <- read.delim(paste0(ROOT_DIR, "raw_data/brandon_liver_ATAC/sample_covariates/regression_covariates_170.txt"))

regression_covariates %>%
  filter(sample_id %in% colnames(fpscore_matrix)) %>%
  arrange(match(sample_id, colnames(fpscore_matrix))) ->
  regression_covariates

cat("Verifying column/row order...\n")
all(colnames(fpscore_matrix) == regression_covariates$sample_id)
all(rownames(fpscore_matrix) == variant_info$variant_id)
all(rownames(genotype_matrix) == variant_info$variant_id)
all(colnames(fpscore_matrix) == colnames(genotype_matrix))


simulate_beta_regression_on_variant <- function(i, link="logit", spike_in_effect=0) {
  x <- genotype_matrix[i,]
  y <- fpscore_matrix_permuted[i,]
  
  if (spike_in_effect != 0) {
    # spike-in genomic effect on FP score
    y <- y + spike_in_effect*x
    y <- pmax(pmin(y, 0.99), 0.01) # clamp values to [0.01,0.99]
  }
  
  model <- betareg(y ~ x + regression_covariates$sex + regression_covariates$batch + 
                     regression_covariates$PC1 + regression_covariates$PC2 + regression_covariates$PC3, link=link)
  
  coef <- summary(model)$coefficients$mean
  beta <- coef[2,1]
  pval <- coef[2,4]
  
  c(beta = beta, pval = pval)
}

simulate_linear_regression_on_variant <- function(i, spike_in_effect=0, link=NULL) {
  x <- genotype_matrix[i,]
  y <- fpscore_matrix_permuted[i,]
  
  if (spike_in_effect != 0) {
    # spike-in genomic effect on FP score
    y <- y + spike_in_effect*x
    y <- pmax(pmin(y, 0.99), 0.01) # clamp values to [0.01,0.99]
  }
  
  model <- lm(y ~ x + regression_covariates$sex + regression_covariates$batch + 
                     regression_covariates$PC1 + regression_covariates$PC2 + regression_covariates$PC3)
  
  coef <- summary(model)$coefficients
  beta <- coef[2,1]
  pval <- coef[2,4]
  
  c(beta = beta, pval = pval)
}

cat("Setting up cluster...\n")
cl <- makeCluster(NCORES, outfile="simulations/cluster_out.txt")
f = file(); sink(file=f) # Silence output
clusterEvalQ(cl, { suppressPackageStartupMessages({library(tidyverse); library(betareg)}) })
sink(); close(f)
clusterExport(cl=cl, varlist=c("genotype_matrix", "regression_covariates"))


# null simulations

N_PERMUTATIONS <- 1000000
sample_i <- sample(1:nrow(fpscore_matrix), N_PERMUTATIONS, replace = TRUE) # Sample SNPs
fpscore_matrix_permuted <- fpscore_matrix[sample_i,] %>% apply(1, sample) %>% t() # Permute FP scores
clusterExport(cl=cl, varlist=c("fpscore_matrix_permuted"))

null_pvals <- list()
for (link in c("logit", "probit", "cloglog")) {
  cat(sprintf("Running null simulations for link=%s\n", link))
  null_pvals[[link]] <- parSapply(cl, 1:N_PERMUTATIONS, simulate_beta_regression_on_variant, link=link)
}
cat(sprintf("Running null simulations for linear model\n"))
null_pvals[["linear"]] <- parSapply(cl, 1:N_PERMUTATIONS, simulate_linear_regression_on_variant)
cat("Saving null pvals...\n")
null_pvals %>% saveRDS("simulations/null_pvals.Rds")


# alt simulations

N_PERMUTATIONS <- 10000
sample_i <- sample(1:nrow(fpscore_matrix), N_PERMUTATIONS, replace = TRUE) # Sample SNPs
fpscore_matrix_permuted <- fpscore_matrix[sample_i,] %>% apply(1, sample) %>% t() # Permute FP scores
clusterExport(cl=cl, varlist=c("fpscore_matrix_permuted"))


alt_pvals <- list()
for (link in c("logit", "probit", "cloglog", "linear")) {
  cat(sprintf("Running spike-in simulations for link=%s\n", link))
  
  regression_function <- ifelse(link == "linear", simulate_linear_regression_on_variant, simulate_beta_regression_on_variant)

  alt_pvals[[link]] <- data.frame()
  for (spike_in_effect in c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12)) {
    cat(sprintf("\tspike-in-effect = %.2f\n", spike_in_effect))

    regression_result <- parSapply(cl, 1:N_PERMUTATIONS, regression_function, link=link, spike_in_effect=spike_in_effect)
    df <- data.frame(p = regression_result["pval", ], 
                     beta = regression_result["beta", ],
                     spike_in_effect=spike_in_effect, 
                     row.names = NULL)
    alt_pvals[[link]] <- rbind(alt_pvals[[link]], df)
  }
}
stopCluster(cl)

cat("Saving alt pvals...\n")
alt_pvals %>% saveRDS("simulations/alt_pvals.Rds")


