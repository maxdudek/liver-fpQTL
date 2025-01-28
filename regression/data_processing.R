suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(preprocessCore)
})

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]
# FP_METHOD = "PRINT_local_cutsites_geq10"
DIR <- paste0("FP_methods/", FP_METHOD)
cat("FP_METHOD = ", FP_METHOD, "\n")

cat("Loading data...\n")
fpscore_matrix <- readRDS(paste0(DIR, "/data/fpscore_matrix.Rds"))
genotype_matrix <- readRDS("genotype_matrix.Rds")
variant_info <- readRDS("variant_info.Rds")

# If fpscore_matrix has less samples than genotype_matrix, correct this
if (ncol(fpscore_matrix) < ncol(genotype_matrix)) {
  cat("Shortening genotype matrix because genotype matrix has",
      ncol(genotype_matrix), "samples while fpscore matrix has",
      ncol(fpscore_matrix), "samples\n")
      genotype_matrix <- genotype_matrix[,colnames(fpscore_matrix)]
}

cat("Number of starting variants:", nrow(variant_info), "\n")

cat("Verifying order of samples in data tables...\n")
all(rownames(fpscore_matrix) == variant_info$variant_id)
all(rownames(fpscore_matrix) == rownames(genotype_matrix))
all(colnames(fpscore_matrix) == colnames(genotype_matrix))

cat("Plotting missing data figure...\n")
# How many samples does each variant have data for?
rowSums(!is.na(fpscore_matrix)) ->
  numSamplesPerVariantWithData

numSamplesPerVariantWithData %>%
  data.frame() %>%
  ggplot(aes(x = .)) +
  geom_bar(color = "black") +
  xlab("Number of samples with data") +
  labs(title = "How many samples does each variant have data for?") +
  theme_classic()

ggsave(paste0(DIR, "/figures/missing_data/num_samples_per_variant_with_data.png"))

# Match missing FP score data in genotype matrix
cat("Finding variants with missing data...\n")
genotype_matrix[is.na(fpscore_matrix)] <- NA

variant_info$N_regression_genotype0 <- rowSums(genotype_matrix == 0, na.rm = TRUE)
variant_info$N_regression_genotype1 <- rowSums(genotype_matrix == 1, na.rm = TRUE)
variant_info$N_regression_genotype2 <- rowSums(genotype_matrix == 2, na.rm = TRUE)

variant_info <- variant_info %>%
  mutate(
    N_regression = N_regression_genotype0 + N_regression_genotype1 + N_regression_genotype2,
    num_genotypes = as.integer(N_regression_genotype0 != 0) + 
                    as.integer(N_regression_genotype1 != 0) +
                    as.integer(N_regression_genotype2 != 0)
  ) %>%
  filter(num_genotypes >= 2 & N_regression >= 10)
cat("Number of variants after filtering missing data:", nrow(variant_info), "\n")

# Filter variants in MHC
MHC <- c(28510120, 33480577)
variant_info <- variant_info %>%
  filter(!(variant_chrom == "chr6" & variant_pos > MHC[1] & variant_pos < MHC[2])) 
cat("Number of variants after filtering MHC:", nrow(variant_info), "\n")


cat("Removing variants with missing data...\n")
fpscore_matrix <- fpscore_matrix[variant_info$variant_id,]
genotype_matrix <- genotype_matrix[variant_info$variant_id,]

cat("Verifying order of variants and samples in data tables...\n")
all(rownames(fpscore_matrix) == variant_info$variant_id)
all(rownames(fpscore_matrix) == rownames(genotype_matrix))
all(colnames(fpscore_matrix) == colnames(genotype_matrix))

cat("5x5 slice of matrices...\n")
fpscore_matrix[1:5, 1:5]
genotype_matrix[1:5, 1:5]

cat("dim(fpscore_matrix) = \n")
dim(fpscore_matrix)

cat("Saving genotype_matrix and variant_info for regression...\n")
fpscore_matrix %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_regression.Rds"))
genotype_matrix %>% saveRDS(paste0(DIR, "/data/genotype_matrix_regression.Rds"))
variant_info %>% saveRDS(paste0(DIR, "/data/variant_info_regression.Rds"))

# cat("Loading data...\n")
# fpscore_matrix <- readRDS(paste0(DIR, "/data/fpscore_matrix_regression.Rds"))
# genotype_matrix <- readRDS(paste0(DIR, "/data/genotype_matrix_regression.Rds"))
# variant_info <- readRDS(paste0(DIR, "/data/variant_info_regression.Rds"))

# Assess missing data
colSums( !is.na( fpscore_matrix )) %>%
  data.frame(num_variants_with_data = .) %>%
  arrange(num_variants_with_data) ->
  sample_completeness

cat("Sample completeness...\n")
sample_completeness

sample_completeness %>%
  ggplot(aes(x = num_variants_with_data)) +
  geom_histogram(bins=15, fill = "black") +
  xlab("Number of variants with data") +
  labs(title = "How many variants does each sample have data for?") +
  theme_classic()
ggsave(paste0(DIR, "/figures/missing_data/num_variants_per_sample_with_data.png"))

rowSums( !is.na( fpscore_matrix )) %>%
  data.frame(
    variant_id = names(.),
    num_samples_with_data = .
  ) %>%
  cbind(variant_info %>% select(-variant_id)) ->
  variant_completeness

cat("Verifying that variant_ids line up in variant_completeness:\n")
all(variant_completeness$variand_id == variant_info$variant_id)

variant_completeness %>%
  ggplot(aes(x = num_samples_with_data)) +
  geom_bar(color = "black") +
  xlab("Number of samples with data") +
  labs(title = "How many samples does each variant have data for?") +
  theme_classic()
ggsave(paste0(DIR, "/figures/missing_data/num_samples_per_variant_with_data2.png"))

# Are we missing data for variants uniformly across the chromosome?
variant_completeness %>%
  filter(variant_chrom == "chr21") %>%
  ggplot(aes(x = variant_pos, y = num_samples_with_data)) +
  geom_point(size = 0.1) +
  theme_classic()
ggsave(paste0(DIR, "/figures/missing_data/chr21_map_of_variant_completeness.png"))

# Plot FP score distributions
cat("Plotting FP score distributions...\n")
fpscore_matrix %>%
  melt() %>%
  na.omit() %>%
  rename(
    variant_id = Var1,
    sample_id = Var2,
    fp_score = value
  ) ->
  fpscore_rows_df

if (grepl("PRINT", FP_METHOD, fixed = TRUE)) {
  fpscore.xlim <- c(0,1)
} else {
  fpscore.xlim <- c(0,0.1)
}



fpscore_rows_df %>%
  ggplot(aes(x = fp_score)) +
  geom_histogram(bins = 200) +
  xlim(fpscore.xlim) +
  theme_classic()

ggsave(paste0(DIR, "/figures/all_fp_score_distribution.png"))


fpscore_rows_df %>%
  ggplot(aes(x = fp_score)) +
  geom_histogram(bins = 100) +
  xlim(fpscore.xlim) +
  facet_wrap(~sample_id) +
  theme_classic()

ggsave(paste0(DIR, "/figures/sample_fp_score_distributions_raw.png"), width = 16, height = 12)

variants_to_display <- sample(variant_info$variant_id, 100)

fpscore_rows_df %>%
  filter(variant_id %in% variants_to_display) %>%
  ggplot(aes(x = fp_score)) +
  geom_histogram(bins = 20) +
  xlim(fpscore.xlim) +
  facet_wrap(~variant_id, nrow = 10) +
  theme_classic()

ggsave(paste0(DIR, "/figures/variant_fp_score_distributions_raw.png"), width = 16, height = 12)

######################
# Normalization      #
######################

cat("Quantile normalize sample distributions...\n")
fpscore_matrix %>%
  normalize.quantiles() ->
  fpscore_matrix_normalized


# Get intermediate figure
fpscore_matrix_normalized %>%
  melt() %>%
  na.omit() %>%
  rename(
    variant_id = Var1,
    sample_id = Var2,
    fp_score = value
  ) %>%
  ggplot(aes(x = fp_score)) +
  geom_histogram(bins = 100) +
  xlim(fpscore.xlim) +
  facet_wrap(~sample_id) +
  theme_classic()

ggsave(paste0(DIR, "/figures/sample_fp_score_distributions_quantile_normalized.png"), width = 16, height = 12)


# Quantile normalize variant distributions to standard normal
cat("Normalize variant distributions...\n")
normalize.quantiles.snorm <- function(x) {
  ((rank(x, na.last = "keep") - 0.5)/sum(!is.na(x))) %>% qnorm()
}

GAUSSIAN_NORMALIZE <- !grepl("no_gaussian", FP_METHOD, fixed = TRUE)

if (GAUSSIAN_NORMALIZE) {
  # Transform variant distributions to standard normal
  cat("Running gaussian transformation...\n")
  fpscore_matrix_normalized %>% 
    apply(1, normalize.quantiles.snorm) %>%
    t() -> 
    fpscore_matrix_normalized
} else {
  cat("Skipping gaussian transformation...\n")
}

rownames(fpscore_matrix_normalized) <- rownames(fpscore_matrix)
colnames(fpscore_matrix_normalized) <- colnames(fpscore_matrix)

cat("Write matrices...\n")
fpscore_matrix_normalized %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_normalized.Rds"))

if (GAUSSIAN_NORMALIZE) {
  cat("Plot variant and sample distributions...\n")
  fpscore_matrix_normalized %>%
    melt() %>%
    na.omit() %>%
    rename(
      variant_id = Var1,
      sample_id = Var2,
      fp_score = value
    ) ->
    fpscore_rows_df_normalized

  fpscore_rows_df_normalized %>%
    ggplot(aes(x = fp_score)) +
    geom_histogram(bins = 100) +
    xlim(c(-3,3)) +
    facet_wrap(~sample_id) +
    theme_classic()

  ggsave(paste0(DIR, "/figures/sample_fp_score_distributions_after_variant_normalization.png"), width = 16, height = 12)

  variants <- rownames(fpscore_matrix_normalized)
  variants_to_display <- sample(variants, 100)

  fpscore_rows_df_normalized %>%
    filter(variant_id %in% variants_to_display) %>%
    ggplot(aes(x = fp_score)) +
    geom_histogram(bins = 20) +
    xlim(c(-3,3)) +
    facet_wrap(~variant_id, nrow = 10) +
    theme_classic()

  ggsave(paste0(DIR, "/figures/variant_fp_score_distributions_standard_normal.png"), width = 16, height = 12)
}

cat("Done\n")

