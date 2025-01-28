
library(tidyverse)
ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

FP_METHOD <- "PRINT_no_gaussian"
DIR <- paste0("../../regression/FP_methods/", FP_METHOD)

cat("Loading data...\n")
fpscore_cov_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))
brandon_peaks <- read.table(paste0(ROOT_DIR, "raw_data/brandon_liver_ATAC/peak_info/genrichAllPeaks_m10_g50_9.15.21.noBLnarrowpeak_chr.bed"))

# Zooming in on the spike in chr 17

cat("Full...\n")
zoom1 <- c(4.53e7, 4.66e7)
fpscore_cov_regression %>%
  filter(variant_chrom == "chr17") %>%
  filter(variant_pos > 0) %>%
  ggplot(aes(x = variant_pos, y = -log10(pval))) +
  geom_vline(xintercept = zoom1[1], color = "red", linetype = "dashed") +
  geom_vline(xintercept = zoom1[2], color = "red", linetype = "dashed") +
  geom_point() + 
  theme_classic()
ggsave("figures/chr17_full.png", width = 11, height = 7)

peaks <- brandon_peaks %>%
  filter(V1 == "chr17") %>%
  filter(V2 > zoom1[1] & V2 < zoom1[2])

cat("Zoom1...\n")
zoom2 <- c(4.62e7, 4.63e7)
fpscore_cov_regression %>%
  filter(variant_chrom == "chr17") %>%
  filter(variant_pos > zoom1[1] & variant_pos < zoom1[2]) %>%
  ggplot(aes(x = variant_pos, y = -log10(pval))) +
  geom_vline(xintercept = zoom2[1], color = "red", linetype = "dashed") +
  geom_vline(xintercept = zoom2[2], color = "red", linetype = "dashed") +
  geom_point() + 
  geom_segment(data = peaks, aes(x = V2, xend = V3, y = -0.5, yend = -0.5), size = 2) + 
  xlim(zoom1) +
  theme_classic()
ggsave("figures/chr17_zoom1.png", width = 11, height = 7)

cat("Zoom2...\n")
zoom3 <- c(4.6259e7, 4.62619e7)
fpscore_cov_regression %>%
  filter(variant_chrom == "chr17") %>%
  filter(variant_pos > zoom2[1] & variant_pos < zoom2[2]) %>%
  ggplot(aes(x = variant_pos, y = -log10(pval))) +
  geom_vline(xintercept = zoom3[1], color = "red", linetype = "dashed") +
  geom_vline(xintercept = zoom3[2], color = "red", linetype = "dashed") +
  geom_point() + 
  geom_segment(data = peaks, aes(x = V2, xend = V3, y = -0.5, yend = -0.5), size = 2) + 
  xlim(zoom2) +
  theme_classic()
ggsave("figures/chr17_zoom2.png", width = 11, height = 7)

cat("Zoom3...\n")
fpscore_cov_regression %>%
  filter(variant_chrom == "chr17") %>%
  filter(variant_pos > zoom3[1] & variant_pos < zoom3[2]) %>%
  ggplot(aes(x = variant_pos, y = -log10(pval))) +
  geom_point() + 
  geom_segment(data = peaks, aes(x = V2, xend = V3, y = -0.5, yend = -0.5), size = 2) + 
  xlim(zoom3) +
  theme_classic()
ggsave("figures/chr17_zoom3.png", width = 11, height = 7)


