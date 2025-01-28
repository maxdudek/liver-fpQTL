suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(GWASTools)
  library(qqman)
  library(ggpubr)
})

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]
DIR <- paste0("FP_methods/", FP_METHOD)
cat("FP_METHOD = ", FP_METHOD, "\n")

cat("Loading data...\n")
fpscore_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_genotype_regression.Rds"))
fpscore_cov_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))
variant_info <- readRDS(paste0(DIR, "/data/variant_info_regression.Rds"))

# For beta regression, join with linear effect sizes
if (grepl("PRINT_beta", FP_METHOD, fixed = TRUE)) {
  linear_regression <- readRDS("FP_methods/PRINT_no_gaussian/regression_results/fp_score_genotype_regression.Rds")
  linear_cov_regression <- readRDS("FP_methods/PRINT_no_gaussian/regression_results/fp_score_covariates_genotype_regression.Rds")
  
  fpscore_regression$linear_beta <- linear_regression$beta
  fpscore_cov_regression$linear_beta <- linear_cov_regression$beta
}

# Save smaller FP score/genotype matrix with only sig variants
# cat("Loading matrices...\n")
# fpscore_matrix <- readRDS(paste0(DIR, "/data/fpscore_matrix_regression.Rds"))
# fpscore_matrix_normalized <- readRDS(paste0(DIR, "/data/fpscore_matrix_normalized.Rds"))
# genotype_matrix <- readRDS(paste0(DIR, "/data/genotype_matrix_regression.Rds"))
# 
# cat("Saving smaller matrices...\n")
# fpQTLs <- fpscore_regression %>% filter(ST_qval < 0.05) %>% pull(variant_id)
# fpscore_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_fpQTLs_no_covariates.Rds"))
# fpscore_matrix_normalized[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_normalized_fpQTLs_no_covariates.Rds"))
# genotype_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/genotype_matrix_fpQTLs_no_covariates.Rds"))
# 
# fpQTLs <- fpscore_cov_regression %>% filter(ST_qval < 0.05) %>% pull(variant_id)
# fpscore_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_fpQTLs_with_covariates.Rds"))
# fpscore_matrix_normalized[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/fpscore_matrix_normalized_fpQTLs_with_covariates.Rds"))
# genotype_matrix[fpQTLs, ] %>% saveRDS(paste0(DIR, "/data/genotype_matrix_fpQTLs_with_covariates.Rds"))
# 
# cat("Freeing memory...\n")
# rm(fpscore_matrix)
# rm(fpscore_matrix_normalized)
# rm(genotype_matrix)

# Calculate p-val thresholds for FDR 5
fpscore_regression %>%
  filter(ST_qval < 0.05) %>%
  pull(pval) %>%
  max() ->
  ST_FDR5_p_thresh_no_cov

fpscore_cov_regression %>%
  filter(ST_qval < 0.05) %>%
  pull(pval) %>%
  max() ->
  ST_FDR5_p_thresh_with_cov

# cat("Writing FDR variant tables...\n")
# get_FDR_fpQTLs <- function(df, fdr) {
#   df %>%
#     filter(ST_qval < fdr) %>%
#     dplyr::select(-c(ST_lfdr)) %>%
#     arrange(ST_qval)
# }
# 
# # Write FDR5 variant table
# fpscore_regression %>%
#   get_FDR_fpQTLs(0.05) %>%
#   write.table(file = paste0(DIR, "/regression_results/fpQTLs_FDR5.txt"),
#               quote = F, row.names = F, sep = "\t")
# 
# fpscore_cov_regression %>%
#   get_FDR_fpQTLs(0.05) %>%
#   write.table(file = paste0(DIR, "/regression_results/fpQTLs_covariates_FDR5.txt"),
#               quote = F, row.names = F, sep = "\t")
# 
# # Write FDR10 variant table
# fpscore_regression %>%
#   get_FDR_fpQTLs(0.1) %>%
#   write.table(file = paste0(DIR, "/regression_results/fpQTLs_FDR10.txt"),
#               quote = F, row.names = F, sep = "\t")
# 
# fpscore_cov_regression %>%
#   get_FDR_fpQTLs(0.1) %>%
#   write.table(file = paste0(DIR, "/regression_results/fpQTLs_covariates_FDR10.txt"),
#               quote = F, row.names = F, sep = "\t")
# 
# # Write FDR1 variant table
# fpscore_regression %>%
#   get_FDR_fpQTLs(0.01) %>%
#   write.table(file = paste0(DIR, "/regression_results/fpQTLs_FDR1.txt"),
#               quote = F, row.names = F, sep = "\t")
# 
# fpscore_cov_regression %>%
#   get_FDR_fpQTLs(0.01) %>%
#   write.table(file = paste0(DIR, "/regression_results/fpQTLs_covariates_FDR1.txt"),
#               quote = F, row.names = F, sep = "\t")
# 
# # Write rsID files
# fpscore_cov_regression %>%
#   filter(ST_qval < 0.05) %>%
#   pull(variant_id) %>%
#   write(paste0(DIR, "/regression_results/fpQTLs_covariates_FDR5_rsIDs.txt"))
# 
# fpscore_cov_regression %>%
#   filter(ST_qval < 0.1) %>%
#   pull(variant_id) %>%
#   write(paste0(DIR, "/regression_results/fpQTLs_covariates_FDR10_rsIDs.txt"))
# 
# fpscore_cov_regression %>%
#   filter(ST_qval < 0.01) %>%
#   pull(variant_id) %>%
#   write(paste0(DIR, "/regression_results/fpQTLs_covariates_FDR1_rsIDs.txt"))
# 
# cat("Plotting pval distributions...\n")
# fpscore_regression %>%
#   ggplot(aes(x = pval)) +
#   xlim(c(0,1)) +
#   geom_histogram(bins = 50) +
#   theme_classic()
# 
# ggsave(paste0(DIR, "/figures/regression_results_no_covariates/regression_pval_distribution.png"), create.dir = TRUE)
# 
# fpscore_cov_regression %>%
#   ggplot(aes(x = pval)) +
#   xlim(c(0,1)) +
#   geom_histogram(bins = 50) +
#   theme_classic()
# 
# ggsave(paste0(DIR, "/figures/regression_results_with_covariates/regression_pval_distribution.png"), create.dir = TRUE)
# 
# fpscore_regression %>%
#   ggplot(aes(x = ST_qval)) +
#   xlim(c(0,1)) +
#   geom_histogram(bins = 50) +
#   theme_classic()
# 
# ggsave(paste0(DIR, "/figures/regression_results_no_covariates/regression_ST_qval_distribution.png"), create.dir = TRUE)
# 
# fpscore_cov_regression %>%
#   ggplot(aes(x = ST_qval)) +
#   xlim(c(0,1)) +
#   geom_histogram(bins = 50) +
#   theme_classic()
# 
# ggsave(paste0(DIR, "/figures/regression_results_with_covariates/regression_ST_qval_distribution.png"), create.dir = TRUE)
# 
# 
# cat("Plotting QQ plots...\n")
# inflation <- function(ps) {
#   chisq <- qchisq(1 - ps, 1)
#   lambda <- median(chisq) / qchisq(0.5, 1)
#   lambda
# }
# # QQ plot
# ggplot() +
#   stat_qq(aes(sample = -log10(fpscore_regression$pval)),
#           color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) +
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
#   ylab(bquote(Observed ~~ -log[10](p))) +
#   xlab(bquote(Expected ~~ -log[10](p))) +
#   theme_classic()
# ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_no_covariates/qqplot.png"),
#        width = 5, height = 5, units = "in", create.dir = TRUE)
# lambda_no_cov <- inflation(fpscore_regression$pval)
# print(lambda_no_cov)
# 
# # QQ plot
# ggplot() +
#   stat_qq(aes(sample = -log10(fpscore_cov_regression$pval)),
#           color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) +
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
#   ylab(bquote(Observed ~~ -log[10](p))) +
#   xlab(bquote(Expected ~~ -log[10](p))) +
#   theme_classic()
# ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_with_covariates/qqplot.png"),
#        width = 5, height = 5, units = "in", create.dir = TRUE)
# lambda_with_cov <- inflation(fpscore_cov_regression$pval)
# print(lambda_with_cov)
# 
# # QQ plot with Genomic Control
# ggplot() +
#   stat_qq(aes(sample = -log10(fpscore_regression$pval / lambda_no_cov)),
#           color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) +
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
#   ylab(bquote(Observed ~~ -log[10](p))) +
#   xlab(bquote(Expected ~~ -log[10](p))) +
#   theme_classic()
# ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_no_covariates/qqplot_GC.png"),
#        width = 5, height = 5, units = "in")
# 
# ggplot() +
#   stat_qq(aes(sample = -log10(fpscore_cov_regression$pval / lambda_with_cov)),
#           color = "black", distribution = stats::qexp, dparams = list(rate = log(10))) +
#   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
#   ylab(bquote(Observed ~~ -log[10](p))) +
#   xlab(bquote(Expected ~~ -log[10](p))) +
#   theme_classic()
# ggsave(paste0("FP_methods/", FP_METHOD, "/figures/regression_results_with_covariates/qqplot_GC.png"),
#        width = 5, height = 5, units = "in")
# 


cat("Plotting Volcano plots...\n")
# Volcano plot
xlim <- max(max(fpscore_regression$beta), -min(fpscore_regression$beta),
            max(fpscore_cov_regression$beta), -min(fpscore_cov_regression$beta)) 
volcano_plot_xlim <- c(-xlim, xlim)

fpscore_regression %>%
  ggplot(aes(x = beta, y = -log10(pval))) +
  geom_point(size = 1) +
  geom_abline(slope = 0, intercept = -log10(ST_FDR5_p_thresh_no_cov), size = 0.3, color = "black") +
  xlim(volcano_plot_xlim) +
  ylab(bquote(-log[10](p))) +
  xlab("Effect size (beta)") +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_no_covariates/volcano_plot.png"),
       width = 3, height = 3, dpi = 600, create.dir = TRUE)

fpscore_cov_regression %>%
  ggplot(aes(x = beta, y = -log10(pval))) +
  geom_point(size = 0.6) +
  geom_abline(slope = 0, intercept = -log10(ST_FDR5_p_thresh_with_cov), size = 0.3, color = "black") +
  xlim(volcano_plot_xlim) +
  ylab(bquote(-log[10](p))) +
  xlab("Effect size (beta)") +
  theme_classic()

ggsave(paste0(DIR, "/figures/regression_results_with_covariates/volcano_plot.png"),
       width = 3, height = 3, dpi = 600, create.dir = TRUE)

xlim <- max(max(fpscore_regression$linear_beta), -min(fpscore_regression$linear_beta),
            max(fpscore_cov_regression$linear_beta), -min(fpscore_cov_regression$linear_beta)) 
volcano_plot_xlim <- c(-xlim, xlim)

# Volcano plots with linear effect
if (grepl("PRINT_beta", FP_METHOD, fixed = TRUE)) {
  fpscore_regression %>%
    ggplot(aes(x = linear_beta, y = -log10(pval))) +
    geom_point(size = 1) +
    geom_abline(slope = 0, intercept = -log10(ST_FDR5_p_thresh_no_cov), size = 0.3, color = "black") +
    xlim(volcano_plot_xlim) +
    ylab(bquote(-log[10](p))) +
    xlab("Linear effect (beta)") +
    theme_classic()
  
  ggsave(paste0(DIR, "/figures/regression_results_no_covariates/volcano_plot_linear_effect.png"),
         width = 3, height = 3, dpi = 600, create.dir = TRUE)
  
  fpscore_cov_regression %>%
    ggplot(aes(x = linear_beta, y = -log10(pval))) +
    geom_point(size = 0.6) +
    geom_abline(slope = 0, intercept = -log10(ST_FDR5_p_thresh_with_cov), size = 0.3, color = "black") +
    xlim(volcano_plot_xlim) +
    ylab(bquote(-log[10](p))) +
    xlab("Linear effect (beta)") +
    theme_classic()
  
  ggsave(paste0(DIR, "/figures/regression_results_with_covariates/volcano_plot_linear_effect.png"),
         width = 3, height = 3, dpi = 600, create.dir = TRUE)
}

# cat("Plotting Manhattan plots...\n")
# # Manhatten plot
# png(file=paste0(DIR, "/figures/regression_results_no_covariates/manhatten_plot.png"),
#     width = 13, height = 6.5, units = "in", res = 600)
# fpscore_regression %>%
#   mutate(variant_chrom = as.integer(gsub("chr", "", variant_chrom))) %>%
#   manhattan(chr = "variant_chrom", bp = "variant_pos",
#             snp = "variant_id", p = "pval",
#             col = c("#dd960f", "#526cb0"), cex = 0.8, cex.axis = 0.9,
#             suggestiveline = FALSE,
#             genomewideline = -log10(ST_FDR5_p_thresh_no_cov))
# dev.off()
# 
# png(file=paste0(DIR, "/figures/regression_results_with_covariates/manhatten_plot.png"),
#     width = 13, height = 6.5, units = "in", res = 600)
# fpscore_cov_regression %>%
#   mutate(variant_chrom = as.integer(gsub("chr", "", variant_chrom))) %>%
#   manhattan(chr = "variant_chrom", bp = "variant_pos",
#             snp = "variant_id", p = "pval",
#             col = c("#dd960f", "#526cb0"), cex = 0.8, cex.axis = 0.9,
#             suggestiveline = FALSE,
#             genomewideline = -log10(ST_FDR5_p_thresh_with_cov))
# dev.off()
# 
# # cat("Comparing regression with and without covariates...\n")
# data.frame(
#   pval_no_cov = -log10(fpscore_regression$pval),
#   pval_with_cov = -log10(fpscore_cov_regression$pval)
# ) %>%
#   ggplot(aes(x = pval_no_cov, y = pval_with_cov)) +
#   geom_point(size = 0.5) +
#   geom_abline(slope = 1, intercept = 0, color = "blue") +
#   stat_cor(
#     aes(label = paste(after_stat(rr.label), ..p.label.., sep = "~`,`~")),
#     size = 3
#   ) +
#   xlab(bquote(-log[10](p)~"without covariates")) +
#   ylab(bquote(-log[10](p)~"WITH covariates")) +
#   theme_classic()
# ggsave(paste0(DIR, "/figures/regression_pvals_with_vs_without_covariates.png"),
#        width = 3, height = 3, dpi = 300, create.dir = TRUE)
# 
# cat("Correlation between pvals with and without covariates...\n")
# cor.test(-log10(fpscore_regression$pval), -log10(fpscore_cov_regression$pval),
#          method = "spearman")
# cor.test(-log10(fpscore_regression$pval), -log10(fpscore_cov_regression$pval),
#          method = "pearson")
# 



