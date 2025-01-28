suppressPackageStartupMessages({
    library(tidyverse)
})

ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

cat("Loading traits...\n")
GWAS_DIR <- paste0(ROOT_DIR, "/raw_data/GWAS/hg38_liftover/")
# sumstats_hg38 <- readRDS(paste0(GWAS_DIR, "liver_sumstats_hg38.Rds"))
# eqtl_all_vars <- readRDS(paste0(ROOT_DIR, "PMACS/bulk_human_liver_footprints/functional_variants/eQTL/GTEx_Analysis_v8_eQTL/GTEx_v8_Liver_all_vars.Rds"))
eqtl_sig_pairs <- read.delim(paste0(ROOT_DIR, "raw_data/GTEx/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt"))
gtex_all_vars_table <- read.delim(paste0(ROOT_DIR, "raw_data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.location.txt"))

GTEx_tested_SNPs <- gtex_all_vars_table %>%
  mutate(variant_location = sprintf("%s:%d", chr, variant_pos)) %>%
  pull(variant_location)

eQTL_SNPs <- eqtl_sig_pairs %>%
  separate(variant_id, into = c("chrom", "pos", NA, NA, NA), sep = "_") %>%
  mutate(variant_location = sprintf("%s:%s", chrom, pos)) %>%
  pull(variant_location)

# caQTL_tests <- readRDS(paste0(ROOT_DIR, "PMACS/bulk_human_liver_footprints/functional_variants/caQTL/caQTL_all_vars.Rds"))
caQTL_Q_threshold <- 0.0009172353
caPeaks_all_tests <- read.delim(paste0(ROOT_DIR, "raw_data/brandon_liver_ATAC/caQTL/fdr5_caqtl_results_all_tests.txt"))
caQTL_SNPs <- caPeaks_all_tests %>% 
  filter(10^Log_10_BH_Qvalue < caQTL_Q_threshold) %>% 
  pull(rs_ID) %>% unique()


FP_METHODS <- c("PRINT_beta_cloglog")

result_df <- data.frame()
for (FP_METHOD in FP_METHODS) {
    for (FDR in c("5")) {
      cat("FP_METHOD = ", FP_METHOD, "\n")
      DIR <- paste0("../../regression/FP_methods/", FP_METHOD)
      
      cat("\tLoading regression results...\n")
      fpscore_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_genotype_regression.Rds"))
      fpscore_cov_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))
  
      if (FDR == "10") {
          Q_THRESHOLD <- 0.10
          FP_METHOD_NAME <- paste0(FP_METHOD, "_FDR10")
      } else {
          Q_THRESHOLD <- 0.05
          FP_METHOD_NAME <- FP_METHOD
      }
  
      fpQTL_data <- data.frame(
          rsID = fpscore_regression$variant_id,
          variant_chrom = fpscore_regression$variant_chrom,
          variant_pos = fpscore_regression$variant_pos,
          fpQTL_no_cov = fpscore_regression$ST_qval < Q_THRESHOLD,
          fpQTL_with_cov = fpscore_cov_regression$ST_qval < Q_THRESHOLD
      ) %>%
        separate(rsID, into = c("rsID", NA), sep = "_", fill = "right") # "rsXXXXXX_2" --> "rsXXXXXX"
        
  
      # cat("\tMaking 2x2 tables for GWAS traits...\n")
      # for (trait in names(sumstats_hg38)) {
      #     cat("\t\tGWAS trait = ", trait, "\n")
      #     sumstats_hg38[[trait]] %>%
      #         mutate(sig_association = P < 5e-8) %>%
      #         inner_join(fpQTL_data, by = join_by(rsID)) %>%
      #         select(rsID, sig_association, fpQTL_no_cov, fpQTL_with_cov) ->
      #         joined_df
      #     
      #     # Without covariates
      #     table2x2 <- table(joined_df$sig_association, joined_df$fpQTL_no_cov)
      #     rownames(table2x2) <- c("no sig", "sig association")
      #     colnames(table2x2) <- c("non-fpQTL", "fpQTL")
      #     ft <- fisher.test(table2x2)
      # 
      #     df <- data.frame(
      #         fp_method = FP_METHOD_NAME,
      #         trait = trait,
      #         covariates = FALSE,
      #         snps_neither = table2x2[1,1],
      #         snps_fpQTL_only = table2x2[1,2],
      #         snps_sig_only = table2x2[2,1],
      #         snps_sig_and_fpQTL = table2x2[2,2],
      #         odds_ratio = ft$estimate,
      #         odds_ratio_95_low = ft$conf.int[1],
      #         odds_ratio_95_high = ft$conf.int[2],
      #         fishers_p = ft$p.value,
      #         row.names = NULL
      #     )
      #     result_df <- rbind(result_df, df)
      # 
      #     # With covariates
      #     table2x2 <- table(joined_df$sig_association, joined_df$fpQTL_with_cov)
      #     rownames(table2x2) <- c("no sig", "sig association")
      #     colnames(table2x2) <- c("non-fpQTL", "fpQTL")
      #     ft <- fisher.test(table2x2)
      # 
      #     df <- data.frame(
      #         fp_method = FP_METHOD_NAME,
      #         trait = trait,
      #         covariates = TRUE,
      #         snps_neither = table2x2[1,1],
      #         snps_fpQTL_only = table2x2[1,2],
      #         snps_sig_only = table2x2[2,1],
      #         snps_sig_and_fpQTL = table2x2[2,2],
      #         odds_ratio = ft$estimate,
      #         odds_ratio_95_low = ft$conf.int[1],
      #         odds_ratio_95_high = ft$conf.int[2],
      #         fishers_p = ft$p.value,
      #         row.names = NULL
      #     )
      #     result_df <- rbind(result_df, df)
      # }
  
      cat("\tMaking 2x2 tables for eQTLs...\n")
      all_variants_eQTL <- fpQTL_data %>%
        mutate(variant_location = sprintf("%s:%d", variant_chrom, variant_pos)) %>%
        filter(variant_location %in% GTEx_tested_SNPs) %>%
        mutate(eQTL = variant_location %in% eQTL_SNPs)
      
      cat(sprintf("\tNumber of variants tested in GTEx = %d/%d\n", nrow(all_variants_eQTL), nrow(fpQTL_data)))
  
      # Without covariates
      table2x2 <- table(all_variants_eQTL$eQTL, all_variants_eQTL$fpQTL_no_cov)
      rownames(table2x2) <- c("no sig", "sig association")
      colnames(table2x2) <- c("non-fpQTL", "fpQTL")
      ft <- fisher.test(table2x2)
  
      df <- data.frame(
          fp_method = FP_METHOD_NAME,
          trait = "liver_eQTL",
          covariates = FALSE,
          snps_neither = table2x2[1,1],
          snps_fpQTL_only = table2x2[1,2],
          snps_sig_only = table2x2[2,1],
          snps_sig_and_fpQTL = table2x2[2,2],
          odds_ratio = ft$estimate,
          odds_ratio_95_low = ft$conf.int[1],
          odds_ratio_95_high = ft$conf.int[2],
          fishers_p = ft$p.value,
          row.names = NULL
      )
      result_df <- rbind(result_df, df)
  
      # With covariates
      table2x2 <- table(all_variants_eQTL$eQTL, all_variants_eQTL$fpQTL_with_cov)
      rownames(table2x2) <- c("no sig", "sig association")
      colnames(table2x2) <- c("non-fpQTL", "fpQTL")
      ft <- fisher.test(table2x2)
  
      df <- data.frame(
          fp_method = FP_METHOD_NAME,
          trait = "liver_eQTL",
          covariates = TRUE,
          snps_neither = table2x2[1,1],
          snps_fpQTL_only = table2x2[1,2],
          snps_sig_only = table2x2[2,1],
          snps_sig_and_fpQTL = table2x2[2,2],
          odds_ratio = ft$estimate,
          odds_ratio_95_low = ft$conf.int[1],
          odds_ratio_95_high = ft$conf.int[2],
          fishers_p = ft$p.value,
          row.names = NULL
      )
      result_df <- rbind(result_df, df)
  
      cat("\tMaking 2x2 tables for caQTLs...\n")
      all_variants_caQTL <- fpQTL_data %>%
        mutate(caQTL = rsID %in% caQTL_SNPs)
  
      # Without covariates
      table2x2 <- table(all_variants_caQTL$caQTL, all_variants_caQTL$fpQTL_no_cov)
      rownames(table2x2) <- c("no sig", "sig association")
      colnames(table2x2) <- c("non-fpQTL", "fpQTL")
      ft <- fisher.test(table2x2)
  
      df <- data.frame(
          fp_method = FP_METHOD_NAME,
          trait = "liver_caQTL",
          covariates = FALSE,
          snps_neither = table2x2[1,1],
          snps_fpQTL_only = table2x2[1,2],
          snps_sig_only = table2x2[2,1],
          snps_sig_and_fpQTL = table2x2[2,2],
          odds_ratio = ft$estimate,
          odds_ratio_95_low = ft$conf.int[1],
          odds_ratio_95_high = ft$conf.int[2],
          fishers_p = ft$p.value,
          row.names = NULL
      )
      result_df <- rbind(result_df, df)
  
      # With covariates
      table2x2 <- table(all_variants_caQTL$caQTL, all_variants_caQTL$fpQTL_with_cov)
      rownames(table2x2) <- c("no sig", "sig association")
      colnames(table2x2) <- c("non-fpQTL", "fpQTL")
      ft <- fisher.test(table2x2)
  
      df <- data.frame(
          fp_method = FP_METHOD_NAME,
          trait = "liver_caQTL",
          covariates = TRUE,
          snps_neither = table2x2[1,1],
          snps_fpQTL_only = table2x2[1,2],
          snps_sig_only = table2x2[2,1],
          snps_sig_and_fpQTL = table2x2[2,2],
          odds_ratio = ft$estimate,
          odds_ratio_95_low = ft$conf.int[1],
          odds_ratio_95_high = ft$conf.int[2],
          fishers_p = ft$p.value,
          row.names = NULL
      )
      result_df <- rbind(result_df, df)
    }
  
  cat("Writing results...\n")
  result_df %>%
    filter(trait == "liver_eQTL" | trait == "liver_caQTL") %>%
    write.table("tables/%s/%s_tables_2x2_QTLs.txt" %>% sprintf(FP_METHOD, FP_METHOD),
                quote = FALSE, row.names = F, sep = "\t")
}

# result_df %>%
#   filter(trait != "liver_eQTL" & trait != "liver_caQTL") %>%
#   write.table("tables/tables_2x2_GWAS_GWsignificant.txt",
#               quote = FALSE, row.names = F, sep = "\t")

