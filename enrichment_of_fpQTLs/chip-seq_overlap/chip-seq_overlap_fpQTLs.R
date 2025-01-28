suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicRanges)
})

ROOT_DIR <- "/mnt/isilon/sfgi/dudekm"

cat("Loading data...\n")
chip_peaks <- read.delim(paste0(ROOT_DIR, "/raw_data/ChIP-seq/tf_clusters_liver_ENCODE/new_hg38/hg38_tf_clusters_liver_labelled.txt"))

chip_peaks %>%
  count(name) %>%
  arrange(-n) %>%
  pull(name) ->
  TF_names

FP_METHODS <- c(
    "PRINT_beta_logit", "PRINT_beta_cloglog"
)


for (FP_METHOD in FP_METHODS) {
    result_df <- data.frame()
    cat("FP_METHOD = ", FP_METHOD, "\n")
    DIR <- paste0("../../regression/FP_methods/", FP_METHOD)
    
    cat("\tLoading regression results...\n")
    fpscore_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_genotype_regression.Rds"))
    fpscore_cov_regression <- readRDS(paste0(DIR, "/regression_results/fp_score_covariates_genotype_regression.Rds"))

    if (grepl("local_cutsites_geq10", FP_METHOD, fixed = TRUE)) {
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
    )

    variant_range <- GRanges(seqnames = fpQTL_data$variant_chrom,
                             ranges = IRanges(start = fpQTL_data$variant_pos, 
                                              end = fpQTL_data$variant_pos))


    cat("\tMaking 2x2 tables for ChIP TFs...\n")
    for (TF_name in TF_names) {
        cat("\t\tTF_name = ", TF_name, "\n")

        chip_peaks %>%
            filter(name == TF_name) ->
            TF_peaks

        TF_peaks_range <- GRanges(seqnames = TF_peaks$chrom,
                                  ranges = IRanges(start = TF_peaks$chromStart, 
                                                   end = TF_peaks$chromEnd))

        overlaps <- GenomicRanges::findOverlaps(variant_range, 
                                                TF_peaks_range,
                                                select = "arbitrary")

        fpQTL_data$chip_overlap <- !is.na(overlaps)

        # If there are no overlaps with any variants, skip
        if (all(!fpQTL_data$chip_overlap)) {
            cat("\t\t\tSkipping ", TF_name, "...\n")
            next
        }

        # Without covariates
        table2x2 <- table(fpQTL_data$chip_overlap, fpQTL_data$fpQTL_no_cov)
        rownames(table2x2) <- c("no chip", "chip overlap")
        colnames(table2x2) <- c("non-fpQTL", "fpQTL")
        ft <- fisher.test(table2x2)

        df <- data.frame(
            fp_method = FP_METHOD_NAME,
            TF_name = TF_name,
            covariates = FALSE,
            snps_neither = table2x2[1,1],
            snps_fpQTL_only = table2x2[1,2],
            snps_chip_only = table2x2[2,1],
            snps_chip_and_fpQTL = table2x2[2,2],
            odds_ratio = ft$estimate,
            odds_ratio_95_low = ft$conf.int[1],
            odds_ratio_95_high = ft$conf.int[2],
            fishers_p = ft$p.value,
            row.names = NULL
        )
        result_df <- rbind(result_df, df)

        # With covariates
        table2x2 <- table(fpQTL_data$chip_overlap, fpQTL_data$fpQTL_with_cov)
        rownames(table2x2) <- c("no chip", "chip overlap")
        colnames(table2x2) <- c("non-fpQTL", "fpQTL")
        ft <- fisher.test(table2x2)

        df <- data.frame(
            fp_method = FP_METHOD_NAME,
            TF_name = TF_name,
            covariates = TRUE,
            snps_neither = table2x2[1,1],
            snps_fpQTL_only = table2x2[1,2],
            snps_chip_only = table2x2[2,1],
            snps_chip_and_fpQTL = table2x2[2,2],
            odds_ratio = ft$estimate,
            odds_ratio_95_low = ft$conf.int[1],
            odds_ratio_95_high = ft$conf.int[2],
            fishers_p = ft$p.value,
            row.names = NULL
        )
        result_df <- rbind(result_df, df)

    }
    
    result_df %>% 
		write.table(
			sprintf("overlap_tables/%s_fpQTL_chip-seq_overlap.txt", FP_METHOD),
			row.names = FALSE, quote = FALSE, sep = "\t"
		)

}
