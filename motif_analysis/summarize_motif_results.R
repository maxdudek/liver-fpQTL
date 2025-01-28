suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})
ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

args = commandArgs(trailingOnly=TRUE)

# MOTIF_THRESHOLD <- "p=5e-4"
MOTIF_THRESHOLD <- args[1]
cat("Running motif summarizer for ", MOTIF_THRESHOLD, "\n")

# Get a list of variants used for motif_results (smaller than the total list)
sample_results <- readRDS(paste0("motif_results/", MOTIF_THRESHOLD, "/MA1489.1_FOXN3.Rds"))

# Get positions of those variants
variant_info <- readRDS("../regression/variant_info.Rds")



# Get list of motifs
motif_result_files <- list.files(paste0("motif_results/", MOTIF_THRESHOLD), "*.Rds", full.names = TRUE)
motif_ids <- str_match(motif_result_files, ".*/(.*?)\\.Rds")[,2]
head(motif_ids)

# Load ChIP peaks
chip_peaks <- read.delim(paste0(ROOT_DIR, "raw_data/ChIP-seq/tf_clusters_liver_ENCODE/new_hg38/hg38_tf_clusters_liver_labelled.txt"))

chip_peaks <- chip_peaks %>%
  mutate(name = toupper(name))

chip_peaks %>%
  count(name) %>%
  arrange(-n) %>%
  pull(name) ->
  ChIP_TF_names

# Create dataframe of motifs, merge with ChIP-seq
motif_info <-
  data.frame(
    motif_id = motif_ids
  ) %>%
  separate(col = "motif_id", into = c("motif_id2", "TF_name"), 
           sep = "_", remove = FALSE) %>%
  mutate(TF_name = toupper(TF_name)) %>%
  separate_rows(TF_name, sep = "::") %>%
  mutate(has_chip_data = TF_name %in% ChIP_TF_names)

motif_info %>%
  write.table("motif_summary/motif_info.txt",
              quote = FALSE, row.names = FALSE, sep = "\t")


FP_METHODS <- c(
  "PRINT_beta_cloglog"
)

for (FP_METHOD in FP_METHODS) {
  cat("FP_METHOD = ", FP_METHOD, "\n", sep = "")
  
  regression_results <- readRDS(paste0("../regression/FP_methods/", FP_METHOD, "/regression_results/fp_score_covariates_genotype_regression.Rds"))
  
  if (grepl("local_cutsites_geq10", FP_METHOD, fixed = TRUE)) {
    Q_THRESHOLD <- 0.10
  } else {
    Q_THRESHOLD <- 0.05
  }
  
  tmp = nrow(regression_results)
  cat(sprintf("nSNPs before filtering = %d\n", tmp))
  regression_results <- regression_results %>%
    # Only consider variants we have results for
    filter(variant_id %in% sample_results$variant_id) %>% 
    mutate(fpQTL = ST_qval < Q_THRESHOLD)
  cat(sprintf("nSNPs after filtering = %d\n", nrow(regression_results)))
  cat(sprintf("nSNPs filtered = %d\n", tmp-nrow(regression_results)))

  # Some regression results have less variants - so we need to filter for them
  motif_result_filter <- (sample_results$variant_id %in% regression_results$variant_id)
  
  # Filter variant_info in the same way
  variant_info_filtered <- variant_info %>%
    filter(variant_id %in% regression_results$variant_id)
  
  variant_range <- GRanges(seqnames = variant_info_filtered$variant_chrom,
                           ranges = IRanges(start = variant_info_filtered$variant_pos, 
                                            end = variant_info_filtered$variant_pos))
    
  all_matches <- data.frame()
  motif_summary <- data.frame()
  
  for (i in seq_along(motif_result_files)) {
    motif_result_file <- motif_result_files[i]
    motif_id <- motif_ids[i]
    cat("\tSummarizing ", motif_id, "\n", sep = "")
    
    
    motif_result <- readRDS(motif_result_file)
    
    motif_result %>%
      filter(motif_result_filter) %>%
      mutate(
        delta_score = alt_score - ref_score,
        match = ref_match | alt_match,
        beta = regression_results$beta,
        fpQTL = regression_results$fpQTL
      ) ->
      motif_result
    
    
    # If we have ChIP data for the motif TF, filter for only variants within a chip peak
    motif_id_ <- motif_id # Fix dplyr namespace issue
    TFs_with_chip_data <- motif_info %>% 
      filter(motif_id == motif_id_, has_chip_data) %>%
      pull(TF_name)
    
    # Filter for variants within a ChIP peak
    if (length(TFs_with_chip_data) != 0) {
      cat("\t\t", motif_id, " has ChIP data, filtering on variants near ChIP peaks...\n", sep = "")
      TF_chip_peaks <- chip_peaks %>%
        filter(name %in% TFs_with_chip_data)
      
      TF_chip_peaks_range <- GRanges(seqnames = TF_chip_peaks$chrom,
                                     ranges = IRanges(start = TF_chip_peaks$chromStart, 
                                                      end = TF_chip_peaks$chromEnd))
      
      overlaps <- GenomicRanges::findOverlaps(variant_range, TF_chip_peaks_range, maxgap = -1L)
      
      motif_result <- motif_result[unique(queryHits(overlaps)), ]
    } else {
    	TFs_with_chip_data <- NA
    }
    
    
    # Save motif matches to all_matches
    motif_result %>%
      filter(match) %>%
      mutate(motif_id = motif_id) %>%
      select(-match) ->
      motif_result_matches
    all_matches <- rbind(all_matches, motif_result_matches)
    
    # Get 2x2 table of matches vs. fpQTLs
    motif_result %>%
      count(match, fpQTL) %>%
      complete(match = c(TRUE, FALSE), fpQTL = c(TRUE, FALSE), fill = list(n = 0L)) %>%
      pull(n) %>%
      matrix(nrow = 2,
             dimnames = list(
               c("non_fpQTL", "fpQTL"),
               c("no_match", "match")
             )) ->
      m
    
    num_match_fpQTL = m["fpQTL", "match"]
    
    # Test - are fpQTLs enriched for motif matches?
    ft <- fisher.test(m)
    
    # How many fpQTL matches are concordant?
    # i.e., the fpQTL effect size and the motif score change are in the same direction
    motif_result %>%
      filter(match & fpQTL) %>%
      filter(sign(beta) == sign(delta_score)) %>%
      nrow() ->
      num_concordant
    
    
    # Test - are fpQTLs more likely to be concordant than at random (p = 0.5)?
    if (num_match_fpQTL > 0) {
      bt <- binom.test(num_concordant, num_match_fpQTL, p = 0.5)
    } else {
      # If there are no fpQTLs that overlap this motif
      bt <- list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
    }
    

    motif_summary_row <- data.frame(
      motif_id = motif_id,
      TFs_with_chip_data = paste(TFs_with_chip_data, collapse = ";"),
      num_matches = sum(m[,"match"]),
      num_fpQTLs = sum(m["fpQTL",]),
      num_match_fpQTL = num_match_fpQTL,
      num_match_fpQTL_concordant = num_concordant,
      fpQTL_match_OR = ft$estimate,
      fpQTL_match_OR_95_low = ft$conf.int[1],
      fpQTL_match_OR_95_high = ft$conf.int[2],
      fpQTL_match_p = ft$p.value,
      concordant_prob = bt$estimate,
      concordant_prob_95_low = bt$conf.int[1],
      concordant_prob_95_high = bt$conf.int[2],
      concordant_p = bt$p.value,
      row.names = NULL
    )
    
    motif_summary <- rbind(motif_summary, motif_summary_row)
  }
  
  all_matches %>% saveRDS(paste0("motif_match_results/", MOTIF_THRESHOLD, "/", FP_METHOD, "_motif_matches.Rds"))
  motif_summary %>%
    write.table(paste0("motif_summary/", MOTIF_THRESHOLD, "/", FP_METHOD, "_motif_summary.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
}




