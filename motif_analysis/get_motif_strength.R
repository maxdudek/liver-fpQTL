suppressPackageStartupMessages({
  library(tidyverse)
  library(motifmatchr)
  library(Biostrings)
  library(TFBSTools)
})
ROOT_DIR <- "/mnt/isilon/sfgi/dudekm/"

args = commandArgs(trailingOnly=TRUE)
#MOTIF_THRESHOLD <- "5e-5"
MOTIF_THRESHOLD <- args[1]

cat("Running motif scanner with p=", MOTIF_THRESHOLD, "\n")

cat("Loading data...\n")
seq.df <- readRDS("variant_seqs.Rds")
radius <- nchar(seq.df$ref_seq[1]) %/% 2
PFMs <- readJASPARMatrix(paste0(ROOT_DIR, "raw_data/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"))

pfm.widths <- PFMs %>% lapply(as.matrix) %>% sapply(ncol)
pfm.widths.unique <- pfm.widths %>% unique() %>% sort()

for (width in pfm.widths.unique) {
  cat("Finding motifs of width ", width, "\n", sep = "")

  # Temporary - skip some widths
  # if (width %in% c(6, 12, 14, 8, 11, 18, 20)) {
  #   cat("Skipping...\n")
  #   next
  # }
  
  pfm.subset <- PFMs[pfm.widths == width]
  
  # Crop sequences to have length 2*width - 1
  ref_seqs <- seq.df$ref_seq %>% str_sub(radius+1 - (width-1), radius+1 + (width-1))
  alt_seqs <- seq.df$alt_seq %>% str_sub(radius+1 - (width-1), radius+1 + (width-1))
  
  
  cat("Ref match...\n")
  ref_motif.match <- matchMotifs(pfm.subset, ref_seqs, 
                                 genome = "hg38", bg = "genome",
                                 out = "matches", p.cutoff = as.numeric(MOTIF_THRESHOLD)) 
  ref_motif.match <- motifMatches(ref_motif.match) %>% as.matrix()
  
  cat("Ref score...\n")
  ref_motif.score <- matchMotifs(pfm.subset, ref_seqs, 
                                 genome = "hg38", bg = "genome",
                                 out = "score", p.cutoff = 1)
  ref_motif.score <- motifScores(ref_motif.score) %>% as.matrix()
  
  cat("Alt match...\n")
  alt_motif.match <- matchMotifs(pfm.subset, alt_seqs, 
                                 genome = "hg38", bg = "genome",
                                 out = "matches", p.cutoff = as.numeric(MOTIF_THRESHOLD)) 
  alt_motif.match <- motifMatches(alt_motif.match) %>% as.matrix()
  
  cat("Alt score...\n")
  alt_motif.score <- matchMotifs(pfm.subset, alt_seqs, 
                                 genome = "hg38", bg = "genome",
                                 out = "score", p.cutoff = 1) 
  alt_motif.score <- motifScores(alt_motif.score) %>% as.matrix()
  
  for (i in 1:length(pfm.subset)) {
    pfm <- pfm.subset[[i]]
    motif.name <- name(pfm)
    motif.id <- ID(pfm)
    
    motif.id <- paste(motif.id, motif.name, sep = "_")
    cat("\tExtracting data for ", motif.id, "...\n", sep = "")
    
    df <- data.frame(
      variant_id = seq.df$variant_id,
      ref_score = ref_motif.score[,i],
      alt_score = alt_motif.score[,i],
      ref_match = ref_motif.match[,i],
      alt_match = alt_motif.match[,i]
    ) # %>%
    #   mutate(delta_score = alt_score - ref_score)
    
    df %>% saveRDS(paste0("motif_results/p=", MOTIF_THRESHOLD, "/", motif.id, ".Rds"))
  }
  
}








