suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(knitr)
})

args = commandArgs(trailingOnly=TRUE)
FP_METHOD <- args[1]
if (grepl("no_gaussian", FP_METHOD, fixed = TRUE) | grepl("beta", FP_METHOD, fixed = TRUE)) {
  FP_SCORE_LIM <- c(0,1)
} else {
  FP_SCORE_LIM <- c(-3, 3)
}

genotype_colors <- c("#d33682", "#6c71c4", "#2aa198")

fpQTLs_FDR5 <- read.delim(paste0("../regression/FP_methods/", FP_METHOD,
                                 "/regression_results/fpQTLs_covariates_FDR5.txt"))
fpscore_matrix <- readRDS(paste0("../regression/FP_methods/", FP_METHOD,
                                 "/data/fpscore_matrix_fpQTLs_with_covariates.Rds"))
fpscore_matrix_normalized <- readRDS(paste0("../regression/FP_methods/", FP_METHOD,
                                            "/data/fpscore_matrix_normalized_fpQTLs_with_covariates.Rds"))

genotype_matrix <- readRDS(paste0("../regression/FP_methods/", FP_METHOD,
                                  "/data/genotype_matrix_fpQTLs_with_covariates.Rds"))

fpscore_residuals_file <- paste0("../regression/FP_methods/", FP_METHOD,
                                 "/data/fpscore_residual_matrix_fpQTLs.Rds")
if (file.exists(fpscore_residuals_file)) {
  fpscore_residuals_matrix <- readRDS(fpscore_residuals_file)
} else {
  # If we don't have residuals, make dummy data
  fpscore_residuals_matrix <- matrix(0, nrow = nrow(fpscore_matrix), 
                                     ncol = ncol(fpscore_matrix),
                                     dimnames = list(rownames(fpscore_matrix), 
                                                     colnames(fpscore_matrix)))
}


# rsIDs <- fpQTLs_FDR5 %>% arrange(pval) %>% head(25) %>% pull(variant_id)
rsIDs <- fpQTLs_FDR5 %>% arrange(pval) %>% pull(variant_id)

# filenames <- paste0("cutsites/", FP_METHOD, "/", rsIDs, ".Rds")
# cutsites <- lapply(filenames, readRDS)
# names(cutsites) <- rsIDs
cutsites <- readRDS(paste0("cutsites/", FP_METHOD, "/fpQTL_cutsites.Rds"))
cutsites <- lapply(cutsites, dplyr::select, sample, position, num_insertions, corrected_insertions)

get_df <- function(rs) {
  data.frame(
    fpscore_normalized = fpscore_matrix_normalized[rs,] %>% na.omit(),
    fpscore = fpscore_matrix[rs,] %>% na.omit(),
    fpscore_residual = fpscore_residuals_matrix[rs,],
    genotype = genotype_matrix[rs,] %>% na.omit()
  ) %>%
    rownames_to_column("sample") %>%
    arrange(sample)
}

get_sample_order <- function(rs) {
  get_df(rs) %>%
    arrange(genotype, fpscore_normalized) %>%
    pull(sample) ->
    sample_order
  
  return(sample_order)
}

get_title <- function(rs) {
  rs.info <- fpQTLs_FDR5 %>% filter(variant_id == rs)
  
  title <- paste0(rs, " (", rs.info$variant_chrom, ":", rs.info$variant_pos, "), beta = ",
                  round(rs.info$beta, 3), ", -log10(p) = ", round(-log10(rs.info$pval), 3))
  
  return(title)
}


plot_boxplot <- function(rs, normalized = TRUE, residuals = FALSE) {
  y.aes = ifelse(normalized, "fpscore_normalized", "fpscore")
  y.lab = ifelse(normalized, "Normalized FP score", "FP score")
  
  y.aes = ifelse(residuals, "fpscore_residual", y.aes)
  y.lab = ifelse(residuals, "FP score residual", y.lab)
  
  get_df(rs) %>%
    ggplot(aes(x = as.factor(genotype), y = .data[[y.aes]])) +
    geom_boxplot(aes(fill = as.factor(genotype))) + 
    geom_jitter(height = 0, width = 0.1) +
    scale_fill_manual(name = "Genotype", values = genotype_colors, labels = 0:2) +
    xlab("Genotype") + 
    ylab(y.lab) +
    ggtitle(get_title(rs)) +
    theme_classic()
}


plot_cutsites_pileup <- function(rs, ncol=7) {
  # Order samples by fp score
  sample_order <- get_sample_order(rs)
  nrow <- length(sample_order) %/% ncol
  
  cutsites[[rs]] %>%
    inner_join(get_df(rs), join_by(sample)) ->
    df
  
  df %>%
    mutate(sample = factor(sample, levels = sample_order)) %>%
    filter(!is.na(sample)) %>% # Remove samples with NA fp score
    arrange(genotype, sample) %>%
    ggplot() +
    geom_bar(aes(x = position, y = num_insertions, fill = as.factor(genotype)),
             stat = "identity") +
    geom_point(x = 0, y = -0.5, size = 0.5, shape = "triangle") +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_text(
      aes(label = paste0("FP score = ", round(fpscore_normalized, 2))),
      x = 0, y = 9, size = 2) +
    xlim(c(-75, 75)) +
    scale_y_continuous(breaks = c(0, 5, 10), limits = c(-1, 10)) +
    ylab("Raw insertions") +
    scale_fill_manual(values = genotype_colors, name = "Genotype") +
    facet_wrap(~as.factor(sample), ncol = ncol) +
    theme_classic() +
    theme(strip.text = element_text(size = 6)) +
    coord_fixed(8) +
    ggtitle(get_title(rs))
}

plot_corrected_cutsites_pileup <- function(rs, ncol=7) {
  # Order samples by fp score
  sample_order <- get_sample_order(rs)
  nrow <- length(sample_order) %/% ncol
  
  cutsites[[rs]] %>%
    inner_join(get_df(rs), join_by(sample)) ->
    df
  
  df %>%
    mutate(sample = factor(sample, levels = sample_order)) %>%
    filter(!is.na(sample)) %>% # Remove samples with NA fp score
    arrange(genotype, sample) %>%
    ggplot() +
    geom_bar(aes(x = position, y = corrected_insertions, fill = as.factor(genotype)),
             stat = "identity") +
    geom_point(x = 0, y = -0.5, size = 0.5, shape = 4) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_text(
      aes(label = paste0("FP score = ", round(fpscore_normalized, 2))),
      x = 0, y = 4, size = 2) +
    xlim(c(-75, 75)) +
    scale_y_continuous(breaks = c(-4, 0, 4), limits = c(-5, 6)) +
    ylab("Corrected insertions") +
    scale_fill_manual(values = genotype_colors, name = "Genotype") +
    facet_wrap(~as.factor(sample), ncol = ncol) +
    theme_classic() +
    theme(strip.text = element_text(size = 6)) +
    coord_fixed(8) +
    ggtitle(get_title(rs))
}

plot_bulk_cutsites_pileup <- function(rs) {
  cutsites[[rs]] %>%
    inner_join(get_df(rs), join_by(sample)) ->
    df
  
  df %>%
    filter(!is.na(sample)) %>% # Remove samples with NA fp score
    group_by(position, genotype) %>%
    summarise(num_insertions = sum(num_insertions)) %>%
    ggplot() +
    geom_bar(aes(x = position, y = num_insertions, fill = as.factor(genotype)),
             stat = "identity") +
    geom_point(x = 0, y = 0, size = 2, shape = "triangle") +
    xlim(c(-75, 75)) +
    ylab("Raw insertions") +
    scale_fill_manual(values = genotype_colors, name = "Genotype") +
    facet_wrap(~as.factor(genotype), ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(get_title(rs))
}

plot_bulk_corrected_cutsites_pileup <- function(rs) {
  cutsites[[rs]] %>%
    inner_join(get_df(rs), join_by(sample)) ->
    df
  
  df %>%
    filter(!is.na(sample)) %>% # Remove samples with NA fp score
    group_by(position, genotype) %>%
    summarise(corrected_insertions = sum(corrected_insertions)) %>%
    ggplot() +
    geom_bar(aes(x = position, y = corrected_insertions, fill = as.factor(genotype)),
             stat = "identity") +
    geom_point(x = 0, y = 0, size = 2, shape = 4) +
    geom_hline(yintercept = 0, size = 0.5) +
    xlim(c(-75, 75)) +
    ylab("Corrected insertions") +
    scale_fill_manual(values = genotype_colors, name = "Genotype") +
    facet_wrap(~as.factor(genotype), ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(get_title(rs))
}


# ComplexHeatmap 
ht_opt$TITLE_PADDING = unit(c(4, 8), "points")

plot_complex_heatmap <- function(rs) {
  sample_order <- get_sample_order(rs)
  
  variant_df <- get_df(rs)
  
  complete_df <- 
    data.frame(
      sample = rep(sample_order, each = length(-50:50)),
      position = rep(-50:50, times = length(sample_order))
    )
  
  cutsites[[rs]] %>%
    select(sample, position, num_insertions) %>%
    filter(sample %in% sample_order) %>% # Remove samples where fp_score = NA
    mutate(position = as.integer(position)) %>%
    right_join(complete_df, by = join_by(sample, position)) %>%
    arrange(sample, position) %>%
    replace(is.na(.), 0) %>%
    pivot_wider(names_from = position, values_from = num_insertions) %>%
    arrange(sample) %>%
    column_to_rownames("sample") %>%
    as.matrix() ->
    mat
  
  mat %>%
    Heatmap(
      name = "Tn5 insertions",
      row_order = order(variant_df$fpscore_normalized),
      show_row_names = FALSE,
      column_labels = c("-50", rep("", 24), "-25",  rep("", 24), "0",
                        rep("", 24), "25",  rep("", 24), "50"),
      column_names_rot = 0,
      column_title_side = "bottom", column_title = "Position",
      column_order = colnames(mat),
      col = colorRamp2(c(0, 6), c("white", "red")),
      row_split = factor(variant_df$genotype, levels = 0:2),
      row_title = gt_render(
        0:(length(unique(variant_df$genotype))-1),
        halign = 0.5, valign = 0.5,
        padding = unit(c(0, 0, 0, 0), "pt")
      ),
      row_title_gp = gpar(fill = genotype_colors, fontsize = 20, lwd = 0),
      row_gap = unit(5, "mm"),
      right_annotation = rowAnnotation(
        "FP score\n(normalized)" = variant_df$fpscore_normalized,
        col = list("FP score\n(normalized)" = colorRamp2(FP_SCORE_LIM, c("white", "dark green")))
      ),
      width = unit(8, "in"), height = unit(7, "in")
    ) %>%
    draw(
      row_title = gt_render("Genotype", gp = gpar(cex = 1.3)),
      column_title = gt_render(get_title(rs), gp = gpar(cex = 1))
    )
}

plot_complex_heatmap_corrected <- function(rs) {
  sample_order <- get_sample_order(rs)
  
  variant_df <- get_df(rs)
  
  complete_df <- 
    data.frame(
      sample = rep(sample_order, each = length(-50:50)),
      position = rep(-50:50, times = length(sample_order))
    )
  
  cutsites[[rs]] %>%
    select(sample, position, corrected_insertions) %>%
    filter(sample %in% sample_order) %>% # Remove samples where fp_score = NA
    mutate(position = as.integer(position)) %>%
    right_join(complete_df, by = join_by(sample, position)) %>%
    arrange(sample, position) %>%
    replace(is.na(.), 0) %>%
    pivot_wider(names_from = position, values_from = corrected_insertions) %>%
    arrange(sample) %>%
    column_to_rownames("sample") %>%
    as.matrix() ->
    mat
  
  mat %>%
    Heatmap(
      name = "Tn5 insertions\nCORRECTED",
      row_order = order(variant_df$fpscore_normalized),
      show_row_names = FALSE,
      column_labels = c("-50", rep("", 24), "-25",  rep("", 24), "0",
                        rep("", 24), "25",  rep("", 24), "50"),
      column_names_rot = 0,
      column_title_side = "bottom", column_title = "Position",
      column_order = colnames(mat),
      col = colorRamp2(c(-5, 0, 6), c("blue", "white", "red")),
      row_split = factor(variant_df$genotype, levels = 0:2),
      row_title = gt_render(
        0:(length(unique(variant_df$genotype))-1),
        halign = 0.5, valign = 0.5,
        padding = unit(c(0, 0, 0, 0), "pt")
      ),
      row_title_gp = gpar(fill = genotype_colors, fontsize = 20, lwd = 0),
      row_gap = unit(5, "mm"),
      right_annotation = rowAnnotation(
        "FP score\n(normalized)" = variant_df$fpscore_normalized,
        col = list("FP score\n(normalized)" = colorRamp2(FP_SCORE_LIM, c("white", "dark green")))
      ),
      width = unit(8, "in"), height = unit(7, "in")
    ) %>%
    draw(
      row_title = gt_render("Genotype", gp = gpar(cex = 1.3)),
      column_title = gt_render(get_title(rs), gp = gpar(cex = 1))
    )
}

dir.create(file.path("figures", FP_METHOD), showWarnings = FALSE)

# Plot cutsite pileups
cat("Plotting pileups...\n")
pdf(width = 10, height = 20, file.path("figures", FP_METHOD, "all_fpQTL_pileups.pdf"))
for (rs in names(cutsites)) {
  print(rs)
  print(plot_cutsites_pileup(rs))
  print(plot_corrected_cutsites_pileup(rs))
}
dev.off()

# Plot bulk cutsite pileups
cat("Plotting bulk pileups...\n")
pdf(width = 7, height = 7, file.path("figures", FP_METHOD, "all_fpQTL_bulk_pileups.pdf"))
for (rs in names(cutsites)) {
  print(rs)
  print(plot_bulk_cutsites_pileup(rs))
  print(plot_bulk_corrected_cutsites_pileup(rs))
}
dev.off()


# # Plot boxplots
# cat("Plotting boxplots...\n")
# pdf(width = 6, height = 4, file.path("figures", FP_METHOD, "all_fpQTL_boxplots.pdf"))
# for (rs in names(cutsites)) {
#   print(rs)
#   print(plot_boxplot(rs, normalized = FALSE))
#   print(plot_boxplot(rs, normalized = TRUE))
#   print(plot_boxplot(rs, residuals = TRUE))
# }
# dev.off()


# Plot complex heatmaps
cat("Plotting heatmaps...\n")
pdf(width = 10, height = 10, file.path("figures", FP_METHOD, "all_fpQTL_heatmaps.pdf"))
for (rs in names(cutsites)) {
  print(rs)
  print(plot_complex_heatmap(rs))
  print(plot_complex_heatmap_corrected(rs))
}
dev.off()


cat("Done\n")

