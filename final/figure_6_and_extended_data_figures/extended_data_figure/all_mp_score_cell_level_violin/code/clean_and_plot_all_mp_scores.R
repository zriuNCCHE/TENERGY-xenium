suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
  library(grid)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  normalizePath("final/extended_figure_6/all_mp_score_cell_level_violin/code/clean_and_plot_all_mp_scores.R")
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

input_path <- "/Users/liuz/Downloads/malignant_cell_metaprograms_k7_10_nMP_10_n_gene_20.csv"
data_dir <- file.path(project_root, "data", "geo_ready", "metaprograms")
out_dir <- file.path(project_root, "final", "extended_figure_6", "all_mp_score_cell_level_violin", "outputs")
result_dir <- file.path(project_root, "final", "extended_figure_6", "all_mp_score_cell_level_violin", "results")
legend_dir <- file.path(project_root, "final", "extended_figure_6", "all_mp_score_cell_level_violin", "legends")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
mp_palette <- unlist(palette_config$malignant_metaprograms)
outcome_cols <- unlist(palette_config$clinical_outcome)

mp_cols <- paste0("MP_", 1:7)
keep_cols <- c(
  "cell_id",
  "revised_cellID",
  "sampleID",
  "sample_timepoint",
  "patientID",
  "cCR",
  "final_cell_type",
  "final_cell_type2",
  mp_cols
)

dt <- fread(input_path)
missing_cols <- setdiff(keep_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

clean_dt <- dt[, ..keep_cols]
clean_dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]
clean_dt[, sample_timepoint := factor(sample_timepoint, levels = c("pre", "postC", "postA"))]

clean_path <- file.path(data_dir, "malignant_cell_metaprogram_scores_k7_nMP10_nGene20_clean.csv")
fwrite(clean_dt, clean_path)

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = plotting_config$font$family) +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.title = element_text(colour = "black", size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.key.height = unit(7, "pt"),
      legend.key.width = unit(7, "pt"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

wilcox_log10_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n_x <- as.numeric(length(x))
  n_y <- as.numeric(length(y))
  n_all <- n_x + n_y
  values <- c(x, y)
  ranks <- rank(values, ties.method = "average")
  rank_sum_x <- sum(ranks[seq_len(n_x)])
  u_x <- rank_sum_x - n_x * (n_x + 1) / 2
  expected <- n_x * n_y / 2
  tie_counts <- data.table(value = values)[, .N, by = value][["N"]]
  tie_sum <- sum(tie_counts^3 - tie_counts)
  variance <- n_x * n_y / 12 * ((n_all + 1) - tie_sum / (n_all * (n_all - 1)))
  z <- (abs(u_x - expected) - 0.5) / sqrt(variance)
  (log(2) + pnorm(-z, log.p = TRUE)) / log(10)
}

format_log10_p <- function(log10_p) {
  if (is.na(log10_p)) return("p = NA")
  if (log10_p < -300) return(paste0("p = 10^", sprintf("%.1f", log10_p)))
  paste0("p = ", sprintf("%.2e", 10^log10_p))
}

format_threshold_p <- function(log10_p) {
  if (is.na(log10_p)) return("n.s.")
  if (log10_p < -300) return("p < 1e-300")
  if (log10_p < -100) return("p < 1e-100")
  if (log10_p < -50) return("p < 1e-50")
  if (log10_p < -30) return("p < 1e-30")
  if (log10_p < -10) return("p < 1e-10")
  if (log10_p < -3) return("p < 0.001")
  if (log10_p < log10(0.05)) return("p < 0.05")
  "n.s."
}

long_dt <- melt(
  clean_dt,
  id.vars = c("cell_id", "revised_cellID", "sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type", "final_cell_type2"),
  measure.vars = mp_cols,
  variable.name = "MP",
  value.name = "mp_score"
)
long_dt[, MP := factor(MP, levels = mp_cols)]

cell_stats <- rbindlist(lapply(mp_cols, function(mp) {
  sub <- long_dt[MP == mp]
  log10_p <- wilcox_log10_p(sub[cCR == "cCR", mp_score], sub[cCR == "non-cCR", mp_score])
  data.table(
    MP = mp,
    n_cCR_cells = sub[cCR == "cCR", .N],
    n_non_cCR_cells = sub[cCR == "non-cCR", .N],
    median_cCR = median(sub[cCR == "cCR", mp_score], na.rm = TRUE),
    median_non_cCR = median(sub[cCR == "non-cCR", mp_score], na.rm = TRUE),
    mean_cCR = mean(sub[cCR == "cCR", mp_score], na.rm = TRUE),
    mean_non_cCR = mean(sub[cCR == "non-cCR", mp_score], na.rm = TRUE),
    delta_median_non_cCR_minus_cCR = median(sub[cCR == "non-cCR", mp_score], na.rm = TRUE) - median(sub[cCR == "cCR", mp_score], na.rm = TRUE),
    log10_p_value_wilcoxon_cell_level = log10_p,
    p_value_label = format_log10_p(log10_p)
  )
}))
cell_stats[, p_adjusted_BH_from_log10 := {
  p_raw <- 10^pmax(log10_p_value_wilcoxon_cell_level, -300)
  p.adjust(p_raw, method = "BH")
}]
cell_stats[, p_adjusted_BH_label := fifelse(
  p_adjusted_BH_from_log10 == 0,
  "BH-adjusted p < 1e-300",
  paste0("BH-adjusted p = ", sprintf("%.2e", p_adjusted_BH_from_log10))
)]
cell_stats[, plot_p_value_label := vapply(log10_p_value_wilcoxon_cell_level, format_threshold_p, character(1))]
fwrite(cell_stats, file.path(result_dir, "all_mp_score_cell_level_response_stats.csv"))
fwrite(long_dt, file.path(result_dir, "all_mp_score_cell_level_long.csv"))

dominant_dt <- copy(clean_dt)
dominant_matrix <- as.matrix(dominant_dt[, ..mp_cols])
dominant_index <- max.col(dominant_matrix, ties.method = "first")
dominant_dt[, dominant_MP := mp_cols[dominant_index]]
dominant_dt[, dominant_MP_score := dominant_matrix[cbind(seq_len(.N), dominant_index)]]
dominant_dt[, second_MP_score := apply(dominant_matrix, 1, function(x) sort(x, decreasing = TRUE)[2])]
dominant_dt[, dominance_margin := dominant_MP_score - second_MP_score]
dominant_dt[, MP7_rank := apply(-dominant_matrix, 1, function(x) match(7, order(x)))]
dominant_dt[, MP7_rank := as.numeric(MP7_rank)]
dominant_dt[, MP7_is_dominant := dominant_MP == "MP_7"]
dominant_dt[, MP7_score := MP_7]
fwrite(dominant_dt, file.path(result_dir, "all_mp_score_dominance_discordance_cell_level.csv"))

dominance_summary <- dominant_dt[, .(
  n_cells = .N,
  median_dominance_margin = median(dominance_margin, na.rm = TRUE),
  MP7_dominant_percent = mean(MP7_is_dominant, na.rm = TRUE) * 100,
  median_MP7_score = median(MP7_score, na.rm = TRUE),
  median_MP7_rank = median(MP7_rank, na.rm = TRUE)
), by = cCR]
fwrite(dominance_summary, file.path(result_dir, "all_mp_score_dominance_discordance_summary.csv"))

score_cor <- cor(clean_dt[, ..mp_cols], use = "pairwise.complete.obs", method = "spearman")
score_cor_dt <- as.data.table(as.table(score_cor))
setnames(score_cor_dt, c("MP_x", "MP_y", "spearman_rho"))
fwrite(score_cor_dt, file.path(result_dir, "all_mp_score_spearman_correlation.csv"))

p_label_dt <- long_dt[, .(
  x = 1,
  y = max(mp_score, na.rm = TRUE)
), by = MP]
p_label_dt <- merge(
  p_label_dt,
  cell_stats[, .(MP, plot_p_value_label)],
  by = "MP",
  all.x = TRUE
)
p_label_dt[, MP := factor(MP, levels = mp_cols)]

all_violin <- ggplot(long_dt, aes(cCR, mp_score, fill = cCR)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_violin(width = 0.72, linewidth = 0.5, colour = "#7FAEC8", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.5, colour = "black", fill = "white") +
  geom_text(
    data = p_label_dt,
    aes(x = x, y = y, label = plot_p_value_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.2,
    size = 2.0,
    family = plotting_config$font$family
  ) +
  facet_wrap(~MP, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = outcome_cols, guide = "none") +
  labs(
    title = "Cell-level metaprogram scores by outcome",
    x = NULL,
    y = "Relative metaprogram score"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(
  file.path(out_dir, "extended_figure6_all_mp_score_cell_level_violin.pdf"),
  all_violin,
  width = 180,
  height = 55,
  units = "mm",
  device = "pdf",
  dpi = 300
)

delta_dt <- copy(cell_stats)
delta_dt[, MP := factor(MP, levels = mp_cols)]
delta_plot <- ggplot(delta_dt, aes(MP, delta_median_non_cCR_minus_cCR, fill = MP)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#666666") +
  geom_col(width = 0.62, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = mp_palette, guide = "none") +
  labs(
    title = "Median score difference by metaprogram",
    x = NULL,
    y = "Median difference: non-cCR - cCR"
  ) +
  theme_nc()
ggsave(
  file.path(out_dir, "extended_figure6_all_mp_score_response_delta.pdf"),
  delta_plot,
  width = 82,
  height = 58,
  units = "mm",
  device = "pdf",
  dpi = 300
)

cor_plot <- ggplot(score_cor_dt, aes(MP_x, MP_y, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.5 * 0.352777778) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.0, family = plotting_config$font$family) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#C52B2F", midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(
    title = "Metaprogram score correlation",
    x = NULL,
    y = NULL
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
ggsave(
  file.path(out_dir, "extended_figure6_all_mp_score_discordance_correlation.pdf"),
  cor_plot,
  width = 82,
  height = 72,
  units = "mm",
  device = "pdf",
  dpi = 300
)

legend_text <- c(
  "Extended Figure 6 all malignant metaprogram score legend draft",
  "",
  "Cell-level relative malignant metaprogram scores by clinical outcome. Violin plots show all pretreatment malignant cells; embedded boxes show median and interquartile range. Scores below zero are expected for relative module scores and indicate expression below the matched background/reference level, not negative expression.",
  "",
  "P values were calculated by two-sided Wilcoxon rank-sum tests comparing cCR and non-cCR cells for each MP, with log-scale p-value storage to avoid numerical underflow. The figure displays threshold-style significance labels only; exact medians, mean scores, log10 p values and BH-adjusted p values across the seven MPs are saved in `results/all_mp_score_cell_level_response_stats.csv`.",
  "",
  "Discordance outputs assign each cell to its highest-scoring MP and summarize the top-score margin, MP7 dominance, MP7 rank and pairwise Spearman correlations among MP scores. These are saved in `results/all_mp_score_dominance_discordance_cell_level.csv`, `results/all_mp_score_dominance_discordance_summary.csv` and `results/all_mp_score_spearman_correlation.csv`."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# All MP Score Cell-Level Violin",
  "",
  "Cleaned continuous malignant-cell MP score table and cell-level response plots.",
  "",
  "## Cleaned Data",
  "",
  paste0("- `", clean_path, "`"),
  "",
  "## Outputs",
  "",
  "- `outputs/extended_figure6_all_mp_score_cell_level_violin.pdf`",
  "- `outputs/extended_figure6_all_mp_score_response_delta.pdf`",
  "- `outputs/extended_figure6_all_mp_score_discordance_correlation.pdf`",
  "- `results/all_mp_score_cell_level_response_stats.csv`",
  "- `results/all_mp_score_dominance_discordance_cell_level.csv`",
  "- `results/all_mp_score_dominance_discordance_summary.csv`",
  "- `results/all_mp_score_spearman_correlation.csv`"
)
writeLines(readme_text, file.path(project_root, "final", "extended_figure_6", "all_mp_score_cell_level_violin", "README.md"))

message("Saved cleaned MP score table to: ", clean_path)
message("Saved all-MP score plots to: ", out_dir)
