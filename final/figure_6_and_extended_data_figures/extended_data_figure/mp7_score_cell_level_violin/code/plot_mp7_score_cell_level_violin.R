suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

score_path <- file.path(project_root, "final", "figure_6", "pretreatment_spatial_signature_boundary", "results", "pretreatment_malignant_cell_signature_scores.csv")
out_dir <- file.path(project_root, "final", "extended_figure_6", "mp7_score_cell_level_violin", "outputs")
legend_dir <- file.path(project_root, "final", "extended_figure_6", "mp7_score_cell_level_violin", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

outcome_cols <- c("cCR" = "#7BB6A4", "non-cCR" = "#F2A38A")

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 7, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

format_scientific <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("0")
  sprintf("%.2e", x)
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
  if (is.na(log10_p)) return("Wilcoxon p = NA")
  if (log10_p < -300) {
    return(paste0("Wilcoxon p = 10^", sprintf("%.1f", log10_p)))
  }
  paste0("Wilcoxon p = ", sprintf("%.2e", 10^log10_p))
}

mm_pdf <- function(filename, width_mm, height_mm) {
  grDevices::pdf(
    file.path(out_dir, filename),
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    useDingbats = FALSE,
    family = "Helvetica"
  )
}

score_dt <- fread(score_path, select = c("sampleID", "patientID", "cCR", "malignant_region", "score__MP7 all"))
setnames(score_dt, "score__MP7 all", "mp7_score")
score_dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

wilcox_res <- wilcox.test(mp7_score ~ cCR, data = score_dt, exact = FALSE)
cell_p <- wilcox_res$p.value
cell_log10_p <- wilcox_log10_p(
  score_dt[cCR == "cCR", mp7_score],
  score_dt[cCR == "non-cCR", mp7_score]
)

stats <- score_dt[, .(
  n_cells = .N,
  median_mp7_score = median(mp7_score, na.rm = TRUE),
  mean_mp7_score = mean(mp7_score, na.rm = TRUE),
  q25_mp7_score = quantile(mp7_score, 0.25, na.rm = TRUE),
  q75_mp7_score = quantile(mp7_score, 0.75, na.rm = TRUE),
  high_score_fraction_percent = mean(mp7_score > 0, na.rm = TRUE) * 100
), by = cCR]
stats[, p_value_wilcoxon_cell_level := cell_p]
stats[, log10_p_value_wilcoxon_cell_level := cell_log10_p]
stats[, p_value_label := format_scientific(cell_p)]
stats[, log10_p_value_label := format_log10_p(cell_log10_p)]

fwrite(score_dt, file.path(out_dir, "extended_figure6_mp7_score_cell_level_values.csv"))
fwrite(stats, file.path(out_dir, "extended_figure6_mp7_score_cell_level_stats.csv"))

p_label <- format_log10_p(cell_log10_p)

cell_violin <- ggplot(score_dt, aes(cCR, mp7_score, fill = cCR)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_violin(width = 0.72, linewidth = 0.5, colour = "#7FAEC8", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.5, colour = "black", fill = "white") +
  scale_fill_manual(values = outcome_cols, guide = "none") +
  labs(
    title = "Cell-level MP7 score by outcome",
    x = NULL,
    y = "MP7 score"
  ) +
  annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.35, size = 2.0, label = p_label) +
  theme_nc()
mm_pdf("extended_figure6_mp7_score_cell_level_violin.pdf", 58, 58)
print(cell_violin)
dev.off()

legend_text <- c(
  "Extended Figure 6 MP7 score cell-level violin legend",
  "",
  paste0("Cell-level MP7 score by clinical outcome. Violin plots show the distribution of continuous MP7 scores across pretreatment malignant cells; embedded boxes show median and interquartile range. P value was calculated by two-sided Wilcoxon rank-sum test using all malignant cells (", p_label, "). This panel is descriptive because cells are not independent biological replicates.")
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# Extended Figure 6 MP7 Score Cell-Level Violin",
  "",
  "Final plot showing cell-level continuous MP7 score by clinical outcome.",
  "",
  "## Output",
  "",
  "- `outputs/extended_figure6_mp7_score_cell_level_violin.pdf`",
  "- `outputs/extended_figure6_mp7_score_cell_level_stats.csv`"
)
writeLines(readme_text, file.path(dirname(out_dir), "README.md"))

message("Saved final MP7 score violin to: ", out_dir)
