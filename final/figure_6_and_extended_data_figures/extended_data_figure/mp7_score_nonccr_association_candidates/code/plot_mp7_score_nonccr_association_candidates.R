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
out_dir <- file.path(project_root, "final", "extended_figure_6", "mp7_score_nonccr_association_candidates", "outputs")
legend_dir <- file.path(project_root, "final", "extended_figure_6", "mp7_score_nonccr_association_candidates", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

outcome_cols <- c("cCR" = "#7BB6A4", "non-cCR" = "#F2A38A")
mp7_col <- "#C52B2F"

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 7, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
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

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
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

score_dt <- fread(
  score_path,
  select = c("sampleID", "patientID", "cCR", "malignant_region", "score__MP7 all")
)
setnames(score_dt, "score__MP7 all", "mp7_score")
score_dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]
score_dt[, malignant_region := factor(malignant_region, levels = c("tumor_inner", "tumor_margin", "isolated_malignant"))]

sample_summary <- score_dt[, .(
  n_cells = .N,
  median_mp7_score = median(mp7_score, na.rm = TRUE),
  mean_mp7_score = mean(mp7_score, na.rm = TRUE),
  q25_mp7_score = quantile(mp7_score, 0.25, na.rm = TRUE),
  q75_mp7_score = quantile(mp7_score, 0.75, na.rm = TRUE),
  high_score_fraction_percent = mean(mp7_score > 0, na.rm = TRUE) * 100
), by = .(sampleID, patientID, cCR)]

region_summary <- score_dt[, .(
  n_cells = .N,
  median_mp7_score = median(mp7_score, na.rm = TRUE),
  mean_mp7_score = mean(mp7_score, na.rm = TRUE),
  high_score_fraction_percent = mean(mp7_score > 0, na.rm = TRUE) * 100
), by = .(sampleID, patientID, cCR, malignant_region)]

sample_p <- tryCatch(wilcox.test(median_mp7_score ~ cCR, data = sample_summary, exact = FALSE)$p.value, error = function(e) NA_real_)
sample_high_p <- tryCatch(wilcox.test(high_score_fraction_percent ~ cCR, data = sample_summary, exact = FALSE)$p.value, error = function(e) NA_real_)
cell_p <- tryCatch(wilcox.test(mp7_score ~ cCR, data = score_dt, exact = FALSE)$p.value, error = function(e) NA_real_)

stats <- data.table(
  level = c("sample_median_score", "sample_fraction_score_gt_0", "cell_level_score_descriptive"),
  test = c("Wilcoxon rank-sum", "Wilcoxon rank-sum", "Wilcoxon rank-sum"),
  n_cCR = c(sum(sample_summary$cCR == "cCR"), sum(sample_summary$cCR == "cCR"), sum(score_dt$cCR == "cCR")),
  n_non_cCR = c(sum(sample_summary$cCR == "non-cCR"), sum(sample_summary$cCR == "non-cCR"), sum(score_dt$cCR == "non-cCR")),
  median_cCR = c(
    median(sample_summary[cCR == "cCR", median_mp7_score], na.rm = TRUE),
    median(sample_summary[cCR == "cCR", high_score_fraction_percent], na.rm = TRUE),
    median(score_dt[cCR == "cCR", mp7_score], na.rm = TRUE)
  ),
  median_non_cCR = c(
    median(sample_summary[cCR == "non-cCR", median_mp7_score], na.rm = TRUE),
    median(sample_summary[cCR == "non-cCR", high_score_fraction_percent], na.rm = TRUE),
    median(score_dt[cCR == "non-cCR", mp7_score], na.rm = TRUE)
  ),
  p_value = c(sample_p, sample_high_p, cell_p),
  note = c("Primary sample-level score comparison", "Sample-level proportion of malignant cells with MP7 score > 0", "Descriptive only; cells are not independent biological replicates")
)

fwrite(score_dt, file.path(out_dir, "extended_figure6_mp7_score_cell_level.csv"))
fwrite(sample_summary, file.path(out_dir, "extended_figure6_mp7_score_sample_summary.csv"))
fwrite(region_summary, file.path(out_dir, "extended_figure6_mp7_score_region_sample_summary.csv"))
fwrite(stats, file.path(out_dir, "extended_figure6_mp7_score_response_stats.csv"))

sample_box <- ggplot(sample_summary, aes(cCR, median_mp7_score)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(aes(colour = cCR), size = 0.9, alpha = 1, position = position_jitter(width = 0.08, height = 0), show.legend = FALSE) +
  scale_colour_manual(values = outcome_cols) +
  labs(
    title = "Pretreatment MP7 score by outcome",
    x = NULL,
    y = "Median MP7 score per sample"
  ) +
  annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.35, size = 2.0, label = paste0("Wilcoxon ", fmt_p(sample_p))) +
  theme_nc()
mm_pdf("extended_figure6_candidate_A_mp7_score_sample_median_boxplot.pdf", 58, 58)
print(sample_box)
dev.off()

high_fraction_box <- ggplot(sample_summary, aes(cCR, high_score_fraction_percent)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(aes(colour = cCR), size = 0.9, alpha = 1, position = position_jitter(width = 0.08, height = 0), show.legend = FALSE) +
  scale_colour_manual(values = outcome_cols) +
  labs(
    title = "MP7-high fraction by outcome",
    x = NULL,
    y = "Cells with MP7 score > 0 (%)"
  ) +
  annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.35, size = 2.0, label = paste0("Wilcoxon ", fmt_p(sample_high_p))) +
  theme_nc()
mm_pdf("extended_figure6_candidate_B_mp7_score_positive_fraction_boxplot.pdf", 58, 58)
print(high_fraction_box)
dev.off()

cell_density <- ggplot(score_dt, aes(mp7_score, colour = cCR)) +
  geom_density(linewidth = 0.55, adjust = 1.1) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  scale_colour_manual(values = outcome_cols, name = "Clinical outcome") +
  labs(
    title = "Cell-level MP7 score distribution",
    x = "MP7 score",
    y = "Density"
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.35, size = 2.0, label = paste0("Cell-level ", fmt_p(cell_p))) +
  theme_nc() +
  theme(legend.position = "right")
mm_pdf("extended_figure6_candidate_C_mp7_score_cell_level_density.pdf", 78, 58)
print(cell_density)
dev.off()

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
  annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.35, size = 2.0, label = paste0("Cell-level ", fmt_p(cell_p))) +
  theme_nc()
mm_pdf("extended_figure6_candidate_C2_mp7_score_cell_level_violin.pdf", 58, 58)
print(cell_violin)
dev.off()

cell_cdf <- score_dt[!is.na(mp7_score), .(
  mp7_score = sort(mp7_score),
  cumulative_fraction = seq_len(.N) / .N
), by = cCR]
cell_cdf_plot <- ggplot(cell_cdf, aes(mp7_score, cumulative_fraction, colour = cCR)) +
  geom_step(linewidth = 0.55) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  scale_colour_manual(values = outcome_cols, name = "Clinical outcome") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Cell-level MP7 score CDF",
    x = "MP7 score",
    y = "Cumulative fraction of malignant cells"
  ) +
  theme_nc() +
  theme(legend.position = "right")
mm_pdf("extended_figure6_candidate_D_mp7_score_cell_level_cdf.pdf", 78, 58)
print(cell_cdf_plot)
dev.off()

rank_dt <- copy(sample_summary)
rank_dt[, sampleID := factor(sampleID, levels = rank_dt[order(cCR, median_mp7_score), sampleID])]
rank_plot <- ggplot(rank_dt, aes(sampleID, median_mp7_score, fill = cCR)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = outcome_cols, name = "Clinical outcome") +
  labs(
    title = "Sample-ranked MP7 score",
    x = NULL,
    y = "Median MP7 score per sample"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
mm_pdf("extended_figure6_candidate_E_mp7_score_sample_rank_barplot.pdf", 105, 58)
print(rank_plot)
dev.off()

region_plot <- ggplot(region_summary[malignant_region != "isolated_malignant"], aes(cCR, median_mp7_score)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(aes(colour = cCR), size = 0.8, alpha = 1, position = position_jitter(width = 0.08, height = 0), show.legend = FALSE) +
  facet_wrap(~ malignant_region, nrow = 1) +
  scale_colour_manual(values = outcome_cols) +
  labs(
    title = "MP7 score by malignant region",
    x = NULL,
    y = "Median MP7 score per sample-region"
  ) +
  theme_nc()
mm_pdf("extended_figure6_candidate_F_mp7_score_region_boxplots.pdf", 100, 58)
print(region_plot)
dev.off()

legend_text <- c(
  "Extended Figure 6 MP7 score non-cCR association candidate legends",
  "",
  "Candidate A. Pretreatment sample-level median MP7 score by clinical outcome. Each dot is one pretreatment sample. Boxes show median and interquartile range. P value was calculated by two-sided Wilcoxon rank-sum test.",
  "",
  "Candidate B. Percentage of malignant cells with MP7 score greater than zero by clinical outcome. Each dot is one pretreatment sample. Boxes show median and interquartile range. P value was calculated by two-sided Wilcoxon rank-sum test.",
  "",
  "Candidate C. Cell-level MP7 score density by clinical outcome. This is descriptive because cells are not independent biological replicates.",
  "",
  "Candidate C2. Cell-level MP7 score violin plot by clinical outcome using all pretreatment malignant cells. This is descriptive because cells are not independent biological replicates.",
  "",
  "Candidate D. Cell-level cumulative distribution of MP7 score by clinical outcome. Right-shifted non-cCR curve indicates higher MP7 scores. This is descriptive because cells are not independent biological replicates.",
  "",
  "Candidate E. Sample-ranked pretreatment median MP7 score.",
  "",
  "Candidate F. Sample-level median MP7 score by malignant region and clinical outcome."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# Extended Figure 6 MP7 Score Non-cCR Association Candidates",
  "",
  "Score-based candidate plots testing whether continuous pretreatment MP7 score is associated with non-cCR clinical outcome.",
  "",
  "This is different from the hard-classification analysis because it uses `score__MP7 all` for every pretreatment malignant cell instead of assigning each cell to a single MP.",
  "",
  "## Outputs",
  "",
  "- `outputs/extended_figure6_candidate_A_mp7_score_sample_median_boxplot.pdf`",
  "- `outputs/extended_figure6_candidate_B_mp7_score_positive_fraction_boxplot.pdf`",
  "- `outputs/extended_figure6_candidate_C_mp7_score_cell_level_density.pdf`",
  "- `outputs/extended_figure6_candidate_C2_mp7_score_cell_level_violin.pdf`",
  "- `outputs/extended_figure6_candidate_D_mp7_score_cell_level_cdf.pdf`",
  "- `outputs/extended_figure6_candidate_E_mp7_score_sample_rank_barplot.pdf`",
  "- `outputs/extended_figure6_candidate_F_mp7_score_region_boxplots.pdf`"
)
writeLines(readme_text, file.path(dirname(out_dir), "README.md"))

message("Saved Extended Figure 6 MP7 score candidate plots to: ", out_dir)
