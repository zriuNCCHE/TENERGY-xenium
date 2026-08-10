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
  normalizePath("exploratory/figure_6/mp_score_margin_association_candidates/code/plot_mp_score_margin_association_candidates.R")
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

distance_path <- "/Users/liuz/Downloads/malignant_distance_info.csv"
mp_score_path <- file.path(project_root, "data", "geo_ready", "metaprograms", "malignant_cell_metaprogram_scores_k7_nMP10_nGene20_clean.csv")
data_dir <- file.path(project_root, "data", "geo_ready", "metaprograms")
out_dir <- file.path(project_root, "exploratory", "figure_6", "mp_score_margin_association_candidates", "outputs")
result_dir <- file.path(project_root, "exploratory", "figure_6", "mp_score_margin_association_candidates", "results")
legend_dir <- file.path(project_root, "exploratory", "figure_6", "mp_score_margin_association_candidates", "legends")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
mp_palette <- unlist(palette_config$malignant_metaprograms)

mp_cols <- paste0("MP_", 1:7)
region_levels <- c("tumor_inner", "tumor_margin")
region_labels <- c("tumor_inner" = "Inner", "tumor_margin" = "Margin")
region_cols <- c("Inner" = "#BFC7D5", "Margin" = "#C52B2F")

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = plotting_config$font$family) +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
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

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

distance_keep <- c(
  "cell_id",
  "revised_cellID",
  "sampleID",
  "sample_timepoint",
  "patientID",
  "cCR",
  "final_cell_type",
  "final_cell_type2",
  "major_cell_type",
  "distance_to_malignant",
  "nearest_malignant_cell",
  "neighbor_count",
  "malignant_neighbor_count",
  "malignant_neighbor_fraction",
  "malignant_region",
  "malignant_proximity",
  "distance_region"
)
distance_dt <- fread(distance_path, select = distance_keep)
mp_dt <- fread(mp_score_path)

distance_dt <- distance_dt[final_cell_type == "malignant" | major_cell_type == "malignant"]
distance_dt <- unique(distance_dt, by = "revised_cellID")
mp_dt <- unique(mp_dt, by = "revised_cellID")

merged <- merge(
  mp_dt[, c("revised_cellID", mp_cols), with = FALSE],
  distance_dt,
  by = "revised_cellID",
  all.x = FALSE,
  all.y = FALSE
)
setcolorder(merged, c(
  "cell_id",
  "revised_cellID",
  "sampleID",
  "sample_timepoint",
  "patientID",
  "cCR",
  "final_cell_type",
  "final_cell_type2",
  "major_cell_type",
  "malignant_region",
  "distance_region",
  "distance_to_malignant",
  "neighbor_count",
  "malignant_neighbor_count",
  "malignant_neighbor_fraction",
  "malignant_proximity",
  "nearest_malignant_cell",
  mp_cols
))

dominant_matrix <- as.matrix(merged[, ..mp_cols])
dominant_index <- max.col(dominant_matrix, ties.method = "first")
merged[, dominant_MP := mp_cols[dominant_index]]
merged[, dominant_MP_score := dominant_matrix[cbind(seq_len(.N), dominant_index)]]
merged[, MP7_is_dominant := dominant_MP == "MP_7"]

clean_path <- file.path(data_dir, "malignant_cell_metaprogram_scores_with_distance_region_clean.csv")
fwrite(merged, clean_path)

analysis_dt <- merged[
  sample_timepoint == "pre" &
    malignant_region %in% region_levels
]
analysis_dt[, malignant_region := factor(malignant_region, levels = region_levels)]

long_dt <- melt(
  analysis_dt,
  id.vars = c("cell_id", "revised_cellID", "sampleID", "patientID", "cCR", "malignant_region"),
  measure.vars = mp_cols,
  variable.name = "MP",
  value.name = "mp_score"
)
long_dt[, MP := factor(MP, levels = mp_cols)]
long_dt[, region_label := factor(region_labels[as.character(malignant_region)], levels = c("Inner", "Margin"))]

sample_region <- long_dt[, .(
  n_cells = .N,
  mean_score = mean(mp_score, na.rm = TRUE),
  median_score = median(mp_score, na.rm = TRUE)
), by = .(sampleID, patientID, cCR, malignant_region, region_label, MP)]
fwrite(sample_region, file.path(result_dir, "mp_score_by_sample_region_long.csv"))

wide <- dcast(
  sample_region,
  sampleID + patientID + cCR + MP ~ malignant_region,
  value.var = "mean_score"
)
wide <- wide[!is.na(tumor_inner) & !is.na(tumor_margin)]
wide[, margin_minus_inner := tumor_margin - tumor_inner]
wide[, margin_over_inner := tumor_margin / tumor_inner]
wide[, MP := factor(MP, levels = mp_cols)]
fwrite(wide, file.path(result_dir, "mp_score_margin_minus_inner_by_sample.csv"))

margin_stats <- rbindlist(lapply(mp_cols, function(mp) {
  sub <- wide[MP == mp]
  p <- tryCatch(wilcox.test(sub$margin_minus_inner, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.table(
    MP = mp,
    n_samples = nrow(sub),
    median_inner = median(sub$tumor_inner, na.rm = TRUE),
    median_margin = median(sub$tumor_margin, na.rm = TRUE),
    median_margin_minus_inner = median(sub$margin_minus_inner, na.rm = TRUE),
    mean_margin_minus_inner = mean(sub$margin_minus_inner, na.rm = TRUE),
    p_value_margin_minus_inner = p,
    p_value_label = fmt_p(p)
  )
}))
margin_stats[, p_adjusted_BH := p.adjust(p_value_margin_minus_inner, method = "BH")]
fwrite(margin_stats, file.path(result_dir, "mp_score_margin_association_stats.csv"))

delta_wide <- dcast(wide, sampleID + patientID + cCR ~ MP, value.var = "margin_minus_inner")
mp7_vs_rows <- lapply(setdiff(mp_cols, "MP_7"), function(mp) {
  sub <- delta_wide[!is.na(MP_7) & !is.na(get(mp))]
  p <- tryCatch(wilcox.test(sub$MP_7 - sub[[mp]], mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.table(
    comparison = paste0("MP_7 delta vs ", mp, " delta"),
    n_samples = nrow(sub),
    median_MP7_delta = median(sub$MP_7, na.rm = TRUE),
    median_other_MP_delta = median(sub[[mp]], na.rm = TRUE),
    median_MP7_minus_other_delta = median(sub$MP_7 - sub[[mp]], na.rm = TRUE),
    p_value_paired_wilcoxon = p,
    p_value_label = fmt_p(p)
  )
})
mp7_vs_stats <- rbindlist(mp7_vs_rows)
mp7_vs_stats[, p_adjusted_BH := p.adjust(p_value_paired_wilcoxon, method = "BH")]
fwrite(mp7_vs_stats, file.path(result_dir, "mp7_margin_delta_vs_other_mp_stats.csv"))

dominant_sample_region <- analysis_dt[, .(
  n_cells = .N,
  MP7_dominant_percent = mean(MP7_is_dominant, na.rm = TRUE) * 100
), by = .(sampleID, patientID, cCR, malignant_region)]
dominant_wide <- dcast(dominant_sample_region, sampleID + patientID + cCR ~ malignant_region, value.var = "MP7_dominant_percent")
dominant_wide[, margin_minus_inner := tumor_margin - tumor_inner]
fwrite(dominant_wide, file.path(result_dir, "mp7_dominant_fraction_margin_by_sample.csv"))

p_label_dt <- margin_stats[, .(MP, p_value_label)]
p_y <- sample_region[, .(y = max(mean_score, na.rm = TRUE)), by = MP]
p_label_dt <- merge(p_label_dt, p_y, by = "MP", all.x = TRUE)
p_label_dt[, MP := factor(MP, levels = mp_cols)]

paired_plot <- ggplot(sample_region, aes(region_label, mean_score, group = sampleID)) +
  geom_line(colour = "#B7B7B7", linewidth = 0.35, alpha = 0.85) +
  geom_point(aes(fill = region_label), shape = 21, colour = "#333333", stroke = 0.25, size = 1.25, alpha = 1) +
  geom_text(
    data = p_label_dt,
    aes(x = 1, y = y, label = p_value_label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.2,
    size = 2.0,
    family = plotting_config$font$family
  ) +
  facet_wrap(~MP, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = region_cols, name = "Region") +
  labs(
    title = "Malignant MP scores by tumor region",
    x = NULL,
    y = "Sample mean MP score"
  ) +
  theme_nc() +
  theme(legend.position = "right")
ggsave(
  file.path(out_dir, "candidate_A_all_mp_inner_margin_paired.pdf"),
  paired_plot,
  width = 180,
  height = 55,
  units = "mm",
  device = "pdf",
  dpi = 300
)

delta_plot <- ggplot(wide, aes(MP, margin_minus_inner, fill = MP)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5, colour = "black") +
  geom_point(size = 0.9, colour = "black", alpha = 1, position = position_jitter(width = 0.08, height = 0)) +
  scale_fill_manual(values = mp_palette, guide = "none") +
  labs(
    title = "Tumor-margin enrichment by MP",
    x = NULL,
    y = "Margin minus inner mean score"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(
  file.path(out_dir, "candidate_B_all_mp_margin_minus_inner_delta.pdf"),
  delta_plot,
  width = 95,
  height = 62,
  units = "mm",
  device = "pdf",
  dpi = 300
)

advantage_long <- melt(
  delta_wide,
  id.vars = c("sampleID", "patientID", "cCR", "MP_7"),
  measure.vars = setdiff(mp_cols, "MP_7"),
  variable.name = "other_MP",
  value.name = "other_delta"
)
advantage_long[, MP7_minus_other_delta := MP_7 - other_delta]
advantage_long[, other_MP := factor(other_MP, levels = setdiff(mp_cols, "MP_7"))]

advantage_plot <- ggplot(advantage_long, aes(other_MP, MP7_minus_other_delta)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(size = 0.9, colour = "#C52B2F", alpha = 1, position = position_jitter(width = 0.08, height = 0)) +
  labs(
    title = "MP7 margin advantage over other MPs",
    x = "Comparator MP",
    y = "MP7 delta minus comparator delta"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(
  file.path(out_dir, "candidate_C_mp7_margin_advantage_vs_other_mps.pdf"),
  advantage_plot,
  width = 95,
  height = 62,
  units = "mm",
  device = "pdf",
  dpi = 300
)

region_mean <- sample_region[, .(mean_of_sample_means = mean(mean_score, na.rm = TRUE)), by = .(MP, region_label)]
heat_plot <- ggplot(region_mean, aes(region_label, MP, fill = mean_of_sample_means)) +
  geom_tile(colour = "white", linewidth = 0.5 * 0.352777778) +
  geom_text(aes(label = sprintf("%.3f", mean_of_sample_means)), size = 2.0, family = plotting_config$font$family) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#C52B2F", midpoint = 0, name = "Mean score") +
  labs(
    title = "Mean MP score by tumor region",
    x = NULL,
    y = NULL
  ) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank())
ggsave(
  file.path(out_dir, "candidate_D_all_mp_region_mean_heatmap.pdf"),
  heat_plot,
  width = 62,
  height = 70,
  units = "mm",
  device = "pdf",
  dpi = 300
)

dominant_plot <- ggplot(dominant_wide, aes("", margin_minus_inner)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_boxplot(width = 0.42, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(size = 0.9, colour = "#C52B2F", alpha = 1, position = position_jitter(width = 0.05, height = 0)) +
  labs(
    title = "MP7-dominant cells at tumor margin",
    x = NULL,
    y = "Margin minus inner MP7-dominant cells (%)"
  ) +
  theme_nc()
ggsave(
  file.path(out_dir, "candidate_E_mp7_dominant_fraction_margin_delta.pdf"),
  dominant_plot,
  width = 55,
  height = 62,
  units = "mm",
  device = "pdf",
  dpi = 300
)

legend_text <- c(
  "Exploratory MP score tumor-margin association candidate legends",
  "",
  "Candidate A. Paired sample-level inner versus margin mean scores for each malignant metaprogram. Each line connects tumor-inner and tumor-margin mean scores from the same pretreatment sample. P values were calculated by two-sided paired Wilcoxon signed-rank tests on sample-level mean scores.",
  "",
  "Candidate B. Sample-level margin-minus-inner mean score for each MP. Values above zero indicate tumor-margin enrichment.",
  "",
  "Candidate C. MP7 margin advantage over other MPs. Each value is the sample-level MP7 margin-minus-inner score minus the corresponding comparator MP margin-minus-inner score.",
  "",
  "Candidate D. Heatmap of mean sample-level MP scores by tumor region.",
  "",
  "Candidate E. Difference in the fraction of MP7-dominant malignant cells between tumor margin and tumor inner regions."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# MP Score Tumor-Margin Association Candidates",
  "",
  "Exploratory candidate plots testing whether MP7 is more associated with tumor margin than other malignant metaprograms.",
  "",
  "## Cleaned Data",
  "",
  paste0("- `", clean_path, "`"),
  "",
  "## Outputs",
  "",
  "- `outputs/candidate_A_all_mp_inner_margin_paired.pdf`",
  "- `outputs/candidate_B_all_mp_margin_minus_inner_delta.pdf`",
  "- `outputs/candidate_C_mp7_margin_advantage_vs_other_mps.pdf`",
  "- `outputs/candidate_D_all_mp_region_mean_heatmap.pdf`",
  "- `outputs/candidate_E_mp7_dominant_fraction_margin_delta.pdf`",
  "- `results/mp_score_margin_association_stats.csv`",
  "- `results/mp7_margin_delta_vs_other_mp_stats.csv`"
)
writeLines(readme_text, file.path(project_root, "exploratory", "figure_6", "mp_score_margin_association_candidates", "README.md"))

message("Saved cleaned merged table to: ", clean_path)
message("Saved MP margin association candidate plots to: ", out_dir)
