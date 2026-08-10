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
  normalizePath("exploratory/figure_6/mp7_subsignature_margin_association/code/plot_mp7_subsignature_margin_association.R")
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

summary_path <- file.path(project_root, "final", "figure_6", "pretreatment_spatial_signature_boundary", "results", "pretreatment_malignant_sample_region_summary_long_from_csv.csv")
delta_path <- file.path(project_root, "final", "figure_6", "pretreatment_spatial_signature_boundary", "results", "pretreatment_malignant_boundary_delta_by_sample_long_from_csv.csv")
stats_path <- file.path(project_root, "final", "figure_6", "pretreatment_spatial_signature_boundary", "results", "pretreatment_malignant_spatial_signature_statistics.csv")

out_dir <- file.path(project_root, "exploratory", "figure_6", "mp7_subsignature_margin_association", "outputs")
result_dir <- file.path(project_root, "exploratory", "figure_6", "mp7_subsignature_margin_association", "results")
legend_dir <- file.path(project_root, "exploratory", "figure_6", "mp7_subsignature_margin_association", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

signature_order <- c(
  "MP7 all",
  "Laminin-332",
  "Adhesion / anchoring",
  "Remodeling / invasion",
  "Stress / adaptation"
)
subsignature_order <- signature_order[-1]
signature_labels <- c(
  "MP7 all" = "MP7 all",
  "Laminin-332" = "Laminin-332",
  "Adhesion / anchoring" = "Adhesion / anchoring",
  "Remodeling / invasion" = "Remodeling / invasion",
  "Stress / adaptation" = "Broader stress"
)
region_cols <- c(
  "Inner" = palette_config$tumor_regions[["inner region"]],
  "Margin" = palette_config$tumor_regions[["margin"]]
)

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", sprintf("%.3f", p))
}

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

summary_dt <- fread(summary_path)
delta_dt <- fread(delta_path)
stats_dt <- fread(stats_path)
stats_dt[, feature_name := sub("^score__", "", feature)]

summary_sig <- summary_dt[
  feature_type == "signature" &
    feature_name %in% subsignature_order &
    region_label %in% c("Inner", "Margin") &
    statistic == "mean"
]
summary_sig[, feature_name := factor(feature_name, levels = subsignature_order)]
summary_sig[, feature_label := factor(signature_labels[as.character(feature_name)], levels = signature_labels[subsignature_order])]
summary_sig[, region_label := factor(region_label, levels = c("Inner", "Margin"))]

p_labels <- stats_dt[
  grepl("^score__", feature) &
    test == "paired_boundary_enrichment" &
    feature_name %in% subsignature_order,
  .(feature_name, p_value)
]
p_labels[, feature_name := factor(feature_name, levels = subsignature_order)]
p_y <- summary_sig[, .(y = max(value, na.rm = TRUE)), by = feature_name]
p_labels <- merge(p_labels, p_y, by = "feature_name", all.x = TRUE)
p_labels[, label := vapply(p_value, fmt_p, character(1))]
p_labels[, feature_label := factor(signature_labels[as.character(feature_name)], levels = signature_labels[subsignature_order])]

paired_plot <- ggplot(summary_sig, aes(region_label, value, group = sampleID)) +
  geom_line(colour = "#B7B7B7", linewidth = 0.35, alpha = 0.9) +
  geom_point(aes(fill = region_label), shape = 21, colour = "#333333", stroke = 0.25, size = 1.3, alpha = 1) +
  geom_text(
    data = p_labels,
    aes(x = 1, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.25,
    size = 2.0,
    family = plotting_config$font$family
  ) +
  facet_wrap(~feature_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = region_cols, name = "Region") +
  labs(
    title = "MP7 sub-signatures at the tumor margin",
    x = NULL,
    y = "Mean relative signature score"
  ) +
  theme_nc() +
  theme(legend.position = "right")

ggsave(
  file.path(out_dir, "exploratory_mp7_subsignature_inner_margin_paired.pdf"),
  paired_plot,
  width = 170,
  height = 55,
  units = "mm",
  device = "pdf",
  dpi = 300
)

delta_sig <- delta_dt[
  feature_type == "signature" &
    feature_name %in% signature_order
]
delta_sig[, feature_name := factor(feature_name, levels = signature_order)]
delta_sig[, feature_label := factor(signature_labels[as.character(feature_name)], levels = signature_labels[signature_order])]
delta_summary <- delta_sig[, .(
  n_samples = .N,
  median_margin_minus_inner = median(margin_minus_inner, na.rm = TRUE),
  mean_margin_minus_inner = mean(margin_minus_inner, na.rm = TRUE),
  q25 = quantile(margin_minus_inner, 0.25, na.rm = TRUE),
  q75 = quantile(margin_minus_inner, 0.75, na.rm = TRUE)
), by = feature_name]
p_delta <- stats_dt[
  grepl("^score__", feature) &
    test == "paired_boundary_enrichment" &
    feature_name %in% signature_order,
  .(feature_name, p_value)
]
p_delta[, feature_name := factor(feature_name, levels = signature_order)]
delta_summary <- merge(delta_summary, p_delta, by = "feature_name", all.x = TRUE)
delta_summary[, p_value_label := vapply(p_value, fmt_p, character(1))]
delta_summary[, feature_label := factor(signature_labels[as.character(feature_name)], levels = signature_labels[signature_order])]
fwrite(delta_sig, file.path(result_dir, "mp7_signature_margin_minus_inner_by_sample.csv"))
fwrite(delta_summary, file.path(result_dir, "mp7_signature_margin_association_summary.csv"))

p_delta_labels <- delta_sig[, .(
  y = max(margin_minus_inner, na.rm = TRUE) +
    0.12 * diff(range(margin_minus_inner, na.rm = TRUE))
), by = feature_label]
p_delta_labels <- merge(
  p_delta_labels,
  delta_summary[, .(feature_label, p_value_label)],
  by = "feature_label",
  all.x = TRUE
)

delta_plot <- ggplot(delta_sig, aes(feature_label, margin_minus_inner)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#777777") +
  geom_boxplot(width = 0.55, outlier.shape = NA, colour = "#7FAEC8", fill = "white", linewidth = 0.5) +
  geom_point(size = 0.9, colour = "#7FAEC8", alpha = 0.85, position = position_jitter(width = 0.08, height = 0)) +
  geom_text(
    data = p_delta_labels,
    aes(x = feature_label, y = y, label = p_value_label),
    inherit.aes = FALSE,
    hjust = 0.5,
    vjust = 0.5,
    size = 2.0,
    family = plotting_config$font$family
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.18))) +
  labs(
    title = "Margin enrichment differs across MP7 components",
    x = NULL,
    y = "Margin minus inner score"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(
  file.path(out_dir, "exploratory_mp7_signature_margin_minus_inner_delta.pdf"),
  delta_plot,
  width = 82,
  height = 82,
  units = "mm",
  device = "pdf",
  dpi = 300
)

legend_text <- c(
  "Exploratory MP7 sub-signature tumor-margin association",
  "",
  "Paired sample-level plots show mean relative signature scores in tumor inner and tumor margin regions for four MP7 sub-signatures. Each line connects the inner and margin values from one pretreatment sample. P values were calculated by two-sided paired Wilcoxon signed-rank tests across samples.",
  "",
  "The margin-minus-inner delta plot includes the integrated MP7 all score alongside the four sub-signatures. Boxes are unfilled and share the same outline color. Stress / adaptation is labeled as Broader stress. This exploratory panel is intended to show why focused MP7 components can be tumor-margin enriched even when the full MP7 score is not uniformly margin-associated."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# MP7 Sub-signature Tumor-Margin Association",
  "",
  "Exploratory plots using existing Figure 6 sample-level tumor inner/margin summaries.",
  "",
  "## Outputs",
  "",
  "- `outputs/exploratory_mp7_subsignature_inner_margin_paired.pdf`",
  "- `outputs/exploratory_mp7_signature_margin_minus_inner_delta.pdf`",
  "- `results/mp7_signature_margin_minus_inner_by_sample.csv`",
  "- `results/mp7_signature_margin_association_summary.csv`"
)
writeLines(readme_text, file.path(project_root, "exploratory", "figure_6", "mp7_subsignature_margin_association", "README.md"))

message("Saved exploratory MP7 sub-signature margin plots to: ", out_dir)
