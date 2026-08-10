suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

data_path <- file.path(project_root, "data", "raw", "combined_final_all_cell_types_obs.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

target_cell_type <- "Macro_CXCL9"
excluded_from_fibro_immune <- c("A2ML1+ epi", "malignant")
group_order <- c("cCR", "non-cCR")
outcome_colors <- unlist(palette_config$clinical_outcome)[group_order]

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.4f", p_value)
}

two_group_wilcox <- function(values, groups) {
  keep <- !is.na(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(groups[keep])
  exact_possible <- length(unique(values)) == length(values)
  result <- suppressWarnings(wilcox.test(values ~ groups, exact = exact_possible))
  list(
    p_value = result$p.value,
    method = ifelse(exact_possible, "exact", "asymptotic due to ties")
  )
}

metadata <- read.csv(data_path, check.names = TRUE)
metadata <- metadata[
  metadata$sample_timepoint == "postC" &
    !metadata$final_cell_type %in% excluded_from_fibro_immune,
  c("sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type")
]
metadata$is_target <- metadata$final_cell_type == target_cell_type

summary_df <- aggregate(
  is_target ~ patientID + sampleID + sample_timepoint + cCR,
  data = metadata,
  FUN = function(x) c(fibro_immune_cells = length(x), target_cells = sum(x))
)
summary_df <- do.call(data.frame, summary_df)
names(summary_df)[names(summary_df) == "is_target.fibro_immune_cells"] <- "fibro_immune_cells"
names(summary_df)[names(summary_df) == "is_target.target_cells"] <- "target_cells"
summary_df$target_cell_type <- target_cell_type
summary_df$target_pct <- summary_df$target_cells / summary_df$fibro_immune_cells * 100
summary_df$cCR <- factor(summary_df$cCR, levels = group_order)

write.csv(
  summary_df,
  file.path(out_dir, "macro_cxcl9_postc_by_response_summary.csv"),
  row.names = FALSE
)

wilcox_result <- two_group_wilcox(summary_df$target_pct, summary_df$cCR)
stats_df <- data.frame(
  target_cell_type = target_cell_type,
  denominator = "fibro-immune cells",
  comparison = "cCR_vs_non-cCR",
  timepoint = "postC",
  n_cCR = sum(summary_df$cCR == "cCR"),
  n_non_cCR = sum(summary_df$cCR == "non-cCR"),
  test = "two-sided Wilcoxon rank-sum test",
  p_value = wilcox_result$p_value,
  p_value_method = wilcox_result$method,
  stringsAsFactors = FALSE
)
write.csv(
  stats_df,
  file.path(out_dir, "macro_cxcl9_postc_by_response_unadjusted_stats.csv"),
  row.names = FALSE
)

y_upper <- ceiling((max(summary_df$target_pct, na.rm = TRUE) * 1.18 + 0.3) / 1) * 1
y_upper <- max(y_upper, 5)
p_label <- data.frame(
  x_start = 1,
  x_end = 2,
  y_start = y_upper * 0.84,
  y_end = y_upper * 0.80,
  label_y = y_upper * 0.89,
  label = format_p_value(stats_df$p_value)
)

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(summary_df, aes(x = cCR, y = target_pct)) +
  geom_boxplot(
    aes(fill = cCR),
    width = 0.55,
    outlier.shape = NA,
    color = "#333333",
    linewidth = line_width,
    alpha = 0.55
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_label,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_label,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_label,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_text(
    data = p_label,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    color = "#D62728",
    size = 2.0
  ) +
  facet_wrap(~ target_cell_type, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = pretty(c(0, y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Macro_CXCL9 after chemo-radio",
    x = NULL,
    y = "Proportion among fibro-immune cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    legend.key.width = unit(0.32, "cm"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.clip = "off",
    strip.text.x = element_text(size = base_size, face = "bold", margin = margin(t = 3)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    plot.margin = margin(4, 4, 8, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.55)))

ggsave(
  filename = file.path(out_dir, "macro_cxcl9_postc_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = 75,
  height = 78,
  units = "mm",
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "macro_cxcl9_postc_by_response.png"),
  plot = plot,
  device = "png",
  width = 75,
  height = 78,
  units = "mm",
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "# Macro_CXCL9 postC by response",
    "",
    "Exploratory Extended Figure 3-style plot. Macro_CXCL9 proportions were calculated among fibro-immune cells at postC, defined by excluding malignant and A2ML1+ epithelial cells.",
    "",
    "Boxes show median and interquartile range. Points denote patient samples. P value is an unadjusted two-sided Wilcoxon rank-sum test."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
