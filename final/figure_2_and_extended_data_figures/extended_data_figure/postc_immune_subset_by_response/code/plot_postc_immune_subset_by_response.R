suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

data_path <- file.path(project_root, "data", "raw", "combined_final_all_cell_types_obs.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

target_cell_types <- c("Macro_CXCL5", "Mono_CDC27", "Mono_SLC2A3", "Neutrophil")
excluded_from_fibro_immune <- c("A2ML1+ epi", "malignant")
group_order <- c("cCR", "non-cCR")
outcome_colors <- unlist(palette_config$clinical_outcome)
outcome_colors <- outcome_colors[group_order]

format_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.4f", p_value)
}

two_group_wilcox <- function(values, groups) {
  exact_possible <- length(unique(values)) == length(values)
  result <- suppressWarnings(wilcox.test(values ~ groups, exact = exact_possible))
  list(
    p_value = result$p.value,
    method = ifelse(exact_possible, "exact", "asymptotic due to ties")
  )
}

metadata <- read.csv(data_path, check.names = TRUE)
metadata <- metadata[
  metadata$sample_timepoint == "postC",
  c("sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type")
]
metadata <- metadata[!metadata$final_cell_type %in% excluded_from_fibro_immune, ]

summary_list <- lapply(target_cell_types, function(target_cell_type) {
  metadata$is_target <- metadata$final_cell_type == target_cell_type
  summary_df <- aggregate(
    is_target ~ patientID + sampleID + cCR,
    data = metadata,
    FUN = function(x) c(fibro_immune_cells = length(x), target_cells = sum(x))
  )
  summary_df <- do.call(data.frame, summary_df)
  names(summary_df)[names(summary_df) == "is_target.fibro_immune_cells"] <- "fibro_immune_cells"
  names(summary_df)[names(summary_df) == "is_target.target_cells"] <- "target_cells"
  summary_df$immune_subset <- target_cell_type
  summary_df$target_pct <- summary_df$target_cells / summary_df$fibro_immune_cells * 100
  summary_df
})
summary_df <- do.call(rbind, summary_list)
summary_df$cCR <- factor(summary_df$cCR, levels = group_order)
summary_df$immune_subset <- factor(summary_df$immune_subset, levels = target_cell_types)

write.csv(
  summary_df,
  file.path(out_dir, "postc_immune_subset_by_response_summary.csv"),
  row.names = FALSE
)

stats_list <- lapply(target_cell_types, function(target_cell_type) {
  subset_df <- summary_df[summary_df$immune_subset == target_cell_type, ]
  wilcox_result <- two_group_wilcox(subset_df$target_pct, subset_df$cCR)
  data.frame(
    immune_subset = target_cell_type,
    comparison = "cCR_vs_non-cCR",
    timepoint = "postC",
    n_cCR = sum(subset_df$cCR == "cCR"),
    n_non_cCR = sum(subset_df$cCR == "non-cCR"),
    test = "two-sided Wilcoxon rank-sum test",
    p_value = wilcox_result$p_value,
    p_value_method = wilcox_result$method,
    stringsAsFactors = FALSE
  )
})
stats_df <- do.call(rbind, stats_list)
write.csv(
  stats_df,
  file.path(out_dir, "postc_immune_subset_by_response_unadjusted_stats.csv"),
  row.names = FALSE
)

p_labels <- stats_df
p_labels$immune_subset <- factor(p_labels$immune_subset, levels = target_cell_types)
p_labels$x_start <- 1
p_labels$x_end <- 2
p_labels$y_start <- 60
p_labels$y_end <- 58.5
p_labels$label_y <- 62
p_labels$label <- paste0("p = ", vapply(p_labels$p_value, format_p_value, character(1)))

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
part_width <- plotting_config$part_sizes_mm$postc_four_subset_outcome$width
part_height <- plotting_config$part_sizes_mm$postc_four_subset_outcome$height

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
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_text(
    data = p_labels,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    color = "#D62728",
    size = 2.0
  ) +
  facet_wrap(~ immune_subset, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(0, 68),
    breaks = c(0, 25, 50),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Proportion among fibro-immune cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
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
    strip.text.x = element_text(size = base_size, face = "bold", angle = 45, hjust = 1, vjust = 1, margin = margin(t = 2, b = 6)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.25, "cm"),
    plot.margin = margin(4, 4, 28, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.55)))

ggsave(
  filename = file.path(out_dir, "postc_immune_subset_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = part_width,
  height = part_height,
  units = "mm",
  dpi = 300
)
