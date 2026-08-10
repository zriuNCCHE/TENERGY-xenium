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

target_cell_type <- "A2ML1+ epi"
timepoint_order <- c("pre", "postC", "postA")
group_order <- c("cCR", "non-cCR")
timepoint_colors <- unlist(palette_config$timepoints)
timepoint_colors <- timepoint_colors[timepoint_order]

format_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.3f", p_value)
}

metadata <- read.csv(data_path, check.names = TRUE)
metadata <- metadata[, c("sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type")]
metadata$is_target <- metadata$final_cell_type == target_cell_type

summary_df <- aggregate(
  is_target ~ patientID + sampleID + sample_timepoint + cCR,
  data = metadata,
  FUN = function(x) c(total_cells = length(x), a2ml1_positive_cells = sum(x))
)
summary_df <- do.call(data.frame, summary_df)
names(summary_df)[names(summary_df) == "is_target.total_cells"] <- "total_cells"
names(summary_df)[names(summary_df) == "is_target.a2ml1_positive_cells"] <- "a2ml1_positive_cells"
summary_df$a2ml1_positive_pct <- summary_df$a2ml1_positive_cells / summary_df$total_cells * 100
summary_df$sample_timepoint <- factor(summary_df$sample_timepoint, levels = timepoint_order)
summary_df$cCR <- factor(summary_df$cCR, levels = group_order)

write.csv(
  summary_df,
  file.path(out_dir, "a2ml1_positive_cell_proportion_summary.csv"),
  row.names = FALSE
)

stats_df <- do.call(
  rbind,
  lapply(group_order, function(group_name) {
    group_df <- summary_df[summary_df$cCR == group_name, ]
    wide_df <- reshape(
      group_df[, c("patientID", "sample_timepoint", "a2ml1_positive_pct")],
      idvar = "patientID",
      timevar = "sample_timepoint",
      direction = "wide"
    )
    p_value <- wilcox.test(
      wide_df$a2ml1_positive_pct.pre,
      wide_df$a2ml1_positive_pct.postC,
      paired = TRUE,
      alternative = "two.sided",
      exact = TRUE
    )$p.value
    data.frame(
      group = group_name,
      comparison = "pre_vs_postC",
      n_patients = nrow(wide_df),
      test = "paired two-sided Wilcoxon signed-rank test",
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })
)
write.csv(
  stats_df,
  file.path(out_dir, "a2ml1_positive_cell_proportion_stats.csv"),
  row.names = FALSE
)

p_labels <- stats_df
p_labels$x <- 1.5
p_labels$x_start <- 1
p_labels$x_end <- 2
p_labels$y_start <- 100
p_labels$y_end <- 99
p_labels$y <- 103
p_labels$label <- vapply(p_labels$p_value, format_p_value, character(1))
p_labels$cCR <- factor(p_labels$group, levels = group_order)

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
part_width <- plotting_config$part_sizes_mm$paired_two_group_boxplot$width
part_height <- plotting_config$part_sizes_mm$paired_two_group_boxplot$height

plot <- ggplot(summary_df, aes(x = sample_timepoint, y = a2ml1_positive_pct)) +
  geom_line(
    aes(group = patientID),
    color = "#BDBDBD",
    linewidth = line_width,
    alpha = 0.50
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(color = sample_timepoint),
    size = 1.25,
    alpha = 0.9
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    linewidth = line_width
  ) +
  geom_text(
    data = p_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.0
  ) +
  facet_wrap(~ cCR, nrow = 1) +
  scale_color_manual(values = timepoint_colors, breaks = timepoint_order, name = "Timepoint") +
  scale_y_continuous(
    limits = c(-5, 140),
    breaks = c(0, 50, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "A2ML1+ epithelial cell proportion",
    x = NULL,
    y = "A2ML1+ epithelial cells (%)"
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
    strip.text = element_text(size = 8, face = "bold", margin = margin(t = 5)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.85, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8)))

ggsave(
  filename = file.path(out_dir, "a2ml1_positive_cell_proportion_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = part_width,
  height = part_height,
  units = "mm",
  dpi = 300
)
