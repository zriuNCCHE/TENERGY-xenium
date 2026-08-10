suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

source_dir <- file.path(project_root, "final", "figure_spatial", "neighborhood_area_by_response", "outputs")
summary_path <- file.path(source_dir, "postc_k7_neighborhood_area_weighted_score_by_response_summary.csv")
stats_path <- file.path(source_dir, "postc_k7_neighborhood_area_weighted_score_by_response_stats.csv")

out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

group_order <- c("cCR", "non-cCR")
selected_neighborhoods <- paste0("neighborhood_", c(1, 2, 4, 7), "_k7")
selected_labels <- sub("_k7$", "", sub("neighborhood_", "Neighborhood ", selected_neighborhoods))
outcome_colors <- unlist(palette_config$clinical_outcome)[group_order]

format_p_value <- function(p_value, prefix = "p") {
  if (is.na(p_value)) {
    return(paste0(prefix, " = NA"))
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("%s = %.6f", prefix, p_value)))
  }
  sprintf("%s = %.3f", prefix, p_value)
}

area_summary <- fread(summary_path)
stats_df <- fread(stats_path)

area_summary <- area_summary[nmf_neighborhood %in% selected_neighborhoods]
stats_df <- stats_df[nmf_neighborhood %in% selected_neighborhoods]
stats_df[, p_adjusted_bh := p.adjust(p_value, method = "BH")]

area_summary[, cCR := factor(cCR, levels = group_order)]
area_summary[, nmf_neighborhood := factor(nmf_neighborhood, levels = selected_neighborhoods)]
area_summary[, neighborhood_label := factor(neighborhood_label, levels = selected_labels)]

stats_df[, nmf_neighborhood := factor(nmf_neighborhood, levels = selected_neighborhoods)]
stats_df[, neighborhood_label := factor(neighborhood_label, levels = selected_labels)]

fwrite(
  area_summary,
  file.path(out_dir, "figure4_ncr_associated_neighborhood_area_weighted_scores.csv")
)
fwrite(
  stats_df,
  file.path(out_dir, "figure4_ncr_associated_neighborhood_stats.csv")
)

p_labels <- merge(
  stats_df,
  area_summary[, .(panel_max = max(area_weighted_score, na.rm = TRUE)), by = .(nmf_neighborhood)],
  by = "nmf_neighborhood",
  all.x = TRUE
)
p_labels[, neighborhood_label := factor(neighborhood_label, levels = selected_labels)]
global_max <- max(area_summary$area_weighted_score, na.rm = TRUE)
label_unit <- max(global_max * 0.055, 0.001)
p_labels[, y_start := panel_max + label_unit * 1.1]
p_labels[, y_end := y_start - label_unit * 0.35]
p_labels[, label_y := y_start + label_unit * 1.05]
p_labels[, label := mapply(
  function(raw_p, adj_p) {
    paste(
      format_p_value(raw_p, "raw p"),
      format_p_value(adj_p, "adj p"),
      sep = "\n"
    )
  },
  p_value,
  p_adjusted_bh,
  USE.NAMES = FALSE
)]
y_upper <- max(c(area_summary$area_weighted_score, p_labels$label_y), na.rm = TRUE) * 1.06

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(area_summary, aes(x = cCR, y = area_weighted_score)) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    color = "#333333",
    linewidth = line_width,
    fill = NA
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.075, height = 0),
    size = 0.9,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = 1, xend = 2, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = 1, xend = 1, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = 2, xend = 2, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_text(
    data = p_labels,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    color = "#D62728",
    size = 1.75,
    lineheight = 0.9
  ) +
  facet_wrap(~ neighborhood_label, nrow = 2, ncol = 2) +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = pretty(c(0, y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Area-weighted neighborhood score"
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
    strip.text.x = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.18, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "figure4_ncr_associated_spatial_neighborhoods.pdf"),
  plot = plot,
  device = "pdf",
  width = 82,
  height = 100,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Figure 4: non-cCR-associated spatial neighborhoods",
    "",
    "Area-weighted neighborhood score was calculated per post-chemoradiotherapy sample from continuous k=7 NMF weights as sum(neighborhood weight x cell area) divided by total annotated cell area in that sample.",
    "Only neighborhoods 1, 2, 4 and 7 are shown because these neighborhoods were higher in non-cCR samples in the full seven-neighborhood screen.",
    "Raw P values are two-sided Wilcoxon rank-sum tests comparing cCR and non-cCR samples. Adjusted P values are Benjamini-Hochberg adjusted across the four neighborhoods shown in this panel."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
