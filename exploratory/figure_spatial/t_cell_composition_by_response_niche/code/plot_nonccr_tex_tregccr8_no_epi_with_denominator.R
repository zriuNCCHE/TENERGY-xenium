suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

source_long <- file.path(
  project_root,
  "final", "extended_figure_5", "t_cell_composition_by_niche", "outputs",
  "extended_figure5_tex_within_cd8_tregccr8_within_treg_by_patient_niche_long.csv"
)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

neighborhood_order <- c(
  "Niche_iCAF_CXCL5",
  "Niche_Mono",
  "Niche_Neutro",
  "Niche_Macro_CXCL9",
  "Niche_iCAF_CXCL6",
  "Niche_T_cell"
)
feature_order <- c("CD8_Tex_PDCD1 within CD8", "CD4_Treg_CCR8 within Treg")

dt <- fread(source_long)
dt <- dt[cCR == "non-cCR" & neighborhood_label != "Niche_Epi"]
dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]
dt[, feature := factor(feature, levels = feature_order)]

summary_dt <- dt[
  ,
  .(
    n_defined_patient_niches = sum(!is.na(composition_percent)),
    n_total_patient_niches = .N,
    n_100_percent = sum(composition_percent == 100, na.rm = TRUE),
    median_percent = median(composition_percent, na.rm = TRUE),
    mean_percent = mean(composition_percent, na.rm = TRUE),
    median_denominator_cells = as.numeric(median(denominator_cells[!is.na(composition_percent)], na.rm = TRUE)),
    mean_denominator_cells = as.numeric(mean(denominator_cells[!is.na(composition_percent)], na.rm = TRUE)),
    max_denominator_cells = as.numeric(max(denominator_cells[!is.na(composition_percent)], na.rm = TRUE))
  ),
  by = .(feature, neighborhood_label)
]

hundred_dt <- dt[
  composition_percent == 100,
  .(patientID, sampleID, neighborhood_label, feature, denominator_cells)
][order(feature, denominator_cells)]

fwrite(summary_dt, file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_summary.csv"))
fwrite(hundred_dt, file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_100_percent_patient_niches.csv"))

order_dt <- dt[
  ,
  .(median_percent = median(composition_percent, na.rm = TRUE)),
  by = .(feature, neighborhood_label)
][order(feature, median_percent, as.character(neighborhood_label))]
order_dt[, feature_neighborhood := paste(as.character(feature), as.character(neighborhood_label), sep = "___")]
feature_neighborhood_levels <- order_dt$feature_neighborhood
feature_neighborhood_labels <- setNames(
  as.character(order_dt$neighborhood_label),
  order_dt$feature_neighborhood
)

plot_dt <- copy(dt)
plot_dt[, feature_neighborhood := factor(
  paste(as.character(feature), as.character(neighborhood_label), sep = "___"),
  levels = feature_neighborhood_levels
)]

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

composition_plot <- ggplot(plot_dt, aes(x = feature_neighborhood, y = composition_percent)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(size = pmin(denominator_cells, 200)),
    position = position_jitter(width = 0.08, height = 0),
    alpha = 0.72,
    color = "#333333",
    stroke = 0
  ) +
  facet_wrap(~ feature, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = feature_neighborhood_labels) +
  scale_size_continuous(range = c(0.55, 1.45), guide = "none") +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "Subset composition (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.35, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

denominator_plot <- ggplot(plot_dt, aes(x = feature_neighborhood, y = denominator_cells)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    alpha = 0.72,
    size = 0.75,
    color = "#333333",
    stroke = 0
  ) +
  facet_wrap(~ feature, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = feature_neighborhood_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = NULL,
    y = "Denominator cells"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.35, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

denominator_zoom_plot <- denominator_plot +
  coord_cartesian(ylim = c(0, 120)) +
  labs(
    x = NULL,
    y = "Denominator cells (zoomed)"
  )

ggsave(
  file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_composition_boxplots.pdf"),
  composition_plot,
  device = "pdf",
  width = 160,
  height = 72,
  units = "mm",
  dpi = 300
)

ggsave(
  file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_denominator_boxplots.pdf"),
  denominator_plot,
  device = "pdf",
  width = 160,
  height = 72,
  units = "mm",
  dpi = 300
)

ggsave(
  file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_denominator_boxplots_zoomed_0_120.pdf"),
  denominator_zoom_plot,
  device = "pdf",
  width = 160,
  height = 72,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Exploratory non-cCR-only within-lineage T-cell composition by niche",
    "",
    "This diagnostic analysis excludes Niche_Epi and includes only non-cCR post-chemoradiotherapy samples.",
    "CD8_Tex_PDCD1 is measured as a percentage of CD8_Teff, CD8_Tex_PDCD1 and CD8_prolif cells within each patient-niche.",
    "CD4_Treg_CCR8 is measured as a percentage of CD4_Treg_CCR8 and CD4_Treg_FOXP3 cells within each patient-niche.",
    "A separate denominator plot shows the number of CD8-lineage or Treg-like cells used for each patient-niche percentage.",
    "A zoomed denominator plot limits the displayed y-axis to 0-120 cells to make low-denominator patient-niches visible; the exported CSV retains the full denominator counts.",
    "Patient-niches with zero denominator cells are undefined and are not plotted in the composition panel."
  ),
  file.path(legend_dir, "legend_nonccr_tex_tregccr8_no_epi.md")
)
