suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

source_path <- file.path(
  project_root,
  "exploratory", "figure_spatial", "t_cell_neighborhood_abundance", "outputs",
  "postc_t_cell_neighborhood_abundance_by_patient.csv"
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
feature_order <- c("CD8_Tex_PDCD1", "CD4_Treg_CCR8")
feature_labels <- c(
  CD8_Tex_PDCD1 = "CD8_Tex_PDCD1 density",
  CD4_Treg_CCR8 = "CD4_Treg_CCR8 density"
)

dt <- fread(source_path)
dt <- dt[
  cCR == "non-cCR" &
    neighborhood_label != "Niche_Epi" &
    feature_type == "Individual T subset" &
    feature %in% feature_order
]
dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]
dt[, feature := factor(feature, levels = feature_order)]
dt[, neighborhood_area_mm2 := total_neighborhood_area / 1e6]
dt[, cells_per_mm2 := fifelse(neighborhood_area_mm2 > 0, subset_cells / neighborhood_area_mm2, NA_real_)]

summary_dt <- dt[
  ,
  .(
    n_patient_niches = .N,
    median_cells = median(subset_cells, na.rm = TRUE),
    mean_cells = mean(subset_cells, na.rm = TRUE),
    median_area_mm2 = median(neighborhood_area_mm2, na.rm = TRUE),
    median_cells_per_mm2 = median(cells_per_mm2, na.rm = TRUE),
    mean_cells_per_mm2 = mean(cells_per_mm2, na.rm = TRUE),
    max_cells_per_mm2 = max(cells_per_mm2, na.rm = TRUE)
  ),
  by = .(feature, neighborhood_label)
]

global_stats <- dt[
  ,
  {
    wide <- dcast(.SD, patientID ~ neighborhood_label, value.var = "cells_per_mm2")
    required_cols <- c("patientID", neighborhood_order)
    wide <- wide[complete.cases(wide[, required_cols, with = FALSE])]
    friedman_p <- if (nrow(wide) >= 3) {
      suppressWarnings(friedman.test(as.matrix(wide[, ..neighborhood_order]))$p.value)
    } else {
      NA_real_
    }
    kruskal_p <- if (length(unique(neighborhood_label)) >= 2) {
      suppressWarnings(kruskal.test(cells_per_mm2 ~ neighborhood_label, data = .SD)$p.value)
    } else {
      NA_real_
    }
    .(
      n_complete_patients_for_friedman = nrow(wide),
      friedman_p = friedman_p,
      kruskal_p = kruskal_p
    )
  },
  by = feature
]
global_stats[, kruskal_p_adjusted_bh := p.adjust(kruskal_p, method = "BH")]

order_dt <- dt[
  ,
  .(median_cells_per_mm2 = median(cells_per_mm2, na.rm = TRUE)),
  by = .(feature, neighborhood_label)
][order(feature, median_cells_per_mm2, as.character(neighborhood_label))]
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
plot_dt[, feature_label := factor(feature_labels[as.character(feature)], levels = unname(feature_labels))]

fwrite(dt, file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_density_by_patient_niche.csv"))
fwrite(summary_dt, file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_density_summary.csv"))
fwrite(global_stats, file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_density_global_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

density_plot <- ggplot(plot_dt, aes(x = feature_neighborhood, y = cells_per_mm2)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 0.75,
    alpha = 0.72,
    color = "#333333",
    stroke = 0
  ) +
  facet_wrap(~ feature_label, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = feature_neighborhood_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = NULL,
    y = expression("Cells per mm"^2)
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

ggsave(
  file.path(out_dir, "postc_nonccr_tex_tregccr8_no_epi_density_per_mm2_boxplots.pdf"),
  density_plot,
  device = "pdf",
  width = 160,
  height = 72,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Exploratory non-cCR-only T-cell subset density by niche",
    "",
    "This analysis excludes Niche_Epi and includes only non-cCR post-chemoradiotherapy samples.",
    "For each patient-niche, density is calculated as target subset cell count divided by total assigned neighborhood area in mm2.",
    "CD8_Tex_PDCD1 and CD4_Treg_CCR8 are plotted as absolute cell density rather than within-lineage percentage, reducing sensitivity to small denominators.",
    "X axes are ordered from low to high median density within each feature."
  ),
  file.path(legend_dir, "legend_nonccr_tex_tregccr8_density_no_epi.md")
)
