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
  "exploratory", "figure_spatial", "cell_type_composition_by_niche", "outputs",
  "postc_patient_niche_cell_type_composition_full.csv"
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
  "Niche_Epi",
  "Niche_Macro_CXCL9",
  "Niche_iCAF_CXCL6",
  "Niche_T_cell"
)

format_p <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < 1e-4) {
    return(paste0("p = ", formatC(p_value, format = "e", digits = 2)))
  }
  if (p_value < 0.001) {
    return(sprintf("p = %.6f", p_value))
  }
  sprintf("p = %.3f", p_value)
}

full_dt <- fread(source_path)
full_dt[, neighborhood_label := as.character(neighborhood_label)]
old_to_canonical <- c(
  "Neighborhood 1" = "Niche_Mono",
  "Neighborhood 2" = "Niche_Neutro",
  "Neighborhood 3" = "Niche_Epi",
  "Neighborhood 4" = "Niche_iCAF_CXCL6",
  "Neighborhood 5" = "Niche_T_cell",
  "Neighborhood 6" = "Niche_Macro_CXCL9",
  "Neighborhood 7" = "Niche_iCAF_CXCL5"
)
full_dt[, neighborhood_label := fifelse(
  neighborhood_label %in% names(old_to_canonical),
  old_to_canonical[neighborhood_label],
  neighborhood_label
)]

all_patient_niches <- unique(full_dt[, .(patientID, sampleID, cCR, neighborhood_label)])
b_dt <- full_dt[
  final_cell_type2 == "B/Plasma",
  .(patientID, sampleID, cCR, neighborhood_label, b_plasma_cells = n_cells)
]
dt <- merge(all_patient_niches, b_dt, by = c("patientID", "sampleID", "cCR", "neighborhood_label"), all.x = TRUE)
dt[is.na(b_plasma_cells), b_plasma_cells := 0L]

total_dt <- unique(full_dt[, .(patientID, sampleID, cCR, neighborhood_label, total_niche_cells)])
dt <- merge(dt, total_dt, by = c("patientID", "sampleID", "cCR", "neighborhood_label"), all.x = TRUE)
dt[, b_plasma_percent := fifelse(total_niche_cells > 0, b_plasma_cells / total_niche_cells * 100, NA_real_)]
dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]
dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

summary_dt <- dt[
  ,
  .(
    n_patients = .N,
    median_cell_percent = median(b_plasma_percent, na.rm = TRUE),
    mean_cell_percent = mean(b_plasma_percent, na.rm = TRUE),
    max_cell_percent = max(b_plasma_percent, na.rm = TRUE),
    median_b_plasma_cells = median(b_plasma_cells, na.rm = TRUE),
    median_total_niche_cells = median(total_niche_cells, na.rm = TRUE)
  ),
  by = neighborhood_label
]

wide <- dcast(dt, patientID ~ neighborhood_label, value.var = "b_plasma_percent")
required_cols <- c("patientID", neighborhood_order)
wide_complete <- wide[complete.cases(wide[, required_cols, with = FALSE])]
friedman_p <- if (nrow(wide_complete) >= 3) {
  suppressWarnings(friedman.test(as.matrix(wide_complete[, neighborhood_order, with = FALSE]))$p.value)
} else {
  NA_real_
}
kruskal_p <- suppressWarnings(kruskal.test(b_plasma_percent ~ neighborhood_label, data = dt)$p.value)
global_stats <- data.table(
  feature = "B/Plasma",
  denominator = "all cells assigned to each niche",
  n_complete_patients_for_friedman = nrow(wide_complete),
  friedman_p = friedman_p,
  kruskal_p = kruskal_p
)

order_dt <- summary_dt[order(median_cell_percent, as.character(neighborhood_label))]
plot_levels <- as.character(order_dt$neighborhood_label)
dt[, neighborhood_plot := factor(as.character(neighborhood_label), levels = plot_levels)]

fwrite(dt, file.path(out_dir, "postc_b_plasma_niche_cell_fraction_by_patient.csv"))
fwrite(summary_dt, file.path(out_dir, "postc_b_plasma_niche_cell_fraction_summary.csv"))
fwrite(global_stats, file.path(out_dir, "postc_b_plasma_niche_cell_fraction_global_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

y_top <- max(dt$b_plasma_percent, na.rm = TRUE) * 1.18
label_dt <- data.table(
  x = ceiling(length(plot_levels) / 2),
  y = y_top * 0.96,
  label = paste0("Friedman ", format_p(friedman_p))
)

plot <- ggplot(dt, aes(x = neighborhood_plot, y = b_plasma_percent)) +
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
  geom_text(
    data = label_dt,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 1.65
  ) +
  scale_y_continuous(
    limits = c(0, y_top),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "B/Plasma cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(out_dir, "postc_b_plasma_niche_cell_fraction_boxplot.pdf"),
  plot,
  device = "pdf",
  width = 85,
  height = 72,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# B/Plasma abundance across spatial niches",
    "",
    "Exploratory focused panel for B/Plasma abundance across post-chemoradiotherapy spatial niches.",
    "For each patient-niche, the plotted value is B/Plasma cells divided by all cells assigned to that niche.",
    "Each point represents one patient-niche. X-axis order follows low-to-high median abundance."
  ),
  file.path(legend_dir, "legend_b_plasma_niche_cell_fraction.md")
)
