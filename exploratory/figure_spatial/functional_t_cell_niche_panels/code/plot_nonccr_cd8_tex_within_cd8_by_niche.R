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
  "final", "extended_figure_5", "t_cell_composition_by_niche", "outputs",
  "extended_figure5_tex_within_cd8_tregccr8_within_treg_by_patient_niche_wide.csv"
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

dt <- fread(source_path)
dt <- dt[cCR == "non-cCR" & neighborhood_label != "Niche_Epi"]
dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]
dt[, cd8_tex_percent_within_cd8 := tex_percent_within_cd8]
dt[, denominator_cd8_cells := total_cd8_cells]
dt[, target_cd8_tex_cells := tex_cells]

plot_dt <- dt[!is.na(cd8_tex_percent_within_cd8)]

summary_dt <- dt[
  ,
  .(
    n_defined_patient_niches = sum(!is.na(cd8_tex_percent_within_cd8)),
    n_total_patient_niches = .N,
    median_percent = median(cd8_tex_percent_within_cd8, na.rm = TRUE),
    mean_percent = mean(cd8_tex_percent_within_cd8, na.rm = TRUE),
    max_percent = max(cd8_tex_percent_within_cd8, na.rm = TRUE),
    median_denominator_cd8_cells = as.numeric(median(denominator_cd8_cells[!is.na(cd8_tex_percent_within_cd8)], na.rm = TRUE)),
    median_target_cd8_tex_cells = as.numeric(median(target_cd8_tex_cells[!is.na(cd8_tex_percent_within_cd8)], na.rm = TRUE))
  ),
  by = neighborhood_label
]

wide <- dcast(plot_dt, patientID ~ neighborhood_label, value.var = "cd8_tex_percent_within_cd8")
required_cols <- c("patientID", neighborhood_order)
wide_complete <- wide[complete.cases(wide[, required_cols, with = FALSE])]
friedman_p <- if (nrow(wide_complete) >= 3) {
  suppressWarnings(friedman.test(as.matrix(wide_complete[, neighborhood_order, with = FALSE]))$p.value)
} else {
  NA_real_
}
kruskal_p <- suppressWarnings(kruskal.test(cd8_tex_percent_within_cd8 ~ neighborhood_label, data = plot_dt)$p.value)
global_stats <- data.table(
  feature = "CD8_Tex_PDCD1 within CD8",
  filter = "non-cCR only",
  n_complete_patients_for_friedman = nrow(wide_complete),
  friedman_p = friedman_p,
  kruskal_p = kruskal_p
)

order_dt <- summary_dt[order(median_percent, as.character(neighborhood_label))]
plot_levels <- as.character(order_dt$neighborhood_label)
plot_dt[, neighborhood_plot := factor(as.character(neighborhood_label), levels = plot_levels)]

fwrite(dt, file.path(out_dir, "postc_nonccr_cd8_tex_within_cd8_no_epi_by_patient_niche.csv"))
fwrite(summary_dt, file.path(out_dir, "postc_nonccr_cd8_tex_within_cd8_no_epi_summary.csv"))
fwrite(global_stats, file.path(out_dir, "postc_nonccr_cd8_tex_within_cd8_no_epi_global_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

y_top <- 108
label_dt <- data.table(
  x = ceiling(length(plot_levels) / 2),
  y = 103,
  label = paste0("Friedman ", format_p(friedman_p))
)

plot <- ggplot(plot_dt, aes(x = neighborhood_plot, y = cd8_tex_percent_within_cd8)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(size = pmin(denominator_cd8_cells, 200)),
    position = position_jitter(width = 0.08, height = 0),
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
  scale_size_continuous(range = c(0.55, 1.45), guide = "none") +
  scale_y_continuous(limits = c(0, y_top), expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "CD8_Tex_PDCD1 within CD8 (%)"
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
  file.path(out_dir, "postc_nonccr_cd8_tex_within_cd8_by_niche_boxplot.pdf"),
  plot,
  device = "pdf",
  width = 90,
  height = 72,
  units = "mm",
  dpi = 300
)

ggsave(
  file.path(out_dir, "postc_nonccr_cd8_tex_within_cd8_no_epi_by_niche_boxplot.pdf"),
  plot,
  device = "pdf",
  width = 85,
  height = 72,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Non-cCR CD8_Tex_PDCD1 within CD8 across spatial niches",
    "",
    "Exploratory plot for post-chemoradiotherapy non-cCR samples only, excluding Niche_Epi.",
    "For each patient-niche, CD8_Tex_PDCD1 is measured as a percentage of CD8-lineage cells: CD8_Teff, CD8_Tex_PDCD1 and CD8_prolif.",
    "Point size is proportional to the CD8-lineage denominator cell count, capped at 200 cells for readability.",
    "Patient-niches with zero CD8-lineage denominator cells are undefined and are not plotted.",
    "X-axis order follows low-to-high median percentage."
  ),
  file.path(legend_dir, "legend_nonccr_cd8_tex_within_cd8_by_niche.md")
)
