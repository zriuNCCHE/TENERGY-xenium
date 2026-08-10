suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

t_path <- file.path(
  project_root,
  "exploratory", "figure_spatial", "t_cell_neighborhood_abundance", "outputs",
  "postc_total_t_cells_neighborhood_abundance_by_patient.csv"
)
macro_path <- file.path(
  project_root,
  "exploratory", "figure_spatial", "macrophage_neighborhood_abundance", "outputs",
  "postc_macrophage_neighborhood_abundance_by_patient.csv"
)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

neighborhood_order <- c(
  "Niche_iCAF_CXCL5",
  "Niche_Mono",
  "Niche_Neutro",
  "Niche_Epi",
  "Niche_Macro_CXCL9",
  "Niche_iCAF_CXCL6",
  "Niche_T_cell"
)
neighborhood_colors <- c(
  "Niche_iCAF_CXCL5" = "#E6550D",
  "Niche_Mono" = "#6A51A3",
  "Niche_Neutro" = "#D95F8D",
  "Niche_Epi" = "#56B4E9",
  "Niche_Macro_CXCL9" = "#2E7D32",
  "Niche_iCAF_CXCL6" = "#FD8D3C",
  "Niche_T_cell" = "#4C78A8"
)
outcome_colors <- unlist(palette_config$clinical_outcome)[c("cCR", "non-cCR")]

format_p <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 1e-4) {
    return(formatC(p_value, format = "e", digits = 2))
  }
  if (p_value < 0.001) {
    return(sprintf("%.6f", p_value))
  }
  sprintf("%.3f", p_value)
}

safe_patient_spearman <- function(dt) {
  keep <- !is.na(dt$t_cell_percent) & !is.na(dt$macro_cxcl9_percent)
  dt <- dt[keep]
  if (nrow(dt) < 3 || sd(dt$t_cell_percent) == 0 || sd(dt$macro_cxcl9_percent) == 0) {
    return(NA_real_)
  }
  suppressWarnings(cor(dt$t_cell_percent, dt$macro_cxcl9_percent, method = "spearman"))
}

t_dt <- fread(t_path)
macro_dt <- fread(macro_path)

t_dt <- t_dt[
  ,
  .(
    patientID,
    sampleID,
    cCR,
    neighborhood_label,
    t_cell_percent = cell_percent,
    t_area_percent = area_percent
  )
]
macro_dt <- macro_dt[
  feature == "Macro_CXCL9",
  .(
    patientID,
    sampleID,
    cCR,
    neighborhood_label,
    macro_cxcl9_percent = cell_percent,
    macro_cxcl9_area_percent = area_percent
  )
]

plot_dt <- merge(
  t_dt,
  macro_dt,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all = FALSE
)
plot_dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]
plot_dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

niche_summary <- plot_dt[
  ,
  .(
    n_patients = uniqueN(patientID),
    median_t_cell_percent = median(t_cell_percent, na.rm = TRUE),
    mean_t_cell_percent = mean(t_cell_percent, na.rm = TRUE),
    median_macro_cxcl9_percent = median(macro_cxcl9_percent, na.rm = TRUE),
    mean_macro_cxcl9_percent = mean(macro_cxcl9_percent, na.rm = TRUE)
  ),
  by = .(neighborhood_label)
]
niche_summary[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]

niche_spearman <- suppressWarnings(cor.test(
  niche_summary$median_t_cell_percent,
  niche_summary$median_macro_cxcl9_percent,
  method = "spearman",
  exact = FALSE
))

patient_rho <- plot_dt[
  ,
  .(
    rho = safe_patient_spearman(.SD),
    n_niches = sum(!is.na(t_cell_percent) & !is.na(macro_cxcl9_percent))
  ),
  by = .(patientID, sampleID, cCR)
]
patient_rho <- patient_rho[!is.na(rho)]
rho_test <- suppressWarnings(wilcox.test(patient_rho$rho, mu = 0, alternative = "greater", exact = FALSE))
rho_summary <- data.table(
  analysis = "patient_level_niche_rank_concordance",
  test = "one-sided Wilcoxon signed-rank test vs rho = 0",
  alternative = "rho greater than 0",
  n_patients = nrow(patient_rho),
  median_rho = median(patient_rho$rho),
  mean_rho = mean(patient_rho$rho),
  p_value = rho_test$p.value
)
niche_stats <- data.table(
  analysis = "niche_median_correlation_visual_summary",
  test = "Spearman correlation across seven niche medians",
  alternative = "two-sided",
  n_niches = nrow(niche_summary),
  rho = unname(niche_spearman$estimate),
  p_value = niche_spearman$p.value
)

fwrite(plot_dt, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_patient_niche_values.csv"))
fwrite(niche_summary, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_niche_median_summary.csv"))
fwrite(patient_rho, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_patient_level_rho.csv"))
fwrite(rho_summary, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_patient_level_rho_stats.csv"))
fwrite(niche_stats, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_niche_median_correlation_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

niche_label_text <- paste0(
  "rho = ", sprintf("%.2f", unname(niche_spearman$estimate)),
  "\np = ", format_p(niche_spearman$p.value)
)

niche_plot <- ggplot(
  niche_summary,
  aes(x = median_t_cell_percent, y = median_macro_cxcl9_percent)
) +
  geom_point(
    aes(color = neighborhood_label),
    size = 2.1,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = neighborhood_label),
    size = 1.8,
    hjust = -0.05,
    vjust = 0.5,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = niche_label_text,
    hjust = 1.05,
    vjust = 1.15,
    size = 1.9
  ) +
  scale_color_manual(values = neighborhood_colors, drop = FALSE, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    x = "Median total T cells in niche (%)",
    y = "Median Macro_CXCL9 cells in niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_total_t_cells_macro_cxcl9_niche_median_scatter.pdf"),
  plot = niche_plot,
  device = "pdf",
  width = 92,
  height = 78,
  units = "mm",
  dpi = 300
)

rho_label_text <- paste0(
  "median rho = ", sprintf("%.2f", rho_summary$median_rho),
  "\np = ", format_p(rho_summary$p_value)
)

rho_plot <- ggplot(patient_rho, aes(x = "Patient-level rho", y = rho)) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = line_width) +
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.07, height = 0),
    size = 1.25,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  annotate(
    "text",
    x = 1,
    y = max(patient_rho$rho, na.rm = TRUE) * 1.05,
    label = rho_label_text,
    size = 1.9,
    lineheight = 0.9
  ) +
  scale_color_manual(values = outcome_colors, drop = FALSE, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(min(-0.1, min(patient_rho$rho, na.rm = TRUE) * 1.08), max(patient_rho$rho, na.rm = TRUE) * 1.18),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Within-patient niche Spearman rho"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = small_size),
    legend.key.height = unit(0.28, "cm"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_total_t_cells_macro_cxcl9_patient_level_rho_boxplot.pdf"),
  plot = rho_plot,
  device = "pdf",
  width = 70,
  height = 78,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Reviewer-safe T-cell and Macro_CXCL9 niche concordance",
    "",
    "This exploratory analysis separates visualization from formal inference.",
    "The niche-median scatter shows seven points, one per spatial niche, using median abundance across patients.",
    "Formal inference uses one Spearman rho per patient: within each patient, total T-cell abundance across the seven niches was correlated with Macro_CXCL9 abundance across the same seven niches.",
    "The distribution of patient-level rho values was tested against zero using a one-sided Wilcoxon signed-rank test.",
    paste0("Niche-median Spearman rho = ", sprintf("%.3f", unname(niche_spearman$estimate)), ", P = ", format_p(niche_spearman$p.value), "."),
    paste0("Patient-level median rho = ", sprintf("%.3f", rho_summary$median_rho), ", one-sided Wilcoxon P = ", format_p(rho_summary$p_value), ".")
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
