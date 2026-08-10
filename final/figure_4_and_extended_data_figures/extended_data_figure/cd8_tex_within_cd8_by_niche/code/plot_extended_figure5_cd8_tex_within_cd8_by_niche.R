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
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

neighborhood_key <- c(
  "Niche_Mono" = "neighborhood_1_k7",
  "Niche_Neutro" = "neighborhood_2_k7",
  "Niche_Epi" = "neighborhood_3_k7",
  "Niche_iCAF_CXCL6" = "neighborhood_4_k7",
  "Niche_T_cell" = "neighborhood_5_k7",
  "Niche_Macro_CXCL9" = "neighborhood_6_k7",
  "Niche_iCAF_CXCL5" = "neighborhood_7_k7"
)
neighborhood_colors <- unlist(palette_config$spatial_niches_functional_overlay[neighborhood_key], use.names = FALSE)
names(neighborhood_colors) <- names(neighborhood_key)

format_p <- function(p_value, prefix = "p") {
  if (is.na(p_value)) {
    return(paste0(prefix, " = NA"))
  }
  if (p_value < 1e-4) {
    return(sprintf("%s = %.1e", prefix, p_value))
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("%s = %.6f", prefix, p_value)))
  }
  sprintf("%s = %.3f", prefix, p_value)
}

dt <- fread(source_path)
dt <- dt[neighborhood_label != "Niche_Epi"]
dt[, tex_percent_within_cd8 := as.numeric(tex_percent_within_cd8)]
dt[, tex_area_percent_within_cd8 := as.numeric(tex_area_percent_within_cd8)]
dt[, total_cd8_cells := as.integer(total_cd8_cells)]
dt[, tex_cells := as.integer(tex_cells)]

plot_dt <- dt[
  total_cd8_cells > 0 & !is.na(tex_percent_within_cd8),
  .(
    patientID,
    sampleID,
    cCR,
    neighborhood_label,
    denominator_cd8_cells = total_cd8_cells,
    cd8_tex_pdcd1_cells = tex_cells,
    cd8_tex_pdcd1_percent = tex_percent_within_cd8,
    cd8_tex_pdcd1_area_percent = tex_area_percent_within_cd8
  )
]

order_dt <- plot_dt[
  ,
  .(
    n_defined_patient_niches = .N,
    median_percent = median(cd8_tex_pdcd1_percent, na.rm = TRUE),
    mean_percent = mean(cd8_tex_pdcd1_percent, na.rm = TRUE),
    max_percent = max(cd8_tex_pdcd1_percent, na.rm = TRUE),
    median_denominator_cd8_cells = median(denominator_cd8_cells, na.rm = TRUE),
    median_cd8_tex_pdcd1_cells = median(cd8_tex_pdcd1_cells, na.rm = TRUE)
  ),
  by = .(neighborhood_label)
][order(median_percent, neighborhood_label)]
order_dt[, x_order := seq_len(.N)]
ordered_niches <- order_dt$neighborhood_label

plot_dt[, neighborhood_label := factor(neighborhood_label, levels = ordered_niches)]

kruskal_p <- if (uniqueN(plot_dt$neighborhood_label) >= 2) {
  suppressWarnings(kruskal.test(cd8_tex_pdcd1_percent ~ neighborhood_label, data = plot_dt)$p.value)
} else {
  NA_real_
}
stats_dt <- data.table(
  feature = "CD8_Tex_PDCD1 within CD8",
  denominator = "all CD8 T cells",
  scope = "non-epithelial niches with defined CD8 denominator",
  test = "Kruskal-Wallis test across non-epithelial niches",
  n_defined_patient_niches = nrow(plot_dt),
  n_patients = uniqueN(plot_dt$patientID),
  p_value = kruskal_p
)

patient_summary <- plot_dt[
  ,
  .(
    n_niches_with_cd8 = .N,
    total_cd8_cells = sum(denominator_cd8_cells, na.rm = TRUE),
    total_cd8_tex_pdcd1_cells = sum(cd8_tex_pdcd1_cells, na.rm = TRUE),
    weighted_cd8_tex_pdcd1_percent = sum(cd8_tex_pdcd1_cells, na.rm = TRUE) / sum(denominator_cd8_cells, na.rm = TRUE) * 100,
    max_niche_cd8_tex_pdcd1_percent = max(cd8_tex_pdcd1_percent, na.rm = TRUE),
    median_niche_cd8_tex_pdcd1_percent = median(cd8_tex_pdcd1_percent, na.rm = TRUE),
    niches_with_at_least_20_tex = sum(cd8_tex_pdcd1_cells >= 20, na.rm = TRUE),
    niches_with_at_least_40_percent = sum(cd8_tex_pdcd1_percent >= 40, na.rm = TRUE)
  ),
  by = .(patientID, sampleID, cCR)
]
patient_summary[, reviewer_image_score := scale(log1p(total_cd8_tex_pdcd1_cells))[, 1] +
  scale(weighted_cd8_tex_pdcd1_percent)[, 1] +
  scale(niches_with_at_least_20_tex)[, 1] +
  0.5 * scale(niches_with_at_least_40_percent)[, 1]]
setorder(patient_summary, -reviewer_image_score)

niche_candidate_summary <- plot_dt[
  ,
  .(
    denominator_cd8_cells = sum(denominator_cd8_cells, na.rm = TRUE),
    cd8_tex_pdcd1_cells = sum(cd8_tex_pdcd1_cells, na.rm = TRUE),
    weighted_cd8_tex_pdcd1_percent = sum(cd8_tex_pdcd1_cells, na.rm = TRUE) / sum(denominator_cd8_cells, na.rm = TRUE) * 100,
    median_patient_percent = median(cd8_tex_pdcd1_percent, na.rm = TRUE)
  ),
  by = .(patientID, sampleID, cCR, neighborhood_label)
]
setorder(niche_candidate_summary, -cd8_tex_pdcd1_cells, -weighted_cd8_tex_pdcd1_percent)

fwrite(plot_dt, file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_by_niche_values_no_epi.csv"))
fwrite(order_dt, file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_by_niche_summary_no_epi.csv"))
fwrite(stats_dt, file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_by_niche_stats_no_epi.csv"))
fwrite(patient_summary, file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_image_candidate_patient_ranking.csv"))
fwrite(niche_candidate_summary, file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_image_candidate_patient_niche_ranking.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
y_top <- ceiling(max(plot_dt$cd8_tex_pdcd1_percent, na.rm = TRUE) / 10) * 10
y_limit <- c(0, max(80, y_top))

plot_obj <- ggplot(plot_dt, aes(x = neighborhood_label, y = cd8_tex_pdcd1_percent)) +
  geom_boxplot(
    fill = NA,
    color = "#4D4D4D",
    width = 0.55,
    outlier.shape = NA,
    linewidth = line_width
  ) +
  geom_point(
    aes(color = neighborhood_label),
    size = 0.85,
    alpha = 0.78,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  annotate(
    "text",
    x = Inf,
    y = y_limit[2],
    label = paste0("Kruskal ", format_p(kruskal_p, "p")),
    hjust = 1.03,
    vjust = 1.15,
    size = 1.9
  ) +
  scale_color_manual(values = neighborhood_colors, drop = FALSE, guide = "none") +
  scale_y_continuous(limits = y_limit, expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "CD8_Tex_PDCD1 cells in CD8 cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = small_size),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_by_niche_no_epi_ordered.pdf"),
  plot = plot_obj,
  device = "pdf",
  width = 92,
  height = 78,
  units = "mm",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "extended_figure5_cd8_tex_within_cd8_by_niche_no_epi_ordered.png"),
  plot = plot_obj,
  device = "png",
  width = 92,
  height = 78,
  units = "mm",
  dpi = 300
)

legend_text <- paste0(
  "# Extended Data Figure 5 Legend Draft - CD8_Tex_PDCD1 within CD8 across spatial niches\n\n",
  "Post-chemoradiotherapy cells were assigned to a unique spatial niche by selecting the highest k = 7 neighborhood score. ",
  "For each patient and non-epithelial niche, the proportion of CD8_Tex_PDCD1 cells was calculated among all CD8 T cells ",
  "(CD8_Teff, CD8_Tex_PDCD1 and CD8_prolif) assigned to the same niche. ",
  "Niche_Epi and patient-niche observations without CD8 cells were excluded. Niches are ordered from lower to higher median proportion. ",
  "The global P value was calculated using a Kruskal-Wallis test across defined patient-niche observations."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))
