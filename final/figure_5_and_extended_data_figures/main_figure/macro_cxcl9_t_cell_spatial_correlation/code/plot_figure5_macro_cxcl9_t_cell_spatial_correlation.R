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

neighborhood_key <- c(
  "Niche_Mono" = "neighborhood_1_k7",
  "Niche_Neutro" = "neighborhood_2_k7",
  "Niche_Epi" = "neighborhood_3_k7",
  "Niche_iCAF_CXCL6" = "neighborhood_4_k7",
  "Niche_T_cell" = "neighborhood_5_k7",
  "Niche_Macro_CXCL9" = "neighborhood_6_k7",
  "Niche_iCAF_CXCL5" = "neighborhood_7_k7"
)
neighborhood_order <- c(
  "Niche_iCAF_CXCL5",
  "Niche_Mono",
  "Niche_Neutro",
  "Niche_Epi",
  "Niche_Macro_CXCL9",
  "Niche_iCAF_CXCL6",
  "Niche_T_cell"
)
neighborhood_colors <- unlist(palette_config$spatial_niches_functional_overlay[neighborhood_key], use.names = FALSE)
names(neighborhood_colors) <- names(neighborhood_key)

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

safe_cor <- function(dt, x_col, y_col, method = "spearman") {
  dt <- dt[!is.na(get(x_col)) & !is.na(get(y_col))]
  if (nrow(dt) < 3 || sd(dt[[x_col]]) == 0 || sd(dt[[y_col]]) == 0) {
    return(list(estimate = NA_real_, p_value = NA_real_, n = nrow(dt)))
  }
  result <- suppressWarnings(cor.test(dt[[x_col]], dt[[y_col]], method = method, exact = FALSE))
  list(
    estimate = unname(result$estimate),
    p_value = result$p.value,
    n = nrow(dt)
  )
}

t_dt <- fread(t_path)[
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
macro_dt <- fread(macro_path)[
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

overall_spearman <- safe_cor(plot_dt, "macro_cxcl9_percent", "t_cell_percent", "spearman")
overall_pearson <- safe_cor(plot_dt, "macro_cxcl9_percent", "t_cell_percent", "pearson")
area_spearman <- safe_cor(plot_dt, "macro_cxcl9_area_percent", "t_area_percent", "spearman")

stats_dt <- data.table(
  comparison = c(
    "Macro_CXCL9_cell_percent_vs_total_T_cell_percent",
    "Macro_CXCL9_cell_percent_vs_total_T_cell_percent",
    "Macro_CXCL9_area_percent_vs_total_T_area_percent"
  ),
  scope = "patient_niche_observations",
  method = c("Spearman", "Pearson", "Spearman"),
  n = c(overall_spearman$n, overall_pearson$n, area_spearman$n),
  estimate = c(overall_spearman$estimate, overall_pearson$estimate, area_spearman$estimate),
  p_value = c(overall_spearman$p_value, overall_pearson$p_value, area_spearman$p_value)
)
stats_dt[, p_adjusted_bh := p.adjust(p_value, method = "BH")]

fwrite(plot_dt, file.path(out_dir, "figure5_macro_cxcl9_t_cell_spatial_correlation_values.csv"))
fwrite(stats_dt, file.path(out_dir, "figure5_macro_cxcl9_t_cell_spatial_correlation_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

label_text <- paste0(
  "Spearman rho = ", sprintf("%.2f", overall_spearman$estimate),
  "\np = ", format_p(overall_spearman$p_value)
)

axis_upper <- ceiling(max(c(plot_dt$macro_cxcl9_percent, plot_dt$t_cell_percent), na.rm = TRUE) / 10) * 10
axis_limits <- c(0, axis_upper)

scatter_plot <- ggplot(plot_dt, aes(x = macro_cxcl9_percent, y = t_cell_percent)) +
  geom_point(
    aes(color = neighborhood_label),
    size = 1.0,
    alpha = 0.88
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#555555",
    fill = "#D9D9D9",
    linewidth = line_width
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_text,
    hjust = 1.05,
    vjust = 1.15,
    size = 1.9
  ) +
  scale_color_manual(values = neighborhood_colors, drop = FALSE, name = "Niche") +
  scale_x_continuous(limits = axis_limits, expand = expansion(mult = c(0, 0.04))) +
  scale_y_continuous(limits = axis_limits, expand = expansion(mult = c(0, 0.04))) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Macro_CXCL9 cells in assigned niche (%)",
    y = "Total T cells in assigned niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = small_size),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.34, "cm"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "figure5_macro_cxcl9_t_cell_spatial_correlation.pdf"),
  plot = scatter_plot,
  device = "pdf",
  width = 108,
  height = 108,
  units = "mm",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "figure5_macro_cxcl9_t_cell_spatial_correlation.png"),
  plot = scatter_plot,
  device = "png",
  width = 108,
  height = 108,
  units = "mm",
  dpi = 300
)

legend_text <- paste0(
  "# Figure 5 Legend Draft - Macro_CXCL9 and T-cell spatial association\n\n",
  "Macro_CXCL9 and total T-cell fractions were calculated for each patient within each post-chemoradiotherapy spatial niche. ",
  "Each point represents one patient-niche observation, with fractions defined as the number of indicated cells divided by the total number of cells assigned to that niche. ",
  "The association between Macro_CXCL9 and total T-cell fractions was evaluated using two-sided Spearman correlation across patient-niche observations ",
  "(n = ", overall_spearman$n, "; rho = ", sprintf('%.2f', overall_spearman$estimate), "; p = ", format_p(overall_spearman$p_value), ")."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))
