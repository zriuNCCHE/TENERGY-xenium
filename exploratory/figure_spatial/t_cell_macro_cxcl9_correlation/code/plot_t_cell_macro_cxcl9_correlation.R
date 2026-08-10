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

overall_spearman <- safe_cor(plot_dt, "t_cell_percent", "macro_cxcl9_percent", "spearman")
overall_pearson <- safe_cor(plot_dt, "t_cell_percent", "macro_cxcl9_percent", "pearson")
overall_area_spearman <- safe_cor(plot_dt, "t_area_percent", "macro_cxcl9_area_percent", "spearman")

overall_stats <- data.table(
  comparison = c(
    "total_T_cell_percent_vs_Macro_CXCL9_cell_percent",
    "total_T_cell_percent_vs_Macro_CXCL9_cell_percent",
    "total_T_area_percent_vs_Macro_CXCL9_area_percent"
  ),
  scope = "all_patient_niche_pairs",
  method = c("Spearman", "Pearson", "Spearman"),
  n = c(overall_spearman$n, overall_pearson$n, overall_area_spearman$n),
  estimate = c(overall_spearman$estimate, overall_pearson$estimate, overall_area_spearman$estimate),
  p_value = c(overall_spearman$p_value, overall_pearson$p_value, overall_area_spearman$p_value)
)

niche_stats <- rbindlist(lapply(neighborhood_order, function(one_niche) {
  niche_dt <- plot_dt[neighborhood_label == one_niche]
  spearman_result <- safe_cor(niche_dt, "t_cell_percent", "macro_cxcl9_percent", "spearman")
  pearson_result <- safe_cor(niche_dt, "t_cell_percent", "macro_cxcl9_percent", "pearson")
  data.table(
    comparison = "total_T_cell_percent_vs_Macro_CXCL9_cell_percent",
    scope = one_niche,
    method = c("Spearman", "Pearson"),
    n = c(spearman_result$n, pearson_result$n),
    estimate = c(spearman_result$estimate, pearson_result$estimate),
    p_value = c(spearman_result$p_value, pearson_result$p_value)
  )
}))
niche_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = method]
overall_stats[, p_adjusted_bh := p_value]

stats_dt <- rbindlist(list(overall_stats, niche_stats), use.names = TRUE)
fwrite(plot_dt, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_correlation_by_patient_niche.csv"))
fwrite(stats_dt, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_correlation_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

label_text <- paste0(
  "Spearman rho = ", sprintf("%.2f", overall_spearman$estimate),
  "\np = ", format_p(overall_spearman$p_value)
)

scatter_plot <- ggplot(plot_dt, aes(x = t_cell_percent, y = macro_cxcl9_percent)) +
  geom_point(
    aes(color = neighborhood_label),
    size = 1.15,
    alpha = 0.82
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#333333",
    fill = "#BDBDBD",
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
  labs(
    x = "Total T cells in assigned niche (%)",
    y = "Macro_CXCL9 cells in assigned niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = small_size),
    legend.key.height = unit(0.28, "cm"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_total_t_cells_macro_cxcl9_correlation_scatter.pdf"),
  plot = scatter_plot,
  device = "pdf",
  width = 95,
  height = 82,
  units = "mm",
  dpi = 300
)

facet_plot <- ggplot(plot_dt, aes(x = t_cell_percent, y = macro_cxcl9_percent)) +
  geom_point(
    color = "#333333",
    size = 0.9,
    alpha = 0.72
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "#2E7D32",
    linewidth = line_width
  ) +
  facet_wrap(~ neighborhood_label, ncol = 4, scales = "free") +
  labs(
    x = "Total T cells in assigned niche (%)",
    y = "Macro_CXCL9 cells in assigned niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.22, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_total_t_cells_macro_cxcl9_correlation_by_niche.pdf"),
  plot = facet_plot,
  device = "pdf",
  width = 150,
  height = 105,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# PostC total T-cell and Macro_CXCL9 niche correlation",
    "",
    "Exploratory analysis based on the highest-score spatial neighborhood assignment.",
    "Each point in the main scatter plot is one patient-niche pair.",
    "Total T-cell abundance and Macro_CXCL9 abundance were both calculated as percentages of all cells assigned to the same niche in the same patient.",
    paste0("Overall Spearman rho = ", sprintf("%.3f", overall_spearman$estimate), ", P = ", format_p(overall_spearman$p_value), "."),
    "Niche-specific statistics are saved in the correlation statistics table."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
