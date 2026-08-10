suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

data_path <- file.path(project_root, "data", "geo_ready", "spatial_neighborhoods", "postC_k7_cell_niche_assignments.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

t_cell_types <- c(
  "T/NK",
  "CD8_Teff",
  "CD8_Tex_PDCD1",
  "CD8_prolif",
  "CD4_Treg_CCR8",
  "CD4_Treg_FOXP3",
  "CD4_CXCL13",
  "CD4_prolif"
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

cells <- fread(
  data_path,
  select = c("cell_id", "sampleID", "patientID", "cCR", "sample_timepoint", "cell_area", "final_cell_type2")
)
cells <- cells[
  sample_timepoint == "postC" &
    !is.na(final_cell_type2) &
    !is.na(cell_area) &
    cell_area > 0
]

sample_summary <- cells[
  ,
  .(
    total_cells = .N,
    total_area = sum(cell_area),
    total_t_cells = sum(final_cell_type2 %in% t_cell_types),
    total_t_area = sum(cell_area[final_cell_type2 %in% t_cell_types]),
    macro_cxcl9_cells = sum(final_cell_type2 == "Macro_CXCL9"),
    macro_cxcl9_area = sum(cell_area[final_cell_type2 == "Macro_CXCL9"])
  ),
  by = .(patientID, sampleID, cCR)
]

sample_summary[, total_t_cell_percent := total_t_cells / total_cells * 100]
sample_summary[, macro_cxcl9_cell_percent := macro_cxcl9_cells / total_cells * 100]
sample_summary[, total_t_area_percent := total_t_area / total_area * 100]
sample_summary[, macro_cxcl9_area_percent := macro_cxcl9_area / total_area * 100]
sample_summary[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

spearman_cell <- safe_cor(sample_summary, "total_t_cell_percent", "macro_cxcl9_cell_percent", "spearman")
pearson_cell <- safe_cor(sample_summary, "total_t_cell_percent", "macro_cxcl9_cell_percent", "pearson")
spearman_area <- safe_cor(sample_summary, "total_t_area_percent", "macro_cxcl9_area_percent", "spearman")
pearson_area <- safe_cor(sample_summary, "total_t_area_percent", "macro_cxcl9_area_percent", "pearson")

stats_dt <- data.table(
  comparison = c(
    "total_T_cell_percent_vs_Macro_CXCL9_cell_percent",
    "total_T_cell_percent_vs_Macro_CXCL9_cell_percent",
    "total_T_area_percent_vs_Macro_CXCL9_area_percent",
    "total_T_area_percent_vs_Macro_CXCL9_area_percent"
  ),
  method = c("Spearman", "Pearson", "Spearman", "Pearson"),
  n = c(spearman_cell$n, pearson_cell$n, spearman_area$n, pearson_area$n),
  estimate = c(spearman_cell$estimate, pearson_cell$estimate, spearman_area$estimate, pearson_area$estimate),
  p_value = c(spearman_cell$p_value, pearson_cell$p_value, spearman_area$p_value, pearson_area$p_value)
)
stats_dt[, p_adjusted_bh := p.adjust(p_value, method = "BH")]

fwrite(sample_summary, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_sample_level_summary.csv"))
fwrite(stats_dt, file.path(out_dir, "postc_total_t_cells_macro_cxcl9_sample_level_correlation_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

label_text <- paste0(
  "Spearman rho = ", sprintf("%.2f", spearman_cell$estimate),
  "\np = ", format_p(spearman_cell$p_value)
)

plot <- ggplot(sample_summary, aes(x = total_t_cell_percent, y = macro_cxcl9_cell_percent)) +
  geom_point(
    aes(color = cCR),
    size = 1.8,
    alpha = 0.9
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
  scale_color_manual(values = outcome_colors, drop = FALSE, name = "Clinical Outcome") +
  labs(
    x = "Total T cells in sample (%)",
    y = "Macro_CXCL9 cells in sample (%)"
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
  filename = file.path(out_dir, "postc_total_t_cells_macro_cxcl9_sample_level_correlation.pdf"),
  plot = plot,
  device = "pdf",
  width = 82,
  height = 72,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# PostC sample-level total T-cell and Macro_CXCL9 correlation",
    "",
    "Exploratory sample-level analysis. Each point is one post-chemoradiotherapy sample/patient.",
    "Total T-cell abundance was calculated as all listed T-cell subsets divided by all annotated cells in the sample.",
    "Macro_CXCL9 abundance was calculated as Macro_CXCL9 cells divided by all annotated cells in the sample.",
    paste0("Cell-fraction Spearman rho = ", sprintf("%.3f", spearman_cell$estimate), ", P = ", format_p(spearman_cell$p_value), ".")
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
