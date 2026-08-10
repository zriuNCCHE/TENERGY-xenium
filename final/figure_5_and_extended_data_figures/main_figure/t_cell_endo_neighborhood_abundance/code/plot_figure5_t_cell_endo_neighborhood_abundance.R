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
endo_path <- file.path(
  project_root,
  "exploratory", "figure_spatial", "vascular_neighborhood_abundance", "outputs",
  "postc_endo_plvap_neighborhood_abundance_by_patient.csv"
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
plot_neighborhoods <- setdiff(names(neighborhood_key), "Niche_Epi")

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

safe_friedman <- function(dt, value_col, ordered_niches) {
  wide <- dcast(dt, patientID ~ neighborhood_label, value.var = value_col, fill = NA_real_)
  complete_cols <- c("patientID", ordered_niches)
  wide <- wide[complete.cases(wide[, ..complete_cols])]
  if (nrow(wide) < 3) {
    return(list(p_value = NA_real_, n_patients = nrow(wide)))
  }
  value_matrix <- as.matrix(wide[, ..ordered_niches])
  if (all(apply(value_matrix, 2, sd) == 0)) {
    return(list(p_value = NA_real_, n_patients = nrow(wide)))
  }
  p_value <- suppressWarnings(friedman.test(value_matrix)$p.value)
  list(p_value = p_value, n_patients = nrow(wide))
}

prepare_feature <- function(path, feature_label, feature_type_label) {
  dt <- fread(path)
  dt <- dt[neighborhood_label %in% plot_neighborhoods]
  dt[
    ,
    .(
      patientID,
      sampleID,
      cCR,
      neighborhood_label,
      cell_percent,
      area_percent,
      feature = feature_label,
      feature_type = feature_type_label
    )
  ]
}

t_dt <- prepare_feature(t_path, "Total T cells", "T-cell abundance")
endo_dt <- prepare_feature(endo_path, "endo_PLVAP", "Endothelial abundance")
plot_dt <- rbindlist(list(t_dt, endo_dt), use.names = TRUE)

ordered_labels <- rbindlist(lapply(unique(plot_dt$feature), function(one_feature) {
  feature_dt <- plot_dt[feature == one_feature]
  ordered <- feature_dt[
    ,
    .(median_cell_percent = median(cell_percent, na.rm = TRUE)),
    by = .(neighborhood_label)
  ][order(median_cell_percent, neighborhood_label)]
  ordered[, feature := one_feature]
  ordered[, x_order := seq_len(.N)]
  ordered[]
}))

plot_dt <- merge(plot_dt, ordered_labels[, .(feature, neighborhood_label, x_order)], by = c("feature", "neighborhood_label"))
plot_dt[, x_label := factor(neighborhood_label, levels = ordered_labels[order(feature, x_order), unique(neighborhood_label)])]

stats_dt <- rbindlist(lapply(unique(plot_dt$feature), function(one_feature) {
  feature_dt <- plot_dt[feature == one_feature]
  ordered_niches <- ordered_labels[feature == one_feature][order(x_order), neighborhood_label]
  friedman_result <- safe_friedman(feature_dt, "cell_percent", ordered_niches)
  data.table(
    feature = one_feature,
    metric = "cell_percent",
    test = "paired Friedman test across six non-epithelial neighborhoods",
    n_patients = friedman_result$n_patients,
    p_value = friedman_result$p_value
  )
}))
stats_dt[, p_adjusted_bh := p.adjust(p_value, method = "BH")]

summary_dt <- plot_dt[
  ,
  .(
    n_patients = uniqueN(patientID),
    median_cell_percent = median(cell_percent, na.rm = TRUE),
    mean_cell_percent = mean(cell_percent, na.rm = TRUE),
    max_cell_percent = max(cell_percent, na.rm = TRUE),
    median_area_percent = median(area_percent, na.rm = TRUE),
    mean_area_percent = mean(area_percent, na.rm = TRUE),
    max_area_percent = max(area_percent, na.rm = TRUE)
  ),
  by = .(feature, neighborhood_label)
]
summary_dt <- merge(summary_dt, ordered_labels[, .(feature, neighborhood_label, x_order)], by = c("feature", "neighborhood_label"))
setorder(summary_dt, feature, x_order)

fwrite(plot_dt, file.path(out_dir, "figure5_t_cell_endo_neighborhood_abundance_values_no_epi.csv"))
fwrite(summary_dt, file.path(out_dir, "figure5_t_cell_endo_neighborhood_abundance_summary_no_epi.csv"))
fwrite(stats_dt, file.path(out_dir, "figure5_t_cell_endo_neighborhood_abundance_friedman_stats_no_epi.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

make_plot <- function(feature_name, y_label, file_stem, y_limit = NULL) {
  one_dt <- copy(plot_dt[feature == feature_name])
  one_order <- ordered_labels[feature == feature_name][order(x_order), neighborhood_label]
  one_dt[, neighborhood_label := factor(neighborhood_label, levels = one_order)]
  one_stats <- stats_dt[feature == feature_name]
  y_top <- max(one_dt$cell_percent, na.rm = TRUE)
  if (is.null(y_limit)) {
    y_limit <- c(0, ceiling((y_top * 1.18) / 5) * 5)
  }
  p_label <- paste0("Friedman ", format_p(one_stats$p_value, "p"))

  p <- ggplot(one_dt, aes(x = neighborhood_label, y = cell_percent)) +
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
      label = p_label,
      hjust = 1.03,
      vjust = 1.15,
      size = 1.9
    ) +
    scale_color_manual(values = neighborhood_colors, drop = FALSE, guide = "none") +
    scale_y_continuous(limits = y_limit, expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = y_label) +
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
    filename = file.path(out_dir, paste0(file_stem, ".pdf")),
    plot = p,
    device = "pdf",
    width = 92,
    height = 78,
    units = "mm",
    dpi = 300
  )
  ggsave(
    filename = file.path(out_dir, paste0(file_stem, ".png")),
    plot = p,
    device = "png",
    width = 92,
    height = 78,
    units = "mm",
    dpi = 300
  )
}

make_plot(
  "Total T cells",
  "Total T cells in assigned niche (%)",
  "figure5_total_t_cells_by_neighborhood_no_epi_ordered"
)
make_plot(
  "endo_PLVAP",
  "endo_PLVAP cells in assigned niche (%)",
  "figure5_endo_plvap_by_neighborhood_no_epi_ordered"
)

legend_text <- paste0(
  "# Figure 5 Legend Draft - T-cell and endothelial abundance across non-epithelial niches\n\n",
  "Post-chemoradiotherapy cells were assigned to a unique spatial niche by selecting the highest k = 7 neighborhood score. ",
  "For each patient and non-epithelial niche, cell fractions were calculated as the number of indicated cells divided by the total number of cells assigned to that niche. ",
  "Niche_Epi was excluded. Niches are ordered independently in each panel from lower to higher median cell fraction. ",
  "Global P values were calculated using paired Friedman tests across the six non-epithelial niches."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))
