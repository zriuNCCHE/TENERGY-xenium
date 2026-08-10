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

neighborhood_cols <- paste0("neighborhood_", 1:7, "_k7")
neighborhood_display_cols <- c(
  "neighborhood_7_k7",
  "neighborhood_1_k7",
  "neighborhood_2_k7",
  "neighborhood_3_k7",
  "neighborhood_6_k7",
  "neighborhood_4_k7",
  "neighborhood_5_k7"
)
neighborhood_label_map <- c(
  "neighborhood_7_k7" = "Niche_iCAF_CXCL5",
  "neighborhood_1_k7" = "Niche_Mono",
  "neighborhood_2_k7" = "Niche_Neutro",
  "neighborhood_3_k7" = "Niche_Epi",
  "neighborhood_6_k7" = "Niche_Macro_CXCL9",
  "neighborhood_4_k7" = "Niche_iCAF_CXCL6",
  "neighborhood_5_k7" = "Niche_T_cell"
)
neighborhood_labels <- unname(neighborhood_label_map[neighborhood_display_cols])
neighborhood_colors <- c(
  "Niche_iCAF_CXCL5" = "#F0E442",
  "Niche_Mono" = "#C43E2F",
  "Niche_Neutro" = "#D86F1D",
  "Niche_Epi" = "#8C8C8C",
  "Niche_Macro_CXCL9" = "#2FBF71",
  "Niche_iCAF_CXCL6" = "#56B4E9",
  "Niche_T_cell" = "#FF8CCB"
)

target_feature <- "Macro_CXCL5"
partner_features <- c("Mono_CDC27", "Mono_SLC2A3", "Neutrophil")
partner_labels <- c(
  "Mono_CDC27" = "Mono_CDC27",
  "Mono_SLC2A3" = "Mono_SLC2A3",
  "Neutrophil" = "Neutrophil",
  "Combined myeloid" = "Mono_CDC27 + Mono_SLC2A3 + Neutrophil"
)
partner_colors <- c(
  "Mono_CDC27" = "#0072B2",
  "Mono_SLC2A3" = "#009E73",
  "Neutrophil" = "#D86F1D",
  "Combined myeloid" = "#C43E2F"
)

format_p <- function(p_value) {
  vapply(p_value, function(one_p) {
    if (is.na(one_p)) {
      return("NA")
    }
    if (one_p < 1e-4) {
      return(formatC(one_p, format = "e", digits = 2))
    }
    if (one_p < 0.001) {
      return(sprintf("%.6f", one_p))
    }
    sprintf("%.3f", one_p)
  }, character(1))
}

safe_cor <- function(dt, x_col, y_col, method = "spearman") {
  dt <- dt[!is.na(get(x_col)) & !is.na(get(y_col))]
  if (nrow(dt) < 3 || sd(dt[[x_col]]) == 0 || sd(dt[[y_col]]) == 0) {
    return(list(estimate = NA_real_, p_value = NA_real_, n = nrow(dt)))
  }
  result <- suppressWarnings(cor.test(dt[[x_col]], dt[[y_col]], method = method, exact = FALSE))
  list(estimate = unname(result$estimate), p_value = result$p.value, n = nrow(dt))
}

cells <- fread(
  data_path,
  select = c(
    "cell_id", "sampleID", "patientID", "cCR", "sample_timepoint",
    "cell_area", "final_cell_type2", neighborhood_cols
  )
)
cells <- cells[
  sample_timepoint == "postC" &
    !is.na(final_cell_type2) &
    !is.na(cell_area) &
    cell_area > 0
]

for (neighborhood in neighborhood_cols) {
  cells[, (neighborhood) := as.numeric(get(neighborhood))]
}

score_matrix <- as.matrix(cells[, ..neighborhood_cols])
score_matrix[is.na(score_matrix)] <- -Inf
max_index <- max.col(score_matrix, ties.method = "first")
all_missing <- apply(!is.finite(score_matrix), 1, all)
cells[, highest_score_neighborhood := neighborhood_cols[max_index]]
cells[all_missing, highest_score_neighborhood := NA_character_]
cells <- cells[!is.na(highest_score_neighborhood)]
cells[, neighborhood_label := factor(
  neighborhood_label_map[highest_score_neighborhood],
  levels = neighborhood_labels
)]

sample_grid <- unique(cells[, .(patientID, sampleID, cCR)])
sample_grid[, grid_key := 1L]
neighborhood_grid <- data.table(
  neighborhood_label = factor(neighborhood_labels, levels = neighborhood_labels),
  grid_key = 1L
)
sample_neighborhood_grid <- merge(sample_grid, neighborhood_grid, by = "grid_key", allow.cartesian = TRUE)
sample_neighborhood_grid[, grid_key := NULL]

sample_neighborhood_totals <- cells[
  ,
  .(
    total_neighborhood_cells = .N,
    total_neighborhood_area = sum(cell_area)
  ),
  by = .(patientID, sampleID, cCR, neighborhood_label)
]

features_for_count <- c(target_feature, partner_features)
feature_counts <- cells[
  final_cell_type2 %in% features_for_count,
  .(
    feature_cells = .N,
    feature_area = sum(cell_area)
  ),
  by = .(patientID, sampleID, cCR, neighborhood_label, feature = final_cell_type2)
]

feature_grid <- sample_neighborhood_grid
feature_grid[, grid_key := 1L]
feature_name_grid <- data.table(feature = features_for_count, grid_key = 1L)
abundance_dt <- merge(feature_grid, feature_name_grid, by = "grid_key", allow.cartesian = TRUE)
abundance_dt[, grid_key := NULL]

abundance_dt <- merge(
  abundance_dt,
  sample_neighborhood_totals,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all.x = TRUE
)
abundance_dt <- merge(
  abundance_dt,
  feature_counts,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label", "feature"),
  all.x = TRUE
)

abundance_dt[is.na(total_neighborhood_cells), total_neighborhood_cells := 0]
abundance_dt[is.na(total_neighborhood_area), total_neighborhood_area := 0]
abundance_dt[is.na(feature_cells), feature_cells := 0]
abundance_dt[is.na(feature_area), feature_area := 0]
abundance_dt[, cell_percent := fifelse(total_neighborhood_cells > 0, feature_cells / total_neighborhood_cells * 100, 0)]
abundance_dt[, area_percent := fifelse(total_neighborhood_area > 0, feature_area / total_neighborhood_area * 100, 0)]

wide_dt <- dcast(
  abundance_dt,
  patientID + sampleID + cCR + neighborhood_label ~ feature,
  value.var = c("cell_percent", "area_percent", "feature_cells", "feature_area", "total_neighborhood_cells", "total_neighborhood_area")
)
wide_dt[, combined_myeloid_cell_percent := cell_percent_Mono_CDC27 + cell_percent_Mono_SLC2A3 + cell_percent_Neutrophil]
wide_dt[, combined_myeloid_area_percent := area_percent_Mono_CDC27 + area_percent_Mono_SLC2A3 + area_percent_Neutrophil]
wide_dt[, combined_myeloid_cells := feature_cells_Mono_CDC27 + feature_cells_Mono_SLC2A3 + feature_cells_Neutrophil]
wide_dt[, neighborhood_label := factor(as.character(neighborhood_label), levels = neighborhood_labels)]
wide_dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

cor_specs <- data.table(
  partner = c(partner_features, "Combined myeloid"),
  x_col = "cell_percent_Macro_CXCL5",
  y_col = c(
    "cell_percent_Mono_CDC27",
    "cell_percent_Mono_SLC2A3",
    "cell_percent_Neutrophil",
    "combined_myeloid_cell_percent"
  ),
  x_area_col = "area_percent_Macro_CXCL5",
  y_area_col = c(
    "area_percent_Mono_CDC27",
    "area_percent_Mono_SLC2A3",
    "area_percent_Neutrophil",
    "combined_myeloid_area_percent"
  )
)

overall_stats <- rbindlist(lapply(seq_len(nrow(cor_specs)), function(i) {
  spec <- cor_specs[i]
  sp <- safe_cor(wide_dt, spec$x_col, spec$y_col, "spearman")
  pe <- safe_cor(wide_dt, spec$x_col, spec$y_col, "pearson")
  area_sp <- safe_cor(wide_dt, spec$x_area_col, spec$y_area_col, "spearman")
  data.table(
    comparison = paste0("Macro_CXCL5_vs_", spec$partner),
    partner = spec$partner,
    scope = "all_patient_niche_pairs",
    metric = c("cell_percent", "cell_percent", "area_percent"),
    method = c("Spearman", "Pearson", "Spearman"),
    n = c(sp$n, pe$n, area_sp$n),
    estimate = c(sp$estimate, pe$estimate, area_sp$estimate),
    p_value = c(sp$p_value, pe$p_value, area_sp$p_value)
  )
}))

niche_stats <- rbindlist(lapply(seq_len(nrow(cor_specs)), function(i) {
  spec <- cor_specs[i]
  rbindlist(lapply(neighborhood_labels, function(one_niche) {
    sub_dt <- wide_dt[neighborhood_label == one_niche]
    sp <- safe_cor(sub_dt, spec$x_col, spec$y_col, "spearman")
    pe <- safe_cor(sub_dt, spec$x_col, spec$y_col, "pearson")
    data.table(
      comparison = paste0("Macro_CXCL5_vs_", spec$partner),
      partner = spec$partner,
      scope = one_niche,
      metric = "cell_percent",
      method = c("Spearman", "Pearson"),
      n = c(sp$n, pe$n),
      estimate = c(sp$estimate, pe$estimate),
      p_value = c(sp$p_value, pe$p_value)
    )
  }))
}))
stats_dt <- rbindlist(list(overall_stats, niche_stats), use.names = TRUE)
stats_dt[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = .(scope, method, metric)]

median_dt <- wide_dt[
  ,
  .(
    macro_cxcl5_median_cell_percent = median(cell_percent_Macro_CXCL5, na.rm = TRUE),
    mono_cdc27_median_cell_percent = median(cell_percent_Mono_CDC27, na.rm = TRUE),
    mono_slc2a3_median_cell_percent = median(cell_percent_Mono_SLC2A3, na.rm = TRUE),
    neutrophil_median_cell_percent = median(cell_percent_Neutrophil, na.rm = TRUE),
    combined_myeloid_median_cell_percent = median(combined_myeloid_cell_percent, na.rm = TRUE)
  ),
  by = .(neighborhood_label)
]
median_long <- melt(
  median_dt,
  id.vars = c("neighborhood_label", "macro_cxcl5_median_cell_percent"),
  measure.vars = c(
    "mono_cdc27_median_cell_percent",
    "mono_slc2a3_median_cell_percent",
    "neutrophil_median_cell_percent",
    "combined_myeloid_median_cell_percent"
  ),
  variable.name = "partner_metric",
  value.name = "partner_median_cell_percent"
)
median_long[, partner := fifelse(
  partner_metric == "mono_cdc27_median_cell_percent", "Mono_CDC27",
  fifelse(
    partner_metric == "mono_slc2a3_median_cell_percent", "Mono_SLC2A3",
    fifelse(partner_metric == "neutrophil_median_cell_percent", "Neutrophil", "Combined myeloid")
  )
)]
median_stats <- rbindlist(lapply(unique(median_long$partner), function(one_partner) {
  sub_dt <- median_long[partner == one_partner]
  sp <- safe_cor(sub_dt, "macro_cxcl5_median_cell_percent", "partner_median_cell_percent", "spearman")
  data.table(
    analysis = "niche_median_correlation_visual_summary",
    comparison = paste0("Macro_CXCL5_vs_", one_partner),
    test = "Spearman correlation across seven niche medians",
    alternative = "two-sided",
    n_niches = sp$n,
    rho = sp$estimate,
    p_value = sp$p_value
  )
}))

plot_long <- rbindlist(lapply(seq_len(nrow(cor_specs)), function(i) {
  spec <- cor_specs[i]
  wide_dt[
    ,
    .(
      patientID,
      sampleID,
      cCR,
      neighborhood_label,
      partner = spec$partner,
      macro_cxcl5_cell_percent = get(spec$x_col),
      partner_cell_percent = get(spec$y_col),
      macro_cxcl5_area_percent = get(spec$x_area_col),
      partner_area_percent = get(spec$y_area_col)
    )
  ]
}))
plot_long[, partner := factor(partner, levels = c(partner_features, "Combined myeloid"))]

fwrite(abundance_dt, file.path(out_dir, "postc_macro_cxcl5_myeloid_abundance_by_patient_niche_long.csv"))
fwrite(wide_dt, file.path(out_dir, "postc_macro_cxcl5_myeloid_abundance_by_patient_niche_wide.csv"))
fwrite(plot_long, file.path(out_dir, "postc_macro_cxcl5_myeloid_correlation_plot_values.csv"))
fwrite(stats_dt, file.path(out_dir, "postc_macro_cxcl5_myeloid_correlation_stats.csv"))
fwrite(median_long, file.path(out_dir, "postc_macro_cxcl5_myeloid_niche_median_values.csv"))
fwrite(median_stats, file.path(out_dir, "postc_macro_cxcl5_myeloid_niche_median_correlation_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

selected_stats <- overall_stats[method == "Spearman" & metric == "cell_percent"]
stat_labels <- selected_stats[
  ,
  .(
    partner,
    label = paste0("rho = ", sprintf("%.2f", estimate), "\np = ", format_p(p_value))
  )
]

scatter_dt <- merge(plot_long, stat_labels, by = "partner", all.x = TRUE)
scatter_plot <- ggplot(scatter_dt, aes(x = macro_cxcl5_cell_percent, y = partner_cell_percent)) +
  geom_point(aes(color = neighborhood_label), size = 1.05, alpha = 0.82) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#333333",
    fill = "#D9D9D9",
    linewidth = line_width
  ) +
  geom_text(
    data = unique(scatter_dt[, .(partner, label)]),
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    size = 1.9
  ) +
  facet_wrap(~ partner, ncol = 2, scales = "free_y", labeller = as_labeller(partner_labels)) +
  scale_color_manual(values = neighborhood_colors, drop = FALSE, name = "Niche") +
  labs(
    x = "Macro_CXCL5 cells in assigned niche (%)",
    y = "Partner cells in assigned niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = small_size),
    legend.key.height = unit(0.28, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.34, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_macro_cxcl5_myeloid_correlation_scatter.pdf"),
  plot = scatter_plot,
  device = "pdf",
  width = 150,
  height = 118,
  units = "mm",
  dpi = 300
)

median_stats_labels <- median_stats[
  ,
  .(
    partner = sub("^Macro_CXCL5_vs_", "", comparison),
    label = paste0("rho = ", sprintf("%.2f", rho), "\np = ", format_p(p_value))
  )
]
median_plot_dt <- merge(median_long, median_stats_labels, by = "partner", all.x = TRUE)
median_plot_dt[, partner := factor(partner, levels = c(partner_features, "Combined myeloid"))]

median_plot <- ggplot(
  median_plot_dt,
  aes(x = macro_cxcl5_median_cell_percent, y = partner_median_cell_percent)
) +
  geom_point(aes(color = neighborhood_label), size = 1.8, alpha = 0.95) +
  geom_text(
    aes(label = neighborhood_label),
    size = 1.8,
    hjust = -0.05,
    vjust = 0.4,
    check_overlap = TRUE
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_text(
    data = unique(median_plot_dt[, .(partner, label)]),
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.15,
    size = 1.9
  ) +
  facet_wrap(~ partner, ncol = 2, scales = "free_y", labeller = as_labeller(partner_labels)) +
  scale_color_manual(values = neighborhood_colors, drop = FALSE, name = "Niche") +
  labs(
    x = "Median Macro_CXCL5 cells in niche (%)",
    y = "Median partner cells in niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = small_size),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.34, "cm"),
    plot.margin = margin(4, 10, 4, 4)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(out_dir, "postc_macro_cxcl5_myeloid_niche_median_concordance.pdf"),
  plot = median_plot,
  device = "pdf",
  width = 150,
  height = 118,
  units = "mm",
  dpi = 300
)

paired_compare_dt <- plot_long[
  partner %in% c("Combined myeloid"),
  .(
    feature = c(target_feature, "Combined myeloid"),
    cell_percent = c(macro_cxcl5_cell_percent, partner_cell_percent)
  ),
  by = .(patientID, sampleID, cCR, neighborhood_label)
]
paired_compare_dt[, feature := factor(feature, levels = c(target_feature, "Combined myeloid"))]
paired_compare_dt[, neighborhood_label := factor(as.character(neighborhood_label), levels = neighborhood_labels)]

pair_plot <- ggplot(paired_compare_dt, aes(x = neighborhood_label, y = cell_percent)) +
  geom_boxplot(
    aes(color = feature),
    fill = NA,
    width = 0.58,
    outlier.shape = NA,
    linewidth = line_width,
    position = position_dodge(width = 0.72)
  ) +
  geom_point(
    aes(color = feature),
    size = 0.75,
    alpha = 0.72,
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.72)
  ) +
  scale_color_manual(values = c("Macro_CXCL5" = "#7A7A7A", "Combined myeloid" = "#C43E2F")) +
  labs(
    x = NULL,
    y = "Cells in assigned niche (%)",
    color = NULL
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    legend.position = "top",
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_macro_cxcl5_combined_myeloid_by_niche_boxplot.pdf"),
  plot = pair_plot,
  device = "pdf",
  width = 112,
  height = 82,
  units = "mm",
  dpi = 300
)

legend_text <- paste0(
  "# Exploratory Legend - Macro_CXCL5 and Myeloid Niches\n\n",
  "Post-chemoradiotherapy cells were assigned to a unique spatial niche by selecting the highest score among the seven k = 7 neighborhood scores. ",
  "For each patient and niche, cell fractions were calculated as the number of cells in the indicated subtype divided by the total number of cells assigned to that niche. ",
  "Correlation plots compare Macro_CXCL5 fraction with Mono_CDC27, Mono_SLC2A3, Neutrophil, or the combined myeloid fraction across patient-niche observations. ",
  "P values are exploratory and were calculated using two-sided Spearman correlation unless otherwise indicated."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))
