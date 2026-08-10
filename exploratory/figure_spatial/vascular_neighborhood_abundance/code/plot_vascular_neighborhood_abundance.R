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
neighborhood_labels <- c(
  unname(neighborhood_label_map[neighborhood_display_cols])
)
vascular_features <- c("endo_PLVAP", "vCAF")

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.3f", p_value)
}

safe_friedman <- function(dt, value_col) {
  wide <- dcast(
    dt,
    patientID ~ neighborhood_label,
    value.var = value_col,
    fill = NA_real_
  )
  complete_cols <- c("patientID", neighborhood_labels)
  wide <- wide[complete.cases(wide[, ..complete_cols])]
  if (nrow(wide) < 3) {
    return(NA_real_)
  }
  value_matrix <- as.matrix(wide[, ..neighborhood_labels])
  if (all(apply(value_matrix, 2, sd) == 0)) {
    return(NA_real_)
  }
  suppressWarnings(friedman.test(value_matrix)$p.value)
}

safe_pairwise_paired_wilcox <- function(dt, value_col, feature_name) {
  pairs <- combn(neighborhood_labels, 2, simplify = FALSE)
  result_list <- lapply(pairs, function(one_pair) {
    one <- dt[neighborhood_label == one_pair[1], .(patientID, value_1 = get(value_col))]
    two <- dt[neighborhood_label == one_pair[2], .(patientID, value_2 = get(value_col))]
    paired <- merge(one, two, by = "patientID")
    paired <- paired[!is.na(value_1) & !is.na(value_2)]
    if (nrow(paired) < 3 || all(paired$value_1 == paired$value_2)) {
      p_value <- NA_real_
    } else {
      p_value <- suppressWarnings(wilcox.test(
        paired$value_2,
        paired$value_1,
        paired = TRUE,
        exact = FALSE
      )$p.value)
    }
    data.table(
      feature = feature_name,
      comparison = paste(one_pair[2], "vs", one_pair[1]),
      neighborhood_1 = one_pair[1],
      neighborhood_2 = one_pair[2],
      n_patients = nrow(paired),
      median_neighborhood_1 = median(paired$value_1, na.rm = TRUE),
      median_neighborhood_2 = median(paired$value_2, na.rm = TRUE),
      median_delta_2_minus_1 = median(paired$value_2 - paired$value_1, na.rm = TRUE),
      p_value = p_value
    )
  })
  result <- rbindlist(result_list)
  result[, p_adjusted_bh := p.adjust(p_value, method = "BH")]
  result[]
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
neighborhood_grid <- data.table(neighborhood_label = factor(neighborhood_labels, levels = neighborhood_labels), grid_key = 1L)
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

feature_counts <- cells[
  final_cell_type2 %in% vascular_features,
  .(
    feature_cells = .N,
    feature_area = sum(cell_area)
  ),
  by = .(patientID, sampleID, cCR, neighborhood_label, feature = final_cell_type2)
]

feature_grid <- sample_neighborhood_grid
feature_grid[, grid_key := 1L]
vascular_feature_grid <- data.table(feature = vascular_features, grid_key = 1L)
plot_dt <- merge(feature_grid, vascular_feature_grid, by = "grid_key", allow.cartesian = TRUE)
plot_dt[, grid_key := NULL]

plot_dt <- merge(
  plot_dt,
  sample_neighborhood_totals,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all.x = TRUE
)
plot_dt <- merge(
  plot_dt,
  feature_counts,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label", "feature"),
  all.x = TRUE
)

plot_dt[is.na(total_neighborhood_cells), total_neighborhood_cells := 0]
plot_dt[is.na(total_neighborhood_area), total_neighborhood_area := 0]
plot_dt[is.na(feature_cells), feature_cells := 0]
plot_dt[is.na(feature_area), feature_area := 0]
plot_dt[, cell_percent := fifelse(total_neighborhood_cells > 0, feature_cells / total_neighborhood_cells * 100, 0)]
plot_dt[, area_percent := fifelse(total_neighborhood_area > 0, feature_area / total_neighborhood_area * 100, 0)]
plot_dt[, feature := factor(feature, levels = vascular_features)]

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
summary_dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]
summary_dt[, feature := factor(feature, levels = vascular_features)]
setorder(summary_dt, feature, neighborhood_label)

global_stats <- rbindlist(lapply(vascular_features, function(feature_name) {
  feature_dt <- plot_dt[feature == feature_name]
  data.table(
    feature = feature_name,
    metric = c("cell_percent", "area_percent"),
    test = "paired Friedman test across seven neighborhoods",
    p_value = c(
      safe_friedman(feature_dt, "cell_percent"),
      safe_friedman(feature_dt, "area_percent")
    )
  )
}))
global_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = metric]

pairwise_cell <- rbindlist(lapply(vascular_features, function(feature_name) {
  safe_pairwise_paired_wilcox(plot_dt[feature == feature_name], "cell_percent", feature_name)
}))
pairwise_area <- rbindlist(lapply(vascular_features, function(feature_name) {
  safe_pairwise_paired_wilcox(plot_dt[feature == feature_name], "area_percent", feature_name)
}))
pairwise_cell[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature]
pairwise_area[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature]

fwrite(plot_dt, file.path(out_dir, "postc_vascular_neighborhood_abundance_by_patient.csv"))
fwrite(summary_dt, file.path(out_dir, "postc_vascular_neighborhood_abundance_summary.csv"))
fwrite(global_stats, file.path(out_dir, "postc_vascular_neighborhood_global_friedman_stats.csv"))
fwrite(pairwise_cell, file.path(out_dir, "postc_vascular_neighborhood_pairwise_cell_percent_pvalues.csv"))
fwrite(pairwise_area, file.path(out_dir, "postc_vascular_neighborhood_pairwise_area_percent_pvalues.csv"))

fwrite(plot_dt[feature == "endo_PLVAP"], file.path(out_dir, "postc_endo_plvap_neighborhood_abundance_by_patient.csv"))
fwrite(summary_dt[feature == "endo_PLVAP"], file.path(out_dir, "postc_endo_plvap_neighborhood_abundance_summary.csv"))
fwrite(global_stats[feature == "endo_PLVAP"], file.path(out_dir, "postc_endo_plvap_neighborhood_global_friedman_stats.csv"))
fwrite(pairwise_cell[feature == "endo_PLVAP"], file.path(out_dir, "postc_endo_plvap_neighborhood_pairwise_cell_percent_pvalues.csv"))
fwrite(pairwise_area[feature == "endo_PLVAP"], file.path(out_dir, "postc_endo_plvap_neighborhood_pairwise_area_percent_pvalues.csv"))

fwrite(plot_dt[feature == "vCAF"], file.path(out_dir, "postc_vcaf_neighborhood_abundance_by_patient.csv"))
fwrite(summary_dt[feature == "vCAF"], file.path(out_dir, "postc_vcaf_neighborhood_abundance_summary.csv"))
fwrite(global_stats[feature == "vCAF"], file.path(out_dir, "postc_vcaf_neighborhood_global_friedman_stats.csv"))
fwrite(pairwise_cell[feature == "vCAF"], file.path(out_dir, "postc_vcaf_neighborhood_pairwise_cell_percent_pvalues.csv"))
fwrite(pairwise_area[feature == "vCAF"], file.path(out_dir, "postc_vcaf_neighborhood_pairwise_area_percent_pvalues.csv"))

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(plot_dt, aes(x = neighborhood_label, y = cell_percent)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    position = position_jitter(width = 0.09, height = 0),
    size = 0.75,
    alpha = 0.7,
    color = "#00A6A6"
  ) +
  facet_wrap(~ feature, nrow = 1, scales = "free_y") +
  scale_y_continuous(
    limits = c(0, max(plot_dt$cell_percent, na.rm = TRUE) * 1.18),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Cells in assigned niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_text(size = plotting_config$font$small_pt, angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_vascular_neighborhood_cell_fraction_boxplots.pdf"),
  plot = plot,
  device = "pdf",
  width = 160,
  height = 70,
  units = "mm",
  dpi = 300
)

make_single_feature_plot <- function(feature_name) {
  feature_dt <- copy(plot_dt[feature == feature_name])
  feature_order <- feature_dt[
    ,
    .(median_cell_percent = median(cell_percent, na.rm = TRUE)),
    by = .(neighborhood_label)
  ][order(median_cell_percent, as.character(neighborhood_label)), as.character(neighborhood_label)]
  feature_dt[, neighborhood_label := factor(as.character(neighborhood_label), levels = feature_order)]
  ggplot(feature_dt, aes(x = neighborhood_label, y = cell_percent)) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      fill = NA,
      color = "#333333",
      linewidth = line_width
    ) +
    geom_point(
      position = position_jitter(width = 0.09, height = 0),
      size = 0.75,
      alpha = 0.7,
      color = "#00A6A6"
    ) +
    scale_y_continuous(
      limits = c(0, max(feature_dt$cell_percent, na.rm = TRUE) * 1.18),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = paste0(feature_name, " cells in assigned niche (%)")
    ) +
    theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
    theme(
      axis.title.y = element_text(size = base_size, face = "bold"),
      axis.text.y = element_text(size = plotting_config$font$small_pt),
      axis.text.x = element_text(size = plotting_config$font$small_pt, angle = 45, hjust = 1, vjust = 1),
      axis.line = element_line(linewidth = line_width),
      axis.ticks = element_line(linewidth = line_width),
      plot.margin = margin(4, 4, 4, 4)
    )
}

for (feature_name in vascular_features) {
  feature_plot <- make_single_feature_plot(feature_name)
  ggsave(
    filename = file.path(out_dir, paste0("postc_", tolower(feature_name), "_neighborhood_cell_fraction_boxplot.pdf")),
    plot = feature_plot,
    device = "pdf",
    width = 105,
    height = 70,
    units = "mm",
    dpi = 300
  )
}

writeLines(
  c(
    "# PostC endothelial and vascular CAF abundance across highest-score spatial neighborhoods",
    "",
    "Exploratory analysis based on the simple highest-score spatial neighborhood assignment.",
    "Each cell was assigned to one unique neighborhood by selecting the largest value among neighborhood_1_k7 through neighborhood_7_k7.",
    "For each patient and neighborhood, endothelial and vascular CAF abundance was calculated as the number of endo_PLVAP or vCAF cells divided by the total number of cells assigned to that neighborhood.",
    "Global P values use a paired Friedman test across the seven neighborhoods. Pairwise P values use paired Wilcoxon signed-rank tests between neighborhood pairs with BH adjustment."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
