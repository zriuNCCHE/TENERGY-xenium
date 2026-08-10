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
t_subset_order <- c(
  "T/NK",
  "CD8_Teff",
  "CD8_Tex_PDCD1",
  "CD8_prolif",
  "CD4_Treg_CCR8",
  "CD4_Treg_FOXP3",
  "CD4_CXCL13",
  "CD4_prolif"
)
treg_like <- c("CD4_Treg_CCR8", "CD4_Treg_FOXP3")
non_treg_t <- setdiff(t_subset_order, treg_like)
combined_order <- c("non-Treg T/NK", "Treg-like")
total_t_order <- "Total T cells"

format_p_value <- function(p_value) {
  vapply(p_value, function(one_p) {
    if (is.na(one_p)) {
      return("NA")
    }
    if (one_p < 0.001) {
      return(sub("\\.?0+$", "", sprintf("%.6f", one_p)))
    }
    sprintf("%.3f", one_p)
  }, character(1))
}

format_plot_p_value <- function(p_value, prefix = "p") {
  vapply(p_value, function(one_p) {
    if (is.na(one_p)) {
      return(paste0(prefix, " = NA"))
    }
    if (one_p < 0.001) {
      return(sub("\\.?0+$", "", sprintf("%s = %.6f", prefix, one_p)))
    }
    sprintf("%s = %.3f", prefix, one_p)
  }, character(1))
}

safe_friedman <- function(dt) {
  wide <- dcast(
    dt,
    patientID ~ neighborhood_label,
    value.var = "cell_percent",
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

safe_pairwise_paired_wilcox <- function(dt, feature_name) {
  pairs <- combn(neighborhood_labels, 2, simplify = FALSE)
  result_list <- lapply(pairs, function(one_pair) {
    one <- dt[neighborhood_label == one_pair[1], .(patientID, value_1 = cell_percent)]
    two <- dt[neighborhood_label == one_pair[2], .(patientID, value_2 = cell_percent)]
    paired <- merge(one, two, by = "patientID")
    paired <- paired[!is.na(value_1) & !is.na(value_2)]
    if (nrow(paired) < 3 || all(paired$value_1 == paired$value_2)) {
      p_value <- NA_real_
      median_delta <- median(paired$value_2 - paired$value_1, na.rm = TRUE)
    } else {
      p_value <- suppressWarnings(wilcox.test(
        paired$value_2,
        paired$value_1,
        paired = TRUE,
        exact = FALSE
      )$p.value)
      median_delta <- median(paired$value_2 - paired$value_1, na.rm = TRUE)
    }
    data.table(
      feature = feature_name,
      neighborhood_1 = one_pair[1],
      neighborhood_2 = one_pair[2],
      n_patients = nrow(paired),
      median_neighborhood_1 = median(paired$value_1, na.rm = TRUE),
      median_neighborhood_2 = median(paired$value_2, na.rm = TRUE),
      median_delta_2_minus_1 = median_delta,
      p_value = p_value
    )
  })
  rbindlist(result_list)
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
cells[, neighborhood_label := neighborhood_label_map[highest_score_neighborhood]]
cells <- cells[!is.na(neighborhood_label)]

sample_neighborhood_totals <- cells[
  ,
  .(total_neighborhood_area = sum(cell_area), total_neighborhood_cells = .N),
  by = .(patientID, sampleID, cCR, neighborhood_label)
]

sample_grid <- unique(cells[, .(patientID, sampleID, cCR)])
sample_grid[, grid_key := 1L]
neighborhood_grid <- data.table(neighborhood_label = neighborhood_labels, grid_key = 1L)
sample_neighborhood_grid <- merge(sample_grid, neighborhood_grid, by = "grid_key", allow.cartesian = TRUE)
sample_neighborhood_grid[, grid_key := NULL]
sample_neighborhood_grid <- merge(
  sample_neighborhood_grid,
  sample_neighborhood_totals,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all.x = TRUE
)
sample_neighborhood_grid[is.na(total_neighborhood_area), total_neighborhood_area := 0]
sample_neighborhood_grid[is.na(total_neighborhood_cells), total_neighborhood_cells := 0]

t_subset_counts <- cells[
  final_cell_type2 %in% t_subset_order,
  .(subset_area = sum(cell_area), subset_cells = .N),
  by = .(patientID, sampleID, cCR, neighborhood_label, feature = final_cell_type2)
]
subset_grid <- sample_neighborhood_grid[, .(patientID, sampleID, cCR, neighborhood_label)]
subset_grid[, grid_key := 1L]
feature_grid <- data.table(feature = t_subset_order, grid_key = 1L)
subset_plot_dt <- merge(subset_grid, feature_grid, by = "grid_key", allow.cartesian = TRUE)
subset_plot_dt[, grid_key := NULL]
subset_plot_dt <- merge(
  subset_plot_dt,
  t_subset_counts,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label", "feature"),
  all.x = TRUE
)
subset_plot_dt <- merge(
  subset_plot_dt,
  sample_neighborhood_grid[, .(patientID, sampleID, cCR, neighborhood_label, total_neighborhood_area, total_neighborhood_cells)],
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all.x = TRUE
)
subset_plot_dt[is.na(subset_area), subset_area := 0]
subset_plot_dt[is.na(subset_cells), subset_cells := 0]
subset_plot_dt[, area_percent := fifelse(total_neighborhood_area > 0, subset_area / total_neighborhood_area * 100, 0)]
subset_plot_dt[, cell_percent := fifelse(total_neighborhood_cells > 0, subset_cells / total_neighborhood_cells * 100, 0)]
subset_plot_dt[, feature_type := "Individual T subset"]

cells[, combined_t_group := fifelse(
  final_cell_type2 %in% treg_like,
  "Treg-like",
  fifelse(final_cell_type2 %in% non_treg_t, "non-Treg T/NK", NA_character_)
)]
combined_counts <- cells[
  !is.na(combined_t_group),
  .(subset_area = sum(cell_area), subset_cells = .N),
  by = .(patientID, sampleID, cCR, neighborhood_label, feature = combined_t_group)
]
combined_grid <- sample_neighborhood_grid[, .(patientID, sampleID, cCR, neighborhood_label)]
combined_grid[, grid_key := 1L]
combined_feature_grid <- data.table(feature = combined_order, grid_key = 1L)
combined_plot_dt <- merge(combined_grid, combined_feature_grid, by = "grid_key", allow.cartesian = TRUE)
combined_plot_dt[, grid_key := NULL]
combined_plot_dt <- merge(
  combined_plot_dt,
  combined_counts,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label", "feature"),
  all.x = TRUE
)
combined_plot_dt <- merge(
  combined_plot_dt,
  sample_neighborhood_grid[, .(patientID, sampleID, cCR, neighborhood_label, total_neighborhood_area, total_neighborhood_cells)],
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all.x = TRUE
)
combined_plot_dt[is.na(subset_area), subset_area := 0]
combined_plot_dt[is.na(subset_cells), subset_cells := 0]
combined_plot_dt[, area_percent := fifelse(total_neighborhood_area > 0, subset_area / total_neighborhood_area * 100, 0)]
combined_plot_dt[, cell_percent := fifelse(total_neighborhood_cells > 0, subset_cells / total_neighborhood_cells * 100, 0)]
combined_plot_dt[, feature_type := "Combined T group"]

total_t_counts <- cells[
  final_cell_type2 %in% t_subset_order,
  .(subset_area = sum(cell_area), subset_cells = .N),
  by = .(patientID, sampleID, cCR, neighborhood_label)
]
total_t_plot_dt <- copy(sample_neighborhood_grid)
total_t_plot_dt <- merge(
  total_t_plot_dt,
  total_t_counts,
  by = c("patientID", "sampleID", "cCR", "neighborhood_label"),
  all.x = TRUE
)
total_t_plot_dt[is.na(subset_area), subset_area := 0]
total_t_plot_dt[is.na(subset_cells), subset_cells := 0]
total_t_plot_dt[, area_percent := fifelse(total_neighborhood_area > 0, subset_area / total_neighborhood_area * 100, 0)]
total_t_plot_dt[, cell_percent := fifelse(total_neighborhood_cells > 0, subset_cells / total_neighborhood_cells * 100, 0)]
total_t_plot_dt[, feature_type := "Total T cells"]
total_t_plot_dt[, feature := total_t_order]

plot_dt <- rbindlist(list(subset_plot_dt, combined_plot_dt, total_t_plot_dt), use.names = TRUE)
plot_dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]
plot_dt[, feature := factor(feature, levels = c(t_subset_order, combined_order, total_t_order))]
subset_plot_dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]
subset_plot_dt[, feature := factor(feature, levels = t_subset_order)]
combined_plot_dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]
combined_plot_dt[, feature := factor(feature, levels = combined_order)]
total_t_plot_dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]
total_t_plot_dt[, feature := factor(feature, levels = total_t_order)]

median_summary <- plot_dt[
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
  by = .(feature_type, feature, neighborhood_label)
]
median_summary[, feature := factor(feature, levels = c(t_subset_order, combined_order, total_t_order))]
median_summary[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]

global_stats <- plot_dt[
  ,
  .(friedman_p_value = safe_friedman(.SD)),
  by = .(feature_type, feature)
]
global_stats[, friedman_p_adjusted_bh := p.adjust(friedman_p_value, method = "BH"), by = feature_type]

pairwise_stats <- rbindlist(lapply(split(plot_dt, plot_dt$feature), function(feature_dt) {
  safe_pairwise_paired_wilcox(feature_dt, as.character(unique(feature_dt$feature)))
}), use.names = TRUE)
feature_type_map <- unique(plot_dt[, .(feature = as.character(feature), feature_type)])
pairwise_stats <- merge(pairwise_stats, feature_type_map, by = "feature", all.x = TRUE)
pairwise_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature]
pairwise_stats[, p_label := format_p_value(p_value)]
pairwise_stats[, p_adjusted_label := format_p_value(p_adjusted_bh)]

fwrite(plot_dt, file.path(out_dir, "postc_t_cell_neighborhood_abundance_by_patient.csv"))
fwrite(median_summary, file.path(out_dir, "postc_t_cell_neighborhood_abundance_median_summary.csv"))
fwrite(global_stats, file.path(out_dir, "postc_t_cell_neighborhood_abundance_global_friedman_stats.csv"))
fwrite(pairwise_stats, file.path(out_dir, "postc_t_cell_neighborhood_abundance_pairwise_pvalues.csv"))
fwrite(total_t_plot_dt, file.path(out_dir, "postc_total_t_cells_neighborhood_abundance_by_patient.csv"))
fwrite(median_summary[feature == total_t_order], file.path(out_dir, "postc_total_t_cells_neighborhood_abundance_summary.csv"))
fwrite(global_stats[feature == total_t_order], file.path(out_dir, "postc_total_t_cells_neighborhood_global_friedman_stats.csv"))
fwrite(pairwise_stats[feature == total_t_order], file.path(out_dir, "postc_total_t_cells_neighborhood_pairwise_pvalues.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

make_p_label_dt <- function(plot_data, stats_data, feature_type_name) {
  label_dt <- plot_data[
    ,
    .(panel_max = max(cell_percent, na.rm = TRUE)),
    by = .(feature)
  ]
  label_dt <- merge(
    label_dt,
    stats_data[feature_type == feature_type_name],
    by = "feature",
    all.x = TRUE
  )
  label_dt[, y_position := fifelse(panel_max > 0, panel_max * 1.12, 0.1)]
  label_dt[, label := paste0(
    "Friedman ",
    format_plot_p_value(friedman_p_value, "p"),
    "\nadj ",
    format_plot_p_value(friedman_p_adjusted_bh, "p")
  )]
  label_dt[]
}

subset_p_labels <- make_p_label_dt(subset_plot_dt, global_stats, "Individual T subset")
subset_p_labels[, feature := factor(feature, levels = t_subset_order)]
subset_order_dt <- subset_plot_dt[
  ,
  .(median_cell_percent = median(cell_percent, na.rm = TRUE)),
  by = .(feature, neighborhood_label)
][order(feature, median_cell_percent, as.character(neighborhood_label))]
subset_order_dt[
  ,
  feature_neighborhood := paste(as.character(neighborhood_label), as.character(feature), sep = "___")
]
subset_x_levels <- subset_order_dt$feature_neighborhood
subset_label_map <- setNames(
  as.character(subset_order_dt$neighborhood_label),
  subset_order_dt$feature_neighborhood
)
subset_plot_ordered_dt <- copy(subset_plot_dt)
subset_plot_ordered_dt[
  ,
  feature_neighborhood := factor(
    paste(as.character(neighborhood_label), as.character(feature), sep = "___"),
    levels = subset_x_levels
  )
]
subset_p_labels[
  ,
  x_position := 4
]
combined_p_labels <- make_p_label_dt(combined_plot_dt, global_stats, "Combined T group")
combined_p_labels[, feature := factor(feature, levels = combined_order)]
total_t_p_labels <- make_p_label_dt(total_t_plot_dt, global_stats, "Total T cells")
total_t_p_labels[, feature := factor(feature, levels = total_t_order)]
total_t_abundance_order <- total_t_plot_dt[
  ,
  .(median_cell_percent = median(cell_percent, na.rm = TRUE)),
  by = .(neighborhood_label)
][order(median_cell_percent, as.character(neighborhood_label)), as.character(neighborhood_label)]
total_t_plot_ordered <- copy(total_t_plot_dt)
total_t_plot_ordered[, neighborhood_label := factor(as.character(neighborhood_label), levels = total_t_abundance_order)]
total_t_p_labels_ordered <- copy(total_t_p_labels)
total_t_p_labels_ordered[, x_position := ceiling(length(total_t_abundance_order) / 2)]

subset_boxplot <- ggplot(subset_plot_ordered_dt, aes(x = feature_neighborhood, y = cell_percent)) +
  geom_boxplot(
    fill = "#C6DBEF",
    color = "#333333",
    width = 0.56,
    outlier.shape = NA,
    linewidth = line_width,
    alpha = 0.62
  ) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 0.48,
    alpha = 0.62,
    color = "#333333"
  ) +
  geom_text(
    data = subset_p_labels,
    aes(x = x_position, y = y_position, label = label),
    inherit.aes = FALSE,
    size = 1.45,
    lineheight = 0.88
  ) +
  facet_wrap(~ feature, ncol = 4, scales = "free") +
  scale_x_discrete(labels = subset_label_map) +
  labs(x = NULL, y = "Cells in assigned niche (%)") +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.22, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )
ggsave(
  filename = file.path(out_dir, "postc_individual_t_subset_neighborhood_boxplots.pdf"),
  plot = subset_boxplot,
  device = "pdf",
  width = 165,
  height = 108,
  units = "mm",
  dpi = 300
)

combined_plot <- ggplot(combined_plot_dt, aes(x = neighborhood_label, y = cell_percent)) +
  geom_boxplot(
    fill = "#BFD3E6",
    color = "#333333",
    width = 0.56,
    outlier.shape = NA,
    linewidth = line_width,
    alpha = 0.62
  ) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 0.62,
    alpha = 0.68,
    color = "#333333"
  ) +
  geom_text(
    data = combined_p_labels,
    aes(x = 4, y = y_position, label = label),
    inherit.aes = FALSE,
    size = 1.65,
    lineheight = 0.88
  ) +
  facet_wrap(~ feature, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "Cells in assigned niche (%)") +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.28, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )
ggsave(
  filename = file.path(out_dir, "postc_combined_t_groups_neighborhood_cell_fraction_boxplots.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 132,
  height = 68,
  units = "mm",
  dpi = 300
)

total_t_plot <- ggplot(total_t_plot_ordered, aes(x = neighborhood_label, y = cell_percent)) +
  geom_boxplot(
    fill = "#BFD3E6",
    color = "#333333",
    width = 0.56,
    outlier.shape = NA,
    linewidth = line_width,
    alpha = 0.62
  ) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 0.62,
    alpha = 0.68,
    color = "#333333"
  ) +
  geom_text(
    data = total_t_p_labels_ordered,
    aes(x = x_position, y = y_position, label = label),
    inherit.aes = FALSE,
    size = 1.65,
    lineheight = 0.88
  ) +
  labs(x = NULL, y = "Total T cells in assigned niche (%)") +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )
ggsave(
  filename = file.path(out_dir, "postc_total_t_cells_neighborhood_cell_fraction_boxplot.pdf"),
  plot = total_t_plot,
  device = "pdf",
  width = 105,
  height = 68,
  units = "mm",
  dpi = 300
)

top_pairwise <- pairwise_stats[
  feature %in% combined_order & !is.na(p_value)
][order(p_value)]
top_pairwise <- top_pairwise[seq_len(min(.N, 12))]
fwrite(top_pairwise, file.path(out_dir, "postc_combined_t_groups_top_pairwise_pvalues.csv"))

writeLines(
  c(
    "# PostC T-cell abundance across highest-score spatial neighborhoods",
    "",
    "Exploratory analysis based on the simple highest-score spatial neighborhood assignment.",
    "Each cell was assigned to one unique neighborhood by selecting the largest value among neighborhood_1_k7 through neighborhood_7_k7.",
    "For each patient and neighborhood, T-cell subset abundance was calculated as the number of cells in that subset divided by the total number of cells assigned to that neighborhood.",
    "",
    "Individual T subsets: T/NK, CD8_Teff, CD8_Tex_PDCD1, CD8_prolif, CD4_Treg_CCR8, CD4_Treg_FOXP3, CD4_CXCL13 and CD4_prolif.",
    "Combined T groups: Treg-like = CD4_Treg_CCR8 + CD4_Treg_FOXP3; non-Treg T/NK = all other listed T/NK, CD8 and CD4 T subsets.",
    "",
    "Global P values use a paired Friedman test across the seven neighborhoods.",
    "Pairwise P values use paired Wilcoxon signed-rank tests between neighborhood pairs, with BH adjustment within each feature.",
    "This is intentionally circular/exploratory because neighborhoods were derived from cell-type neighborhood information."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
