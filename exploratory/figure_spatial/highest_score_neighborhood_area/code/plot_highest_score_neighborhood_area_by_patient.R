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
neighborhood_labels <- sub("_k7$", "", sub("neighborhood_", "Neighborhood ", neighborhood_cols))
neighborhood_colors <- c(
  "Neighborhood 1" = "#4E79A7",
  "Neighborhood 2" = "#F28E2B",
  "Neighborhood 3" = "#59A14F",
  "Neighborhood 4" = "#E15759",
  "Neighborhood 5" = "#B07AA1",
  "Neighborhood 6" = "#76B7B2",
  "Neighborhood 7" = "#EDC948"
)

cells <- fread(
  data_path,
  select = c(
    "cell_id", "sampleID", "patientID", "cCR", "sample_timepoint",
    "cell_area", neighborhood_cols
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
cells[, highest_score_neighborhood_label := sub(
  "_k7$",
  "",
  sub("neighborhood_", "Neighborhood ", highest_score_neighborhood)
)]

cell_assignments <- cells[, .(
  cell_id,
  sampleID,
  patientID,
  cCR,
  cell_area,
  highest_score_neighborhood,
  highest_score_neighborhood_label
)]

area_summary <- cells[
  !is.na(highest_score_neighborhood),
  .(assigned_area = sum(cell_area), n_cells = .N),
  by = .(patientID, sampleID, cCR, highest_score_neighborhood, highest_score_neighborhood_label)
]
sample_totals <- cells[, .(total_area = sum(cell_area), total_cells = .N), by = .(patientID, sampleID, cCR)]

sample_grid <- sample_totals[, .(patientID, sampleID, cCR)]
sample_grid[, grid_key := 1L]
neighborhood_grid <- data.table(
  highest_score_neighborhood = neighborhood_cols,
  highest_score_neighborhood_label = neighborhood_labels,
  grid_key = 1L
)
plot_dt <- merge(sample_grid, neighborhood_grid, by = "grid_key", allow.cartesian = TRUE)
plot_dt[, grid_key := NULL]
plot_dt <- merge(
  plot_dt,
  area_summary,
  by = c("patientID", "sampleID", "cCR", "highest_score_neighborhood", "highest_score_neighborhood_label"),
  all.x = TRUE
)
plot_dt <- merge(plot_dt, sample_totals, by = c("patientID", "sampleID", "cCR"), all.x = TRUE)
plot_dt[is.na(assigned_area), assigned_area := 0]
plot_dt[is.na(n_cells), n_cells := 0]
plot_dt[, area_fraction := assigned_area / total_area]
plot_dt[, area_percent := area_fraction * 100]

patient_order <- sample_totals[
  order(factor(cCR, levels = c("cCR", "non-cCR")), patientID),
  patientID
]
plot_dt[, patientID := factor(patientID, levels = patient_order)]
plot_dt[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]
plot_dt[, highest_score_neighborhood_label := factor(
  highest_score_neighborhood_label,
  levels = neighborhood_labels
)]

fwrite(cell_assignments, file.path(out_dir, "postc_highest_score_neighborhood_cell_assignments.csv"))
fwrite(plot_dt, file.path(out_dir, "postc_highest_score_neighborhood_area_by_patient.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(plot_dt, aes(x = patientID, y = area_percent, fill = highest_score_neighborhood_label)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.08) +
  facet_grid(. ~ cCR, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = neighborhood_colors, drop = FALSE, name = "Highest-score neighborhood") +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = NULL,
    y = "Assigned cell area (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = small_size),
    legend.key.size = unit(0.28, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.35, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "postc_highest_score_neighborhood_area_by_patient_stacked_bar.pdf"),
  plot = plot,
  device = "pdf",
  width = 180,
  height = 86,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# PostC highest-score spatial neighborhood area by patient",
    "",
    "Exploratory stacked bar plot based on postC spatial neighborhood analysis.",
    "Each cell was assigned to one unique neighborhood by selecting the largest value among neighborhood_1_k7 through neighborhood_7_k7.",
    "No score threshold was applied. Ties are assigned to the first matching neighborhood by column order.",
    "Bars show the proportion of total annotated cell area in each patient assigned to each highest-score neighborhood."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
