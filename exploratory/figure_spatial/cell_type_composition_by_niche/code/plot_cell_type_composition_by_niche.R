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

cell_type_colors <- c(
  "A2ML1+ epi" = "#56B4E9",
  "malignant" = "#CC79A7",
  "Macro_CXCL9" = "#2E7D32",
  "Macro_CXCL5" = "#8DAA22",
  "Mono_CDC27" = "#6A51A3",
  "Mono_SLC2A3" = "#807DBA",
  "Neutrophil" = "#D95F8D",
  "cDC3_LAMP3" = "#1B9E77",
  "B/Plasma" = "#377EB8",
  "T/NK" = "#4C78A8",
  "CD8_Teff" = "#1F78B4",
  "CD8_Tex_PDCD1" = "#6BAED6",
  "CD8_prolif" = "#9ECAE1",
  "CD4_Treg_CCR8" = "#756BB1",
  "CD4_Treg_FOXP3" = "#9E9AC8",
  "CD4_CXCL13" = "#BCBDDC",
  "CD4_prolif" = "#DADAEB",
  "iCAF_CXCL5" = "#E6550D",
  "iCAF_CXCL6" = "#FD8D3C",
  "iCAF_CXCL12" = "#FDD0A2",
  "iCAF_TNC" = "#A63603",
  "myCAF_MMP11" = "#6B4C1F",
  "myCAF_DUX4" = "#8C6D31",
  "myCAF_TNFRSF21" = "#A6761D",
  "myCAF_COL10A1" = "#B5A642",
  "endo_PLVAP" = "#00A6A6",
  "vCAF" = "#5AB4AC",
  "Other" = "#BDBDBD"
)

plot_order <- c(
  "A2ML1+ epi", "malignant",
  "iCAF_CXCL5", "iCAF_CXCL6", "iCAF_CXCL12", "iCAF_TNC",
  "myCAF_MMP11", "myCAF_DUX4", "myCAF_TNFRSF21", "myCAF_COL10A1",
  "vCAF", "endo_PLVAP",
  "Macro_CXCL9", "Macro_CXCL5", "Mono_CDC27", "Mono_SLC2A3",
  "Neutrophil", "cDC3_LAMP3",
  "B/Plasma", "T/NK", "CD8_Teff", "CD8_Tex_PDCD1", "CD8_prolif",
  "CD4_Treg_CCR8", "CD4_Treg_FOXP3", "CD4_CXCL13", "CD4_prolif",
  "Other"
)

cells <- fread(
  data_path,
  select = c(
    "cell_id", "sampleID", "patientID", "cCR", "sample_timepoint",
    "final_cell_type2", neighborhood_cols
  )
)
cells <- cells[sample_timepoint == "postC" & !is.na(final_cell_type2)]

for (neighborhood in neighborhood_cols) {
  cells[, (neighborhood) := as.numeric(get(neighborhood))]
}

score_matrix <- as.matrix(cells[, ..neighborhood_cols])
score_matrix[is.na(score_matrix)] <- -Inf
max_index <- max.col(score_matrix, ties.method = "first")
all_missing <- apply(!is.finite(score_matrix), 1, all)
cells[, highest_score_neighborhood := neighborhood_cols[max_index]]
cells[all_missing, highest_score_neighborhood := NA_character_]
cells[, neighborhood_label := sub(
  "_k7$",
  "",
  sub("neighborhood_", "Neighborhood ", highest_score_neighborhood)
)]
cells <- cells[!is.na(neighborhood_label)]
cells[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_labels)]
cells[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]

full_combined <- cells[, .(n_cells = .N), by = .(neighborhood_label, final_cell_type2)]
full_combined[, total_niche_cells := sum(n_cells), by = neighborhood_label]
full_combined[, proportion := n_cells / total_niche_cells]
full_combined[, percent := proportion * 100]

full_patient <- cells[, .(n_cells = .N), by = .(patientID, sampleID, cCR, neighborhood_label, final_cell_type2)]
full_patient[, total_niche_cells := sum(n_cells), by = .(patientID, sampleID, neighborhood_label)]
full_patient[, proportion := n_cells / total_niche_cells]
full_patient[, percent := proportion * 100]

max_props <- full_combined[, .(max_percent = max(percent), total_cells = sum(n_cells)), by = final_cell_type2]
kept_cell_types <- max_props[max_percent >= 2 | total_cells >= 1500, final_cell_type2]

group_cell_types <- function(dt, denominator_by) {
  grouped <- copy(dt)
  grouped[, plot_cell_type := fifelse(final_cell_type2 %in% kept_cell_types, final_cell_type2, "Other")]
  grouped <- grouped[
    ,
    .(n_cells = sum(n_cells)),
    by = setdiff(names(grouped), c("final_cell_type2", "n_cells", "total_niche_cells", "proportion", "percent"))
  ]
  grouped[, total_niche_cells := sum(n_cells), by = denominator_by]
  grouped[, proportion := n_cells / total_niche_cells]
  grouped[, percent := proportion * 100]
  grouped
}

combined_plot_dt <- group_cell_types(full_combined, denominator_by = c("neighborhood_label"))
patient_plot_dt <- group_cell_types(
  full_patient,
  denominator_by = c("patientID", "sampleID", "cCR", "neighborhood_label")
)

present_order <- plot_order[plot_order %in% unique(c(combined_plot_dt$plot_cell_type, patient_plot_dt$plot_cell_type))]
combined_plot_dt[, plot_cell_type := factor(plot_cell_type, levels = present_order)]
patient_plot_dt[, plot_cell_type := factor(plot_cell_type, levels = present_order)]

patient_order <- unique(cells[order(cCR, patientID), patientID])
patient_plot_dt[, patientID := factor(patientID, levels = patient_order)]

fwrite(cells[, .(cell_id, sampleID, patientID, cCR, final_cell_type2, highest_score_neighborhood, neighborhood_label)],
       file.path(out_dir, "postc_highest_score_niche_cell_type_assignments.csv"))
fwrite(full_combined, file.path(out_dir, "postc_all_patient_niche_cell_type_composition_full.csv"))
fwrite(full_patient, file.path(out_dir, "postc_patient_niche_cell_type_composition_full.csv"))
fwrite(combined_plot_dt, file.path(out_dir, "postc_all_patient_niche_cell_type_composition_grouped.csv"))
fwrite(patient_plot_dt, file.path(out_dir, "postc_patient_niche_cell_type_composition_grouped.csv"))
fwrite(
  data.table(cell_type = names(cell_type_colors), color = unname(cell_type_colors))[cell_type %in% present_order],
  file.path(out_dir, "postc_niche_cell_type_composition_palette.csv")
)

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

common_theme <- theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = 5.4),
    legend.key.size = unit(0.24, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = small_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.18, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

combined_plot <- ggplot(combined_plot_dt, aes(x = neighborhood_label, y = percent, fill = plot_cell_type)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.07) +
  scale_fill_manual(values = cell_type_colors[present_order], breaks = present_order, drop = FALSE, name = "Cell type") +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Cells in niche (%)", title = "All patients combined") +
  common_theme +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5)
  )

patient_plot <- ggplot(patient_plot_dt, aes(x = neighborhood_label, y = percent, fill = plot_cell_type)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.05) +
  facet_wrap(~ patientID + cCR, ncol = 6) +
  scale_fill_manual(values = cell_type_colors[present_order], breaks = present_order, drop = FALSE, name = "Cell type") +
  scale_y_continuous(
    breaks = c(0, 50, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Cells in niche (%)", title = "Individual patients") +
  common_theme +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5)
  )

pdf(
  file = file.path(out_dir, "postc_cell_type_composition_by_highest_score_niche_stacked_bars.pdf"),
  width = 180 / 25.4,
  height = 130 / 25.4,
  family = plotting_config$font$family
)
print(combined_plot)
print(patient_plot)
dev.off()

writeLines(
  c(
    "# PostC cell-type composition by highest-score spatial niche",
    "",
    "Exploratory stacked bar plots based on simple highest-score neighborhood assignment.",
    "Each cell was assigned to one unique niche by selecting the largest value among neighborhood_1_k7 through neighborhood_7_k7.",
    "The first PDF page shows all patients combined. The second PDF page shows patient-level composition.",
    "Values are cell-count proportions: number of cells of a given final_cell_type2 within a niche divided by total cells assigned to that niche.",
    "Cell types with maximum all-patient niche proportion <2% and total cell count <1500 were grouped as Other for readability; full ungrouped tables are also saved."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
