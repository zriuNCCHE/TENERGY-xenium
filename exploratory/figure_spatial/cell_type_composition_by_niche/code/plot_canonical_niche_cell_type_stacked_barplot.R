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
  "exploratory", "figure_spatial", "cell_type_composition_by_niche", "outputs",
  "postc_all_patient_niche_cell_type_composition_grouped.csv"
)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

old_to_canonical <- c(
  "Neighborhood 1" = "Niche_Mono",
  "Neighborhood 2" = "Niche_Neutro",
  "Neighborhood 3" = "Niche_Epi",
  "Neighborhood 4" = "Niche_iCAF_CXCL6",
  "Neighborhood 5" = "Niche_T_cell",
  "Neighborhood 6" = "Niche_Macro_CXCL9",
  "Neighborhood 7" = "Niche_iCAF_CXCL5"
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

dt <- fread(source_path)
dt[, original_neighborhood_label := as.character(neighborhood_label)]
dt[, neighborhood_label := fifelse(
  original_neighborhood_label %in% names(old_to_canonical),
  old_to_canonical[original_neighborhood_label],
  original_neighborhood_label
)]
dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]

present_order <- plot_order[plot_order %in% unique(dt$plot_cell_type)]
dt[, plot_cell_type := factor(plot_cell_type, levels = present_order)]

fwrite(
  dt[, .(neighborhood_label, plot_cell_type, n_cells, total_niche_cells, proportion, percent)],
  file.path(out_dir, "postc_canonical_niche_cell_type_composition_grouped.csv")
)

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(dt, aes(x = neighborhood_label, y = percent, fill = plot_cell_type)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.06) +
  scale_fill_manual(
    values = cell_type_colors[present_order],
    breaks = present_order,
    drop = FALSE,
    name = "Cell type"
  ) +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = NULL,
    y = "Cells in niche (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = 5.3),
    legend.key.size = unit(0.24, "cm"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(out_dir, "postc_canonical_niche_cell_type_composition_stacked_barplot.pdf"),
  plot,
  device = "pdf",
  width = 150,
  height = 92,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Canonical spatial niche cell-type composition",
    "",
    "Exploratory stacked barplot showing the proportion of cell types within each post-chemoradiotherapy spatial niche.",
    "Each cell was assigned to a unique spatial niche using the highest-score neighborhood assignment.",
    "The y-axis is the percentage of all cells assigned to each niche.",
    "Small cell types were grouped as Other using the existing grouped composition table."
  ),
  file.path(legend_dir, "legend_canonical_niche_cell_type_composition.md")
)
