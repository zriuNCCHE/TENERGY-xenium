suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
source(file.path(project_root, "scripts", "plotting", "umap_helpers.R"))
source(file.path(project_root, "scripts", "plotting", "subset_umap_helpers.R"))

out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

immune_order <- c("T/NK", "B/Plasma", "cDC3_LAMP3", "Macro_CXCL5", "Macro_CXCL9", "Mono_CDC27", "Mono_SLC2A3", "Neutrophil")
immune_colors <- c(
  "T/NK" = "#D95F9F",
  "B/Plasma" = "#2C7FB8",
  "cDC3_LAMP3" = "#8C8C8C",
  "Macro_CXCL5" = "#F28E2B",
  "Macro_CXCL9" = "#31A354",
  "Mono_CDC27" = "#E15759",
  "Mono_SLC2A3" = "#8C6BB1",
  "Neutrophil" = "#8D6E63"
)

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "immune_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, final_cell_type, cCR, sample_timepoint, X_umap_1, X_umap_2)
)) %>%
  mutate(final_cell_type = factor(final_cell_type, levels = immune_order))

write_csv(count(umap, final_cell_type, name = "n_cells"), file.path(out_dir, "immune_umap_subtype_counts.csv"))
write_csv(data.frame(cell_type = names(immune_colors), color = unname(immune_colors)), file.path(out_dir, "immune_umap_subtype_palette.csv"))
save_palette_labels(immune_colors, file.path(out_dir, "immune_umap_subtype_color_labels"), plotting_config$font$family, plotting_config$font$standard_pt, width_mm = 75)

plot <- make_subset_umap(
  umap,
  color_col = "final_cell_type",
  colors = immune_colors,
  title = "Immune cell states",
  point_size = 0.18,
  base_family = plotting_config$font$family,
  base_size = plotting_config$font$standard_pt,
  show_labels = FALSE
)

save_umap_raster(plot, file.path(out_dir, "immune_umap_subtypes"), width_mm = 105, height_mm = 90, dpi = 600)
