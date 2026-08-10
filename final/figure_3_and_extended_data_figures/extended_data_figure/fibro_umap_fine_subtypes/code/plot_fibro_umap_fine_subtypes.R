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

fibro_order <- c("myCAF_MMP11", "myCAF_TNFRSF21", "myCAF_DUX4", "myCAF_COL10A1", "iCAF_CXCL5", "iCAF_CXCL6", "iCAF_CXCL12", "iCAF_TNC", "endo_PLVAP", "vCAF")
fibro_colors <- c(
  "myCAF_MMP11" = "#31A354",
  "myCAF_TNFRSF21" = "#E78AC3",
  "myCAF_DUX4" = "#8D6E63",
  "myCAF_COL10A1" = "#BDB76B",
  "iCAF_CXCL5" = "#F28E2B",
  "iCAF_CXCL6" = "#8C6BB1",
  "iCAF_CXCL12" = "#B3A369",
  "iCAF_TNC" = "#4DBBD5",
  "endo_PLVAP" = "#2C7FB8",
  "vCAF" = "#C44E52"
)

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "fibro_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, final_cell_type2, cCR, sample_timepoint, X_umap_1, X_umap_2)
)) %>%
  mutate(final_cell_type2 = factor(final_cell_type2, levels = fibro_order))

write_csv(count(umap, final_cell_type2, name = "n_cells"), file.path(out_dir, "fibro_umap_fine_subtype_counts.csv"))
write_csv(data.frame(cell_type = names(fibro_colors), color = unname(fibro_colors)), file.path(out_dir, "fibro_umap_fine_subtype_palette.csv"))
save_palette_labels(fibro_colors, file.path(out_dir, "fibro_umap_fine_subtype_color_labels"), plotting_config$font$family, plotting_config$font$standard_pt, width_mm = 82)

plot <- make_subset_umap(
  umap,
  color_col = "final_cell_type2",
  colors = fibro_colors,
  title = "Fibroblast and endothelial subtypes",
  point_size = 0.19,
  base_family = plotting_config$font$family,
  base_size = plotting_config$font$standard_pt,
  show_labels = FALSE
)

save_umap_raster(plot, file.path(out_dir, "fibro_umap_fine_subtypes"), width_mm = 115, height_mm = 95, dpi = 600)
