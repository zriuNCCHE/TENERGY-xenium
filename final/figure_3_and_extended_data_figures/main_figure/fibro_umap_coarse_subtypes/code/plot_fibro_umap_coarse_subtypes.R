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

fibro_order <- c("myCAF", "iCAF", "endothelial", "vCAF")
fibro_colors <- c("myCAF" = "#0072B2", "iCAF" = "#D55E00", "endothelial" = "#009E73", "vCAF" = "#E69F00")

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "fibro_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, final_cell_type, cCR, sample_timepoint, X_umap_1, X_umap_2)
)) %>%
  mutate(final_cell_type = factor(final_cell_type, levels = fibro_order))

write_csv(count(umap, final_cell_type, name = "n_cells"), file.path(out_dir, "fibro_umap_coarse_subtype_counts.csv"))
write_csv(data.frame(cell_type = names(fibro_colors), color = unname(fibro_colors)), file.path(out_dir, "fibro_umap_coarse_subtype_palette.csv"))

plot <- make_subset_umap(
  umap,
  color_col = "final_cell_type",
  colors = fibro_colors,
  title = "Fibroblast and endothelial states",
  point_size = 0.09,
  base_family = plotting_config$font$family,
  base_size = plotting_config$font$standard_pt
)

save_umap_raster(plot, file.path(out_dir, "fibro_umap_coarse_subtypes"), width_mm = 95, height_mm = 85, dpi = 600)
