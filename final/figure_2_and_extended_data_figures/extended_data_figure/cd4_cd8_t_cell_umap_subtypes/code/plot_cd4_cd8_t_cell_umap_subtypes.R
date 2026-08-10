suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
source(file.path(project_root, "scripts", "plotting", "umap_helpers.R"))
source(file.path(project_root, "scripts", "plotting", "subset_umap_helpers.R"))

out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

cd4_order <- c("CD4_Treg_FOXP3", "CD4_Treg_CCR8", "CD4_CXCL13", "CD4_prolif")
cd4_colors <- c(
  "CD4_Treg_FOXP3" = "#7B3294",
  "CD4_Treg_CCR8" = "#C51B7D",
  "CD4_CXCL13" = "#008837",
  "CD4_prolif" = "#E6AB02"
)

cd8_order <- c("CD8_Teff", "CD8_Tex_PDCD1", "CD8_prolif")
cd8_colors <- c(
  "CD8_Teff" = "#1F78B4",
  "CD8_Tex_PDCD1" = "#D95F02",
  "CD8_prolif" = "#7570B3"
)

read_t_cell_umap <- function(path, cell_type_order) {
  suppressMessages(read_csv(
    path,
    show_col_types = FALSE,
    col_select = c(
      cell_id,
      final_cell_type,
      final_cell_type_2,
      cCR,
      sampleID,
      sample_timepoint,
      patientID,
      X_umap_1,
      X_umap_2
    )
  )) %>%
    rename(t_cell_state = final_cell_type_2) %>%
    mutate(t_cell_state = factor(t_cell_state, levels = cell_type_order))
}

cd4_umap <- read_t_cell_umap(
  file.path(project_root, "data", "geo_ready", "cd4_adata_obs_with_all_obsm.csv"),
  cd4_order
)
cd8_umap <- read_t_cell_umap(
  file.path(project_root, "data", "geo_ready", "cd8_adata_obs_with_all_obsm.csv"),
  cd8_order
)

write_csv(count(cd4_umap, t_cell_state, name = "n_cells"), file.path(out_dir, "cd4_t_cell_umap_state_counts.csv"))
write_csv(count(cd8_umap, t_cell_state, name = "n_cells"), file.path(out_dir, "cd8_t_cell_umap_state_counts.csv"))
write_csv(data.frame(cell_type = names(cd4_colors), color = unname(cd4_colors)), file.path(out_dir, "cd4_t_cell_umap_state_palette.csv"))
write_csv(data.frame(cell_type = names(cd8_colors), color = unname(cd8_colors)), file.path(out_dir, "cd8_t_cell_umap_state_palette.csv"))

save_palette_labels(
  cd4_colors,
  file.path(out_dir, "cd4_t_cell_umap_state_color_labels"),
  plotting_config$font$family,
  plotting_config$font$standard_pt,
  width_mm = 70
)
save_palette_labels(
  cd8_colors,
  file.path(out_dir, "cd8_t_cell_umap_state_color_labels"),
  plotting_config$font$family,
  plotting_config$font$standard_pt,
  width_mm = 70
)

cd4_plot <- make_subset_umap(
  cd4_umap,
  color_col = "t_cell_state",
  colors = cd4_colors,
  title = "CD4 T cell states",
  point_size = 0.34,
  base_family = plotting_config$font$family,
  base_size = plotting_config$font$standard_pt,
  show_labels = FALSE
)

cd8_plot <- make_subset_umap(
  cd8_umap,
  color_col = "t_cell_state",
  colors = cd8_colors,
  title = "CD8 T cell states",
  point_size = 0.34,
  base_family = plotting_config$font$family,
  base_size = plotting_config$font$standard_pt,
  show_labels = FALSE
)

save_umap_raster(cd4_plot, file.path(out_dir, "cd4_t_cell_umap_subtypes"), width_mm = 82, height_mm = 72, dpi = 600)
save_umap_raster(cd8_plot, file.path(out_dir, "cd8_t_cell_umap_subtypes"), width_mm = 82, height_mm = 72, dpi = 600)

combined_plot <- cd4_plot + cd8_plot +
  plot_layout(nrow = 1, widths = c(1, 1)) &
  theme(plot.margin = margin(2, 2, 2, 2))

ggsave(
  filename = file.path(out_dir, "cd4_cd8_t_cell_umap_subtypes.png"),
  plot = combined_plot,
  width = 170,
  height = 76,
  units = "mm",
  dpi = 600,
  bg = "white"
)
ggsave(
  filename = file.path(out_dir, "cd4_cd8_t_cell_umap_subtypes.tiff"),
  plot = combined_plot,
  width = 170,
  height = 76,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

writeLines(
  c(
    "# CD4 and CD8 T cell UMAPs",
    "",
    "Extended Data Figure 3 draft part. UMAPs show CD4 and CD8 T cell states colored by final_cell_type_2 from the exported AnnData observation table. The source final_cell_type column contains only the parent T/NK label in these subset exports.",
    "",
    "Points represent individual cells. Colors are fixed in the accompanying palette CSV files and should be reused for downstream CD4/CD8 T cell state plots."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
