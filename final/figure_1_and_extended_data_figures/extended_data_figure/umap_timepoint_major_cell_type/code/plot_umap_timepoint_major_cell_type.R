suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(normalizePath(file.path(dirname(script_path), "..", "..", "..", "..")), "scripts", "plotting", "umap_helpers.R"))

project_root <- get_project_root(script_path)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
major_colors <- unlist(palette_config$figure2_major_cell_types)
major_order <- names(major_colors)
timepoint_order <- c("pre", "postC", "postA")

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "combined_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, sample_timepoint, major_cell_type, X_umap_1, X_umap_2)
)) %>%
  mutate(
    major_cell_type = recode_major_cell_type(major_cell_type),
    major_cell_type = factor(major_cell_type, levels = major_order),
    sample_timepoint = factor(sample_timepoint, levels = timepoint_order)
  )

write_csv(
  count(umap, sample_timepoint, major_cell_type, name = "n_cells"),
  file.path(out_dir, "umap_timepoint_major_cell_type_counts.csv")
)

limits <- get_umap_limits(umap)
plot <- ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = major_cell_type)) +
  geom_point(size = 0.075, alpha = 1, stroke = 0) +
  facet_wrap(~ sample_timepoint, nrow = 1) +
  scale_color_manual(values = major_colors, drop = FALSE) +
  coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  labs(title = NULL) +
  theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = plotting_config$font$standard_pt, face = "bold"),
    strip.background = element_blank()
  )

save_umap_raster(
  plot,
  file.path(out_dir, "umap_timepoint_major_cell_type"),
  width_mm = 180,
  height_mm = 72,
  dpi = 600
)
