suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(normalizePath(file.path(dirname(script_path), "..", "..", "..", "..")), "scripts", "plotting", "umap_helpers.R"))

project_root <- get_project_root(script_path)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "combined_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, sampleID, patientID, X_umap_1, X_umap_2)
))
patient_order <- sort(unique(umap$patientID))
patient_colors <- setNames(
  hue_pal(l = 45, c = 110, h = c(15, 375))(length(patient_order)),
  patient_order
)
write_csv(
  data.frame(patientID = patient_order, color = unname(patient_colors)),
  file.path(out_dir, "umap_patient_id_palette.csv")
)

limits <- get_umap_limits(umap)
plot <- ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = patientID)) +
  geom_point(size = 0.075, alpha = 1, stroke = 0) +
  scale_color_manual(values = patient_colors, drop = FALSE, name = "Patient ID") +
  coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), nrow = 3, byrow = TRUE)) +
  labs(title = "Patient ID") +
  theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = plotting_config$font$standard_pt),
    legend.text = element_text(size = plotting_config$font$small_pt),
    legend.box.margin = margin(2, 6, 2, 6),
    plot.margin = margin(4, 8, 6, 8)
  )

save_umap_raster(
  plot,
  file.path(out_dir, "umap_sample_id"),
  width_mm = 132,
  height_mm = 112,
  dpi = 600
)

save_umap_raster(
  plot,
  file.path(out_dir, "umap_patient_id"),
  width_mm = 132,
  height_mm = 112,
  dpi = 600
)
