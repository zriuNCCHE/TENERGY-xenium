suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
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

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "fibro_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, sampleID, X_umap_1, X_umap_2)
)) %>%
  mutate(sampleID = factor(sampleID, levels = sort(unique(sampleID))))

sample_colors <- setNames(hue_pal(l = 35, c = 120, h = c(10, 370))(length(levels(umap$sampleID))), levels(umap$sampleID))
write_csv(data.frame(sampleID = names(sample_colors), color = unname(sample_colors)), file.path(out_dir, "fibro_umap_sample_id_palette.csv"))

limits <- get_umap_limits(umap)
plot <- ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = sampleID)) +
  geom_point(size = 0.15, alpha = 1, stroke = 0) +
  scale_color_manual(values = sample_colors, drop = FALSE, name = "Sample ID") +
  coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1), nrow = 6, byrow = TRUE)) +
  labs(title = "Fibroblast and endothelial cells by sample ID") +
  theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = plotting_config$font$standard_pt),
    legend.text = element_text(size = plotting_config$font$small_pt),
    legend.box.margin = margin(2, 6, 2, 6),
    plot.margin = margin(4, 8, 6, 8)
  )

save_palette_labels(sample_colors, file.path(out_dir, "fibro_umap_sample_id_color_labels"), plotting_config$font$family, plotting_config$font$small_pt, width_mm = 105, height_mm = 150)
save_umap_raster(plot, file.path(out_dir, "fibro_umap_sample_id"), width_mm = 180, height_mm = 135, dpi = 600)
