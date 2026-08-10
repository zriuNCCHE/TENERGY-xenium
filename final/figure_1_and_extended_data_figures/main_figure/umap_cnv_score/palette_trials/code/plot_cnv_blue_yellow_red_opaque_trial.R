suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", "..", ".."))
source(file.path(project_root, "scripts", "plotting", "umap_helpers.R"))

out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "combined_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, cnv_score, X_umap_1, X_umap_2)
)) %>%
  arrange(cnv_score)

limits <- get_umap_limits(umap)

plot <- ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = cnv_score)) +
  geom_point(size = 0.08, alpha = 1, stroke = 0) +
  scale_color_gradientn(
    colors = c("#2166AC", "#FFFFBF", "#B2182B"),
    name = "CNV score"
  ) +
  coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
  labs(title = "Blue-yellow-red, opaque") +
  theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
  theme(legend.title = element_text(size = plotting_config$font$standard_pt))

ggsave(
  filename = file.path(out_dir, "trial_cnv_umap_blue_yellow_red_opaque.png"),
  plot = plot,
  device = "png",
  width = 92,
  height = 78,
  units = "mm",
  dpi = 400,
  bg = "white"
)
