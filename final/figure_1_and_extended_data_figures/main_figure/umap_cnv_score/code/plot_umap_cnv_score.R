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

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "combined_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, cnv_score, X_umap_1, X_umap_2)
)) %>%
  arrange(cnv_score)

write_csv(
  data.frame(
    metric = c("min", "median", "mean", "max"),
    cnv_score = c(min(umap$cnv_score, na.rm = TRUE), median(umap$cnv_score, na.rm = TRUE), mean(umap$cnv_score, na.rm = TRUE), max(umap$cnv_score, na.rm = TRUE))
  ),
  file.path(out_dir, "umap_cnv_score_summary.csv")
)

limits <- get_umap_limits(umap)
cnv_range <- range(umap$cnv_score, na.rm = TRUE)
cnv_quantiles <- as.numeric(quantile(
  umap$cnv_score,
  probs = c(0.50, 0.90),
  na.rm = TRUE
))
cnv_scale_values <- (c(cnv_range[1], cnv_quantiles[1], cnv_quantiles[2], cnv_range[2]) - cnv_range[1]) / diff(cnv_range)

plot <- ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = cnv_score)) +
  geom_point(size = 0.08, alpha = 1, stroke = 0) +
  scale_color_gradientn(
    colors = c("#2166AC", "#2166AC", "#F4C430", "#B2182B"),
    values = cnv_scale_values,
    limits = cnv_range,
    name = "CNV score"
  ) +
  coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
  labs(title = "CNV score") +
  theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
  theme(legend.title = element_text(size = plotting_config$font$standard_pt))

save_umap_raster(
  plot,
  file.path(out_dir, "umap_cnv_score"),
  width_mm = 92,
  height_mm = 78,
  dpi = 600
)
