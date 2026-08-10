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

cnv_range <- range(umap$cnv_score, na.rm = TRUE)
cnv_quantiles <- as.numeric(quantile(
  umap$cnv_score,
  probs = c(0, 0.50, 0.65, 0.75, 0.90, 1),
  na.rm = TRUE
))

rescale_values <- function(x) {
  (x - cnv_range[1]) / diff(cnv_range)
}

trials <- list(
  blue_until_median = list(
    label = "Blue held to median",
    colors = c("#2166AC", "#2166AC", "#F4C430", "#B2182B"),
    values = rescale_values(c(cnv_range[1], cnv_quantiles[2], cnv_quantiles[5], cnv_range[2]))
  ),
  blue_until_65th = list(
    label = "Blue held to 65th percentile",
    colors = c("#2166AC", "#2166AC", "#F4C430", "#B2182B"),
    values = rescale_values(c(cnv_range[1], cnv_quantiles[3], cnv_quantiles[5], cnv_range[2]))
  ),
  blue_until_75th = list(
    label = "Blue held to 75th percentile",
    colors = c("#2166AC", "#2166AC", "#F4C430", "#B2182B"),
    values = rescale_values(c(cnv_range[1], cnv_quantiles[4], cnv_quantiles[5], cnv_range[2]))
  ),
  smoother_blue_bias = list(
    label = "Smoother blue-biased scale",
    colors = c("#08306B", "#2166AC", "#6BAED6", "#F4C430", "#B2182B"),
    values = c(0, 0.30, 0.62, 0.84, 1)
  )
)

make_plot <- function(trial) {
  ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = cnv_score)) +
    geom_point(size = 0.08, alpha = 1, stroke = 0) +
    scale_color_gradientn(
      colors = trial$colors,
      values = trial$values,
      limits = cnv_range,
      name = "CNV score"
    ) +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = trial$label) +
    theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
    theme(legend.title = element_text(size = plotting_config$font$standard_pt))
}

for (trial_name in names(trials)) {
  ggsave(
    filename = file.path(out_dir, paste0("trial_cnv_umap_blue_threshold_", trial_name, ".png")),
    plot = make_plot(trials[[trial_name]]),
    device = "png",
    width = 92,
    height = 78,
    units = "mm",
    dpi = 400,
    bg = "white"
  )
}

write_csv(
  tibble(
    trial = names(trials),
    label = vapply(trials, `[[`, character(1), "label"),
    value_positions = vapply(trials, function(x) paste(round(x$values, 4), collapse = ";"), character(1))
  ),
  file.path(out_dir, "cnv_blue_threshold_trials.csv")
)
