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

contrast_trials <- list(
  stronger = list(
    label = "Blue-yellow-red, stronger",
    colors = c("#1F4E99", "#F4C430", "#B2182B"),
    alpha = 0.75,
    size = 0.08,
    limits = range(umap$cnv_score, na.rm = TRUE),
    oob = scales::squish
  ),
  darker = list(
    label = "Blue-gold-red, darker",
    colors = c("#08306B", "#FDBF11", "#99000D"),
    alpha = 0.85,
    size = 0.08,
    limits = range(umap$cnv_score, na.rm = TRUE),
    oob = scales::squish
  ),
  clipped = list(
    label = "Blue-yellow-red, clipped",
    colors = c("#08306B", "#FEC44F", "#A50F15"),
    alpha = 0.85,
    size = 0.08,
    limits = as.numeric(quantile(umap$cnv_score, probs = c(0.01, 0.99), na.rm = TRUE)),
    oob = scales::squish
  )
)

make_plot <- function(trial) {
  ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = cnv_score)) +
    geom_point(size = trial$size, alpha = trial$alpha, stroke = 0) +
    scale_color_gradientn(
      colors = trial$colors,
      limits = trial$limits,
      oob = trial$oob,
      name = "CNV score"
    ) +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = trial$label) +
    theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
    theme(legend.title = element_text(size = plotting_config$font$standard_pt))
}

for (trial_name in names(contrast_trials)) {
  ggsave(
    filename = file.path(out_dir, paste0("trial_cnv_umap_blue_yellow_red_", trial_name, ".png")),
    plot = make_plot(contrast_trials[[trial_name]]),
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
    trial = names(contrast_trials),
    label = vapply(contrast_trials, `[[`, character(1), "label"),
    low_color = vapply(contrast_trials, function(x) x$colors[1], character(1)),
    mid_color = vapply(contrast_trials, function(x) x$colors[2], character(1)),
    high_color = vapply(contrast_trials, function(x) x$colors[3], character(1)),
    alpha = vapply(contrast_trials, `[[`, numeric(1), "alpha"),
    lower_limit = vapply(contrast_trials, function(x) x$limits[1], numeric(1)),
    upper_limit = vapply(contrast_trials, function(x) x$limits[2], numeric(1))
  ),
  file.path(out_dir, "cnv_blue_yellow_red_contrast_trials.csv")
)
