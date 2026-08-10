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

palettes <- list(
  blue_yellow_red = c("#2166AC", "#FFFFBF", "#B2182B"),
  grey_purple = c("#F2F2F2", "#B2ABD2", "#542788"),
  grey_orange_red = c("#F7F7F7", "#FDAE61", "#B2182B"),
  black_green = c("#000000", "#006837", "#66BD63", "#D9F0A3")
)

palette_labels <- c(
  blue_yellow_red = "Blue-yellow-red",
  grey_purple = "Grey-purple",
  grey_orange_red = "Grey-orange-red",
  black_green = "Black-green"
)

make_plot <- function(colors, title) {
  ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = cnv_score)) +
    geom_point(size = 0.08, alpha = 0.55, stroke = 0) +
    scale_color_gradientn(colors = colors, name = "CNV score") +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = title) +
    theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
    theme(legend.title = element_text(size = plotting_config$font$standard_pt))
}

for (palette_name in names(palettes)) {
  ggsave(
    filename = file.path(out_dir, paste0("trial_cnv_umap_", palette_name, ".png")),
    plot = make_plot(palettes[[palette_name]], palette_labels[[palette_name]]),
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
    palette = rep(names(palettes), lengths(palettes)),
    label = rep(unname(palette_labels), lengths(palettes)),
    color = unlist(palettes, use.names = FALSE)
  ),
  file.path(out_dir, "cnv_palette_trials.csv")
)
