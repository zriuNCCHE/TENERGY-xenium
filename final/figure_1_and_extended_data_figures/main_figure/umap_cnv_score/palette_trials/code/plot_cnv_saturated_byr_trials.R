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

trials <- list(
  saturated = list(
    label = "Saturated blue-yellow-red",
    colors = c("#0057B8", "#FFD400", "#D7191C")
  ),
  deep = list(
    label = "Deep blue-gold-red",
    colors = c("#003F8C", "#FFB000", "#C00000")
  ),
  vivid = list(
    label = "Vivid blue-yellow-crimson",
    colors = c("#0047FF", "#FFC400", "#E0002A")
  ),
  balanced = list(
    label = "Balanced publication blue-gold-red",
    colors = c("#005A9C", "#F2B701", "#B30021")
  )
)

make_plot <- function(trial) {
  ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = cnv_score)) +
    geom_point(size = 0.08, alpha = 1, stroke = 0) +
    scale_color_gradientn(colors = trial$colors, name = "CNV score") +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = trial$label) +
    theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
    theme(legend.title = element_text(size = plotting_config$font$standard_pt))
}

for (trial_name in names(trials)) {
  ggsave(
    filename = file.path(out_dir, paste0("trial_cnv_umap_saturated_byr_", trial_name, ".png")),
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
    low_blue = vapply(trials, function(x) x$colors[1], character(1)),
    mid_yellow = vapply(trials, function(x) x$colors[2], character(1)),
    high_red = vapply(trials, function(x) x$colors[3], character(1)),
    alpha = 1
  ),
  file.path(out_dir, "cnv_saturated_byr_trials.csv")
)
