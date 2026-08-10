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

cell_type_order <- c(
  "A2ML1+ epithelial cells",
  "fibroblast & endothelial cells",
  "immune cells",
  "malignant cells"
)

palettes <- list(
  okabe_ito_no_he = c(
    "A2ML1+ epithelial cells" = "#E69F00",
    "fibroblast & endothelial cells" = "#999999",
    "immune cells" = "#009E73",
    "malignant cells" = "#D55E00"
  ),
  colorbrewer_dark2_no_he = c(
    "A2ML1+ epithelial cells" = "#E6AB02",
    "fibroblast & endothelial cells" = "#666666",
    "immune cells" = "#1B9E77",
    "malignant cells" = "#D95F02"
  ),
  tableau_muted_no_he = c(
    "A2ML1+ epithelial cells" = "#EDC948",
    "fibroblast & endothelial cells" = "#79706E",
    "immune cells" = "#59A14F",
    "malignant cells" = "#F28E2B"
  ),
  nature_muted_no_he = c(
    "A2ML1+ epithelial cells" = "#C49A00",
    "fibroblast & endothelial cells" = "#6F6F6F",
    "immune cells" = "#007A3D",
    "malignant cells" = "#C43C00"
  )
)

palette_labels <- c(
  okabe_ito_no_he = "Okabe-Ito without H&E-like colors",
  colorbrewer_dark2_no_he = "ColorBrewer Dark2-like without H&E colors",
  tableau_muted_no_he = "Tableau muted without H&E-like colors",
  nature_muted_no_he = "Muted publication palette without H&E-like colors"
)

umap <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "combined_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, major_cell_type, X_umap_1, X_umap_2)
)) %>%
  mutate(
    major_cell_type = recode_major_cell_type(major_cell_type),
    major_cell_type = factor(major_cell_type, levels = cell_type_order)
  )

limits <- get_umap_limits(umap)

make_plot <- function(colors, title) {
  ggplot(umap, aes(x = X_umap_1, y = X_umap_2, color = major_cell_type)) +
    geom_point(size = 0.08, alpha = 0.55, stroke = 0) +
    scale_color_manual(values = colors, drop = FALSE) +
    coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    labs(title = title) +
    theme_umap(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt)
}

for (palette_name in names(palettes)) {
  ggsave(
    filename = file.path(out_dir, paste0("trial_umap_major_cell_type_literature_", palette_name, ".png")),
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
  bind_rows(lapply(names(palettes), function(palette_name) {
    tibble(
      palette = palette_name,
      label = palette_labels[[palette_name]],
      cell_type = names(palettes[[palette_name]]),
      color = unname(palettes[[palette_name]])
    )
  })),
  file.path(out_dir, "major_cell_type_literature_palette_trials.csv")
)
