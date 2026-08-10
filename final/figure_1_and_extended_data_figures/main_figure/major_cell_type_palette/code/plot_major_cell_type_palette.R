suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

major_colors <- unlist(palette_config$figure2_major_cell_types)
major_order <- names(major_colors)

palette_df <- data.frame(
  major_cell_type = factor(major_order, levels = major_order),
  color = unname(major_colors)
)
write.csv(palette_df, file.path(out_dir, "figure2_major_cell_type_palette.csv"), row.names = FALSE)

plot <- ggplot(palette_df, aes(x = major_cell_type, y = 1, fill = major_cell_type)) +
  geom_tile(width = 0.92, height = 0.55, color = "#333333", linewidth = 0.2) +
  geom_text(aes(label = major_cell_type), y = 0.28, size = 2.0, family = plotting_config$font$family, color = "#222222") +
  scale_fill_manual(values = major_colors, guide = "none") +
  coord_cartesian(ylim = c(0.05, 1.35), clip = "off") +
  labs(x = NULL, y = NULL, title = "Figure 2 major cell type palette") +
  theme_void(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
    plot.margin = margin(4, 4, 12, 4)
  )

ggsave(
  file.path(out_dir, "figure2_major_cell_type_palette.pdf"),
  plot = plot,
  device = "pdf",
  width = 120,
  height = 35,
  units = "mm",
  dpi = 300,
  bg = "white"
)
