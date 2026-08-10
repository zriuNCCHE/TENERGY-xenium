suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

data_path <- file.path(project_root, "data", "raw", "combined_final_all_cell_types_obs.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

timepoint_order <- c("pre", "postC", "postA")
group_order <- c("cCR", "non-cCR")
compartment_order <- c("malignant", "A2ML1+ epithelial", "immune", "fibro/stromal")
compartment_colors <- unlist(palette_config$broad_compartments)[compartment_order]

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
part_width <- plotting_config$part_sizes_mm$broad_compartment_composition$width
part_height <- plotting_config$part_sizes_mm$broad_compartment_composition$height

immune_types <- c("B/Plasma", "cDC3_LAMP3", "Macro_CXCL5", "Macro_CXCL9", "Mono_CDC27", "Mono_SLC2A3", "Neutrophil", "T/NK")
fibro_stromal_types <- c("myCAF", "iCAF", "vCAF", "endothelial")

metadata <- suppressMessages(
  read_csv(
    data_path,
    show_col_types = FALSE,
    col_select = c(sampleID, sample_timepoint, patientID, cCR, final_cell_type)
  )
) %>%
  mutate(
    compartment = case_when(
      final_cell_type == "malignant" ~ "malignant",
      final_cell_type == "A2ML1+ epi" ~ "A2ML1+ epithelial",
      final_cell_type %in% immune_types ~ "immune",
      final_cell_type %in% fibro_stromal_types ~ "fibro/stromal",
      TRUE ~ "fibro/stromal"
    ),
    sample_timepoint = factor(sample_timepoint, levels = timepoint_order),
    cCR = factor(cCR, levels = group_order),
    compartment = factor(compartment, levels = compartment_order)
  )

composition <- metadata %>%
  count(patientID, cCR, sampleID, sample_timepoint, compartment, name = "n_cells") %>%
  complete(nesting(patientID, cCR, sampleID, sample_timepoint), compartment, fill = list(n_cells = 0)) %>%
  group_by(patientID, cCR, sampleID, sample_timepoint) %>%
  mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells * 100) %>%
  ungroup() %>%
  arrange(cCR, suppressWarnings(as.integer(sub("^p", "", patientID))), patientID, sample_timepoint) %>%
  mutate(sample_label = factor(sampleID, levels = unique(sampleID)))

write_csv(composition, file.path(out_dir, "broad_compartment_composition_summary.csv"))

plot <- ggplot(composition, aes(x = sample_label, y = proportion, fill = compartment)) +
  geom_col(width = 0.82, color = "white", linewidth = line_width * 0.4) +
  facet_grid(. ~ cCR, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = compartment_colors, name = "Compartment") +
  scale_y_continuous(limits = c(0, 100.5), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Proportion of all cells (%)") +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_text(size = plotting_config$font$small_pt, angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.35, "cm"),
    plot.margin = margin(4, 4, 16, 4)
  )

ggsave(
  filename = file.path(out_dir, "broad_compartment_composition.pdf"),
  plot = plot,
  device = "pdf",
  width = part_width,
  height = part_height,
  units = "mm",
  dpi = 300
)
