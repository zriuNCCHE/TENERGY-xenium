suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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
timepoint_colors <- unlist(palette_config$timepoints)[timepoint_order]

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
part_width <- plotting_config$part_sizes_mm$supplementary_qc_overview$width
part_height <- plotting_config$part_sizes_mm$supplementary_qc_overview$height

metadata <- suppressMessages(
  read_csv(
    data_path,
    show_col_types = FALSE,
    col_select = c(
      sampleID, sample_timepoint, patientID, cCR,
      transcript_counts, n_genes, cell_area, nucleus_area, nucleus_ratio,
      low_quality_200, gene_low_quality_200
    )
  )
)

sample_summary <- metadata %>%
  mutate(
    sample_timepoint = factor(sample_timepoint, levels = timepoint_order),
    cCR = factor(cCR, levels = group_order),
    low_quality_200 = low_quality_200 == "low_quality",
    gene_low_quality_200 = gene_low_quality_200 == "low_quality"
  ) %>%
  group_by(patientID, cCR, sampleID, sample_timepoint) %>%
  summarise(
    total_cells = n(),
    median_transcripts = median(transcript_counts, na.rm = TRUE),
    median_genes = median(n_genes, na.rm = TRUE),
    median_cell_area = median(cell_area, na.rm = TRUE),
    median_nucleus_area = median(nucleus_area, na.rm = TRUE),
    median_nucleus_ratio = median(nucleus_ratio, na.rm = TRUE),
    low_quality_cell_fraction = mean(low_quality_200, na.rm = TRUE),
    low_gene_quality_fraction = mean(gene_low_quality_200, na.rm = TRUE),
    segmented_cell_area_mm2 = sum(cell_area, na.rm = TRUE) / 1e6,
    cells_per_segmented_mm2 = total_cells / segmented_cell_area_mm2,
    .groups = "drop"
  ) %>%
  arrange(cCR, suppressWarnings(as.integer(sub("^p", "", patientID))), patientID, sample_timepoint) %>%
  mutate(sample_label = factor(sampleID, levels = sampleID))

write_csv(sample_summary, file.path(out_dir, "qc_overview_sample_summary.csv"))

qc_long <- sample_summary %>%
  select(sampleID, sample_timepoint, cCR, median_transcripts, median_genes, median_cell_area, median_nucleus_ratio) %>%
  pivot_longer(
    cols = c(median_transcripts, median_genes, median_cell_area, median_nucleus_ratio),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("median_transcripts", "median_genes", "median_cell_area", "median_nucleus_ratio"),
      labels = c("Median transcripts", "Median genes", "Median cell area", "Median nucleus ratio")
    )
  )

area_long <- sample_summary %>%
  select(patientID, sampleID, sample_timepoint, cCR, segmented_cell_area_mm2, cells_per_segmented_mm2) %>%
  pivot_longer(
    cols = c(segmented_cell_area_mm2, cells_per_segmented_mm2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("segmented_cell_area_mm2", "cells_per_segmented_mm2"),
      labels = c("Summed segmented cell area (mm2)", "Cells per segmented mm2")
    )
  )

theme_pub <- theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = plotting_config$font$small_pt),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width)
  )

p_cells <- ggplot(sample_summary, aes(x = sample_label, y = total_cells, fill = sample_timepoint)) +
  geom_col(width = 0.72, color = "#333333", linewidth = line_width) +
  facet_grid(. ~ cCR, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = timepoint_colors, name = "Timepoint") +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = NULL, y = "Cells per sample") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p_qc <- ggplot(qc_long, aes(x = sample_timepoint, y = value, color = sample_timepoint)) +
  geom_boxplot(aes(fill = sample_timepoint), width = 0.55, outlier.shape = NA, alpha = 0.25, linewidth = line_width) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 0.8, alpha = 0.85, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = timepoint_colors, guide = "none") +
  scale_color_manual(values = timepoint_colors, guide = "none") +
  labs(x = NULL, y = "Sample median") +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_area <- ggplot(area_long, aes(x = sample_timepoint, y = value, color = sample_timepoint)) +
  geom_line(aes(group = patientID), color = "#BDBDBD", linewidth = line_width, alpha = 0.55) +
  geom_point(size = 1.0, alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = timepoint_colors, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

plot <- (p_cells / p_qc / p_area) +
  plot_layout(heights = c(1.2, 1.0, 0.95), guides = "collect") &
  theme(legend.position = "top")

ggsave(
  filename = file.path(out_dir, "qc_overview.pdf"),
  plot = plot,
  device = "pdf",
  width = part_width,
  height = part_height,
  units = "mm",
  dpi = 300
)
