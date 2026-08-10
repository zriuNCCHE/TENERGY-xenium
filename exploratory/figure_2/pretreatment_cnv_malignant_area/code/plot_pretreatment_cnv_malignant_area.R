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
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
outcome_colors <- unlist(palette_config$clinical_outcome)
outcome_order <- c("cCR", "non-cCR")

obs <- suppressMessages(read_csv(
  file.path(project_root, "data", "geo_ready", "combined_adata_obs_with_all_obsm.csv"),
  show_col_types = FALSE,
  col_select = c(cell_id, sampleID, patientID, sample_timepoint, cCR, major_cell_type, CNV_high_status, cell_area)
)) %>%
  filter(sample_timepoint == "pre") %>%
  mutate(
    cCR = factor(cCR, levels = outcome_order),
    cnv_high_area = if_else(CNV_high_status == "CNV_high", cell_area, 0),
    malignant_area = if_else(major_cell_type == "malignant", cell_area, 0)
  )

sample_summary <- obs %>%
  group_by(sampleID, patientID, cCR) %>%
  summarise(
    total_cell_area_um2 = sum(cell_area, na.rm = TRUE),
    cnv_high_area_um2 = sum(cnv_high_area, na.rm = TRUE),
    malignant_area_um2 = sum(malignant_area, na.rm = TRUE),
    n_cells = n(),
    n_cnv_high = sum(CNV_high_status == "CNV_high", na.rm = TRUE),
    n_malignant = sum(major_cell_type == "malignant", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cnv_high_area_mm2 = cnv_high_area_um2 / 1e6,
    malignant_area_mm2 = malignant_area_um2 / 1e6,
    cnv_high_area_fraction = cnv_high_area_um2 / total_cell_area_um2,
    malignant_area_fraction = malignant_area_um2 / total_cell_area_um2
  )

plot_df <- sample_summary %>%
  select(sampleID, patientID, cCR, cnv_high_area_mm2, malignant_area_mm2, cnv_high_area_fraction, malignant_area_fraction) %>%
  pivot_longer(
    cols = c(cnv_high_area_mm2, malignant_area_mm2, cnv_high_area_fraction, malignant_area_fraction),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    compartment = case_when(
      grepl("^cnv_high", metric) ~ "CNV-high cells",
      grepl("^malignant", metric) ~ "Malignant cells",
      TRUE ~ metric
    ),
    scale = case_when(
      grepl("mm2$", metric) ~ "Area (mm2)",
      grepl("fraction$", metric) ~ "Area fraction (%)",
      TRUE ~ metric
    ),
    value_plot = if_else(scale == "Area fraction (%)", value * 100, value)
  )

stats <- plot_df %>%
  group_by(compartment, scale) %>%
  summarise(
    n_cCR = sum(cCR == "cCR"),
    n_non_cCR = sum(cCR == "non-cCR"),
    median_cCR = median(value_plot[cCR == "cCR"], na.rm = TRUE),
    median_non_cCR = median(value_plot[cCR == "non-cCR"], na.rm = TRUE),
    p_value = wilcox.test(value_plot ~ cCR, exact = FALSE)$p.value,
    .groups = "drop"
  )

write_csv(sample_summary, file.path(out_dir, "pretreatment_cnv_malignant_area_sample_summary.csv"))
write_csv(stats, file.path(out_dir, "pretreatment_cnv_malignant_area_stats.csv"))

make_plot <- function(scale_name, y_label, output_name) {
  dat <- filter(plot_df, scale == scale_name)
  ymax <- dat %>%
    group_by(compartment) %>%
    summarise(y = max(value_plot, na.rm = TRUE) * 1.18, .groups = "drop")
  stat_labels <- stats %>%
    filter(scale == scale_name) %>%
    mutate(label = paste0("p = ", signif(p_value, 3))) %>%
    left_join(ymax, by = "compartment")

  plot <- ggplot(dat, aes(x = cCR, y = value_plot, fill = cCR, color = cCR)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5, alpha = 0.25) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.6, alpha = 1) +
    facet_wrap(~ compartment, scales = "free_y", nrow = 1) +
    geom_text(
      data = stat_labels,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = plotting_config$font$standard_pt / 2.85,
      family = plotting_config$font$family
    ) +
    scale_fill_manual(values = outcome_colors, drop = FALSE) +
    scale_color_manual(values = outcome_colors, drop = FALSE) +
    labs(x = NULL, y = y_label, title = "Pretreatment CNV-high and malignant cell area") +
    theme_classic(base_family = plotting_config$font$family, base_size = plotting_config$font$standard_pt) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = plotting_config$font$standard_pt, face = "bold"),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5)
    )

  ggsave(file.path(out_dir, paste0(output_name, ".pdf")), plot, width = 86, height = 62, units = "mm", device = "pdf", bg = "white")
  ggsave(file.path(out_dir, paste0(output_name, ".png")), plot, width = 86, height = 62, units = "mm", dpi = 600, device = "png", bg = "white")
}

make_plot("Area (mm2)", "Summed cell area (mm2)", "pretreatment_cnv_malignant_absolute_area")
make_plot("Area fraction (%)", "Fraction of total cell area (%)", "pretreatment_cnv_malignant_area_fraction")
