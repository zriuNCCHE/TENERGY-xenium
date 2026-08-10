suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(yaml)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

mp2_data_dir <- file.path(project_root, "data", "geo_ready", "figure_7_mp_cd4cxcl13_distance")
mp7_data_dir <- file.path(project_root, "data", "geo_ready", "figure_7_mp7_myeloid_distance")
out_dir <- file.path(project_root, "final", "extended_figure_7", "mp_spatial_proximity_advantage", "outputs")
legend_dir <- file.path(project_root, "final", "extended_figure_7", "mp_spatial_proximity_advantage", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

mp_levels <- paste0("MP_", 1:7)
mp_palette <- unlist(palette_config$malignant_metaprograms)
mp_palette <- mp_palette[mp_levels]

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 7, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "p NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
}

mm_pdf <- function(filename, width_mm, height_mm) {
  grDevices::pdf(
    file.path(out_dir, filename),
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    useDingbats = FALSE,
    family = "Helvetica"
  )
}

make_delta_table <- function(sample_summary, reference_mp) {
  rbindlist(lapply(setdiff(mp_levels, reference_mp), function(other_mp) {
    ref <- sample_summary[MP == reference_mp, .(sampleID, ref_distance = median_distance)]
    other <- sample_summary[MP == other_mp, .(sampleID, other_distance = median_distance)]
    paired <- merge(ref, other, by = "sampleID")
    paired[, other_MP := other_mp]
    paired[, reference_MP := reference_mp]
    paired[, delta_other_minus_reference := other_distance - ref_distance]
    paired
  }))
}

make_delta_plot <- function(delta_dt, stats_dt, reference_mp, title, filename, highlight_color) {
  delta_dt[, other_MP := factor(other_MP, levels = setdiff(mp_levels, reference_mp))]
  stats_dt <- stats_dt[other_MP %in% setdiff(mp_levels, reference_mp)]
  stats_dt[, other_MP := factor(other_MP, levels = setdiff(mp_levels, reference_mp))]
  stats_dt[, label := fmt_p(p_adjusted_BH)]
  y_top <- max(delta_dt$delta_other_minus_reference, na.rm = TRUE)
  y_min <- min(delta_dt$delta_other_minus_reference, na.rm = TRUE)
  y_pad <- max(10, 0.10 * (y_top - y_min))

  p <- ggplot(delta_dt, aes(other_MP, delta_other_minus_reference)) +
    geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#555555") +
    geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.5, fill = "white", colour = "#7FAEC8") +
    geom_point(size = 0.9, alpha = 1, position = position_jitter(width = 0.08, height = 0), colour = "#4D4D4D") +
    geom_text(
      data = stats_dt,
      aes(x = other_MP, y = y_top + y_pad, label = label),
      inherit.aes = FALSE,
      size = 1.9,
      vjust = 0
    ) +
    labs(
      title = title,
      x = "Comparison metaprogram",
      y = paste0("Median distance(other MP) - median distance(", reference_mp, ")")
    ) +
    theme_nc()
  mm_pdf(filename, 92, 68)
  print(p)
  dev.off()
}

mp2_sample <- fread(file.path(mp2_data_dir, "MP_to_CD4_CXCL13_nearest_distance_sample_summary.csv"))
mp2_sample[, MP := factor(MP, levels = mp_levels)]
mp2_delta <- make_delta_table(mp2_sample, "MP_2")
mp2_stats <- fread(file.path(project_root, "exploratory", "figure_7", "mp2_cd4cxcl13_distance_candidates", "outputs", "figure7_recomputed_sample_level_pairwise_stats.csv"))
setnames(mp2_stats, "comparison", "comparison_original")
mp2_stats[, other_MP := factor(other_MP, levels = setdiff(mp_levels, "MP_2"))]
fwrite(mp2_delta, file.path(out_dir, "extended_figure7_mp2_cd4cxcl13_delta_by_sample.csv"))
make_delta_plot(
  mp2_delta,
  mp2_stats,
  reference_mp = "MP_2",
  title = "MP2 proximity advantage to CD4_CXCL13",
  filename = "extended_figure7_mp2_cd4cxcl13_proximity_advantage.pdf",
  highlight_color = mp_palette["MP_2"]
)

mp7_sample <- fread(file.path(mp7_data_dir, "MP_to_monocyte_neutrophil_nearest_distance_sample_summary.csv"))
setnames(mp7_sample, "query_label", "MP")
mp7_sample <- mp7_sample[target_label == "monocyte_neutrophil"]
mp7_sample[, MP := factor(MP, levels = mp_levels)]
mp7_delta <- make_delta_table(mp7_sample, "MP_7")
mp7_stats <- fread(file.path(mp7_data_dir, "MP_to_monocyte_neutrophil_nearest_distance_monocyte_neutrophil_sample_level_stats.csv"))
mp7_stats[, other_MP := factor(other_MP, levels = setdiff(mp_levels, "MP_7"))]
fwrite(mp7_delta, file.path(out_dir, "extended_figure7_mp7_myeloid_delta_by_sample.csv"))
make_delta_plot(
  mp7_delta,
  mp7_stats,
  reference_mp = "MP_7",
  title = "MP7 proximity advantage to myeloid cells",
  filename = "extended_figure7_mp7_myeloid_proximity_advantage.pdf",
  highlight_color = mp_palette["MP_7"]
)

legend_text <- c(
  "Extended Figure 7 proximity-advantage legends",
  "",
  "MP2 proximity advantage. For each sample, the MP2 sample-level median distance to CD4_CXCL13 T cells was subtracted from the sample-level median distance for each comparison metaprogram. Values above zero indicate that MP2 is closer to CD4_CXCL13 T cells than the comparison metaprogram. P values are Benjamini-Hochberg-adjusted paired Wilcoxon P values.",
  "",
  "MP7 proximity advantage. For each sample, the MP7 sample-level median distance to pooled myeloid cells was subtracted from the sample-level median distance for each comparison metaprogram. Values above zero indicate that MP7 is closer to pooled myeloid cells than the comparison metaprogram. P values are Benjamini-Hochberg-adjusted paired Wilcoxon P values."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# Extended Figure 7 MP Spatial Proximity Advantage",
  "",
  "Supplementary proximity-advantage plots for the Figure 7 MP2-CD4_CXCL13 and MP7-myeloid spatial proximity story.",
  "",
  "## Outputs",
  "",
  "- `outputs/extended_figure7_mp2_cd4cxcl13_proximity_advantage.pdf`",
  "- `outputs/extended_figure7_mp7_myeloid_proximity_advantage.pdf`"
)
writeLines(readme_text, file.path(dirname(out_dir), "README.md"))

message("Saved Extended Figure 7 proximity-advantage panels to: ", out_dir)
