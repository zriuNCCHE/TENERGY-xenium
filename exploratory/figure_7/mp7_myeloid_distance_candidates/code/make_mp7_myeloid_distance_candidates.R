suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))
data_dir <- file.path(project_root, "data", "geo_ready", "figure_7_mp7_myeloid_distance")
out_dir <- file.path(project_root, "exploratory", "figure_7", "mp7_myeloid_distance_candidates", "outputs")
legend_dir <- file.path(project_root, "exploratory", "figure_7", "mp7_myeloid_distance_candidates", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.text = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 7, face = "bold"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
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

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  )
}

fmt_fdr <- function(p) {
  ifelse(is.na(p), "FDR NA",
    ifelse(p < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p))
  )
}

target_pretty <- c(
  monocyte_neutrophil = "Pooled myeloid",
  Neutrophil = "Neutrophil",
  Mono_CDC27 = "Mono_CDC27",
  Mono_SLC2A3 = "Mono_SLC2A3"
)
target_levels <- names(target_pretty)
mp_levels <- paste0("MP_", 1:7)

mp_colours <- c(
  MP_1 = "#8F8F8F",
  MP_2 = "#6B8BA4",
  MP_3 = "#9A9A9A",
  MP_4 = "#B0B0B0",
  MP_5 = "#7A6FA8",
  MP_6 = "#4E9A7A",
  MP_7 = "#D55E00"
)
mp_fill <- setNames(rep("#B8B8B8", 7), mp_levels)
mp_fill["MP_7"] <- "#D55E00"
mp_point <- setNames(rep("#4D4D4D", 7), mp_levels)
mp_point["MP_7"] <- "#B24700"

sample_summary <- fread(file.path(data_dir, "MP_to_monocyte_neutrophil_nearest_distance_sample_summary.csv"))
cell_level <- fread(file.path(data_dir, "MP_to_monocyte_neutrophil_nearest_distance_cell_level.csv"))
sample_stats <- fread(file.path(data_dir, "MP_to_monocyte_neutrophil_nearest_distance_sample_level_stats.csv"))
global_stats <- fread(file.path(data_dir, "MP_to_monocyte_neutrophil_nearest_distance_global_stats.csv"))

setnames(sample_summary, "query_label", "MP")
sample_summary <- sample_summary[target_label %in% target_levels]
sample_summary[, MP := factor(MP, levels = mp_levels)]
sample_summary[, target_label := factor(target_label, levels = target_levels, labels = target_pretty)]
sample_stats <- sample_stats[target_label %in% target_levels]
sample_stats[, target_label := factor(target_label, levels = target_levels, labels = target_pretty)]
sample_stats[, other_MP := factor(other_MP, levels = setdiff(mp_levels, "MP_7"))]
global_stats <- global_stats[target_label %in% target_levels]
global_stats[, target_label := factor(target_label, levels = target_levels, labels = target_pretty)]

target_file_stub <- c(
  monocyte_neutrophil = "pooled_myeloid",
  Neutrophil = "neutrophil",
  Mono_CDC27 = "mono_cdc27",
  Mono_SLC2A3 = "mono_slc2a3"
)

for (target_key in target_levels) {
  target_name <- unname(target_pretty[target_key])
  file_stub <- unname(target_file_stub[target_key])
  sample_dt <- sample_summary[target_label == target_name]
  global_p <- global_stats[target_label == target_name, p_value][1]

  box_plot <- ggplot(sample_dt, aes(MP, median_distance, fill = MP)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.5, colour = "black") +
    geom_point(aes(colour = MP), size = 1.35, alpha = 1, position = position_jitter(width = 0.09, height = 0), show.legend = FALSE) +
    scale_fill_manual(values = mp_fill, guide = "none") +
    scale_colour_manual(values = mp_point, guide = "none") +
    labs(
      title = paste0(target_name, " proximity by malignant MP"),
      x = "Malignant-cell metaprogram",
      y = paste0("Median nearest distance to ", target_name)
    ) +
    annotate("text", x = 1, y = Inf, vjust = 1.4, hjust = 0, size = 2.1,
             label = paste0("Friedman ", fmt_p(global_p))) +
    theme_nc()
  mm_pdf(paste0("figure7_mp7_boxplot_", file_stub, ".pdf"), 86, 70)
  print(box_plot)
  dev.off()

  cell_dt <- cell_level[target_label == target_key]
  cell_dt[, MP := factor(MP, levels = mp_levels)]
  cell_cdf <- cell_dt[!is.na(nearest_target_distance), .(
    distance = sort(nearest_target_distance),
    cumulative_fraction = seq_len(.N) / .N,
    n_cells = .N
  ), by = MP]
  cell_cdf[, MP := factor(MP, levels = mp_levels)]
  fwrite(cell_cdf, file.path(out_dir, paste0("figure7_mp7_cell_cdf_", file_stub, ".csv")))

  cell_cdf_plot <- ggplot(cell_cdf, aes(distance, cumulative_fraction, colour = MP, linewidth = MP)) +
    geom_step(lineend = "butt") +
    scale_colour_manual(values = mp_colours) +
    scale_linewidth_manual(values = c(MP_1 = 0.35, MP_2 = 0.35, MP_3 = 0.35, MP_4 = 0.35, MP_5 = 0.35, MP_6 = 0.35, MP_7 = 0.75), guide = "none") +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(0, quantile(cell_dt$nearest_target_distance, 0.985, na.rm = TRUE))) +
    labs(
      title = paste0("Cell-level ", target_name, " proximity"),
      x = paste0("Nearest distance to ", target_name),
      y = "Cumulative fraction of cells",
      colour = "Malignant MP"
    ) +
    theme_nc() +
    theme(legend.position = "right")
  mm_pdf(paste0("figure7_mp7_cell_cdf_", file_stub, ".pdf"), 96, 70)
  print(cell_cdf_plot)
  dev.off()
}

legend_text <- c(
  "Figure 7 MP7-myeloid exploratory legends",
  "",
  "Boxplot panels. Sample-level median nearest-neighbor distance from each malignant-cell metaprogram to the indicated myeloid target. Each dot is one pretreatment sample. Boxes show median and interquartile range. Pooled myeloid is defined as neutrophils, Mono_CDC27 and Mono_SLC2A3 cells.",
  "",
  "Cell-level CDF panels. Cumulative distribution of nearest-neighbor distances from malignant-cell metaprogram cells to the indicated myeloid target. A left-shifted curve indicates shorter cell-level distances. These panels are descriptive; manuscript statistics should prioritize the sample-level boxplots."
)
writeLines(legend_text, file.path(legend_dir, "figure7_mp7_myeloid_distance_candidate_legends.txt"))

readme_text <- c(
  "# Figure 7 MP7-Myeloid Distance Candidates",
  "",
  "This exploratory folder uses server-exported nearest-neighbor distance metrics to test whether MP7 malignant cells are spatially closer to neutrophil/monocyte populations.",
  "",
  "## Source Data",
  "",
  "`data/geo_ready/figure_7_mp7_myeloid_distance/`",
  "",
  "Only CSV files from the uploaded tarball were preserved; server-generated PDFs were not copied.",
  "",
  "## Preferred Interpretation",
  "",
  "Use the pooled myeloid analysis as the simple main story and the individual neutrophil / Mono_CDC27 / Mono_SLC2A3 analyses as sensitivity support.",
  "",
  "The sample-level analysis should be prioritized for manuscript statistics because each sample is the biological replicate.",
  "",
  "## Recommended Exploratory Outputs",
  "",
  "`outputs/figure7_mp7_boxplot_pooled_myeloid.pdf`",
  "",
  "`outputs/figure7_mp7_boxplot_neutrophil.pdf`",
  "",
  "`outputs/figure7_mp7_boxplot_mono_cdc27.pdf`",
  "",
  "`outputs/figure7_mp7_boxplot_mono_slc2a3.pdf`",
  "",
  "`outputs/figure7_mp7_cell_cdf_pooled_myeloid.pdf`",
  "",
  "`outputs/figure7_mp7_cell_cdf_neutrophil.pdf`",
  "",
  "`outputs/figure7_mp7_cell_cdf_mono_cdc27.pdf`",
  "",
  "`outputs/figure7_mp7_cell_cdf_mono_slc2a3.pdf`"
)
writeLines(readme_text, file.path(dirname(out_dir), "README.md"))

message("Saved MP7-myeloid candidate PDFs to: ", out_dir)
