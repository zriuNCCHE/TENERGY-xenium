suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
  library(yaml)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

mp2_data_dir <- file.path(project_root, "data", "geo_ready", "figure_7_mp_cd4cxcl13_distance")
mp7_data_dir <- file.path(project_root, "data", "geo_ready", "figure_7_mp7_myeloid_distance")
out_dir <- file.path(project_root, "final", "figure_7", "mp_spatial_proximity", "outputs")
legend_dir <- file.path(project_root, "final", "figure_7", "mp_spatial_proximity", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

mp_levels <- paste0("MP_", 1:7)
mp_palette <- unlist(palette_config$malignant_metaprograms)
mp_palette <- mp_palette[mp_levels]
design_blue <- "#7FAEC8"
other_mp_gray <- "#777777"
mp_point_mp2 <- setNames(rep(other_mp_gray, 7), mp_levels)
mp_point_mp2["MP_2"] <- mp_palette["MP_2"]
mp_point_mp7 <- setNames(rep(other_mp_gray, 7), mp_levels)
mp_point_mp7["MP_7"] <- mp_palette["MP_7"]

fwrite(
  data.table(MP = names(mp_palette), color = unname(mp_palette)),
  file.path(out_dir, "figure7_malignant_metaprogram_palette.csv")
)

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
      legend.key.height = unit(7, "pt"),
      legend.key.width = unit(7, "pt"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
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

make_sample_boxplot <- function(dt, highlight_mp, point_values, title, y_title, global_p, filename) {
  p <- ggplot(dt, aes(MP, median_distance, colour = MP)) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      linewidth = 0.5,
      fill = NA
    ) +
    geom_point(
      size = 0.7,
      alpha = 0.9,
      position = position_jitter(width = 0.08, height = 0),
      show.legend = FALSE
    ) +
    scale_colour_manual(values = point_values, guide = "none") +
    labs(title = title, x = NULL, y = y_title) +
    annotate("text", x = 1, y = Inf, vjust = 1.35, hjust = 0, size = 2.0,
             label = paste0("Friedman ", fmt_p(global_p))) +
    theme_nc()
  mm_pdf(filename, 62, 55)
  print(p)
  dev.off()
}

make_cell_cdf <- function(dt, highlight_mp, title, x_title, filename, legend_position = "right") {
  cdf_dt <- dt[!is.na(distance), .(
    distance = sort(distance),
    cumulative_fraction = seq_len(.N) / .N,
    n_cells = .N
  ), by = MP]
  cdf_dt[, MP := factor(MP, levels = mp_levels)]
  cdf_colours <- c(
    MP_1 = "#9A8FC8",
    MP_2 = "#A7C9DD",
    MP_3 = "#EDB370",
    MP_4 = "#F5D5AC",
    MP_5 = "#86C986",
    MP_6 = "#C5E8BF",
    MP_7 = "#DD7A7D"
  )
  cdf_colours[highlight_mp] <- mp_palette[highlight_mp]
  cdf_linewidth <- setNames(rep(0.45, 7), mp_levels)
  cdf_linewidth[highlight_mp] <- 0.9
  fwrite(cdf_dt, file.path(out_dir, sub("\\.pdf$", ".csv", filename)))

  p <- ggplot(cdf_dt, aes(distance, cumulative_fraction, colour = MP, linewidth = MP)) +
    geom_step(lineend = "butt") +
    scale_colour_manual(values = cdf_colours) +
    scale_linewidth_manual(values = cdf_linewidth, guide = "none") +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(0, quantile(dt$distance, 0.985, na.rm = TRUE))) +
    labs(title = title, x = x_title, y = "Cumulative fraction of cells", colour = "Malignant MP") +
    theme_nc() +
    theme(legend.position = legend_position)
  mm_pdf(filename, 82, 52)
  print(p)
  dev.off()
}

mp2_sample <- fread(file.path(mp2_data_dir, "MP_to_CD4_CXCL13_nearest_distance_sample_summary.csv"))
mp2_sample[, MP := factor(MP, levels = mp_levels)]
mp2_wide <- dcast(mp2_sample, sampleID ~ MP, value.var = "median_distance")
mp2_complete <- mp2_wide[complete.cases(mp2_wide[, ..mp_levels])]
mp2_global_p <- friedman.test(as.matrix(mp2_complete[, ..mp_levels]))$p.value

mp2_cell <- fread(file.path(mp2_data_dir, "MP_to_CD4_CXCL13_nearest_distance_cell_level.csv"))
setnames(mp2_cell, "nearest_CD4_CXCL13_distance", "distance")
mp2_cell[, MP := factor(MP, levels = mp_levels)]

mp7_sample <- fread(file.path(mp7_data_dir, "MP_to_monocyte_neutrophil_nearest_distance_sample_summary.csv"))
setnames(mp7_sample, "query_label", "MP")
mp7_sample <- mp7_sample[target_label == "monocyte_neutrophil"]
mp7_sample[, MP := factor(MP, levels = mp_levels)]
mp7_global <- fread(file.path(mp7_data_dir, "MP_to_monocyte_neutrophil_nearest_distance_global_stats.csv"))
mp7_global_p <- mp7_global[target_label == "monocyte_neutrophil", p_value][1]

mp7_cell <- fread(file.path(mp7_data_dir, "MP_to_monocyte_neutrophil_nearest_distance_cell_level.csv"))
mp7_cell <- mp7_cell[target_label == "monocyte_neutrophil"]
setnames(mp7_cell, "nearest_target_distance", "distance")
mp7_cell[, MP := factor(MP, levels = mp_levels)]

make_sample_boxplot(
  mp2_sample,
  highlight_mp = "MP_2",
  point_values = mp_point_mp2,
  title = "MP2-CD4_CXCL13 proximity",
  y_title = "Median distance to CD4_CXCL13",
  global_p = mp2_global_p,
  filename = "figure7_main_B_mp2_cd4cxcl13_boxplot.pdf"
)

make_sample_boxplot(
  mp7_sample,
  highlight_mp = "MP_7",
  point_values = mp_point_mp7,
  title = "MP7-myeloid proximity",
  y_title = "Median distance to myeloid cells",
  global_p = mp7_global_p,
  filename = "figure7_main_C_mp7_myeloid_boxplot.pdf"
)

make_cell_cdf(
  mp2_cell,
  highlight_mp = "MP_2",
  title = "Cell-level CD4_CXCL13 proximity",
  x_title = "Nearest distance to CD4_CXCL13",
  filename = "figure7_main_D_mp2_cd4cxcl13_cell_cdf.pdf"
)

make_cell_cdf(
  mp7_cell,
  highlight_mp = "MP_7",
  title = "Cell-level myeloid proximity",
  x_title = "Nearest distance to myeloid cells",
  filename = "figure7_main_E_mp7_myeloid_cell_cdf.pdf"
)

legend_text <- c(
  "Figure 7 main proximity panel legends",
  "",
  "MP2-CD4_CXCL13 boxplot. Sample-level median nearest-neighbor distances from malignant-cell metaprogram cells to CD4_CXCL13 T cells. Each dot is one pretreatment sample. Boxes show median and interquartile range. Global P value was calculated by Friedman test across malignant metaprograms.",
  "",
  "MP7-myeloid boxplot. Sample-level median nearest-neighbor distances from malignant-cell metaprogram cells to pooled myeloid cells, defined as neutrophils, Mono_CDC27 cells and Mono_SLC2A3 cells. Each dot is one pretreatment sample. Boxes show median and interquartile range. Global P value was calculated by Friedman test across malignant metaprograms.",
  "",
  "CD4_CXCL13 cell-level CDF. Cumulative distribution of nearest-neighbor distances from malignant-cell metaprogram cells to CD4_CXCL13 T cells. A left-shifted curve indicates shorter cell-level distances.",
  "",
  "Myeloid cell-level CDF. Cumulative distribution of nearest-neighbor distances from malignant-cell metaprogram cells to pooled myeloid cells. A left-shifted curve indicates shorter cell-level distances.",
  "",
  "The sample-level boxplots should be prioritized for manuscript statistics because each sample is the biological replicate. Cell-level CDFs are descriptive visual support."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# Figure 7 MP Spatial Proximity",
  "",
  "Final candidate right-side panels for Figure 7, designed to sit next to the pretreatment co-occurrence heatmap.",
  "",
  "The shared malignant metaprogram palette is saved in `config/palettes.yaml` under `malignant_metaprograms`, and exported as `outputs/figure7_malignant_metaprogram_palette.csv`.",
  "",
  "## Main Panels",
  "",
  "- `outputs/figure7_main_B_mp2_cd4cxcl13_boxplot.pdf`",
  "- `outputs/figure7_main_C_mp7_myeloid_boxplot.pdf`",
  "- `outputs/figure7_main_D_mp2_cd4cxcl13_cell_cdf.pdf`",
  "- `outputs/figure7_main_E_mp7_myeloid_cell_cdf.pdf`"
)
writeLines(readme_text, file.path(dirname(out_dir), "README.md"))

message("Saved Figure 7 main proximity panels to: ", out_dir)
