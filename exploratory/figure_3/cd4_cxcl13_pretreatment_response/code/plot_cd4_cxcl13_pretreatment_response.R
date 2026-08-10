suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

data_path <- file.path(project_root, "data", "raw", "combined_final_all_cell_types_obs.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

group_order <- c("cCR", "non-cCR")
outcome_colors <- unlist(palette_config$clinical_outcome)[group_order]
target_cell_type <- "CD4_CXCL13"
t_cell_types <- c(
  "T/NK",
  "CD8_Teff",
  "CD8_Tex_PDCD1",
  "CD8_prolif",
  "CD4_Treg_CCR8",
  "CD4_Treg_FOXP3",
  "CD4_CXCL13",
  "CD4_prolif"
)

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.3f", p_value)
}

safe_wilcox <- function(values, groups) {
  keep <- !is.na(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(groups[keep])
  if (length(unique(groups)) != 2 || any(table(groups) < 2)) {
    return(NA_real_)
  }
  suppressWarnings(wilcox.test(values ~ groups, exact = FALSE)$p.value)
}

cells <- fread(
  data_path,
  select = c("cell_id", "sampleID", "patientID", "sample_timepoint", "cCR", "final_cell_type2")
)
cells <- cells[sample_timepoint == "pre" & !is.na(final_cell_type2)]

sample_info <- unique(cells[, .(patientID, sampleID, cCR)])
target_counts <- cells[
  final_cell_type2 == target_cell_type,
  .(target_cells = .N),
  by = .(patientID, sampleID)
]
non_epi_non_malig_denominator <- cells[
  !(final_cell_type2 %in% c("A2ML1+ epi", "malignant")),
  .(denominator_cells = .N),
  by = .(patientID, sampleID)
]
t_cell_denominator <- cells[
  final_cell_type2 %in% t_cell_types,
  .(denominator_cells = .N),
  by = .(patientID, sampleID)
]

make_summary <- function(denominator_dt, denominator_label) {
  summary_dt <- merge(sample_info, denominator_dt, by = c("patientID", "sampleID"), all.x = TRUE)
  summary_dt <- merge(summary_dt, target_counts, by = c("patientID", "sampleID"), all.x = TRUE)
  summary_dt[is.na(target_cells), target_cells := 0]
  summary_dt[is.na(denominator_cells), denominator_cells := 0]
  summary_dt[, proportion_percent := fifelse(
    denominator_cells > 0,
    target_cells / denominator_cells * 100,
    NA_real_
  )]
  summary_dt[, denominator := denominator_label]
  summary_dt
}

plot_dt <- rbindlist(list(
  make_summary(non_epi_non_malig_denominator, "All non-epithelial, non-malignant cells"),
  make_summary(t_cell_denominator, "All T cells")
), use.names = TRUE)
plot_dt[, cCR := factor(cCR, levels = group_order)]
plot_dt[, denominator := factor(
  denominator,
  levels = c("All non-epithelial, non-malignant cells", "All T cells")
)]

stats_dt <- plot_dt[
  ,
  .(
    n_cCR = sum(cCR == "cCR"),
    n_non_cCR = sum(cCR == "non-cCR"),
    median_cCR = median(proportion_percent[cCR == "cCR"], na.rm = TRUE),
    median_non_cCR = median(proportion_percent[cCR == "non-cCR"], na.rm = TRUE),
    p_value = safe_wilcox(proportion_percent, cCR)
  ),
  by = denominator
]
stats_dt[, p_label := vapply(p_value, format_p_value, character(1))]

fwrite(plot_dt, file.path(out_dir, "pretreatment_cd4_cxcl13_response_summary.csv"))
fwrite(stats_dt, file.path(out_dir, "pretreatment_cd4_cxcl13_response_stats.csv"))

label_dt <- merge(
  stats_dt,
  plot_dt[, .(label_y = max(proportion_percent, na.rm = TRUE) * 1.14 + 0.02), by = denominator],
  by = "denominator",
  all.x = TRUE
)

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(plot_dt, aes(x = cCR, y = proportion_percent)) +
  geom_boxplot(
    aes(fill = cCR),
    width = 0.55,
    outlier.shape = NA,
    color = "#333333",
    linewidth = line_width,
    alpha = 0.6
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.2,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_dt,
    aes(x = 1.5, y = label_y, label = p_label),
    inherit.aes = FALSE,
    size = 2.0
  ) +
  facet_wrap(~ denominator, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  labs(
    x = NULL,
    y = "CD4_CXCL13 cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.28, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, "pretreatment_cd4_cxcl13_response_exploratory.pdf"),
  plot = plot,
  device = "pdf",
  width = 120,
  height = 78,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# Pretreatment CD4_CXCL13 abundance by response",
    "",
    "Exploratory plot comparing pretreatment CD4_CXCL13 cell abundance between cCR and non-cCR patients.",
    "Each point is one pretreatment sample/patient.",
    "The left panel uses all non-epithelial, non-malignant cells as the denominator; the right panel uses all T cells as the denominator.",
    "P values are two-sided Wilcoxon rank-sum tests."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
