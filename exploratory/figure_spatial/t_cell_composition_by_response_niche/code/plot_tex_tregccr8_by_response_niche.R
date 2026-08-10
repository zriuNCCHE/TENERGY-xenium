suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

source_long <- file.path(
  project_root,
  "final", "extended_figure_5", "t_cell_composition_by_niche", "outputs",
  "extended_figure5_tex_within_cd8_tregccr8_within_treg_by_patient_niche_long.csv"
)
source_wide <- file.path(
  project_root,
  "final", "extended_figure_5", "t_cell_composition_by_niche", "outputs",
  "extended_figure5_tex_within_cd8_tregccr8_within_treg_by_patient_niche_wide.csv"
)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

neighborhood_order <- c(
  "Niche_iCAF_CXCL5",
  "Niche_Mono",
  "Niche_Neutro",
  "Niche_Epi",
  "Niche_Macro_CXCL9",
  "Niche_iCAF_CXCL6",
  "Niche_T_cell"
)
feature_order <- c("CD8_Tex_PDCD1 within CD8", "CD4_Treg_CCR8 within Treg")
outcome_order <- c("cCR", "non-cCR")

format_p <- function(p_value) {
  vapply(p_value, function(one_p) {
    if (is.na(one_p)) {
      return("p = NA")
    }
    if (one_p < 1e-4) {
      return(paste0("p = ", formatC(one_p, format = "e", digits = 2)))
    }
    if (one_p < 0.001) {
      return(sprintf("p = %.6f", one_p))
    }
    sprintf("p = %.3f", one_p)
  }, character(1))
}

dt <- fread(source_long)
wide_dt <- fread(source_wide)
dt[, neighborhood_label := factor(neighborhood_label, levels = neighborhood_order)]
dt[, feature := factor(feature, levels = feature_order)]
dt[, cCR := factor(cCR, levels = outcome_order)]

summary_dt <- dt[
  ,
  .(
    n_defined_patient_niches = sum(!is.na(composition_percent)),
    n_total_patient_niches = .N,
    n_100_percent = sum(composition_percent == 100, na.rm = TRUE),
    median_percent = median(composition_percent, na.rm = TRUE),
    mean_percent = mean(composition_percent, na.rm = TRUE),
    median_denominator_cells = as.numeric(median(denominator_cells[!is.na(composition_percent)], na.rm = TRUE)),
    min_denominator_cells_at_100 = as.numeric(suppressWarnings(min(denominator_cells[composition_percent == 100], na.rm = TRUE))),
    median_denominator_cells_at_100 = as.numeric(median(denominator_cells[composition_percent == 100], na.rm = TRUE))
  ),
  by = .(feature, cCR, neighborhood_label)
]
summary_dt[is.infinite(min_denominator_cells_at_100), min_denominator_cells_at_100 := NA_real_]

stats_dt <- dt[
  !is.na(composition_percent),
  {
    if (length(unique(cCR)) < 2 || .N < 3) {
      p_value <- NA_real_
    } else {
      p_value <- suppressWarnings(wilcox.test(
        composition_percent ~ cCR,
        exact = FALSE
      )$p.value)
    }
    .(
      n_defined = .N,
      n_cCR = sum(cCR == "cCR"),
      n_non_cCR = sum(cCR == "non-cCR"),
      median_cCR = median(composition_percent[cCR == "cCR"], na.rm = TRUE),
      median_non_cCR = median(composition_percent[cCR == "non-cCR"], na.rm = TRUE),
      p_value = p_value
    )
  },
  by = .(feature, neighborhood_label)
]
stats_dt[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature]

hundred_dt <- dt[
  composition_percent == 100,
  .(patientID, sampleID, cCR, neighborhood_label, feature, denominator_cells)
][order(feature, cCR, denominator_cells)]

fwrite(summary_dt, file.path(out_dir, "postc_tex_tregccr8_by_response_niche_summary.csv"))
fwrite(stats_dt, file.path(out_dir, "postc_tex_tregccr8_cCR_vs_non_cCR_by_niche_stats.csv"))
fwrite(hundred_dt, file.path(out_dir, "postc_tex_tregccr8_100_percent_patient_niches.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
outcome_colors <- unlist(palette_config$clinical_outcome)

order_dt <- dt[
  ,
  .(median_percent = median(composition_percent, na.rm = TRUE)),
  by = .(feature, neighborhood_label)
][order(feature, median_percent, as.character(neighborhood_label))]
order_dt[, feature_neighborhood := paste(as.character(feature), as.character(neighborhood_label), sep = "___")]
feature_neighborhood_levels <- order_dt$feature_neighborhood
feature_neighborhood_labels <- setNames(
  as.character(order_dt$neighborhood_label),
  order_dt$feature_neighborhood
)

plot_dt <- copy(dt)
plot_dt[, feature_neighborhood := factor(
  paste(as.character(feature), as.character(neighborhood_label), sep = "___"),
  levels = feature_neighborhood_levels
)]

plot <- ggplot(plot_dt, aes(x = feature_neighborhood, y = composition_percent, color = cCR)) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = NA,
    linewidth = line_width,
    position = position_dodge(width = 0.68)
  ) +
  geom_point(
    aes(size = pmin(denominator_cells, 200)),
    position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.68),
    alpha = 0.78,
    stroke = 0
  ) +
  facet_wrap(~ feature + cCR, nrow = 2, scales = "free_x") +
  scale_color_manual(values = outcome_colors, drop = FALSE) +
  scale_size_continuous(range = c(0.55, 1.35), guide = "none") +
  scale_x_discrete(labels = feature_neighborhood_labels) +
  scale_y_continuous(limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "Subset composition (%)",
    color = "Clinical outcome"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = small_size),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing = unit(0.35, "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(out_dir, "postc_tex_tregccr8_by_response_niche_boxplots.pdf"),
  plot,
  device = "pdf",
  width = 180,
  height = 115,
  units = "mm",
  dpi = 300
)

denominator_plot_dt <- dt[
  composition_percent == 100,
  .(
    n_100_patient_niches = .N,
    median_denominator_cells = median(denominator_cells),
    max_denominator_cells = max(denominator_cells)
  ),
  by = .(feature, cCR)
]

denominator_plot <- ggplot(
  hundred_dt,
  aes(x = cCR, y = denominator_cells, color = cCR)
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    fill = NA,
    linewidth = line_width
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 0.9,
    alpha = 0.78
  ) +
  facet_wrap(~ feature, nrow = 1, scales = "free_y") +
  scale_color_manual(values = outcome_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Denominator cells among 100% points",
    color = "Clinical outcome"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    legend.position = "bottom",
    legend.title = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(out_dir, "postc_tex_tregccr8_100_percent_denominator_diagnostic.pdf"),
  denominator_plot,
  device = "pdf",
  width = 120,
  height = 70,
  units = "mm",
  dpi = 300
)

fwrite(denominator_plot_dt, file.path(out_dir, "postc_tex_tregccr8_100_percent_denominator_summary.csv"))

writeLines(
  c(
    "# Exploratory cCR/non-cCR split for within-lineage T-cell composition by niche",
    "",
    "This exploratory plot asks whether 100% composition values occur preferentially in cCR or non-cCR patient-niches and whether they are driven by small denominators.",
    "CD8_Tex_PDCD1 is measured as a percentage of CD8_Teff, CD8_Tex_PDCD1 and CD8_prolif cells within each patient-niche.",
    "CD4_Treg_CCR8 is measured as a percentage of CD4_Treg_CCR8 and CD4_Treg_FOXP3 cells within each patient-niche.",
    "Point size is proportional to denominator cell count, capped at 200 cells for readability.",
    "Patient-niches with zero denominator cells are undefined and are not plotted."
  ),
  file.path(legend_dir, "legend_draft.md")
)
