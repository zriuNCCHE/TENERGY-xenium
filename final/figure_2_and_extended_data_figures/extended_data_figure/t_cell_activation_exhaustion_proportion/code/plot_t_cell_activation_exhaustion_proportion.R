suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
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

t_cell_types <- c(
  "CD8_Teff",
  "CD8_Tex_PDCD1",
  "CD4_Treg_FOXP3",
  "CD4_Treg_CCR8",
  "CD8_prolif",
  "CD4_CXCL13",
  "CD4_prolif"
)
excluded_from_fibro_immune <- c("A2ML1+ epi", "malignant")
group_order <- c("cCR", "non-cCR")
timepoint_order <- c("pre", "postC", "postA")
outcome_colors <- unlist(palette_config$clinical_outcome)[group_order]
timepoint_colors <- unlist(palette_config$timepoints)[timepoint_order]

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.4f", p_value)
}

two_group_wilcox <- function(values, groups) {
  keep <- !is.na(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(groups[keep])
  exact_possible <- length(unique(values)) == length(values)
  result <- suppressWarnings(wilcox.test(values ~ groups, exact = exact_possible))
  list(
    p_value = result$p.value,
    method = ifelse(exact_possible, "exact", "asymptotic due to ties")
  )
}

paired_wilcox <- function(value_1, value_2) {
  keep <- !is.na(value_1) & !is.na(value_2)
  value_1 <- value_1[keep]
  value_2 <- value_2[keep]
  if (length(value_1) < 2) {
    return(list(p_value = NA_real_, method = "not tested: fewer than two paired patients"))
  }

  differences <- value_1 - value_2
  if (all(differences == 0)) {
    return(list(p_value = 1, method = "all paired differences are zero"))
  }

  nonzero_abs_differences <- abs(differences[differences != 0])
  exact_possible <- all(differences != 0) &&
    length(nonzero_abs_differences) > 0 &&
    length(unique(nonzero_abs_differences)) == length(nonzero_abs_differences)
  result <- suppressWarnings(wilcox.test(
    value_1,
    value_2,
    paired = TRUE,
    alternative = "two.sided",
    exact = exact_possible
  ))
  list(
    p_value = result$p.value,
    method = ifelse(exact_possible, "exact", "asymptotic")
  )
}

summarise_target <- function(metadata, target_cell_type2) {
  metadata$is_target <- metadata$final_cell_type2 == target_cell_type2
  summary_df <- aggregate(
    is_target ~ patientID + sampleID + sample_timepoint + cCR,
    data = metadata,
    FUN = function(x) c(t_cells = length(x), target_cells = sum(x))
  )
  summary_df <- do.call(data.frame, summary_df)
  names(summary_df)[names(summary_df) == "is_target.t_cells"] <- "denominator_cells"
  names(summary_df)[names(summary_df) == "is_target.target_cells"] <- "target_cells"
  summary_df$target_cell_type <- target_cell_type2
  summary_df$target_pct <- summary_df$target_cells / summary_df$denominator_cells * 100
  summary_df$sample_timepoint <- factor(summary_df$sample_timepoint, levels = timepoint_order)
  summary_df$cCR <- factor(summary_df$cCR, levels = group_order)
  summary_df
}

metadata <- read.csv(data_path, check.names = TRUE)
metadata <- metadata[, c("sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type", "final_cell_type2")]
metadata_all <- metadata
metadata <- metadata_all[metadata_all$final_cell_type2 %in% t_cell_types, ]
metadata_fibro_immune <- metadata_all[!metadata_all$final_cell_type %in% excluded_from_fibro_immune, ]

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

## Panel 1: postC CD4_Treg_CCR8 by clinical response.
treg_summary <- summarise_target(metadata, "CD4_Treg_CCR8")
treg_postc <- treg_summary[treg_summary$sample_timepoint == "postC", ]

write.csv(
  treg_postc,
  file.path(out_dir, "cd4_treg_ccr8_postc_by_response_summary.csv"),
  row.names = FALSE
)

treg_test <- two_group_wilcox(treg_postc$target_pct, treg_postc$cCR)
treg_stats <- data.frame(
  target_cell_type = "CD4_Treg_CCR8",
  denominator = "T cells",
  comparison = "cCR_vs_non-cCR",
  timepoint = "postC",
  n_cCR = sum(treg_postc$cCR == "cCR"),
  n_non_cCR = sum(treg_postc$cCR == "non-cCR"),
  test = "two-sided Wilcoxon rank-sum test",
  p_value = treg_test$p_value,
  p_value_method = treg_test$method,
  stringsAsFactors = FALSE
)
write.csv(
  treg_stats,
  file.path(out_dir, "cd4_treg_ccr8_postc_by_response_unadjusted_stats.csv"),
  row.names = FALSE
)

treg_y_upper <- ceiling((max(treg_postc$target_pct, na.rm = TRUE) * 1.20 + 0.3) / 1) * 1
treg_y_upper <- max(treg_y_upper, 4)
treg_label <- data.frame(
  x_start = 1,
  x_end = 2,
  y_start = treg_y_upper * 0.84,
  y_end = treg_y_upper * 0.80,
  label_y = treg_y_upper * 0.89,
  label = format_p_value(treg_stats$p_value)
)

treg_plot <- ggplot(treg_postc, aes(x = cCR, y = target_pct)) +
  geom_boxplot(
    aes(fill = cCR),
    width = 0.36,
    outlier.shape = NA,
    color = "#333333",
    linewidth = line_width,
    alpha = 0.55
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_segment(
    data = treg_label,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = treg_label,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = treg_label,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_text(
    data = treg_label,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    color = "#D62728",
    size = 2.0
  ) +
  facet_wrap(~ target_cell_type, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(0, treg_y_upper),
    breaks = pretty(c(0, treg_y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "CCR8+ Treg after chemo-radio",
    x = NULL,
    y = "Proportion of T cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    legend.key.width = unit(0.32, "cm"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.clip = "off",
    strip.text.x = element_text(size = base_size, face = "bold", margin = margin(t = 3)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    plot.margin = margin(4, 4, 8, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.55)))

ggsave(
  filename = file.path(out_dir, "cd4_treg_ccr8_postc_by_response.pdf"),
  plot = treg_plot,
  device = "pdf",
  width = 88,
  height = 72,
  units = "mm",
  dpi = 300
)

## Panel 2: CD8_Tex_PDCD1/PDCD1-positive cell state over time by clinical response.
pdcd1_summary <- summarise_target(metadata, "CD8_Tex_PDCD1")
write.csv(
  pdcd1_summary,
  file.path(out_dir, "pdcd1_positive_timepoint_by_response_summary.csv"),
  row.names = FALSE
)

comparisons <- data.frame(
  comparison = c("pre_vs_postC", "postC_vs_postA"),
  timepoint_1 = c("pre", "postC"),
  timepoint_2 = c("postC", "postA"),
  x_start = c(1, 2),
  x_end = c(2, 3),
  color = c("#D62728", "#000000"),
  stringsAsFactors = FALSE
)

stats_list <- list()
for (group_name in group_order) {
  group_df <- pdcd1_summary[pdcd1_summary$cCR == group_name, ]
  wide_df <- reshape(
    group_df[, c("patientID", "sample_timepoint", "target_pct")],
    idvar = "patientID",
    timevar = "sample_timepoint",
    direction = "wide"
  )
  for (i in seq_len(nrow(comparisons))) {
    comp <- comparisons[i, ]
    wilcox_result <- paired_wilcox(
      wide_df[[paste0("target_pct.", comp$timepoint_1)]],
      wide_df[[paste0("target_pct.", comp$timepoint_2)]]
    )
    stats_list[[length(stats_list) + 1]] <- data.frame(
      target_cell_type = "CD8_Tex_PDCD1",
      denominator = "T cells",
      group = group_name,
      comparison = comp$comparison,
      timepoint_1 = comp$timepoint_1,
      timepoint_2 = comp$timepoint_2,
      n_patients = nrow(wide_df),
      test = "paired two-sided Wilcoxon signed-rank test",
      p_value = wilcox_result$p_value,
      p_value_method = wilcox_result$method,
      stringsAsFactors = FALSE
    )
  }
}
pdcd1_stats <- do.call(rbind, stats_list)
write.csv(
  pdcd1_stats,
  file.path(out_dir, "pdcd1_positive_timepoint_by_response_unadjusted_stats.csv"),
  row.names = FALSE
)

p_labels <- merge(pdcd1_stats, comparisons, by = "comparison")
p_labels$label <- vapply(p_labels$p_value, format_p_value, character(1))
p_labels$cCR <- factor(p_labels$group, levels = group_order)

panel_max <- aggregate(target_pct ~ cCR, data = pdcd1_summary, FUN = max)
names(panel_max)[names(panel_max) == "target_pct"] <- "panel_max"
p_labels <- merge(p_labels, panel_max, by = "cCR", all.x = TRUE)
p_labels$y_start <- p_labels$panel_max + pmax(0.9, p_labels$panel_max * 0.10) +
  ifelse(p_labels$comparison == "postC_vs_postA", pmax(0.9, p_labels$panel_max * 0.12), 0)
p_labels$y_end <- p_labels$y_start - pmax(0.25, p_labels$panel_max * 0.025)
p_labels$label_y <- p_labels$y_start + pmax(0.35, p_labels$panel_max * 0.04)

pdcd1_y_upper <- ceiling((max(c(pdcd1_summary$target_pct, p_labels$label_y), na.rm = TRUE) + 0.4) / 1) * 1
pdcd1_y_upper <- max(pdcd1_y_upper, 5)

pdcd1_plot <- ggplot(pdcd1_summary, aes(x = sample_timepoint, y = target_pct)) +
  geom_line(
    aes(group = patientID),
    color = "#BDBDBD",
    linewidth = line_width,
    alpha = 0.55
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(color = sample_timepoint),
    size = 1.0,
    alpha = 0.95
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start, color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end, color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end, color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_text(
    data = p_labels,
    aes(x = (x_start + x_end) / 2, y = label_y, label = label, color = comparison),
    inherit.aes = FALSE,
    size = 2.0,
    show.legend = FALSE
  ) +
  facet_grid(. ~ cCR) +
  scale_color_manual(
    values = c(timepoint_colors, pre_vs_postC = "#D62728", postC_vs_postA = "#000000"),
    breaks = timepoint_order,
    name = "Timepoint"
  ) +
  scale_y_continuous(
    limits = c(0, pdcd1_y_upper),
    breaks = pretty(c(0, pdcd1_y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "PDCD1-positive CD8 T cells after chemo-radio",
    x = "CD8_Tex_PDCD1",
    y = "Proportion of T cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.title.x = element_text(size = base_size, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_text(size = base_size, face = "bold", angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    legend.key.width = unit(0.32, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8, face = "bold", margin = margin(t = 4)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.65, "cm"),
    plot.margin = margin(4, 4, 8, 4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))

ggsave(
  filename = file.path(out_dir, "pdcd1_positive_timepoint_by_response.pdf"),
  plot = pdcd1_plot,
  device = "pdf",
  width = 88,
  height = 82,
  units = "mm",
  dpi = 300
)

combined_plot <- treg_plot / pdcd1_plot +
  plot_layout(heights = c(0.9, 1.15))

ggsave(
  filename = file.path(out_dir, "t_cell_activation_exhaustion_proportion.pdf"),
  plot = combined_plot,
  device = "pdf",
  width = 95,
  height = 150,
  units = "mm",
  dpi = 300
)

## Denominator check: exclude A2ML1+ epithelial and malignant cells.
## This asks whether the T-cell-state signal is consistent when scaled to the
## fibro-immune compartment rather than only the T-cell compartment.
treg_fi_summary <- summarise_target(metadata_fibro_immune, "CD4_Treg_CCR8")
treg_fi_postc <- treg_fi_summary[treg_fi_summary$sample_timepoint == "postC", ]
write.csv(
  treg_fi_postc,
  file.path(out_dir, "cd4_treg_ccr8_postc_by_response_fibro_immune_denominator_summary.csv"),
  row.names = FALSE
)

treg_fi_test <- two_group_wilcox(treg_fi_postc$target_pct, treg_fi_postc$cCR)
treg_fi_stats <- data.frame(
  target_cell_type = "CD4_Treg_CCR8",
  denominator = "fibro-immune cells",
  comparison = "cCR_vs_non-cCR",
  timepoint = "postC",
  n_cCR = sum(treg_fi_postc$cCR == "cCR"),
  n_non_cCR = sum(treg_fi_postc$cCR == "non-cCR"),
  test = "two-sided Wilcoxon rank-sum test",
  p_value = treg_fi_test$p_value,
  p_value_method = treg_fi_test$method,
  stringsAsFactors = FALSE
)
write.csv(
  treg_fi_stats,
  file.path(out_dir, "cd4_treg_ccr8_postc_by_response_fibro_immune_denominator_unadjusted_stats.csv"),
  row.names = FALSE
)

treg_fi_y_upper <- ceiling((max(treg_fi_postc$target_pct, na.rm = TRUE) * 1.22 + 0.2) / 1) * 1
treg_fi_y_upper <- max(treg_fi_y_upper, 4)
treg_fi_label <- data.frame(
  x_start = 1,
  x_end = 2,
  y_start = treg_fi_y_upper * 0.84,
  y_end = treg_fi_y_upper * 0.80,
  label_y = treg_fi_y_upper * 0.89,
  label = format_p_value(treg_fi_stats$p_value)
)

treg_fi_plot <- ggplot(treg_fi_postc, aes(x = cCR, y = target_pct)) +
  geom_boxplot(
    aes(fill = cCR),
    width = 0.36,
    outlier.shape = NA,
    color = "#333333",
    linewidth = line_width,
    alpha = 0.55
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_segment(
    data = treg_fi_label,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = treg_fi_label,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = treg_fi_label,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_text(
    data = treg_fi_label,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    color = "#D62728",
    size = 2.0
  ) +
  facet_wrap(~ target_cell_type, nrow = 1, strip.position = "bottom") +
  scale_fill_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(0, treg_fi_y_upper),
    breaks = pretty(c(0, treg_fi_y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "CCR8+ Treg after chemo-radio",
    x = NULL,
    y = "Proportion among fibro-immune cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    legend.key.width = unit(0.32, "cm"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.clip = "off",
    strip.text.x = element_text(size = base_size, face = "bold", margin = margin(t = 3)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    plot.margin = margin(4, 4, 8, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.55)))

ggsave(
  filename = file.path(out_dir, "cd4_treg_ccr8_postc_by_response_fibro_immune_denominator.pdf"),
  plot = treg_fi_plot,
  device = "pdf",
  width = 88,
  height = 72,
  units = "mm",
  dpi = 300
)

pdcd1_fi_summary <- summarise_target(metadata_fibro_immune, "CD8_Tex_PDCD1")
write.csv(
  pdcd1_fi_summary,
  file.path(out_dir, "pdcd1_positive_timepoint_by_response_fibro_immune_denominator_summary.csv"),
  row.names = FALSE
)

fi_stats_list <- list()
for (group_name in group_order) {
  group_df <- pdcd1_fi_summary[pdcd1_fi_summary$cCR == group_name, ]
  wide_df <- reshape(
    group_df[, c("patientID", "sample_timepoint", "target_pct")],
    idvar = "patientID",
    timevar = "sample_timepoint",
    direction = "wide"
  )
  for (i in seq_len(nrow(comparisons))) {
    comp <- comparisons[i, ]
    wilcox_result <- paired_wilcox(
      wide_df[[paste0("target_pct.", comp$timepoint_1)]],
      wide_df[[paste0("target_pct.", comp$timepoint_2)]]
    )
    fi_stats_list[[length(fi_stats_list) + 1]] <- data.frame(
      target_cell_type = "CD8_Tex_PDCD1",
      denominator = "fibro-immune cells",
      group = group_name,
      comparison = comp$comparison,
      timepoint_1 = comp$timepoint_1,
      timepoint_2 = comp$timepoint_2,
      n_patients = nrow(wide_df),
      test = "paired two-sided Wilcoxon signed-rank test",
      p_value = wilcox_result$p_value,
      p_value_method = wilcox_result$method,
      stringsAsFactors = FALSE
    )
  }
}
pdcd1_fi_stats <- do.call(rbind, fi_stats_list)
write.csv(
  pdcd1_fi_stats,
  file.path(out_dir, "pdcd1_positive_timepoint_by_response_fibro_immune_denominator_unadjusted_stats.csv"),
  row.names = FALSE
)

fi_p_labels <- merge(pdcd1_fi_stats, comparisons, by = "comparison")
fi_p_labels$label <- vapply(fi_p_labels$p_value, format_p_value, character(1))
fi_p_labels$cCR <- factor(fi_p_labels$group, levels = group_order)

fi_panel_max <- aggregate(target_pct ~ cCR, data = pdcd1_fi_summary, FUN = max)
names(fi_panel_max)[names(fi_panel_max) == "target_pct"] <- "panel_max"
fi_p_labels <- merge(fi_p_labels, fi_panel_max, by = "cCR", all.x = TRUE)
fi_p_labels$y_start <- fi_p_labels$panel_max + pmax(0.25, fi_p_labels$panel_max * 0.10) +
  ifelse(fi_p_labels$comparison == "postC_vs_postA", pmax(0.25, fi_p_labels$panel_max * 0.12), 0)
fi_p_labels$y_end <- fi_p_labels$y_start - pmax(0.08, fi_p_labels$panel_max * 0.025)
fi_p_labels$label_y <- fi_p_labels$y_start + pmax(0.12, fi_p_labels$panel_max * 0.04)

pdcd1_fi_y_upper <- ceiling((max(c(pdcd1_fi_summary$target_pct, fi_p_labels$label_y), na.rm = TRUE) + 0.2) / 1) * 1
pdcd1_fi_y_upper <- max(pdcd1_fi_y_upper, 5)

pdcd1_fi_plot <- ggplot(pdcd1_fi_summary, aes(x = sample_timepoint, y = target_pct)) +
  geom_line(
    aes(group = patientID),
    color = "#BDBDBD",
    linewidth = line_width,
    alpha = 0.55
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(color = sample_timepoint),
    size = 1.0,
    alpha = 0.95
  ) +
  geom_segment(
    data = fi_p_labels,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start, color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_segment(
    data = fi_p_labels,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end, color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_segment(
    data = fi_p_labels,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end, color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_text(
    data = fi_p_labels,
    aes(x = (x_start + x_end) / 2, y = label_y, label = label, color = comparison),
    inherit.aes = FALSE,
    size = 2.0,
    show.legend = FALSE
  ) +
  facet_grid(. ~ cCR) +
  scale_color_manual(
    values = c(timepoint_colors, pre_vs_postC = "#D62728", postC_vs_postA = "#000000"),
    breaks = timepoint_order,
    name = "Timepoint"
  ) +
  scale_y_continuous(
    limits = c(0, pdcd1_fi_y_upper),
    breaks = pretty(c(0, pdcd1_fi_y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "PDCD1-positive CD8 T cells after chemo-radio",
    x = "CD8_Tex_PDCD1",
    y = "Proportion among fibro-immune cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = plotting_config$font$title_pt, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.title.x = element_text(size = base_size, face = "bold", margin = margin(t = 6)),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_text(size = base_size, face = "bold", angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    legend.key.width = unit(0.32, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8, face = "bold", margin = margin(t = 4)),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.65, "cm"),
    plot.margin = margin(4, 4, 8, 4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))

ggsave(
  filename = file.path(out_dir, "pdcd1_positive_timepoint_by_response_fibro_immune_denominator.pdf"),
  plot = pdcd1_fi_plot,
  device = "pdf",
  width = 88,
  height = 82,
  units = "mm",
  dpi = 300
)

combined_fi_plot <- treg_fi_plot / pdcd1_fi_plot +
  plot_layout(heights = c(0.9, 1.15))

ggsave(
  filename = file.path(out_dir, "t_cell_activation_exhaustion_proportion_fibro_immune_denominator.pdf"),
  plot = combined_fi_plot,
  device = "pdf",
  width = 95,
  height = 150,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# T cell activation and exhaustion proportions",
    "",
    "Extended Data Figure 3 draft part. CD4_Treg_CCR8 and CD8_Tex_PDCD1 proportions were calculated among T cells.",
    "Additional denominator-check PDFs calculate the same target populations among fibro-immune cells, defined by excluding A2ML1+ epithelial cells and malignant cells.",
    "",
    "Boxes show median and interquartile range. Points denote patient samples; grey lines connect paired patient samples in the longitudinal PDCD1-positive plot. P values are unadjusted Wilcoxon tests."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
