suppressPackageStartupMessages({
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

target_cell_type <- "Macro_CXCL9"
excluded_from_fibro_immune <- c("A2ML1+ epi", "malignant")
timepoint_order <- c("pre", "postC", "postA")
group_order <- c("cCR", "non-cCR")
timepoint_colors <- unlist(palette_config$timepoints)[timepoint_order]

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("p = NA")
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("p = %.6f", p_value)))
  }
  sprintf("p = %.3f", p_value)
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

metadata <- read.csv(data_path, check.names = TRUE)
metadata <- metadata[
  !metadata$final_cell_type %in% excluded_from_fibro_immune,
  c("sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type")
]
metadata$is_target <- metadata$final_cell_type == target_cell_type

summary_df <- aggregate(
  is_target ~ patientID + sampleID + sample_timepoint + cCR,
  data = metadata,
  FUN = function(x) c(fibro_immune_cells = length(x), target_cells = sum(x))
)
summary_df <- do.call(data.frame, summary_df)
names(summary_df)[names(summary_df) == "is_target.fibro_immune_cells"] <- "fibro_immune_cells"
names(summary_df)[names(summary_df) == "is_target.target_cells"] <- "target_cells"
summary_df$target_cell_type <- target_cell_type
summary_df$target_pct <- summary_df$target_cells / summary_df$fibro_immune_cells * 100
summary_df$sample_timepoint <- factor(summary_df$sample_timepoint, levels = timepoint_order)
summary_df$cCR <- factor(summary_df$cCR, levels = group_order)

write.csv(
  summary_df,
  file.path(out_dir, "macro_cxcl9_timepoint_by_response_summary.csv"),
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
  group_df <- summary_df[summary_df$cCR == group_name, ]
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
      target_cell_type = target_cell_type,
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
stats_df <- do.call(rbind, stats_list)
write.csv(
  stats_df,
  file.path(out_dir, "macro_cxcl9_timepoint_by_response_unadjusted_stats.csv"),
  row.names = FALSE
)

p_labels <- merge(stats_df, comparisons, by = "comparison")
p_labels$label <- vapply(p_labels$p_value, format_p_value, character(1))
p_labels$cCR <- factor(p_labels$group, levels = group_order)

panel_max <- aggregate(target_pct ~ cCR, data = summary_df, FUN = max)
names(panel_max)[names(panel_max) == "target_pct"] <- "panel_max"
p_labels <- merge(p_labels, panel_max, by = "cCR", all.x = TRUE)
p_labels$y_start <- p_labels$panel_max + pmax(1.5, p_labels$panel_max * 0.09) +
  ifelse(p_labels$comparison == "postC_vs_postA", pmax(1.6, p_labels$panel_max * 0.11), 0)
p_labels$y_end <- p_labels$y_start - pmax(0.4, p_labels$panel_max * 0.02)
p_labels$label_y <- p_labels$y_start + pmax(0.6, p_labels$panel_max * 0.035)

y_upper <- ceiling((max(c(summary_df$target_pct, p_labels$label_y), na.rm = TRUE) + 1) / 5) * 5
y_upper <- max(y_upper, 10)

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(summary_df, aes(x = sample_timepoint, y = target_pct)) +
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
    limits = c(0, y_upper),
    breaks = pretty(c(0, y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Macro_CXCL9 across treatment",
    x = target_cell_type,
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
  filename = file.path(out_dir, "macro_cxcl9_timepoint_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = 88,
  height = 82,
  units = "mm",
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "macro_cxcl9_timepoint_by_response.png"),
  plot = plot,
  device = "png",
  width = 88,
  height = 82,
  units = "mm",
  dpi = 300,
  bg = "white"
)

writeLines(
  c(
    "# Macro_CXCL9 longitudinal by response",
    "",
    "Exploratory Extended Figure 3-style plot. Macro_CXCL9 proportions were calculated among fibro-immune cells across pre, postC and postA timepoints, stratified by cCR and non-cCR.",
    "",
    "Boxes show median and interquartile range. Points denote patient samples; grey lines connect paired samples. P values are unadjusted paired two-sided Wilcoxon signed-rank tests."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
