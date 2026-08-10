suppressPackageStartupMessages({
  library(ggplot2)
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

target_cell_types <- c("Macro_CXCL5", "Mono_CDC27", "Mono_SLC2A3", "Neutrophil")
excluded_from_fibro_immune <- c("A2ML1+ epi", "malignant")
timepoint_order <- c("pre", "postC", "postA")
group_order <- c("cCR", "non-cCR")
timepoint_colors <- unlist(palette_config$timepoints)
timepoint_colors <- timepoint_colors[timepoint_order]

format_p_value <- function(p_value, prefix = "p") {
  if (is.na(p_value)) {
    return(paste0(prefix, " = NA"))
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("%s = %.6f", prefix, p_value)))
  }
  sprintf("%s = %.3f", prefix, p_value)
}

paired_wilcox <- function(value_1, value_2) {
  differences <- value_1 - value_2
  nonzero_abs_differences <- abs(differences[differences != 0])
  exact_possible <- all(differences != 0) &&
    length(nonzero_abs_differences) > 0 &&
    length(unique(nonzero_abs_differences)) == length(nonzero_abs_differences)
  result <- wilcox.test(
    value_1,
    value_2,
    paired = TRUE,
    alternative = "two.sided",
    exact = exact_possible
  )
  list(
    p_value = result$p.value,
    method = ifelse(exact_possible, "exact", "asymptotic")
  )
}

metadata <- read.csv(data_path, check.names = TRUE)
metadata <- metadata[, c("sampleID", "sample_timepoint", "patientID", "cCR", "final_cell_type")]
metadata <- metadata[!metadata$final_cell_type %in% excluded_from_fibro_immune, ]

summary_list <- lapply(target_cell_types, function(target_cell_type) {
  metadata$is_target <- metadata$final_cell_type == target_cell_type
  summary_df <- aggregate(
    is_target ~ patientID + sampleID + sample_timepoint + cCR,
    data = metadata,
    FUN = function(x) c(fibro_immune_cells = length(x), target_cells = sum(x))
  )
  summary_df <- do.call(data.frame, summary_df)
  names(summary_df)[names(summary_df) == "is_target.fibro_immune_cells"] <- "fibro_immune_cells"
  names(summary_df)[names(summary_df) == "is_target.target_cells"] <- "target_cells"
  summary_df$immune_subset <- target_cell_type
  summary_df$target_pct <- summary_df$target_cells / summary_df$fibro_immune_cells * 100
  summary_df
})
summary_df <- do.call(rbind, summary_list)
summary_df$sample_timepoint <- factor(summary_df$sample_timepoint, levels = timepoint_order)
summary_df$cCR <- factor(summary_df$cCR, levels = group_order)
summary_df$immune_subset <- factor(summary_df$immune_subset, levels = target_cell_types)

write.csv(
  summary_df,
  file.path(out_dir, "immune_subset_proportion_summary.csv"),
  row.names = FALSE
)

comparisons <- data.frame(
  comparison = "pre_vs_postC",
  timepoint_1 = "pre",
  timepoint_2 = "postC",
  x_start = 1,
  x_end = 2,
  color = "#D62728",
  stringsAsFactors = FALSE
)

stats_list <- list()
for (target_cell_type in target_cell_types) {
  for (group_name in group_order) {
    group_df <- summary_df[
      summary_df$immune_subset == target_cell_type & summary_df$cCR == group_name,
    ]
    wide_df <- reshape(
      group_df[, c("patientID", "sample_timepoint", "target_pct")],
      idvar = "patientID",
      timevar = "sample_timepoint",
      direction = "wide"
    )
    for (i in seq_len(nrow(comparisons))) {
      comp <- comparisons[i, ]
      value_1 <- wide_df[[paste0("target_pct.", comp$timepoint_1)]]
      value_2 <- wide_df[[paste0("target_pct.", comp$timepoint_2)]]
      wilcox_result <- paired_wilcox(value_1, value_2)
      stats_list[[length(stats_list) + 1]] <- data.frame(
        immune_subset = target_cell_type,
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
}
stats_df <- do.call(rbind, stats_list)
stats_df$p_adjusted_bh <- p.adjust(stats_df$p_value, method = "BH")
write.csv(
  stats_df,
  file.path(out_dir, "immune_subset_proportion_unadjusted_stats.csv"),
  row.names = FALSE
)
write.csv(
  stats_df,
  file.path(out_dir, "immune_subset_proportion_pre_postc_raw_adjusted_stats.csv"),
  row.names = FALSE
)

p_labels <- merge(stats_df, comparisons, by = "comparison")
p_labels$label <- paste(
  vapply(p_labels$p_value, format_p_value, character(1), prefix = "raw p"),
  vapply(p_labels$p_adjusted_bh, format_p_value, character(1), prefix = "adj p"),
  sep = "\n"
)
p_labels$cCR <- factor(p_labels$group, levels = group_order)
p_labels$immune_subset <- factor(p_labels$immune_subset, levels = target_cell_types)

panel_max <- aggregate(
  target_pct ~ cCR + immune_subset,
  data = summary_df,
  FUN = max
)
names(panel_max)[names(panel_max) == "target_pct"] <- "panel_max"
p_labels <- merge(p_labels, panel_max, by = c("cCR", "immune_subset"), all.x = TRUE)
p_labels$y_start <- p_labels$panel_max + pmax(2.0, p_labels$panel_max * 0.08)
p_labels$y_end <- p_labels$y_start - pmax(0.6, p_labels$panel_max * 0.02)
p_labels$label_y <- p_labels$y_start + pmax(2.5, p_labels$panel_max * 0.05)

y_upper <- ceiling((max(c(summary_df$target_pct, p_labels$label_y), na.rm = TRUE) + 2) / 5) * 5
y_breaks <- pretty(c(0, y_upper), n = 4)

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm
part_width <- 105
part_height <- 118
part_width <- plotting_config$part_sizes_mm$immune_four_subset_grid$width

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
    size = 0.95,
    alpha = 0.9
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_end, y = y_start, yend = y_start, group = interaction(comparison, immune_subset, cCR), color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_start, xend = x_start, y = y_start, yend = y_end, group = interaction(comparison, immune_subset, cCR), color = comparison),
    inherit.aes = FALSE,
    linewidth = line_width,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = x_end, xend = x_end, y = y_start, yend = y_end, group = interaction(comparison, immune_subset, cCR), color = comparison),
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
  facet_grid(cCR ~ immune_subset, switch = "x") +
  scale_color_manual(
    values = c(timepoint_colors, pre_vs_postC = "#D62728"),
    breaks = timepoint_order,
    name = "Timepoint"
  ) +
  scale_y_continuous(
    limits = c(-1, y_upper),
    breaks = y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Proportion among fibro-immune cells (%)"
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.title.y = element_text(size = base_size, face = "bold"),
    axis.text.y = element_text(size = plotting_config$font$small_pt),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    legend.key.width = unit(0.30, "cm"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.clip = "off",
    strip.text.x = element_text(size = base_size, face = "bold", angle = 45, hjust = 1, vjust = 1, margin = margin(t = 2, b = 6)),
    strip.text.y = element_text(size = 8, face = "bold", angle = 270),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.15, "cm"),
    panel.spacing.y = unit(0.12, "cm"),
    plot.margin = margin(4, 4, 30, 4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))

ggsave(
  filename = file.path(out_dir, "immune_subset_proportion_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = part_width,
  height = part_height,
  units = "mm",
  dpi = 300
)
