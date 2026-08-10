suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

data_path <- file.path(project_root, "data", "geo_ready", "spatial_neighborhoods", "postC_k7_cell_niche_assignments.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

group_order <- c("cCR", "non-cCR")
neighborhood_order <- paste0("neighborhood_", 1:7, "_k7")
outcome_colors <- unlist(palette_config$clinical_outcome)[group_order]

format_p_value <- function(p_value, prefix = "p") {
  if (is.na(p_value)) {
    return(paste0(prefix, " = NA"))
  }
  if (p_value < 0.001) {
    return(sub("\\.?0+$", "", sprintf("%s = %.6f", prefix, p_value)))
  }
  sprintf("%s = %.3f", prefix, p_value)
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

cells <- fread(
  data_path,
  select = c("cell_id", "sampleID", "patientID", "cCR", "sample_timepoint", "cell_area", neighborhood_order)
)
cells <- cells[
  sample_timepoint == "postC" &
    !is.na(cell_area) &
    cell_area > 0
]

for (neighborhood in neighborhood_order) {
  cells[, (neighborhood) := as.numeric(get(neighborhood))]
}

weighted_summary <- cells[
  ,
  {
    weighted_scores <- lapply(neighborhood_order, function(neighborhood) {
      neighborhood_weight <- get(neighborhood)
      keep <- !is.na(neighborhood_weight) & !is.na(cell_area) & cell_area > 0
      if (!any(keep)) {
        return(NA_real_)
      }
      sum(neighborhood_weight[keep] * cell_area[keep]) / sum(cell_area[keep])
    })
    names(weighted_scores) <- neighborhood_order
    as.data.table(weighted_scores)
  },
  by = .(patientID, sampleID, cCR)
]

area_summary <- melt(
  weighted_summary,
  id.vars = c("patientID", "sampleID", "cCR"),
  measure.vars = neighborhood_order,
  variable.name = "nmf_neighborhood",
  value.name = "area_weighted_score"
)
area_summary[, cCR := factor(cCR, levels = group_order)]
area_summary[, nmf_neighborhood := factor(nmf_neighborhood, levels = neighborhood_order)]
area_summary[, neighborhood_label := factor(
  sub("_k7$", "", sub("neighborhood_", "Neighborhood ", as.character(nmf_neighborhood))),
  levels = sub("_k7$", "", sub("neighborhood_", "Neighborhood ", neighborhood_order))
)]

fwrite(area_summary, file.path(out_dir, "postc_k7_neighborhood_area_by_response_summary.csv"))
fwrite(area_summary, file.path(out_dir, "postc_k7_neighborhood_area_weighted_score_by_response_summary.csv"))

stats_list <- list()
for (neighborhood in neighborhood_order) {
  neighborhood_df <- area_summary[nmf_neighborhood == neighborhood]
  wilcox_result <- two_group_wilcox(neighborhood_df$area_weighted_score, neighborhood_df$cCR)
  stats_list[[length(stats_list) + 1]] <- data.table(
    nmf_neighborhood = neighborhood,
    neighborhood_label = sub("_k7$", "", sub("neighborhood_", "Neighborhood ", neighborhood)),
    comparison = "cCR_vs_non-cCR",
    timepoint = "postC",
    n_cCR = sum(neighborhood_df$cCR == "cCR"),
    n_non_cCR = sum(neighborhood_df$cCR == "non-cCR"),
    test = "two-sided Wilcoxon rank-sum test",
    p_value = wilcox_result$p_value,
    p_value_method = wilcox_result$method
  )
}
stats_df <- rbindlist(stats_list)
stats_df[, p_adjusted_bh := p.adjust(p_value, method = "BH")]
fwrite(stats_df, file.path(out_dir, "postc_k7_neighborhood_area_by_response_stats.csv"))
fwrite(stats_df, file.path(out_dir, "postc_k7_neighborhood_area_weighted_score_by_response_stats.csv"))

p_labels <- merge(
  stats_df,
  area_summary[, .(panel_max = max(area_weighted_score, na.rm = TRUE)), by = .(nmf_neighborhood)],
  by = "nmf_neighborhood",
  all.x = TRUE
)
p_labels[, neighborhood_label := factor(neighborhood_label, levels = levels(area_summary$neighborhood_label))]
global_max <- max(area_summary$area_weighted_score, na.rm = TRUE)
label_unit <- max(global_max * 0.045, 0.001)
p_labels[, y_start := panel_max + label_unit * 1.4]
p_labels[, y_end := y_start - label_unit * 0.45]
p_labels[, label_y := y_start + label_unit * 0.95]
p_labels[, label := mapply(
  function(raw_p, adj_p) {
    paste(
      format_p_value(raw_p, "raw p"),
      format_p_value(adj_p, "adj p"),
      sep = "\n"
    )
  },
  p_value,
  p_adjusted_bh,
  USE.NAMES = FALSE
)]
y_upper <- max(c(area_summary$area_weighted_score, p_labels$label_y), na.rm = TRUE) * 1.06
if (!is.finite(y_upper) || y_upper <= 0) {
  y_upper <- 1
}

base_size <- plotting_config$font$standard_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot <- ggplot(area_summary, aes(x = cCR, y = area_weighted_score)) +
  geom_boxplot(
    aes(fill = cCR),
    width = 0.55,
    outlier.shape = NA,
    color = "#333333",
    linewidth = line_width,
    alpha = 0.58
  ) +
  geom_point(
    aes(color = cCR),
    position = position_jitter(width = 0.08, height = 0),
    size = 0.95,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_segment(
    data = p_labels,
    aes(x = 1, xend = 2, y = y_start, yend = y_start),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = 1, xend = 1, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_segment(
    data = p_labels,
    aes(x = 2, xend = 2, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "#D62728",
    linewidth = line_width
  ) +
  geom_text(
    data = p_labels,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    color = "#D62728",
    size = 1.75,
    lineheight = 0.9
  ) +
  facet_wrap(~ neighborhood_label, nrow = 1) +
  scale_fill_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  scale_y_continuous(
    limits = c(0, y_upper),
    breaks = pretty(c(0, y_upper), n = 4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Area-weighted neighborhood score"
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
    legend.key.width = unit(0.32, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = base_size, face = "bold"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    panel.spacing.x = unit(0.16, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.58)))

ggsave(
  filename = file.path(out_dir, "postc_k7_neighborhood_area_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = 180,
  height = 82,
  units = "mm",
  dpi = 300
)
ggsave(
  filename = file.path(out_dir, "postc_k7_neighborhood_area_weighted_score_by_response.pdf"),
  plot = plot,
  device = "pdf",
  width = 180,
  height = 82,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# PostC spatial neighborhood score by response",
    "",
    "Area-weighted neighborhood score was calculated per sample from continuous k=7 NMF weights as sum(neighborhood weight x cell area) divided by total annotated cell area in that sample.",
    "Raw P values are two-sided Wilcoxon rank-sum tests comparing cCR and non-cCR samples at postC. Adjusted P values are Benjamini-Hochberg adjusted across the seven neighborhoods."
  ),
  con = file.path(legend_dir, "legend_draft.md")
)
