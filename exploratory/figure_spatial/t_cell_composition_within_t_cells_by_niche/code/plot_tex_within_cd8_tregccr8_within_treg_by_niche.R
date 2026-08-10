suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

source_path <- file.path(
  project_root,
  "exploratory", "figure_spatial", "t_cell_neighborhood_abundance", "outputs",
  "postc_t_cell_neighborhood_abundance_by_patient.csv"
)
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

neighborhood_order <- c(
  "Niche_iCAF_CXCL5",
  "Niche_Mono",
  "Niche_Neutro",
  "Niche_Epi",
  "Niche_Macro_CXCL9",
  "Niche_iCAF_CXCL6",
  "Niche_T_cell"
)
cd8_types <- c("CD8_Teff", "CD8_Tex_PDCD1", "CD8_prolif")
treg_types <- c("CD4_Treg_CCR8", "CD4_Treg_FOXP3")
needed_types <- c(cd8_types, treg_types)
composition_order <- c("CD8_Tex_PDCD1 within CD8", "CD4_Treg_CCR8 within Treg")

format_plot_p_value <- function(p_value, prefix = "p") {
  vapply(p_value, function(one_p) {
    if (is.na(one_p)) {
      return(paste0(prefix, " = NA"))
    }
    if (one_p < 1e-4) {
      return(paste0(prefix, " = ", formatC(one_p, format = "e", digits = 2)))
    }
    if (one_p < 0.001) {
      return(sprintf("%s = %.6f", prefix, one_p))
    }
    sprintf("%s = %.3f", prefix, one_p)
  }, character(1))
}

safe_friedman <- function(dt) {
  wide <- dcast(
    dt,
    patientID ~ neighborhood_label,
    value.var = "composition_percent",
    fill = NA_real_
  )
  complete_cols <- c("patientID", neighborhood_order)
  wide <- wide[complete.cases(wide[, ..complete_cols])]
  if (nrow(wide) < 3) {
    return(list(p_value = NA_real_, n_patients = nrow(wide)))
  }
  value_matrix <- as.matrix(wide[, ..neighborhood_order])
  if (all(apply(value_matrix, 2, sd) == 0)) {
    return(list(p_value = NA_real_, n_patients = nrow(wide)))
  }
  list(
    p_value = suppressWarnings(friedman.test(value_matrix)$p.value),
    n_patients = nrow(wide)
  )
}

safe_pairwise_paired_wilcox <- function(dt, feature_name) {
  pairs <- combn(neighborhood_order, 2, simplify = FALSE)
  result_list <- lapply(pairs, function(one_pair) {
    one <- dt[neighborhood_label == one_pair[1], .(patientID, value_1 = composition_percent)]
    two <- dt[neighborhood_label == one_pair[2], .(patientID, value_2 = composition_percent)]
    paired <- merge(one, two, by = "patientID")
    paired <- paired[!is.na(value_1) & !is.na(value_2)]
    if (nrow(paired) < 3 || all(paired$value_1 == paired$value_2)) {
      p_value <- NA_real_
    } else {
      p_value <- suppressWarnings(wilcox.test(
        paired$value_2,
        paired$value_1,
        paired = TRUE,
        exact = FALSE
      )$p.value)
    }
    data.table(
      feature = feature_name,
      neighborhood_1 = one_pair[1],
      neighborhood_2 = one_pair[2],
      n_paired_patients = nrow(paired),
      median_neighborhood_1 = median(paired$value_1, na.rm = TRUE),
      median_neighborhood_2 = median(paired$value_2, na.rm = TRUE),
      median_delta_2_minus_1 = median(paired$value_2 - paired$value_1, na.rm = TRUE),
      p_value = p_value
    )
  })
  result <- rbindlist(result_list)
  result[, p_adjusted_bh := p.adjust(p_value, method = "BH")]
  result[]
}

dt <- fread(source_path)
dt <- dt[feature_type == "Individual T subset" & feature %in% needed_types]
dt[, neighborhood_label := as.character(neighborhood_label)]
dt[, feature := as.character(feature)]

wide_counts <- dt[
  ,
  .(subset_cells = sum(subset_cells), subset_area = sum(subset_area)),
  by = .(patientID, sampleID, cCR, neighborhood_label, feature)
]

composition_dt <- wide_counts[
  ,
  .(
    total_cd8_cells = sum(subset_cells[feature %in% cd8_types]),
    total_cd8_area = sum(subset_area[feature %in% cd8_types]),
    tex_cells = sum(subset_cells[feature == "CD8_Tex_PDCD1"]),
    tex_area = sum(subset_area[feature == "CD8_Tex_PDCD1"]),
    total_treg_cells = sum(subset_cells[feature %in% treg_types]),
    total_treg_area = sum(subset_area[feature %in% treg_types]),
    ccr8_treg_cells = sum(subset_cells[feature == "CD4_Treg_CCR8"]),
    ccr8_treg_area = sum(subset_area[feature == "CD4_Treg_CCR8"])
  ),
  by = .(patientID, sampleID, cCR, neighborhood_label)
]
composition_dt[, tex_percent_within_cd8 := fifelse(
  total_cd8_cells > 0,
  tex_cells / total_cd8_cells * 100,
  NA_real_
)]
composition_dt[, ccr8_percent_within_treg := fifelse(
  total_treg_cells > 0,
  ccr8_treg_cells / total_treg_cells * 100,
  NA_real_
)]
composition_dt[, tex_area_percent_within_cd8 := fifelse(
  total_cd8_area > 0,
  tex_area / total_cd8_area * 100,
  NA_real_
)]
composition_dt[, ccr8_area_percent_within_treg := fifelse(
  total_treg_area > 0,
  ccr8_treg_area / total_treg_area * 100,
  NA_real_
)]

plot_dt <- rbindlist(list(
  composition_dt[
    ,
    .(
      patientID, sampleID, cCR, neighborhood_label,
      feature = "CD8_Tex_PDCD1 within CD8",
      composition_percent = tex_percent_within_cd8,
      area_composition_percent = tex_area_percent_within_cd8,
      denominator_cells = total_cd8_cells
    )
  ],
  composition_dt[
    ,
    .(
      patientID, sampleID, cCR, neighborhood_label,
      feature = "CD4_Treg_CCR8 within Treg",
      composition_percent = ccr8_percent_within_treg,
      area_composition_percent = ccr8_area_percent_within_treg,
      denominator_cells = total_treg_cells
    )
  ]
))
plot_dt[, feature := factor(feature, levels = composition_order)]

summary_dt <- plot_dt[
  ,
  .(
    n_defined_patient_niches = sum(!is.na(composition_percent)),
    n_total_patient_niches = .N,
    median_percent = median(composition_percent, na.rm = TRUE),
    mean_percent = mean(composition_percent, na.rm = TRUE),
    max_percent = max(composition_percent, na.rm = TRUE),
    median_area_percent = median(area_composition_percent, na.rm = TRUE),
    mean_area_percent = mean(area_composition_percent, na.rm = TRUE)
  ),
  by = .(feature, neighborhood_label)
]

global_stats <- rbindlist(lapply(composition_order, function(feature_name) {
  one_result <- safe_friedman(plot_dt[feature == feature_name])
  feature_dt <- plot_dt[feature == feature_name & !is.na(composition_percent)]
  kruskal_p <- if (length(unique(feature_dt$neighborhood_label)) >= 2) {
    suppressWarnings(kruskal.test(composition_percent ~ neighborhood_label, data = feature_dt)$p.value)
  } else {
    NA_real_
  }
  data.table(
    feature = feature_name,
    test = c(
      "paired Friedman test across seven niches",
      "exploratory Kruskal-Wallis test across defined patient-niches"
    ),
    n_complete_patients = one_result$n_patients,
    n_defined_patient_niches = c(NA_integer_, nrow(feature_dt)),
    p_value = c(one_result$p_value, kruskal_p)
  )
}))
global_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = test]

pairwise_stats <- rbindlist(lapply(composition_order, function(feature_name) {
  safe_pairwise_paired_wilcox(plot_dt[feature == feature_name], feature_name)
}))

prefix <- "postc_tex_within_cd8_tregccr8_within_treg"
fwrite(composition_dt, file.path(out_dir, paste0(prefix, "_by_patient_niche_wide.csv")))
fwrite(plot_dt, file.path(out_dir, paste0(prefix, "_by_patient_niche_long.csv")))
fwrite(summary_dt, file.path(out_dir, paste0(prefix, "_by_niche_summary.csv")))
fwrite(global_stats, file.path(out_dir, paste0(prefix, "_global_stats.csv")))
fwrite(pairwise_stats, file.path(out_dir, paste0(prefix, "_pairwise_pvalues.csv")))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

order_dt <- plot_dt[
  ,
  .(median_percent = median(composition_percent, na.rm = TRUE)),
  by = .(feature, neighborhood_label)
][order(feature, median_percent, as.character(neighborhood_label))]
order_dt[, feature_neighborhood := paste(neighborhood_label, feature, sep = "___")]
feature_x_levels <- order_dt$feature_neighborhood
label_map <- setNames(order_dt$neighborhood_label, order_dt$feature_neighborhood)

plot_ordered_dt <- copy(plot_dt)
plot_ordered_dt[, feature_neighborhood := factor(
  paste(neighborhood_label, as.character(feature), sep = "___"),
  levels = feature_x_levels
)]

p_labels <- plot_ordered_dt[
  ,
  .(panel_max = max(composition_percent, na.rm = TRUE)),
  by = .(feature)
]
p_labels <- merge(
  p_labels,
  global_stats[test == "exploratory Kruskal-Wallis test across defined patient-niches"],
  by = "feature",
  all.x = TRUE
)
p_labels[, x_position := 4]
p_labels[, y_position := fifelse(is.finite(panel_max) & panel_max > 0, min(panel_max * 1.08, 105), 5)]
p_labels[, label := paste0(
  "Kruskal ",
  format_plot_p_value(p_value, "p"),
  "\nadj ",
  format_plot_p_value(p_adjusted_bh, "p")
)]
p_labels[, feature := factor(feature, levels = composition_order)]

plot <- ggplot(plot_ordered_dt, aes(x = feature_neighborhood, y = composition_percent)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = NA,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 0.72,
    alpha = 0.72,
    color = "#333333"
  ) +
  geom_text(
    data = p_labels,
    aes(x = x_position, y = y_position, label = label),
    inherit.aes = FALSE,
    size = 1.65,
    lineheight = 0.88
  ) +
  facet_wrap(~ feature, nrow = 1, scales = "free_x") +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(
    limits = c(0, 110),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Subset composition (%)"
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
    panel.spacing = unit(0.28, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  filename = file.path(out_dir, paste0(prefix, "_by_niche.pdf")),
  plot = plot,
  device = "pdf",
  width = 145,
  height = 72,
  units = "mm",
  dpi = 300
)

writeLines(
  c(
    "# CD8_Tex_PDCD1 within CD8 and CD4_Treg_CCR8 within Treg by niche",
    "",
    "Exploratory analysis based on the highest-score spatial niche assignment.",
    "For each patient-niche, CD8_Tex_PDCD1 is expressed as a percentage of all CD8-lineage T cells: CD8_Teff, CD8_Tex_PDCD1 and CD8_prolif.",
    "For each patient-niche, CD4_Treg_CCR8 is expressed as a percentage of all Treg-like cells: CD4_Treg_CCR8 and CD4_Treg_FOXP3.",
    "Patient-niches with zero denominator cells are treated as undefined rather than zero.",
    "X axes are ordered from low to high median composition within each feature.",
    "Strict paired Friedman tests require defined composition values in all seven niches and may be underpowered when many patient-niches contain no denominator cells.",
    "On-panel global P values use exploratory Kruskal-Wallis tests across defined patient-niche values."
  ),
  con = file.path(legend_dir, "legend_tex_within_cd8_tregccr8_within_treg.md")
)
