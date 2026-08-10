suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
  library(survival)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

spatial_path <- file.path(project_root, "data", "geo_ready", "spatial_neighborhoods", "postC_k7_cell_niche_assignments.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
legend_dir <- normalizePath(file.path(dirname(script_path), "..", "legends"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

group_order <- c("cCR", "non-cCR")
outcome_colors <- unlist(palette_config$clinical_outcome)[group_order]
event_colors <- c("No PD recorded" = "#7BB6A4", "PD recorded" = "#D95F02")
feature_group_colors <- c(
  "Spatial neighborhood" = "#4C78A8",
  "T immune cell type" = "#756BB1",
  "Fibroblast/endothelial cell type" = "#E6550D"
)

neighborhood_order <- paste0("neighborhood_", 1:7, "_k7")
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
fibro_cell_types <- c(
  "iCAF_CXCL5",
  "iCAF_CXCL6",
  "iCAF_CXCL12",
  "iCAF_TNC",
  "myCAF_MMP11",
  "myCAF_DUX4",
  "myCAF_TNFRSF21",
  "myCAF_COL10A1",
  "vCAF",
  "endo_PLVAP"
)

format_p_value <- function(p_value) {
  vapply(p_value, function(one_p) {
    if (is.na(one_p)) {
      return("NA")
    }
    if (one_p < 0.001) {
      return(sub("\\.?0+$", "", sprintf("%.6f", one_p)))
    }
    sprintf("%.3f", one_p)
  }, character(1))
}

safe_wilcox <- function(value, group) {
  keep <- !is.na(value) & !is.na(group)
  value <- value[keep]
  group <- droplevels(factor(group[keep]))
  if (length(unique(group)) != 2 || any(table(group) < 2)) {
    return(list(p_value = NA_real_, statistic = NA_real_))
  }
  result <- suppressWarnings(wilcox.test(value ~ group, exact = FALSE))
  list(p_value = result$p.value, statistic = unname(result$statistic))
}

safe_spearman <- function(value, time) {
  keep <- !is.na(value) & !is.na(time)
  value <- value[keep]
  time <- time[keep]
  if (length(value) < 4 || length(unique(value)) < 3 || length(unique(time)) < 3) {
    return(list(rho = NA_real_, p_value = NA_real_, n = length(value)))
  }
  result <- suppressWarnings(cor.test(value, time, method = "spearman", exact = FALSE))
  list(rho = unname(result$estimate), p_value = result$p.value, n = length(value))
}

safe_logrank <- function(time, event, group) {
  keep <- !is.na(time) & !is.na(event) & !is.na(group)
  time <- time[keep]
  event <- event[keep]
  group <- droplevels(factor(group[keep]))
  if (length(unique(group)) != 2 || any(table(group) < 2) || length(unique(event)) < 2) {
    return(list(p_value = NA_real_, chisq = NA_real_))
  }
  result <- survdiff(Surv(time, event) ~ group)
  p_value <- pchisq(result$chisq, df = length(result$n) - 1, lower.tail = FALSE)
  list(p_value = p_value, chisq = result$chisq)
}

safe_cox_continuous <- function(time, event, value) {
  keep <- !is.na(time) & !is.na(event) & !is.na(value)
  time <- time[keep]
  event <- event[keep]
  value <- value[keep]
  if (length(value) < 8 || sum(event) < 3 || sd(value) == 0) {
    return(list(
      n = length(value), events = sum(event), hr = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, p_value = NA_real_, zph_p_value = NA_real_,
      coef = NA_real_, se = NA_real_
    ))
  }
  z_value <- as.numeric(scale(value))
  model_df <- data.frame(time = time, event = event, z_value = z_value)
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(time, event) ~ z_value, data = model_df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(
      n = length(value), events = sum(event), hr = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, p_value = NA_real_, zph_p_value = NA_real_,
      coef = NA_real_, se = NA_real_
    ))
  }
  fit_summary <- summary(fit)
  zph_p <- tryCatch(
    suppressWarnings(cox.zph(fit)$table["z_value", "p"]),
    error = function(e) NA_real_
  )
  list(
    n = length(value),
    events = sum(event),
    hr = fit_summary$coefficients["z_value", "exp(coef)"],
    ci_low = fit_summary$conf.int["z_value", "lower .95"],
    ci_high = fit_summary$conf.int["z_value", "upper .95"],
    p_value = fit_summary$coefficients["z_value", "Pr(>|z|)"],
    zph_p_value = zph_p,
    coef = fit_summary$coefficients["z_value", "coef"],
    se = fit_summary$coefficients["z_value", "se(coef)"]
  )
}

safe_cox_median_split <- function(time, event, value) {
  keep <- !is.na(time) & !is.na(event) & !is.na(value)
  time <- time[keep]
  event <- event[keep]
  value <- value[keep]
  if (length(value) < 8 || sum(event) < 3 || length(unique(value)) < 3) {
    return(list(
      median_cutpoint = NA_real_, n_low = NA_integer_, n_high = NA_integer_,
      events_low = NA_integer_, events_high = NA_integer_, hr_high_vs_low = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, zph_p_value = NA_real_
    ))
  }
  median_cutpoint <- median(value, na.rm = TRUE)
  split_group <- factor(ifelse(value > median_cutpoint, "High", "Low"), levels = c("Low", "High"))
  if (length(unique(split_group)) != 2 || any(table(split_group) < 2)) {
    return(list(
      median_cutpoint = median_cutpoint, n_low = sum(split_group == "Low"),
      n_high = sum(split_group == "High"), events_low = sum(event[split_group == "Low"]),
      events_high = sum(event[split_group == "High"]), hr_high_vs_low = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, zph_p_value = NA_real_
    ))
  }
  model_df <- data.frame(time = time, event = event, split_group = split_group)
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(time, event) ~ split_group, data = model_df)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(
      median_cutpoint = median_cutpoint, n_low = sum(split_group == "Low"),
      n_high = sum(split_group == "High"), events_low = sum(event[split_group == "Low"]),
      events_high = sum(event[split_group == "High"]), hr_high_vs_low = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, zph_p_value = NA_real_
    ))
  }
  fit_summary <- summary(fit)
  coef_name <- "split_groupHigh"
  zph_p <- tryCatch(
    suppressWarnings(cox.zph(fit)$table[coef_name, "p"]),
    error = function(e) NA_real_
  )
  list(
    median_cutpoint = median_cutpoint,
    n_low = sum(split_group == "Low"),
    n_high = sum(split_group == "High"),
    events_low = sum(event[split_group == "Low"]),
    events_high = sum(event[split_group == "High"]),
    hr_high_vs_low = fit_summary$coefficients[coef_name, "exp(coef)"],
    ci_low = fit_summary$conf.int[coef_name, "lower .95"],
    ci_high = fit_summary$conf.int[coef_name, "upper .95"],
    p_value = fit_summary$coefficients[coef_name, "Pr(>|z|)"],
    zph_p_value = zph_p
  )
}

cells <- fread(
  spatial_path,
  select = c(
    "cell_id", "sampleID", "patientID", "cCR", "sample_timepoint",
    "cell_area", "time_to_PD_month", "time_to_CR_month",
    "major_cell_type", "final_cell_type2", neighborhood_order
  )
)
cells <- cells[sample_timepoint == "postC"]
for (neighborhood in neighborhood_order) {
  cells[, (neighborhood) := as.numeric(get(neighborhood))]
}

survival_df <- unique(cells[, .(patientID, sampleID, cCR, time_to_PD_month, time_to_CR_month)])
survival_df[, pd_event := !is.na(time_to_PD_month)]
survival_df[, pd_event_label := fifelse(pd_event, "PD recorded", "No PD recorded")]
max_pd_time <- max(survival_df$time_to_PD_month, na.rm = TRUE)
survival_df[, pfs_time_proxy_month := fifelse(pd_event, time_to_PD_month, max_pd_time)]
survival_df[, cCR := factor(cCR, levels = group_order)]
survival_df[, pd_event_label := factor(pd_event_label, levels = names(event_colors))]

spatial_scores <- cells[
  !is.na(cell_area) & cell_area > 0,
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
  by = .(patientID, sampleID)
]
spatial_long <- melt(
  spatial_scores,
  id.vars = c("patientID", "sampleID"),
  measure.vars = neighborhood_order,
  variable.name = "feature",
  value.name = "value"
)
spatial_long[, feature_label := sub("_k7$", "", sub("neighborhood_", "Neighborhood ", feature))]
spatial_long[, feature_group := "Spatial neighborhood"]

cell_denominator <- cells[
  !(major_cell_type %in% c("A2ML1+ epi", "malignant")) &
    !is.na(final_cell_type2),
  .(denominator_cells = .N),
  by = .(patientID, sampleID)
]
cell_counts <- cells[
  !(major_cell_type %in% c("A2ML1+ epi", "malignant")) &
    final_cell_type2 %in% c(t_cell_types, fibro_cell_types),
  .(n_cells = .N),
  by = .(patientID, sampleID, final_cell_type2)
]
cell_grid <- CJ(
  patientID = unique(cells$patientID),
  sampleID = unique(cells$sampleID),
  final_cell_type2 = c(t_cell_types, fibro_cell_types),
  unique = TRUE
)
cell_grid <- merge(cell_grid, survival_df[, .(patientID, sampleID)], by = c("patientID", "sampleID"))
cell_props <- merge(cell_grid, cell_counts, by = c("patientID", "sampleID", "final_cell_type2"), all.x = TRUE)
cell_props <- merge(cell_props, cell_denominator, by = c("patientID", "sampleID"), all.x = TRUE)
cell_props[is.na(n_cells), n_cells := 0]
cell_props[, value := n_cells / denominator_cells * 100]
cell_props[, feature := final_cell_type2]
cell_props[, feature_label := final_cell_type2]
cell_props[, feature_group := fifelse(
  final_cell_type2 %in% t_cell_types,
  "T immune cell type",
  "Fibroblast/endothelial cell type"
)]
cell_long <- cell_props[, .(patientID, sampleID, feature, feature_label, feature_group, value)]

feature_long <- rbindlist(
  list(
    spatial_long[, .(patientID, sampleID, feature, feature_label, feature_group, value)],
    cell_long
  ),
  use.names = TRUE
)
feature_long <- merge(feature_long, survival_df, by = c("patientID", "sampleID"), all.x = TRUE)
feature_long[, feature_group := factor(feature_group, levels = names(feature_group_colors))]

wide_features <- dcast(
  feature_long,
  patientID + sampleID + cCR + time_to_PD_month + time_to_CR_month + pd_event + pd_event_label + pfs_time_proxy_month ~ feature,
  value.var = "value"
)

fwrite(survival_df, file.path(out_dir, "postc_patient_survival_metadata.csv"))
fwrite(feature_long, file.path(out_dir, "postc_survival_feature_long.csv"))
fwrite(wide_features, file.path(out_dir, "postc_survival_feature_wide.csv"))

event_stats <- feature_long[
  ,
  {
    wilcox <- safe_wilcox(value, pd_event_label)
    no_pd_values <- value[pd_event_label == "No PD recorded"]
    pd_values <- value[pd_event_label == "PD recorded"]
    .(
      n_no_pd = sum(!is.na(no_pd_values)),
      n_pd = sum(!is.na(pd_values)),
      median_no_pd = median(no_pd_values, na.rm = TRUE),
      median_pd = median(pd_values, na.rm = TRUE),
      delta_median_pd_minus_no_pd = median(pd_values, na.rm = TRUE) - median(no_pd_values, na.rm = TRUE),
      wilcox_statistic = wilcox$statistic,
      p_value = wilcox$p_value
    )
  },
  by = .(feature_group, feature, feature_label)
]
event_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature_group]
setorder(event_stats, p_value)
fwrite(event_stats, file.path(out_dir, "postc_pd_event_association_stats.csv"))

time_stats <- feature_long[
  pd_event == TRUE,
  {
    spearman <- safe_spearman(value, time_to_PD_month)
    .(
      n_progressed = spearman$n,
      spearman_rho = spearman$rho,
      p_value = spearman$p_value
    )
  },
  by = .(feature_group, feature, feature_label)
]
time_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature_group]
setorder(time_stats, p_value)
fwrite(time_stats, file.path(out_dir, "postc_time_to_pd_spearman_stats_progressed_only.csv"))

proxy_stats <- feature_long[
  ,
  {
    median_value <- median(value, na.rm = TRUE)
    split_group <- factor(ifelse(value > median_value, "High", "Low"), levels = c("Low", "High"))
    logrank <- safe_logrank(pfs_time_proxy_month, pd_event, split_group)
    .(
      median_cutpoint = median_value,
      n_low = sum(split_group == "Low", na.rm = TRUE),
      n_high = sum(split_group == "High", na.rm = TRUE),
      logrank_chisq = logrank$chisq,
      p_value = logrank$p_value
    )
  },
  by = .(feature_group, feature, feature_label)
]
proxy_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature_group]
setorder(proxy_stats, p_value)
fwrite(proxy_stats, file.path(out_dir, "postc_proxy_survival_median_split_logrank_stats.csv"))

cox_continuous_stats <- feature_long[
  ,
  {
    cox_result <- safe_cox_continuous(pfs_time_proxy_month, pd_event, value)
    .(
      endpoint = "proxy_progression_free_survival",
      model = "univariable Cox PH",
      predictor_scale = "continuous, per 1 SD increase",
      n = cox_result$n,
      events = cox_result$events,
      hazard_ratio = cox_result$hr,
      ci_low_95 = cox_result$ci_low,
      ci_high_95 = cox_result$ci_high,
      cox_coef = cox_result$coef,
      cox_se = cox_result$se,
      p_value = cox_result$p_value,
      proportional_hazards_test_p = cox_result$zph_p_value
    )
  },
  by = .(feature_group, feature, feature_label)
]
cox_continuous_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature_group]
setorder(cox_continuous_stats, p_value)
fwrite(cox_continuous_stats, file.path(out_dir, "postc_proxy_pfs_cox_continuous_stats.csv"))

cox_median_stats <- feature_long[
  ,
  {
    cox_result <- safe_cox_median_split(pfs_time_proxy_month, pd_event, value)
    .(
      endpoint = "proxy_progression_free_survival",
      model = "univariable Cox PH",
      predictor_scale = "median split, High versus Low",
      median_cutpoint = cox_result$median_cutpoint,
      n_low = cox_result$n_low,
      n_high = cox_result$n_high,
      events_low = cox_result$events_low,
      events_high = cox_result$events_high,
      hazard_ratio_high_vs_low = cox_result$hr_high_vs_low,
      ci_low_95 = cox_result$ci_low,
      ci_high_95 = cox_result$ci_high,
      p_value = cox_result$p_value,
      proportional_hazards_test_p = cox_result$zph_p_value
    )
  },
  by = .(feature_group, feature, feature_label)
]
cox_median_stats[, p_adjusted_bh := p.adjust(p_value, method = "BH"), by = feature_group]
setorder(cox_median_stats, p_value)
fwrite(cox_median_stats, file.path(out_dir, "postc_proxy_pfs_cox_median_split_stats.csv"))

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

top_event_features <- event_stats[!is.na(p_value)][order(p_value)][seq_len(min(.N, 12)), .(feature)]
event_plot_data <- feature_long[feature %in% top_event_features$feature]
event_plot_data[, feature_label := factor(feature_label, levels = event_stats[feature %in% top_event_features$feature][order(p_value), feature_label])]
event_label_data <- event_stats[feature %in% top_event_features$feature]
event_label_data[, feature_label := factor(feature_label, levels = levels(event_plot_data$feature_label))]
event_label_data[, label := paste0(
  "p=", format_p_value(p_value),
  "\nBH=", format_p_value(p_adjusted_bh)
)]
event_label_data <- merge(
  event_label_data,
  event_plot_data[, .(label_y = max(value, na.rm = TRUE) * 1.12 + 0.01), by = feature_label],
  by = "feature_label",
  all.x = TRUE
)

event_plot <- ggplot(event_plot_data, aes(x = pd_event_label, y = value)) +
  geom_boxplot(
    aes(fill = pd_event_label),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.65,
    color = "#333333",
    linewidth = line_width
  ) +
  geom_point(
    aes(color = pd_event_label),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.0,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  geom_text(
    data = event_label_data,
    aes(x = 1.5, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 1.75,
    lineheight = 0.9
  ) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = event_colors, drop = FALSE, name = NULL) +
  scale_color_manual(values = event_colors, drop = FALSE, name = NULL) +
  labs(x = NULL, y = "Feature value") +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = small_size),
    axis.title.y = element_text(size = base_size, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = base_size),
    strip.background = element_blank(),
    strip.text = element_text(size = small_size, face = "bold"),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.spacing = unit(0.18, "cm")
  )
ggsave(
  file.path(out_dir, "postc_top_features_by_pd_event.pdf"),
  event_plot,
  width = 180,
  height = 135,
  units = "mm",
  dpi = 300
)

top_time_features <- time_stats[!is.na(p_value)][order(p_value)][seq_len(min(.N, 12)), .(feature)]
time_plot_data <- feature_long[pd_event == TRUE & feature %in% top_time_features$feature]
time_plot_data[, feature_label := factor(feature_label, levels = time_stats[feature %in% top_time_features$feature][order(p_value), feature_label])]
time_label_data <- time_stats[feature %in% top_time_features$feature]
time_label_data[, feature_label := factor(feature_label, levels = levels(time_plot_data$feature_label))]
time_label_data[, label := paste0(
  "rho=", sprintf("%.2f", spearman_rho),
  "\np=", format_p_value(p_value)
)]
time_label_data <- merge(
  time_label_data,
  time_plot_data[, .(
    label_x = min(time_to_PD_month, na.rm = TRUE),
    label_y = max(value, na.rm = TRUE) * 1.05 + 0.01
  ), by = feature_label],
  by = "feature_label",
  all.x = TRUE
)

time_plot <- ggplot(time_plot_data, aes(x = time_to_PD_month, y = value)) +
  geom_point(aes(color = cCR), size = 1.4, alpha = 0.95) +
  geom_smooth(method = "lm", se = FALSE, color = "#333333", linewidth = line_width) +
  geom_text(
    data = time_label_data,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 1.75,
    lineheight = 0.9
  ) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = outcome_colors, breaks = group_order, name = "Clinical Outcome") +
  labs(x = "Time to PD (months; progressed patients only)", y = "Feature value") +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    axis.text = element_text(size = small_size),
    axis.title = element_text(size = base_size, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    strip.background = element_blank(),
    strip.text = element_text(size = small_size, face = "bold"),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.spacing = unit(0.18, "cm")
  )
ggsave(
  file.path(out_dir, "postc_top_features_vs_time_to_pd_progressed_only.pdf"),
  time_plot,
  width = 180,
  height = 135,
  units = "mm",
  dpi = 300
)

top_proxy_features <- proxy_stats[!is.na(p_value)][order(p_value)][seq_len(min(.N, 6)), .(feature)]
km_data <- feature_long[feature %in% top_proxy_features$feature]
km_data <- merge(km_data, proxy_stats[, .(feature, median_cutpoint)], by = "feature", all.x = TRUE)
km_data[, split_group := factor(ifelse(value > median_cutpoint, "High", "Low"), levels = c("Low", "High"))]
km_data[, feature_label := factor(feature_label, levels = proxy_stats[feature %in% top_proxy_features$feature][order(p_value), feature_label])]

km_curves <- rbindlist(lapply(split(km_data, km_data$feature_label), function(sub_df) {
  fit <- survfit(Surv(pfs_time_proxy_month, pd_event) ~ split_group, data = sub_df)
  fit_summary <- summary(fit)
  if (length(fit_summary$time) == 0) {
    return(NULL)
  }
  data.table(
    feature_label = as.character(unique(sub_df$feature_label)),
    time = fit_summary$time,
    survival = fit_summary$surv,
    strata = sub("^split_group=", "", fit_summary$strata)
  )
}), fill = TRUE)

if (nrow(km_curves) > 0) {
  km_curves[, feature_label := factor(feature_label, levels = levels(km_data$feature_label))]
  proxy_plot <- ggplot(km_curves, aes(x = time, y = survival, color = strata)) +
    geom_step(linewidth = line_width * 1.5) +
    facet_wrap(~ feature_label, ncol = 3) +
    scale_color_manual(values = c("Low" = "#4C78A8", "High" = "#D95F02"), name = "Median split") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.01, 0.02))) +
    labs(
      x = "Proxy time to PD (months)",
      y = "PD-free fraction"
    ) +
    theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
    theme(
      axis.text = element_text(size = small_size),
      axis.title = element_text(size = base_size, face = "bold"),
      legend.position = "top",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      strip.background = element_blank(),
      strip.text = element_text(size = small_size, face = "bold"),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
      axis.line = element_line(linewidth = line_width),
      axis.ticks = element_line(linewidth = line_width),
      panel.spacing = unit(0.22, "cm")
    )
  ggsave(
    file.path(out_dir, "postc_top_features_proxy_pd_free_curves.pdf"),
    proxy_plot,
    width = 180,
    height = 115,
    units = "mm",
    dpi = 300
  )
}

top_summary <- rbindlist(
  list(
    event_stats[1:min(.N, 10)][, analysis := "PD event association"],
    time_stats[1:min(.N, 10)][, analysis := "Time to PD Spearman among progressed"],
    proxy_stats[1:min(.N, 10)][, analysis := "Proxy median-split log-rank"]
  ),
  fill = TRUE
)
fwrite(top_summary, file.path(out_dir, "top_ranked_feature_summary.csv"))

writeLines(
  c(
    "# Exploratory postC survival association analysis",
    "",
    "This exploratory analysis tests whether postC spatial neighborhood scores and T/fibroblast/endothelial cell-type proportions are associated with progression-related outcomes.",
    "",
    "Survival metadata found in the old OBS/spatial metadata: time_to_PD_month and time_to_CR_month. No explicit censoring/follow-up time column was found.",
    "",
    "Spatial neighborhood features were calculated per patient as area-weighted continuous NMF scores: sum(neighborhood weight x cell area) divided by total annotated cell area.",
    "",
    "Cell-type features were calculated as percent of all non-epithelial, non-malignant cells at postC.",
    "",
    "Primary exploratory outputs:",
    "- PD-event association: Wilcoxon rank-sum test comparing patients with versus without recorded PD.",
    "- Time-to-PD association: Spearman correlation among patients with recorded PD only.",
    "- Proxy PD-free curves: median split by feature value; patients without recorded PD were censored at the maximum observed PD time. This is only a rough visual screen because true censoring times are not available.",
    "- Cox PH continuous models: univariable Cox proportional hazards models using each postC feature scaled to mean 0 and SD 1; hazard ratios therefore correspond to a 1 SD increase in the feature.",
    "- Cox PH median-split models: univariable Cox proportional hazards models comparing high versus low feature groups split at the median; this is the same grouping used for the proxy PD-free curves.",
    "- Proportional hazards diagnostics: Schoenfeld residual global p values were calculated for each univariable Cox model when possible.",
    "",
    "These results should not be treated as final survival statistics until an explicit event indicator and censor/follow-up time are provided."
  ),
  con = file.path(legend_dir, "analysis_notes.md")
)
