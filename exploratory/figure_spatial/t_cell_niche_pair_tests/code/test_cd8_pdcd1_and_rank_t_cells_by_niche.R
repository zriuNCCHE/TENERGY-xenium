suppressPackageStartupMessages({
  library(data.table)
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

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value < 1e-4) {
    return(formatC(p_value, format = "e", digits = 2))
  }
  if (p_value < 0.001) {
    return(sprintf("%.6f", p_value))
  }
  sprintf("%.3f", p_value)
}

target_feature <- "CD8_Tex_PDCD1"
niche_1 <- "Niche_iCAF_CXCL5"
niche_2 <- "Niche_Macro_CXCL9"

dt <- fread(source_path)
dt <- dt[feature_type == "Individual T subset"]
dt[, feature := as.character(feature)]
dt[, neighborhood_label := as.character(neighborhood_label)]

target_dt <- dt[
  feature == target_feature &
    neighborhood_label %in% c(niche_1, niche_2)
]

paired_dt <- dcast(
  target_dt,
  patientID + sampleID + cCR ~ neighborhood_label,
  value.var = "cell_percent"
)
paired_dt <- paired_dt[
  !is.na(get(niche_1)) &
    !is.na(get(niche_2))
]
paired_dt[, delta_iCAF_CXCL5_minus_Macro_CXCL9 := get(niche_1) - get(niche_2)]

paired_test <- suppressWarnings(wilcox.test(
  paired_dt[[niche_1]],
  paired_dt[[niche_2]],
  paired = TRUE,
  alternative = "two.sided",
  exact = FALSE
))
paired_greater_test <- suppressWarnings(wilcox.test(
  paired_dt[[niche_1]],
  paired_dt[[niche_2]],
  paired = TRUE,
  alternative = "greater",
  exact = FALSE
))

unpaired_test <- suppressWarnings(wilcox.test(
  cell_percent ~ neighborhood_label,
  data = target_dt,
  alternative = "two.sided",
  exact = FALSE
))
unpaired_greater_test <- suppressWarnings(wilcox.test(
  target_dt[neighborhood_label == niche_1, cell_percent],
  target_dt[neighborhood_label == niche_2, cell_percent],
  alternative = "greater",
  exact = FALSE
))

stats_dt <- data.table(
  feature = target_feature,
  comparison = paste(niche_1, "vs", niche_2),
  metric = "cell_percent",
  test = c(
    "paired Wilcoxon signed-rank",
    "paired Wilcoxon signed-rank",
    "unpaired Wilcoxon rank-sum",
    "unpaired Wilcoxon rank-sum"
  ),
  alternative = c("two-sided", paste(niche_1, "greater"), "two-sided", paste(niche_1, "greater")),
  n_iCAF_CXCL5 = sum(target_dt$neighborhood_label == niche_1),
  n_Macro_CXCL9 = sum(target_dt$neighborhood_label == niche_2),
  n_paired_patients = nrow(paired_dt),
  median_iCAF_CXCL5 = median(target_dt[neighborhood_label == niche_1, cell_percent], na.rm = TRUE),
  median_Macro_CXCL9 = median(target_dt[neighborhood_label == niche_2, cell_percent], na.rm = TRUE),
  median_delta_iCAF_CXCL5_minus_Macro_CXCL9 = median(paired_dt$delta_iCAF_CXCL5_minus_Macro_CXCL9, na.rm = TRUE),
  p_value = c(
    paired_test$p.value,
    paired_greater_test$p.value,
    unpaired_test$p.value,
    unpaired_greater_test$p.value
  )
)

rank_dt <- dt[
  ,
  .(
    n_patients = uniqueN(patientID),
    median_cell_percent = median(cell_percent, na.rm = TRUE),
    mean_cell_percent = mean(cell_percent, na.rm = TRUE),
    max_cell_percent = max(cell_percent, na.rm = TRUE),
    median_area_percent = median(area_percent, na.rm = TRUE),
    mean_area_percent = mean(area_percent, na.rm = TRUE),
    max_area_percent = max(area_percent, na.rm = TRUE)
  ),
  by = .(neighborhood_label, feature)
]
rank_dt <- rank_dt[
  order(neighborhood_label, -median_cell_percent, -mean_cell_percent, feature)
]
rank_dt[, rank_by_median_cell_percent := seq_len(.N), by = neighborhood_label]
setcolorder(
  rank_dt,
  c(
    "neighborhood_label", "rank_by_median_cell_percent", "feature",
    "n_patients", "median_cell_percent", "mean_cell_percent", "max_cell_percent",
    "median_area_percent", "mean_area_percent", "max_area_percent"
  )
)

fwrite(paired_dt, file.path(out_dir, "cd8_tex_pdcd1_icaf_cxcl5_vs_macro_cxcl9_paired_values.csv"))
fwrite(stats_dt, file.path(out_dir, "cd8_tex_pdcd1_icaf_cxcl5_vs_macro_cxcl9_tests.csv"))
fwrite(rank_dt, file.path(out_dir, "t_cell_subset_rank_by_abundance_within_each_niche.csv"))

writeLines(
  c(
    "# CD8_Tex_PDCD1 niche comparison and T-cell abundance ranking",
    "",
    paste0("Target comparison: `", target_feature, "` in `", niche_1, "` versus `", niche_2, "`."),
    paste0("Median ", niche_1, " = ", sprintf("%.3f", stats_dt$median_iCAF_CXCL5[1]), "%; median ", niche_2, " = ", sprintf("%.3f", stats_dt$median_Macro_CXCL9[1]), "%."),
    paste0("Paired two-sided Wilcoxon P = ", format_p_value(stats_dt[test == "paired Wilcoxon signed-rank" & alternative == "two-sided", p_value]), "."),
    paste0("Paired one-sided Wilcoxon P for ", niche_1, " greater = ", format_p_value(stats_dt[test == "paired Wilcoxon signed-rank" & alternative != "two-sided", p_value]), "."),
    "",
    "T-cell subsets were ranked within each niche by median cell fraction across patients."
  ),
  con = file.path(legend_dir, "analysis_notes.md")
)
