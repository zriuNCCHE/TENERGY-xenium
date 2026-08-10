suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

mp_cell_path <- file.path(project_root, "data", "geo_ready", "figure_7_mp7_myeloid_distance", "MP_to_monocyte_neutrophil_nearest_distance_cell_level.csv")
metadata_path <- file.path(project_root, "data", "raw", "combined_final_all_cell_types_obs.csv")
out_dir <- file.path(project_root, "final", "extended_figure_6", "mp7_nonccr_association_candidates", "outputs")
legend_dir <- file.path(project_root, "final", "extended_figure_6", "mp7_nonccr_association_candidates", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

mp_levels <- paste0("MP_", 1:7)
mp_palette <- c(
  MP_1 = "#1F77B4",
  MP_2 = "#6FA4C8",
  MP_3 = "#E1812C",
  MP_4 = "#F0C08A",
  MP_5 = "#2CA02C",
  MP_6 = "#A1D99B",
  MP_7 = "#C52B2F"
)
outcome_cols <- c("cCR" = "#7BB6A4", "non-cCR" = "#F2A38A")

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
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
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

mp_cells <- fread(mp_cell_path, select = c("sampleID", "cell_id", "MP", "target_label"))
mp_cells <- unique(mp_cells[target_label == "monocyte_neutrophil", .(sampleID, cell_id, MP)])
mp_cells[, MP := factor(MP, levels = mp_levels)]

metadata <- fread(metadata_path, select = c("sampleID", "patientID", "cCR", "sample_timepoint"))
sample_meta <- unique(metadata[sample_timepoint == "pre", .(sampleID, patientID, cCR)])

sample_mp_counts <- mp_cells[, .N, by = .(sampleID, MP)]
sample_total <- mp_cells[, .(total_malignant_mp_cells = uniqueN(cell_id)), by = sampleID]
sample_mp <- merge(sample_mp_counts, sample_total, by = "sampleID")
sample_mp[, mp_percent := N / total_malignant_mp_cells * 100]
sample_mp <- merge(sample_mp, sample_meta, by = "sampleID", all.x = TRUE)
sample_mp[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]
sample_mp[, MP := factor(MP, levels = mp_levels)]

all_combos <- CJ(sampleID = unique(mp_cells$sampleID), MP = factor(mp_levels, levels = mp_levels))
sample_mp <- merge(all_combos, sample_mp, by = c("sampleID", "MP"), all.x = TRUE)
sample_mp <- merge(sample_mp, sample_total, by = "sampleID", all.x = TRUE, suffixes = c("", "_total"))
sample_mp <- merge(sample_mp, sample_meta, by = "sampleID", all.x = TRUE, suffixes = c("", "_meta"))
sample_mp[is.na(N), N := 0]
sample_mp[is.na(mp_percent), mp_percent := 0]
sample_mp[, cCR := fifelse(!is.na(cCR), as.character(cCR), as.character(cCR_meta))]
sample_mp[, patientID := fifelse(!is.na(patientID), as.character(patientID), as.character(patientID_meta))]
sample_mp[, cCR := factor(cCR, levels = c("cCR", "non-cCR"))]
sample_mp[, MP := factor(MP, levels = mp_levels)]
sample_mp <- sample_mp[, .(sampleID, patientID, cCR, MP, mp_cells = N, total_malignant_mp_cells, mp_percent)]

mp7 <- sample_mp[MP == "MP_7"]
wilcox_p <- tryCatch(
  wilcox.test(mp_percent ~ cCR, data = mp7, exact = FALSE)$p.value,
  error = function(e) NA_real_
)
count_p <- tryCatch(
  wilcox.test(mp_cells ~ cCR, data = mp7, exact = FALSE)$p.value,
  error = function(e) NA_real_
)

mp_stats <- rbindlist(lapply(mp_levels, function(mp) {
  sub <- sample_mp[MP == mp]
  p <- tryCatch(wilcox.test(mp_percent ~ cCR, data = sub, exact = FALSE)$p.value, error = function(e) NA_real_)
  data.table(
    MP = mp,
    n_cCR = sum(sub$cCR == "cCR", na.rm = TRUE),
    n_non_cCR = sum(sub$cCR == "non-cCR", na.rm = TRUE),
    median_cCR = median(sub[cCR == "cCR", mp_percent], na.rm = TRUE),
    median_non_cCR = median(sub[cCR == "non-cCR", mp_percent], na.rm = TRUE),
    delta_non_cCR_minus_cCR = median(sub[cCR == "non-cCR", mp_percent], na.rm = TRUE) - median(sub[cCR == "cCR", mp_percent], na.rm = TRUE),
    p_value = p
  )
}))
mp_stats[, p_adjusted_BH := p.adjust(p_value, method = "BH")]

fwrite(sample_mp, file.path(out_dir, "extended_figure6_mp_abundance_by_sample.csv"))
fwrite(mp7, file.path(out_dir, "extended_figure6_mp7_abundance_by_sample.csv"))
fwrite(mp_stats, file.path(out_dir, "extended_figure6_all_mp_response_stats.csv"))

box_plot <- ggplot(mp7, aes(cCR, mp_percent)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(aes(colour = cCR), size = 0.9, alpha = 1, position = position_jitter(width = 0.08, height = 0), show.legend = FALSE) +
  scale_colour_manual(values = outcome_cols) +
  labs(
    title = "MP7 abundance by clinical outcome",
    x = NULL,
    y = "MP7 cells among malignant MPs (%)"
  ) +
  annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.35, size = 2.0, label = paste0("Wilcoxon ", fmt_p(wilcox_p))) +
  theme_nc()
mm_pdf("extended_figure6_candidate_A_mp7_percent_boxplot.pdf", 55, 58)
print(box_plot)
dev.off()

count_plot <- ggplot(mp7, aes(cCR, mp_cells)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white") +
  geom_point(aes(colour = cCR), size = 0.9, alpha = 1, position = position_jitter(width = 0.08, height = 0), show.legend = FALSE) +
  scale_colour_manual(values = outcome_cols) +
  labs(
    title = "MP7 cell count by clinical outcome",
    x = NULL,
    y = "MP7 cells per sample"
  ) +
  annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.35, size = 2.0, label = paste0("Wilcoxon ", fmt_p(count_p))) +
  theme_nc()
mm_pdf("extended_figure6_candidate_B_mp7_count_boxplot.pdf", 55, 58)
print(count_plot)
dev.off()

rank_dt <- copy(mp7)
rank_dt[, sampleID := factor(sampleID, levels = rank_dt[order(cCR, mp_percent), sampleID])]
rank_plot <- ggplot(rank_dt, aes(sampleID, mp_percent, fill = cCR)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = outcome_cols, name = "Clinical outcome") +
  labs(
    title = "Pretreatment MP7 fraction by sample",
    x = NULL,
    y = "MP7 cells among malignant MPs (%)"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
mm_pdf("extended_figure6_candidate_C_mp7_sample_rank_barplot.pdf", 105, 58)
print(rank_plot)
dev.off()

composition_dt <- copy(sample_mp)
sample_order <- mp7[order(cCR, mp_percent), sampleID]
composition_dt[, sampleID := factor(sampleID, levels = sample_order)]
composition_dt[, MP := factor(MP, levels = rev(mp_levels))]
composition_plot <- ggplot(composition_dt, aes(sampleID, mp_percent, fill = MP)) +
  geom_col(width = 0.74, colour = NA) +
  scale_fill_manual(values = mp_palette, name = "Malignant MP") +
  labs(
    title = "Malignant MP composition by sample",
    x = NULL,
    y = "Cells among malignant MPs (%)"
  ) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
mm_pdf("extended_figure6_candidate_D_all_mp_composition_stacked_barplot.pdf", 120, 62)
print(composition_plot)
dev.off()

all_mp_plot <- ggplot(sample_mp, aes(MP, mp_percent, colour = cCR)) +
  geom_boxplot(aes(group = interaction(MP, cCR)), width = 0.58, outlier.shape = NA, linewidth = 0.5, colour = "#7FAEC8", fill = "white", position = position_dodge(width = 0.68)) +
  geom_point(size = 0.75, alpha = 1, position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.68), show.legend = FALSE) +
  scale_colour_manual(values = outcome_cols, name = "Clinical outcome") +
  labs(
    title = "All malignant MP fractions by outcome",
    x = "Malignant-cell metaprogram",
    y = "Cells among malignant MPs (%)"
  ) +
  theme_nc() +
  theme(legend.position = "right")
mm_pdf("extended_figure6_candidate_E_all_mp_response_boxplots.pdf", 104, 62)
print(all_mp_plot)
dev.off()

delta_plot_dt <- copy(mp_stats)
delta_plot_dt[, MP := factor(MP, levels = mp_levels)]
delta_plot <- ggplot(delta_plot_dt, aes(MP, delta_non_cCR_minus_cCR, fill = MP)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#555555") +
  geom_col(width = 0.62, colour = "black", linewidth = 0.35) +
  scale_fill_manual(values = mp_palette, guide = "none") +
  labs(
    title = "Response-associated MP abundance shift",
    x = "Malignant-cell metaprogram",
    y = "Median difference: non-cCR - cCR (%)"
  ) +
  theme_nc()
mm_pdf("extended_figure6_candidate_F_all_mp_nonccr_minus_ccr_delta.pdf", 82, 58)
print(delta_plot)
dev.off()

legend_text <- c(
  "Extended Figure 6 MP7 non-cCR association candidate legends",
  "",
  "Candidate A. Pretreatment MP7 abundance by clinical outcome. Each dot is one pretreatment sample. The y-axis shows the percentage of MP7 cells among malignant metaprogram-assigned cells. Boxes show median and interquartile range. P value was calculated by two-sided Wilcoxon rank-sum test.",
  "",
  "Candidate B. Pretreatment MP7 cell counts by clinical outcome. Each dot is one pretreatment sample. Boxes show median and interquartile range. P value was calculated by two-sided Wilcoxon rank-sum test.",
  "",
  "Candidate C. Pretreatment MP7 fraction ranked by sample and clinical outcome.",
  "",
  "Candidate D. Pretreatment malignant metaprogram composition per sample. Samples are ordered by clinical outcome and MP7 fraction.",
  "",
  "Candidate E. Pretreatment abundance of all malignant metaprograms by clinical outcome. Each dot is one pretreatment sample.",
  "",
  "Candidate F. Difference in median pretreatment metaprogram fraction between non-cCR and cCR samples."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# Extended Figure 6 MP7 Non-cCR Association Candidates",
  "",
  "Candidate plots testing whether pretreatment MP7 abundance is associated with non-cCR clinical outcome.",
  "",
  "Source data are malignant-cell MP assignments from `data/geo_ready/figure_7_mp7_myeloid_distance/MP_to_monocyte_neutrophil_nearest_distance_cell_level.csv`, using the pooled target rows to avoid duplicate cells.",
  "",
  "## Outputs",
  "",
  "- `outputs/extended_figure6_candidate_A_mp7_percent_boxplot.pdf`",
  "- `outputs/extended_figure6_candidate_B_mp7_count_boxplot.pdf`",
  "- `outputs/extended_figure6_candidate_C_mp7_sample_rank_barplot.pdf`",
  "- `outputs/extended_figure6_candidate_D_all_mp_composition_stacked_barplot.pdf`",
  "- `outputs/extended_figure6_candidate_E_all_mp_response_boxplots.pdf`",
  "- `outputs/extended_figure6_candidate_F_all_mp_nonccr_minus_ccr_delta.pdf`"
)
writeLines(readme_text, file.path(dirname(out_dir), "README.md"))

message("Saved Extended Figure 6 MP7 non-cCR candidate plots to: ", out_dir)
