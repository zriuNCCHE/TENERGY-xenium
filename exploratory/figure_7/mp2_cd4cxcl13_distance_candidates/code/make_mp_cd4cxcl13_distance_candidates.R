suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else getwd()
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))
data_dir <- file.path(project_root, "data", "geo_ready", "figure_7_mp_cd4cxcl13_distance")
out_dir <- file.path(project_root, "exploratory", "figure_7", "mp2_cd4cxcl13_distance_candidates", "outputs")
legend_dir <- file.path(project_root, "exploratory", "figure_7", "mp2_cd4cxcl13_distance_candidates", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  )
}

fmt_fdr <- function(p) {
  ifelse(is.na(p), "FDR NA",
    ifelse(p < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", p))
  )
}

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
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
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

mp_levels <- paste0("MP_", 1:7)
mp_fill <- setNames(rep("#B8B8B8", 7), mp_levels)
mp_fill["MP_2"] <- "#C44E52"
mp_point <- setNames(rep("#4D4D4D", 7), mp_levels)
mp_point["MP_2"] <- "#9F2F3A"

sample_summary <- fread(file.path(data_dir, "MP_to_CD4_CXCL13_nearest_distance_sample_summary.csv"))
cell_summary <- fread(file.path(data_dir, "MP_to_CD4_CXCL13_nearest_distance_cell_summary.csv"))
cell_level <- fread(file.path(data_dir, "MP_to_CD4_CXCL13_nearest_distance_cell_level.csv"))

sample_summary[, MP := factor(MP, levels = mp_levels)]
cell_summary[, MP := factor(MP, levels = mp_levels)]
cell_level[, MP := factor(MP, levels = mp_levels)]

wide <- dcast(sample_summary, sampleID ~ MP, value.var = "median_distance")
complete_wide <- wide[complete.cases(wide[, ..mp_levels])]
friedman_p <- tryCatch(
  friedman.test(as.matrix(complete_wide[, ..mp_levels]))$p.value,
  error = function(e) NA_real_
)

pairwise <- rbindlist(lapply(setdiff(mp_levels, "MP_2"), function(other_mp) {
  paired <- merge(
    sample_summary[MP == "MP_2", .(sampleID, median_MP2 = median_distance)],
    sample_summary[MP == other_mp, .(sampleID, median_other = median_distance)],
    by = "sampleID"
  )
  p <- tryCatch(
    wilcox.test(paired$median_MP2, paired$median_other, paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  data.table(
    comparison = paste("MP_2 vs", other_mp),
    other_MP = other_mp,
    n_samples = nrow(paired),
    median_MP2 = median(paired$median_MP2, na.rm = TRUE),
    median_other = median(paired$median_other, na.rm = TRUE),
    median_delta_other_minus_MP2 = median(paired$median_other - paired$median_MP2, na.rm = TRUE),
    p_value = p
  )
}))
pairwise[, p_adjusted_BH := p.adjust(p_value, method = "BH")]
fwrite(pairwise, file.path(out_dir, "figure7_recomputed_sample_level_pairwise_stats.csv"))

sample_rank <- sample_summary[, .(
  n_samples = .N,
  median_of_sample_medians = median(median_distance, na.rm = TRUE),
  q25 = quantile(median_distance, 0.25, na.rm = TRUE),
  q75 = quantile(median_distance, 0.75, na.rm = TRUE)
), by = MP][order(median_of_sample_medians)]
fwrite(sample_rank, file.path(out_dir, "figure7_sample_level_mp_distance_rank.csv"))

delta_df <- rbindlist(lapply(setdiff(mp_levels, "MP_2"), function(other_mp) {
  paired <- merge(
    sample_summary[MP == "MP_2", .(sampleID, median_MP2 = median_distance)],
    sample_summary[MP == other_mp, .(sampleID, median_other = median_distance)],
    by = "sampleID"
  )
  paired[, other_MP := other_mp]
  paired[, delta_other_minus_MP2 := median_other - median_MP2]
  paired
}))
delta_df[, other_MP := factor(other_MP, levels = setdiff(mp_levels, "MP_2"))]
fwrite(delta_df, file.path(out_dir, "figure7_sample_level_mp2_distance_delta_by_sample.csv"))

main_box <- ggplot(sample_summary, aes(MP, median_distance, fill = MP)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.5, colour = "black") +
  geom_point(aes(colour = MP), size = 1.35, alpha = 1, position = position_jitter(width = 0.09, height = 0), show.legend = FALSE) +
  scale_fill_manual(values = mp_fill, guide = "none") +
  scale_colour_manual(values = mp_point, guide = "none") +
  labs(
    title = "MP2 is closest to CD4_CXCL13 T cells",
    x = "Malignant-cell metaprogram",
    y = "Median nearest distance to CD4_CXCL13 cell"
  ) +
  annotate("text", x = 1, y = Inf, vjust = 1.4, hjust = 0, size = 2.1,
           label = paste0("Friedman ", fmt_p(friedman_p))) +
  theme_nc()
mm_pdf("figure7_main_candidate_A_sample_median_distance_boxplot.pdf", 86, 70)
print(main_box)
dev.off()

main_rank <- ggplot(sample_rank, aes(reorder(MP, median_of_sample_medians), median_of_sample_medians)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.12, linewidth = 0.5, colour = "black") +
  geom_point(aes(fill = MP), shape = 21, size = 2.2, stroke = 0.5, colour = "black") +
  scale_fill_manual(values = mp_fill, guide = "none") +
  coord_flip() +
  labs(
    title = "Ranked proximity to CD4_CXCL13 T cells",
    x = NULL,
    y = "Median distance across samples"
  ) +
  annotate("text", x = Inf, y = Inf, vjust = 1.4, hjust = 1, size = 2.1,
           label = "Lower values indicate closer proximity") +
  theme_nc()
mm_pdf("figure7_main_candidate_B_ranked_sample_median_distance.pdf", 75, 62)
print(main_rank)
dev.off()

delta_labels <- pairwise[, .(other_MP = factor(other_MP, levels = setdiff(mp_levels, "MP_2")),
                             label = fmt_fdr(p_adjusted_BH))]
delta_y <- max(delta_df$delta_other_minus_MP2, na.rm = TRUE)
delta_min <- min(delta_df$delta_other_minus_MP2, na.rm = TRUE)
delta_pad <- max(10, 0.10 * (delta_y - delta_min))

main_delta <- ggplot(delta_df, aes(other_MP, delta_other_minus_MP2)) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed", colour = "#555555") +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.5, fill = "#D8D8D8", colour = "black") +
  geom_point(size = 1.25, alpha = 1, position = position_jitter(width = 0.09, height = 0), colour = "#333333") +
  geom_text(data = delta_labels, aes(x = other_MP, y = delta_y + delta_pad, label = label),
            inherit.aes = FALSE, size = 2.0, vjust = 0) +
  labs(
    title = "Sample-paired MP2 proximity advantage",
    x = "Comparison metaprogram",
    y = "Median distance(other MP) - median distance(MP2)"
  ) +
  theme_nc()
mm_pdf("figure7_main_candidate_C_mp2_closeness_delta.pdf", 92, 68)
print(main_delta)
dev.off()

cell_plot_df <- copy(cell_level)
y_cap <- quantile(cell_plot_df$nearest_CD4_CXCL13_distance, 0.985, na.rm = TRUE)
cell_box <- ggplot(cell_plot_df, aes(MP, nearest_CD4_CXCL13_distance, fill = MP)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.5, colour = "black") +
  coord_cartesian(ylim = c(0, y_cap)) +
  scale_fill_manual(values = mp_fill, guide = "none") +
  labs(
    title = "Cell-level distance distribution",
    x = "Malignant-cell metaprogram",
    y = "Nearest distance to CD4_CXCL13 cell"
  ) +
  annotate("text", x = 1, y = y_cap, vjust = 1.4, hjust = 0, size = 2.1,
           label = "Display truncated at 98.5th percentile") +
  theme_nc()
mm_pdf("figure7_supp_candidate_A_cell_level_distance_boxplot_zoomed.pdf", 92, 68)
print(cell_box)
dev.off()

heat_df <- copy(sample_summary)
sample_order <- sample_summary[MP == "MP_2"][order(median_distance), sampleID]
heat_df[, sampleID := factor(sampleID, levels = sample_order)]
heat_df[, MP := factor(MP, levels = mp_levels)]
heat_cap <- quantile(heat_df$median_distance, 0.95, na.rm = TRUE)
heat <- ggplot(heat_df, aes(MP, sampleID, fill = pmin(median_distance, heat_cap))) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colours = c("#204A87", "#F7F7F7", "#B2182B"),
    name = "Median\ndistance"
  ) +
  labs(
    title = "Sample-level MP proximity profile",
    x = "Malignant-cell metaprogram",
    y = "Sample"
  ) +
  theme_nc() +
  theme(panel.border = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
mm_pdf("figure7_supp_candidate_B_sample_level_distance_heatmap.pdf", 86, 92)
print(heat)
dev.off()

cdf_palette <- c(
  MP_1 = "#8F8F8F",
  MP_2 = "#C44E52",
  MP_3 = "#6B8BA4",
  MP_4 = "#E0A458",
  MP_5 = "#7A6FA8",
  MP_6 = "#4E9A7A",
  MP_7 = "#A66A5A"
)

cell_cdf <- cell_level[!is.na(nearest_CD4_CXCL13_distance), .(
  distance = sort(nearest_CD4_CXCL13_distance),
  cumulative_fraction = seq_len(.N) / .N,
  n_cells = .N
), by = MP]
cell_cdf[, MP := factor(MP, levels = mp_levels)]
fwrite(cell_cdf, file.path(out_dir, "figure7_cell_level_cumulative_distribution.csv"))

cell_cdf_plot <- ggplot(cell_cdf, aes(distance, cumulative_fraction, colour = MP, linewidth = MP)) +
  geom_step(lineend = "butt") +
  scale_colour_manual(values = cdf_palette) +
  scale_linewidth_manual(values = c(MP_1 = 0.35, MP_2 = 0.75, MP_3 = 0.35, MP_4 = 0.35, MP_5 = 0.35, MP_6 = 0.35, MP_7 = 0.35), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, quantile(cell_level$nearest_CD4_CXCL13_distance, 0.985, na.rm = TRUE))) +
  labs(
    title = "Cell-level cumulative distance distribution",
    x = "Nearest distance to CD4_CXCL13 cell",
    y = "Cumulative fraction of cells",
    colour = "Malignant MP"
  ) +
  theme_nc() +
  theme(legend.position = "right")
mm_pdf("figure7_supp_candidate_C_cell_level_cumulative_distribution.pdf", 96, 70)
print(cell_cdf_plot)
dev.off()

legend_text <- c(
  "Figure 7 candidate legends",
  "",
  "Main candidate A. Distance from malignant-cell metaprograms to CD4_CXCL13 T cells. Each dot is one pretreatment sample, summarized by the median nearest-neighbor distance from cells in the indicated malignant-cell metaprogram to the closest CD4_CXCL13 T cell. Boxes show median and interquartile range. Global P value was calculated by Friedman test across metaprograms.",
  "",
  "Main candidate B. Ranked sample-level proximity of malignant-cell metaprograms to CD4_CXCL13 T cells. Points show the median of sample-level median distances and whiskers show interquartile range across samples. Lower values indicate closer proximity.",
  "",
  "Main candidate C. Sample-paired MP2 proximity advantage. For each sample, the median distance for MP2 was subtracted from the median distance for each other malignant-cell metaprogram. Values above zero indicate that MP2 is closer to CD4_CXCL13 T cells than the comparison metaprogram. FDR values are Benjamini-Hochberg-adjusted paired Wilcoxon P values.",
  "",
  "Supplementary candidate A. Cell-level nearest-neighbor distance distributions. Boxes summarize all malignant cells by metaprogram; the y-axis is truncated at the 98.5th percentile for readability. Statistical interpretation should prioritize the sample-level analysis.",
  "",
  "Supplementary candidate B. Heatmap of sample-level median nearest-neighbor distances from each malignant-cell metaprogram to CD4_CXCL13 T cells. Samples are ordered by MP2 distance.",
  "",
  "Supplementary candidate C. Cell-level cumulative distribution of nearest-neighbor distances. This is useful as a visual reference, but statistical interpretation should prioritize sample-level analyses."
)
writeLines(legend_text, file.path(legend_dir, "figure7_mp_cd4cxcl13_distance_candidate_legends.txt"))

message("Saved candidate PDFs to: ", out_dir)
message("Saved recomputed sample-level stats to: ", file.path(out_dir, "figure7_recomputed_sample_level_pairwise_stats.csv"))
