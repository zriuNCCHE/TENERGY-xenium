suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

input_path <- file.path(project_root, "data", "raw", "figure_6_deg", "non-cCR_postC_vs_pre.csv")
out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
results_dir <- normalizePath(file.path(dirname(script_path), "..", "results"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

lfc_cutoff <- 0.25
padj_cutoff <- 0.05
label_genes <- c("COL17A1", "MMP14", "LAMA3", "LAMB2", "LAMC3", "LAMB3")

deg <- read.csv(input_path, check.names = FALSE)
deg <- deg[deg$group == "postC", c("names", "scores", "logfoldchanges", "pvals", "pvals_adj")]
names(deg) <- c("gene", "score", "log2fc_postc_vs_pre", "p_value", "p_adj")
deg$p_adj_plot <- pmax(deg$p_adj, 1e-300)
deg$neg_log10_p_adj <- -log10(deg$p_adj_plot)
deg$direction <- "Not significant"
deg$direction[deg$p_adj < padj_cutoff & deg$log2fc_postc_vs_pre >= lfc_cutoff] <- "Higher postC"
deg$direction[deg$p_adj < padj_cutoff & deg$log2fc_postc_vs_pre <= -lfc_cutoff] <- "Higher pre"
deg$is_candidate <- deg$gene %in% label_genes
deg$label <- ifelse(deg$is_candidate, deg$gene, "")

write.csv(deg, file.path(results_dir, "nonccr_postc_vs_pre_volcano_table.csv"), row.names = FALSE)

summary_df <- data.frame(
  contrast = "non-cCR postC vs pre",
  canonical_group = "postC",
  positive_log2fc = "Higher in post-chemo-radiotherapy non-cCR malignant cells",
  lfc_cutoff = lfc_cutoff,
  padj_cutoff = padj_cutoff,
  n_genes = nrow(deg),
  n_higher_postc = sum(deg$direction == "Higher postC"),
  n_higher_pre = sum(deg$direction == "Higher pre"),
  n_candidate_detected = sum(deg$is_candidate),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(results_dir, "nonccr_postc_vs_pre_summary.csv"), row.names = FALSE)

candidate_df <- deg[deg$is_candidate, c("gene", "log2fc_postc_vs_pre", "p_value", "p_adj", "direction")]
write.csv(candidate_df, file.path(results_dir, "nonccr_postc_vs_pre_candidate_genes.csv"), row.names = FALSE)

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
title_size <- plotting_config$font$title_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

timepoint_colors <- unlist(palette_config$timepoints)
plot_colors <- c(
  "Higher pre" = timepoint_colors[["pre"]],
  "Higher postC" = timepoint_colors[["postC"]],
  "Not significant" = "#C7C7C7"
)

x_limit <- ceiling(max(abs(deg$log2fc_postc_vs_pre), na.rm = TRUE) * 1.08)
y_limit <- min(300, ceiling(max(deg$neg_log10_p_adj, na.rm = TRUE) * 1.04))
label_df <- deg[deg$label != "", ]
label_df <- label_df[order(-label_df$log2fc_postc_vs_pre), ]
label_df$label_x <- pmin(label_df$log2fc_postc_vs_pre + 0.45, x_limit - 0.4)
label_df$label_y <- y_limit - seq(10, by = 16, length.out = nrow(label_df))

plot <- ggplot(deg, aes(x = log2fc_postc_vs_pre, y = neg_log10_p_adj)) +
  geom_point(aes(color = direction), size = 0.55, alpha = 0.70, stroke = 0) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "#666666", linewidth = line_width) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "#666666", linewidth = line_width) +
  geom_point(
    data = deg[deg$is_candidate, ],
    aes(x = log2fc_postc_vs_pre, y = neg_log10_p_adj),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "#222222",
    size = 1.6,
    stroke = line_width
  ) +
  geom_segment(
    data = label_df,
    aes(x = log2fc_postc_vs_pre, y = neg_log10_p_adj, xend = label_x, yend = label_y),
    inherit.aes = FALSE,
    color = "#444444",
    linewidth = line_width
  ) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE,
    size = 2.0,
    family = plotting_config$font$family,
    color = "#222222",
    hjust = 0
  ) +
  scale_color_manual(values = plot_colors, breaks = c("Higher postC", "Higher pre", "Not significant")) +
  scale_x_continuous(limits = c(-x_limit, x_limit), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Longitudinal non-cCR malignant-cell DEG",
    x = expression(log[2]~"fold change (postC / pre)"),
    y = expression(-log[10]~"adjusted P"),
    color = NULL
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = small_size),
    legend.position = "top",
    legend.text = element_text(size = small_size),
    legend.key.width = unit(0.32, "cm"),
    axis.line = element_line(linewidth = line_width),
    axis.ticks = element_line(linewidth = line_width),
    panel.border = element_rect(color = "#333333", fill = NA, linewidth = line_width),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.4, alpha = 0.9)))

ggsave(
  filename = file.path(out_dir, "figure_6c_nonccr_postc_vs_pre_volcano.pdf"),
  plot = plot,
  device = "pdf",
  width = 85,
  height = 82,
  units = "mm",
  dpi = 300
)
