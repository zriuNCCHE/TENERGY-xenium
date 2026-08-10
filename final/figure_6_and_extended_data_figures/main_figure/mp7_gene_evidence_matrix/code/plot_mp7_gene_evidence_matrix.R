suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", "..", ".."))

out_dir <- normalizePath(file.path(dirname(script_path), "..", "outputs"), mustWork = FALSE)
results_dir <- normalizePath(file.path(dirname(script_path), "..", "results"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))

mp_genes <- read.csv(
  file.path(project_root, "final", "figure_6", "metaprogram_similarity_heatmap", "results", "xenium_metaprogram_genes.csv"),
  check.names = FALSE
)
mp_anno <- read.csv(
  file.path(project_root, "final", "figure_6", "metaprogram_similarity_heatmap", "results", "xenium_metaprogram_annotations.csv"),
  check.names = FALSE
)

pre_deg <- read.csv(
  file.path(project_root, "data", "raw", "figure_6_deg", "pre-treatment-cCR_vs_non-cCR.csv"),
  check.names = FALSE
)
post_deg <- read.csv(
  file.path(project_root, "data", "raw", "figure_6_deg", "non-cCR_postC_vs_pre.csv"),
  check.names = FALSE
)

lfc_cutoff <- 0.25
padj_cutoff <- 0.05

pre <- pre_deg[pre_deg$group == "non-cCR", c("names", "scores", "logfoldchanges", "pvals", "pvals_adj")]
names(pre) <- c("gene", "pre_score", "pre_log2fc_nonccr_vs_ccr", "pre_p_value", "pre_p_adj")
pre$pre_direction <- "not_significant"
pre$pre_direction[pre$pre_p_adj < padj_cutoff & pre$pre_log2fc_nonccr_vs_ccr >= lfc_cutoff] <- "higher_non_cCR"
pre$pre_direction[pre$pre_p_adj < padj_cutoff & pre$pre_log2fc_nonccr_vs_ccr <= -lfc_cutoff] <- "higher_cCR"

post <- post_deg[post_deg$group == "postC", c("names", "scores", "logfoldchanges", "pvals", "pvals_adj")]
names(post) <- c("gene", "post_score", "post_log2fc_postC_vs_pre", "post_p_value", "post_p_adj")
post$post_direction <- "not_significant"
post$post_direction[post$post_p_adj < padj_cutoff & post$post_log2fc_postC_vs_pre >= lfc_cutoff] <- "higher_postC"
post$post_direction[post$post_p_adj < padj_cutoff & post$post_log2fc_postC_vs_pre <= -lfc_cutoff] <- "higher_pre"

evidence <- merge(mp_genes, pre, by = "gene", all.x = TRUE)
evidence <- merge(evidence, post, by = "gene", all.x = TRUE)
evidence <- merge(evidence, mp_anno[, c("meta_program", "suggested_label")], by = "meta_program", all.x = TRUE)
evidence <- evidence[order(evidence$meta_program, evidence$rank), ]

evidence$pre_neg_log10_p_adj <- -log10(pmax(evidence$pre_p_adj, 1e-300))
evidence$post_neg_log10_p_adj <- -log10(pmax(evidence$post_p_adj, 1e-300))
evidence$shared_signal <- evidence$pre_direction == "higher_non_cCR" & evidence$post_direction == "higher_postC"
evidence$pretreatment_nonccr_feature <- evidence$pre_direction == "higher_non_cCR"
evidence$post_treatment_enriched <- evidence$post_direction == "higher_postC"

write.csv(
  evidence,
  file.path(results_dir, "metaprogram_deg_gene_evidence.csv"),
  row.names = FALSE
)

summary_df <- do.call(
  rbind,
  lapply(split(evidence, evidence$meta_program), function(df) {
    data.frame(
      meta_program = df$meta_program[1],
      suggested_label = df$suggested_label[1],
      n_metaprogram_genes = nrow(df),
      n_pre_higher_nonccr = sum(df$pre_direction == "higher_non_cCR", na.rm = TRUE),
      n_pre_higher_ccr = sum(df$pre_direction == "higher_cCR", na.rm = TRUE),
      n_post_higher_postC = sum(df$post_direction == "higher_postC", na.rm = TRUE),
      n_post_higher_pre = sum(df$post_direction == "higher_pre", na.rm = TRUE),
      n_shared_pre_nonccr_and_postC = sum(df$shared_signal, na.rm = TRUE),
      shared_genes = paste(df$gene[df$shared_signal], collapse = ", "),
      pre_nonccr_genes = paste(df$gene[df$pre_direction == "higher_non_cCR"], collapse = ", "),
      postC_genes = paste(df$gene[df$post_direction == "higher_postC"], collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
)
summary_df <- summary_df[order(-summary_df$n_shared_pre_nonccr_and_postC, -summary_df$n_pre_higher_nonccr), ]
write.csv(
  summary_df,
  file.path(results_dir, "metaprogram_deg_overlap_summary.csv"),
  row.names = FALSE
)

mp7 <- evidence[evidence$meta_program == "MP_7", ]
mp7 <- mp7[
  mp7$pre_direction %in% c("higher_non_cCR", "higher_cCR") |
    mp7$post_direction %in% c("higher_postC", "higher_pre") |
    mp7$rank <= 15,
]
mp7 <- mp7[order(mp7$rank), ]
write.csv(
  mp7,
  file.path(results_dir, "mp7_matrix_remodeling_gene_evidence.csv"),
  row.names = FALSE
)

top_shared <- evidence[evidence$shared_signal, ]
top_shared <- top_shared[order(top_shared$meta_program, top_shared$rank), ]
write.csv(
  top_shared,
  file.path(results_dir, "shared_pre_nonccr_postc_metaprogram_genes.csv"),
  row.names = FALSE
)

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
title_size <- plotting_config$font$title_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

plot_summary <- summary_df
plot_summary$meta_label <- paste0(plot_summary$meta_program, "\n", plot_summary$suggested_label)
plot_summary$meta_label <- factor(plot_summary$meta_label, levels = rev(plot_summary$meta_label))
plot_df <- rbind(
  data.frame(meta_label = plot_summary$meta_label, contrast = "Pre non-cCR high", count = plot_summary$n_pre_higher_nonccr),
  data.frame(meta_label = plot_summary$meta_label, contrast = "PostC high", count = plot_summary$n_post_higher_postC),
  data.frame(meta_label = plot_summary$meta_label, contrast = "Shared", count = plot_summary$n_shared_pre_nonccr_and_postC)
)
plot_df$contrast <- factor(plot_df$contrast, levels = c("Pre non-cCR high", "PostC high", "Shared"))

overlap_plot <- ggplot(plot_df, aes(x = contrast, y = meta_label, fill = count)) +
  geom_tile(color = "white", linewidth = line_width) +
  geom_text(aes(label = count), size = 2.2, family = plotting_config$font$family) +
  scale_fill_gradient(low = "#F2F2F2", high = "#7A3E9D", name = "Genes") +
  labs(
    title = "Metaprogram-DEG overlap",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = small_size, angle = 35, hjust = 1),
    axis.text.y = element_text(size = small_size),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(out_dir, "exploratory_metaprogram_deg_overlap_heatmap.pdf"),
  overlap_plot,
  width = 90,
  height = 80,
  units = "mm",
  device = "pdf",
  dpi = 300
)

candidate_gene_order <- c(
  "LAMA3", "LAMB3", "LAMC2",
  "COL17A1", "ITGB4", "ITGA3",
  "MMP14", "PLAU", "SERPINE1", "TNC", "INHBA",
  "NDRG1", "MYH9", "SLC2A1", "SLC7A5", "CDKN1A"
)
candidate_genes <- unique(candidate_gene_order)
candidate <- evidence[evidence$meta_program == "MP_7" & evidence$gene %in% candidate_genes, ]
candidate <- candidate[match(candidate_gene_order[candidate_gene_order %in% candidate$gene], candidate$gene), ]
candidate_long <- rbind(
  data.frame(
    gene = candidate$gene,
    rank = candidate$rank,
    contrast = "Pre non-cCR vs cCR",
    log2fc = candidate$pre_log2fc_nonccr_vs_ccr,
    neg_log10_p_adj = candidate$pre_neg_log10_p_adj,
    direction = candidate$pre_direction
  ),
  data.frame(
    gene = candidate$gene,
    rank = candidate$rank,
    contrast = "non-cCR postC vs pre",
    log2fc = candidate$post_log2fc_postC_vs_pre,
    neg_log10_p_adj = candidate$post_neg_log10_p_adj,
    direction = candidate$post_direction
  )
)
candidate_long$gene <- factor(candidate_long$gene, levels = candidate_gene_order[candidate_gene_order %in% candidate$gene])
candidate_long$contrast <- factor(candidate_long$contrast, levels = rev(c("Pre non-cCR vs cCR", "non-cCR postC vs pre")))

evidence_plot <- ggplot(candidate_long, aes(x = gene, y = contrast)) +
  geom_point(aes(size = pmin(neg_log10_p_adj, 50), fill = log2fc), shape = 21, color = "#333333", stroke = line_width) +
  scale_fill_gradient2(low = "#5BB3C8", mid = "white", high = "#E76F51", midpoint = 0, name = expression(log[2]~FC)) +
  scale_size_continuous(range = c(1.2, 4.4), name = expression(-log[10]~adj.~P), breaks = c(5, 25, 50)) +
  labs(
    title = "MP_7 matrix-remodeling genes",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = small_size, angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = base_size, face = "bold"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    plot.margin = margin(4, 4, 12, 4)
  )

ggsave(
  file.path(out_dir, "figure_6c_mp7_gene_evidence_matrix.pdf"),
  evidence_plot,
  width = 135,
  height = 52,
  units = "mm",
  device = "pdf",
  dpi = 300
)
