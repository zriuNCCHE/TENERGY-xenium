suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  normalizePath("final/figure_6/metaprogram_top50_deg_overlap/code/analyze_top50_deg_overlap.R")
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."))

out_dir <- file.path(project_root, "final", "figure_6", "metaprogram_top50_deg_overlap", "outputs")
results_dir <- file.path(project_root, "final", "figure_6", "metaprogram_top50_deg_overlap", "results")
legend_dir <- file.path(project_root, "final", "figure_6", "metaprogram_top50_deg_overlap", "legends")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(legend_dir, recursive = TRUE, showWarnings = FALSE)

plotting_config <- yaml::read_yaml(file.path(project_root, "config", "plotting.yaml"))
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))
mp_palette <- unlist(palette_config$malignant_metaprograms)

mp_genes <- read.csv(
  file.path(project_root, "final", "figure_6", "metaprogram_similarity_heatmap", "results", "xenium_metaprogram_genes.csv"),
  check.names = FALSE
)
mp_anno <- read.csv(
  file.path(project_root, "final", "figure_6", "metaprogram_similarity_heatmap", "results", "xenium_metaprogram_annotations.csv"),
  check.names = FALSE
)

deg_dir <- file.path(project_root, "data", "raw", "figure_5_deg")
pre_deg <- read.csv(file.path(deg_dir, "pre-treatment-cCR_vs_non-cCR.csv"), check.names = FALSE)
post_deg <- read.csv(file.path(deg_dir, "non-cCR_postC_vs_pre.csv"), check.names = FALSE)

top_n <- 50
mp_levels <- paste0("MP_", 1:7)
mp_genes$meta_program <- factor(mp_genes$meta_program, levels = mp_levels)

get_top <- function(deg, group_name, top_n = 50) {
  x <- deg[deg[["group"]] == group_name, ]
  x <- x[order(-x[["scores"]], x[["pvals_adj"]]), ]
  x <- x[seq_len(min(nrow(x), top_n)), ]
  x$rank_in_deg <- seq_len(nrow(x))
  x
}

top_sets <- list(
  pre_non_cCR_top50 = get_top(pre_deg, "non-cCR", top_n),
  pre_cCR_top50 = get_top(pre_deg, "cCR", top_n),
  noncCR_postC_top50 = get_top(post_deg, "postC", top_n),
  noncCR_pre_top50 = get_top(post_deg, "pre", top_n)
)

top_gene_sets <- lapply(top_sets, function(x) unique(x$names))
resistance_shared_top50 <- intersect(top_gene_sets$pre_non_cCR_top50, top_gene_sets$noncCR_postC_top50)

contrast_table <- data.frame(
  contrast = names(top_gene_sets),
  comparison_label = c(
    "Pre non-cCR high",
    "Pre cCR high",
    "non-cCR postC high",
    "non-cCR pre high"
  ),
  stringsAsFactors = FALSE
)

overlap_rows <- list()
gene_rows <- list()
for (mp in mp_levels) {
  genes <- unique(mp_genes$gene[mp_genes$meta_program == mp])
  label <- mp_anno$suggested_label[match(mp, mp_anno$meta_program)]
  for (contrast in names(top_gene_sets)) {
    overlap <- intersect(genes, top_gene_sets[[contrast]])
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      meta_program = mp,
      suggested_label = label,
      contrast = contrast,
      n_metaprogram_genes = length(genes),
      n_top_deg = length(top_gene_sets[[contrast]]),
      n_overlap = length(overlap),
      overlap_fraction_of_mp = length(overlap) / length(genes),
      overlap_genes = paste(overlap, collapse = ", "),
      stringsAsFactors = FALSE
    )
    if (length(overlap) > 0) {
      deg_source <- top_sets[[contrast]]
      gene_rows[[length(gene_rows) + 1]] <- data.frame(
        meta_program = mp,
        suggested_label = label,
        contrast = contrast,
        gene = overlap,
        mp_rank = mp_genes$rank[match(overlap, mp_genes$gene[mp_genes$meta_program == mp])],
        deg_rank = deg_source$rank_in_deg[match(overlap, deg_source$names)],
        deg_score = deg_source$scores[match(overlap, deg_source$names)],
        logfoldchanges = deg_source$logfoldchanges[match(overlap, deg_source$names)],
        pvals_adj = deg_source$pvals_adj[match(overlap, deg_source$names)],
        stringsAsFactors = FALSE
      )
    }
  }
  shared_overlap <- intersect(genes, resistance_shared_top50)
  overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
    meta_program = mp,
    suggested_label = label,
    contrast = "shared_pre_non_cCR_and_postC_top50",
    n_metaprogram_genes = length(genes),
    n_top_deg = length(resistance_shared_top50),
    n_overlap = length(shared_overlap),
    overlap_fraction_of_mp = length(shared_overlap) / length(genes),
    overlap_genes = paste(shared_overlap, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

overlap_summary <- do.call(rbind, overlap_rows)
gene_evidence <- if (length(gene_rows) > 0) do.call(rbind, gene_rows) else data.frame()

plot_summary <- overlap_summary[overlap_summary$contrast %in% c(
  "pre_non_cCR_top50",
  "noncCR_postC_top50",
  "shared_pre_non_cCR_and_postC_top50"
), ]
plot_summary$contrast_label <- factor(
  plot_summary$contrast,
  levels = c("pre_non_cCR_top50", "noncCR_postC_top50", "shared_pre_non_cCR_and_postC_top50"),
  labels = c("Pre non-cCR\ntop 50", "non-cCR postC\ntop 50", "Shared\nboth top 50")
)
plot_summary$mp_label <- paste0(plot_summary$meta_program, "\n", plot_summary$suggested_label)
plot_summary$mp_label <- factor(
  plot_summary$mp_label,
  levels = rev(unique(plot_summary$mp_label[order(plot_summary$meta_program)]))
)

write.csv(overlap_summary, file.path(results_dir, "metaprogram_top50_deg_overlap_summary.csv"), row.names = FALSE)
write.csv(gene_evidence, file.path(results_dir, "metaprogram_top50_deg_overlap_gene_evidence.csv"), row.names = FALSE)
write.csv(
  data.frame(gene = resistance_shared_top50, stringsAsFactors = FALSE),
  file.path(results_dir, "shared_pre_nonccr_and_postc_top50_genes.csv"),
  row.names = FALSE
)

base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
title_size <- plotting_config$font$title_pt
pt_to_mm <- 0.352777778
line_width <- plotting_config$line$standard_pt * pt_to_mm

overlap_plot <- ggplot(plot_summary, aes(contrast_label, mp_label)) +
  geom_tile(aes(fill = n_overlap), colour = "white", linewidth = line_width) +
  geom_text(aes(label = n_overlap), size = 2.2, family = plotting_config$font$family) +
  scale_fill_gradient(low = "#F5F5F5", high = mp_palette[["MP_7"]], name = "Genes") +
  labs(
    title = "Top-50 DEG overlap with malignant metaprograms",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, face = "bold", hjust = 0),
    axis.text.x = element_text(size = small_size),
    axis.text.y = element_text(size = small_size),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

ggsave(
  file.path(out_dir, "figure_6_metaprogram_top50_deg_overlap.pdf"),
  overlap_plot,
  width = 95,
  height = 72,
  units = "mm",
  device = "pdf",
  dpi = 300
)

mp7_gene_plot_df <- gene_evidence[
  gene_evidence$meta_program == "MP_7" &
    gene_evidence$contrast %in% c("pre_non_cCR_top50", "noncCR_postC_top50"),
]
mp7_gene_plot_df$contrast_label <- factor(
  mp7_gene_plot_df$contrast,
  levels = c("pre_non_cCR_top50", "noncCR_postC_top50"),
  labels = c("Pre non-cCR top 50", "non-cCR postC top 50")
)
mp7_gene_order <- unique(mp_genes$gene[mp_genes$meta_program == "MP_7"])
mp7_gene_plot_df$gene <- factor(mp7_gene_plot_df$gene, levels = rev(mp7_gene_order))

mp7_gene_plot <- ggplot(mp7_gene_plot_df, aes(contrast_label, gene)) +
  geom_point(aes(size = pmin(-log10(pmax(pvals_adj, 1e-300)), 50), fill = logfoldchanges), shape = 21, colour = "#333333", stroke = line_width) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#C52B2F", midpoint = 0, name = expression(log[2]~FC)) +
  scale_size_continuous(range = c(1.1, 4.2), name = expression(-log[10]~adj.~P), breaks = c(5, 25, 50)) +
  labs(
    title = "MP7 genes in top-50 DEG lists",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_family = plotting_config$font$family, base_size = base_size) +
  theme(
    plot.title = element_text(size = title_size, face = "bold", hjust = 0),
    axis.text.x = element_text(size = small_size),
    axis.text.y = element_text(size = small_size),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = small_size),
    legend.text = element_text(size = small_size),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

ggsave(
  file.path(out_dir, "figure_6_mp7_top50_deg_gene_evidence.pdf"),
  mp7_gene_plot,
  width = 82,
  height = 88,
  units = "mm",
  device = "pdf",
  dpi = 300
)

legend_text <- c(
  "Figure 6 top-50 DEG overlap legend draft",
  "",
  paste0(
    "Top-50 DEG overlap between malignant-cell metaprograms and DEG lists. ",
    "Top genes were selected separately from pretreatment non-cCR-enriched genes and non-cCR post-chemoradiotherapy-enriched genes by absolute Scanpy rank_genes_groups score. ",
    "Numbers indicate genes shared between each metaprogram gene set and the indicated top-50 DEG list."
  ),
  "",
  "The shared-both-top-50 column counts metaprogram genes present in both the pretreatment non-cCR top-50 list and the non-cCR postC top-50 list."
)
writeLines(legend_text, file.path(legend_dir, "legend_draft.md"))

readme_text <- c(
  "# Metaprogram Top-50 DEG Overlap",
  "",
  "Checks whether MP7 is preferentially represented among top-ranked differentially expressed genes.",
  "",
  "## Inputs",
  "",
  "- `final/figure_6/metaprogram_similarity_heatmap/results/xenium_metaprogram_genes.csv`",
  "- `data/raw/figure_5_deg/pre-treatment-cCR_vs_non-cCR.csv`",
  "- `data/raw/figure_5_deg/non-cCR_postC_vs_pre.csv`",
  "",
  "## Outputs",
  "",
  "- `outputs/figure_6_metaprogram_top50_deg_overlap.pdf`",
  "- `outputs/figure_6_mp7_top50_deg_gene_evidence.pdf`",
  "- `results/metaprogram_top50_deg_overlap_summary.csv`",
  "- `results/metaprogram_top50_deg_overlap_gene_evidence.csv`"
)
writeLines(readme_text, file.path(project_root, "final", "figure_6", "metaprogram_top50_deg_overlap", "README.md"))

message("Saved top-50 DEG overlap analysis to: ", results_dir)
