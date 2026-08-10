suppressPackageStartupMessages({
  library(grid)
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
palette_config <- yaml::read_yaml(file.path(project_root, "config", "palettes.yaml"))

mp7 <- read.csv(
  file.path(project_root, "final", "figure_6", "metaprogram_deg_overlap", "results", "mp7_matrix_remodeling_gene_evidence.csv"),
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

set_mp7 <- unique(mp7$gene)
pre <- pre_deg[pre_deg$group == "non-cCR", ]
post <- post_deg[post_deg$group == "postC", ]
set_pre <- unique(pre$names[pre$logfoldchanges >= lfc_cutoff & pre$pvals_adj < padj_cutoff])
set_post <- unique(post$names[post$logfoldchanges >= lfc_cutoff & post$pvals_adj < padj_cutoff])

abc <- Reduce(intersect, list(set_mp7, set_pre, set_post))
ab_only <- setdiff(intersect(set_mp7, set_pre), set_post)
ac_only <- setdiff(intersect(set_mp7, set_post), set_pre)
bc_only <- setdiff(intersect(set_pre, set_post), set_mp7)
a_only <- setdiff(set_mp7, union(set_pre, set_post))
b_only <- setdiff(set_pre, union(set_mp7, set_post))
c_only <- setdiff(set_post, union(set_mp7, set_pre))

venn_counts <- data.frame(
  region = c("MP7 only", "Pre non-cCR only", "PostC only", "MP7 + pre only", "MP7 + postC only", "Pre + postC only", "MP7 + pre + postC"),
  n_genes = c(length(a_only), length(b_only), length(c_only), length(ab_only), length(ac_only), length(bc_only), length(abc)),
  genes = c(
    paste(a_only, collapse = ", "),
    paste(b_only, collapse = ", "),
    paste(c_only, collapse = ", "),
    paste(ab_only, collapse = ", "),
    paste(ac_only, collapse = ", "),
    paste(bc_only, collapse = ", "),
    paste(abc, collapse = ", ")
  ),
  stringsAsFactors = FALSE
)
write.csv(venn_counts, file.path(results_dir, "mp7_deg_overlap_venn_counts.csv"), row.names = FALSE)

shared_groups <- data.frame(
  module = c(
    rep("Laminin-332", 3),
    rep("Adhesion / anchoring", 3),
    rep("Remodeling / invasion", 5),
    rep("Broader stress / adaptation", 5)
  ),
  gene = c(
    "LAMA3", "LAMB3", "LAMC2",
    "COL17A1", "ITGB4", "ITGA3",
    "MMP14", "PLAU", "SERPINE1", "TNC", "INHBA",
    "NDRG1", "MYH9", "SLC2A1", "SLC7A5", "CDKN1A"
  ),
  stringsAsFactors = FALSE
)
shared_groups <- shared_groups[shared_groups$gene %in% abc, ]
write.csv(shared_groups, file.path(results_dir, "mp7_deg_overlap_venn_shared_groups.csv"), row.names = FALSE)

font_family <- plotting_config$font$family
base_size <- plotting_config$font$standard_pt
small_size <- plotting_config$font$small_pt
title_size <- plotting_config$font$title_pt
line_width <- plotting_config$line$standard_pt

col_mp7 <- "#7A3E9D"
col_pre <- "#E69F00"
col_post <- "#56B4E9"
col_text <- "#2B2B2B"
col_muted <- "#666666"

pdf(
  file = file.path(out_dir, "figure_6b_mp7_deg_overlap_venn.pdf"),
  width = 125 / 25.4,
  height = 90 / 25.4,
  family = font_family,
  useDingbats = FALSE
)

grid.newpage()
grid.rect(gp = gpar(fill = "white", col = NA))

grid.text(
  "B",
  x = unit(0.035, "npc"),
  y = unit(0.955, "npc"),
  gp = gpar(fontfamily = font_family, fontsize = 10, fontface = "bold", col = col_text)
)
grid.text(
  "MP7 overlaps non-cCR and post-treatment DEG",
  x = unit(0.52, "npc"),
  y = unit(0.955, "npc"),
  gp = gpar(fontfamily = font_family, fontsize = title_size, fontface = "bold", col = col_text)
)

draw_circle <- function(x, y, r, fill, col) {
  grid.circle(
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    r = unit(r, "npc"),
    gp = gpar(fill = adjustcolor(fill, alpha.f = 0.22), col = col, lwd = line_width * 1.25)
  )
}

draw_label <- function(x, y, label, size = base_size, face = "plain", col = col_text) {
  grid.text(
    label,
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    gp = gpar(fontfamily = font_family, fontsize = size, fontface = face, col = col, lineheight = 0.95)
  )
}

draw_circle(0.52, 0.71, 0.225, col_mp7, col_mp7)
draw_circle(0.40, 0.52, 0.225, col_pre, col_pre)
draw_circle(0.64, 0.52, 0.225, col_post, col_post)

draw_label(0.52, 0.92, "MP7 genes", base_size, "bold", col_mp7)
draw_label(0.23, 0.44, "Pre non-cCR-high DEG", base_size, "bold", col_pre)
draw_label(0.81, 0.44, "PostC-high DEG", base_size, "bold", col_post)

grid.roundrect(
  x = unit(0.52, "npc"),
  y = unit(0.585, "npc"),
  width = unit(0.18, "npc"),
  height = unit(0.12, "npc"),
  r = unit(1.4, "mm"),
  gp = gpar(fill = "white", col = "#333333", lwd = line_width)
)
draw_label(0.52, 0.610, "shared genes", small_size, "bold", col_muted)
draw_label(0.52, 0.565, paste0("n = ", length(abc)), 15, "bold", col_text)

grid.text(
  "Selected shared genes",
  x = unit(0.50, "npc"),
  y = unit(0.285, "npc"),
  gp = gpar(fontfamily = font_family, fontsize = base_size, fontface = "bold", col = col_text)
)

draw_module <- function(x, y, w, h, module, genes, col) {
  grid.roundrect(
    x = unit(x, "npc"),
    y = unit(y, "npc"),
    width = unit(w, "npc"),
    height = unit(h, "npc"),
    r = unit(1.2, "mm"),
    gp = gpar(fill = adjustcolor(col, alpha.f = 0.12), col = col, lwd = line_width)
  )
  grid.text(
    module,
    x = unit(x, "npc"),
    y = unit(y + h * 0.25, "npc"),
    gp = gpar(fontfamily = font_family, fontsize = small_size, fontface = "bold", col = col)
  )
  grid.text(
    paste(genes, collapse = "  "),
    x = unit(x, "npc"),
    y = unit(y - h * 0.14, "npc"),
    gp = gpar(fontfamily = font_family, fontsize = small_size, fontface = "bold", col = col_text)
  )
}

module_cols <- c(
  "Laminin-332" = "#D55E00",
  "Adhesion / anchoring" = "#0072B2",
  "Remodeling / invasion" = "#009E73",
  "Broader stress / adaptation" = "#6F6F6F"
)
module_genes <- split(shared_groups$gene, shared_groups$module)
draw_module(0.23, 0.190, 0.32, 0.10, "Laminin-332", module_genes[["Laminin-332"]], module_cols[["Laminin-332"]])
draw_module(0.60, 0.190, 0.38, 0.10, "Adhesion / anchoring", module_genes[["Adhesion / anchoring"]], module_cols[["Adhesion / anchoring"]])
draw_module(0.31, 0.085, 0.47, 0.10, "Remodeling / invasion", module_genes[["Remodeling / invasion"]], module_cols[["Remodeling / invasion"]])
draw_module(0.75, 0.085, 0.39, 0.10, "Broader stress", module_genes[["Broader stress / adaptation"]], module_cols[["Broader stress / adaptation"]])

grid.text(
  "Areas are not proportional. DEG threshold: adjusted P < 0.05 and log2FC >= 0.25.",
  x = unit(0.5, "npc"),
  y = unit(0.025, "npc"),
  gp = gpar(fontfamily = font_family, fontsize = small_size, col = col_muted)
)

dev.off()
