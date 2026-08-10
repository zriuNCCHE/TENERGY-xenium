# Metaprogram DEG Overlap

Exploratory Figure 6 analysis integrating malignant-cell metaprogram genes with the two DEG contrasts.

## Inputs

- `final/figure_6/metaprogram_similarity_heatmap/results/xenium_metaprogram_genes.csv`
- `final/figure_6/metaprogram_similarity_heatmap/results/xenium_metaprogram_annotations.csv`
- `data/raw/figure_6_deg/pre-treatment-cCR_vs_non-cCR.csv`
- `data/raw/figure_6_deg/non-cCR_postC_vs_pre.csv`

## Contrast Direction

- Pretreatment DEG uses the `non-cCR` rows, so positive log2FC means higher in pretreatment non-cCR malignant cells than cCR malignant cells.
- Longitudinal non-cCR DEG uses the `postC` rows, so positive log2FC means higher in post-chemo-radiotherapy non-cCR malignant cells than pretreatment non-cCR malignant cells.

## Outputs

- `results/metaprogram_deg_gene_evidence.csv`: all metaprogram genes joined to both DEG contrasts.
- `results/metaprogram_deg_overlap_summary.csv`: overlap counts per metaprogram.
- `results/mp7_matrix_remodeling_gene_evidence.csv`: focused evidence table for MP_7.
- `results/shared_pre_nonccr_postc_metaprogram_genes.csv`: metaprogram genes higher in pretreatment non-cCR and higher after chemo-radiotherapy in non-cCR.
- `outputs/exploratory_metaprogram_deg_overlap_heatmap.pdf`: metaprogram-level overlap count heatmap.
- `outputs/exploratory_mp7_gene_evidence_matrix.pdf`: candidate MP_7 gene evidence matrix.
