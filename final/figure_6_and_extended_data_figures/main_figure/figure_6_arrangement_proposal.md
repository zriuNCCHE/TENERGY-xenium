# Figure 6 Arrangement Proposal

## Working Title

Pretreatment malignant-cell matrix-remodeling program marks non-cCR tumors and is reinforced after chemo-radiotherapy.

## Current Data-Driven Message

The clearest integrated signal is `MP_7`, annotated as basal / invasive matrix remodeling. Among 27 MP_7 genes:

- 20 are higher in pretreatment non-cCR malignant cells versus cCR.
- 20 are higher in postC non-cCR malignant cells versus pretreatment non-cCR.
- 16 are shared across both directions.

Key shared genes include `LAMB3`, `LAMC2`, `COL17A1`, `LAMA3`, `ITGB4`, `TNC`, `PLAU`, `INHBA`, `ITGA3`, `NDRG1`, `MYH9`, `SLC2A1`, `SERPINE1`, `SLC7A5`, `MMP14`, and `CDKN1A`.

This supports a coherent laminin-332 / hemidesmosome / invasive ECM-remodeling axis:

- Laminin-332 chains: `LAMA3`, `LAMB3`, `LAMC2`
- Epithelial adhesion and anchoring: `COL17A1`, `ITGA6`, `ITGB4`, `PLEC`
- Invasion and remodeling: `MMP14`, `PLAU`, `SERPINE1`, `TNC`, `INHBA`
- Stress and resistant-cell state: `NDRG1`, `SLC2A1`, `CDKN1A`

## Recommended Main Figure 6 Layout

### Figure 6A

Use the existing metaprogram similarity/correlation heatmap:

`metaprogram_similarity_heatmap/outputs/figure_6a_metaprogram_similarity_heatmap.pdf`

Purpose: introduce malignant-cell metaprograms and show that MP_7 is a distinct program.

### Figure 6B

Use a metaprogram-DEG overlap heatmap rather than a volcano plot:

`metaprogram_deg_overlap/outputs/exploratory_metaprogram_deg_overlap_heatmap.pdf`

Purpose: show that MP_7 is the strongest overlap between metaprograms, pretreatment non-cCR DEG, and non-cCR postC DEG.

### Figure 6C

Use an MP_7 candidate-gene evidence matrix:

`metaprogram_deg_overlap/outputs/exploratory_mp7_gene_evidence_matrix.pdf`

Purpose: highlight the interpretable genes driving the story: `LAMA3`, `LAMB3`, `LAMC2`, `COL17A1`, `MMP14`, plus related adhesion/remodeling genes.

### Figure 6D

Recommended next panel: patient-level pseudobulk or sample-level expression plot for selected MP_7 genes.

Suggested genes:

`COL17A1`, `LAMA3`, `LAMB3`, `LAMC2`, `MMP14`, `ITGA6`, `ITGB4`, `PLAU`, `SERPINE1`, `TNC`

Preferred design:

- X-axis: `pre cCR`, `pre non-cCR`, `postC non-cCR`
- Y-axis: aggregated expression or module score per patient/sample
- Point: patient/sample
- Summary: median and interquartile range

Purpose: avoid over-reliance on cell-level P values and show patient-level consistency.

### Figure 6E

Recommended next panel: spatial validation map or representative Xenium transcript/gene-score map.

Options:

- MP_7 score over malignant cells in representative cCR and non-cCR samples.
- Transcript maps for `COL17A1`, `LAMA3`/`LAMB3`/`LAMC2`, and `MMP14`.
- Overlay with tumor boundary or malignant-cell compartment.

Purpose: make the Xenium-specific value visible; this is stronger than showing only DEG statistics.

## Where Volcano Plots Fit

The current volcano plots are useful as screening/supporting panels, but they should probably not be the main evidence because cell-level DEG yields extremely small P values. If used, place them in an extended figure or supplementary figure.

## Current Exploratory Outputs

- `metaprogram_deg_overlap/results/metaprogram_deg_gene_evidence.csv`
- `metaprogram_deg_overlap/results/metaprogram_deg_overlap_summary.csv`
- `metaprogram_deg_overlap/results/mp7_matrix_remodeling_gene_evidence.csv`
- `metaprogram_deg_overlap/results/shared_pre_nonccr_postc_metaprogram_genes.csv`
