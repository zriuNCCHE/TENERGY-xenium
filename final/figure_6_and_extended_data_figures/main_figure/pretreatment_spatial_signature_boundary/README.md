# Pretreatment Spatial Signature Boundary Analysis

Purpose: test whether pretreatment malignant-cell MP7-related signatures are enriched at the tumor margin versus tumor inner compartment, and whether that spatial pattern differs between cCR and non-cCR patients.

Primary script:

- `code/analyze_pretreatment_spatial_signature_boundary.py`

Upstream server-side annotation script:

- `code/run_pretreatment_malignant_distance_annotation.py`

This upstream script loads the full combined AnnData object, filters to
pretreatment samples, runs `xenium_malignant_distance`, and exports
`malignant_distance_info.csv` with malignant-region labels such as tumor inner,
tumor margin, isolated malignant, and fibro/immune compartments. The exported
region annotations are then used by the main Figure 6 boundary/signature
analysis scripts in this folder.

Default server workflow: write compact CSV summaries only, then download the `results/` folder for final plotting in Codex.

```python
from pathlib import Path
from analyze_pretreatment_spatial_signature_boundary import run_analysis

res = run_analysis(
    pre_treat_epi_adata,
    output_root=Path("final/figure_6/pretreatment_spatial_signature_boundary"),
    layer=None,
    min_cells_per_region=20,
    make_plots=False,
)
```

Expected AnnData fields:

- `obs["major_cell_type"]`: malignant-cell subset label.
- `obs["sample_timepoint"]`: pretreatment subset label.
- `obs["malignant_region"]`: expected values include `tumor_inner`, `tumor_margin`, and optionally `isolated_malignant`.
- `obs["cCR"]`: clinical response, expected values compatible with `cCR` and `non-cCR`.
- `obs["sampleID"]`: biological replicate unit used for statistics.
- `obs["patientID"]`: retained in summary tables when available.

Main outputs:

- `outputs/figure_6d_signature_region_by_response_boxplot.pdf`: optional sample-level signature plot when `make_plots=True`.
- `outputs/figure_6e_signature_boundary_enrichment_by_response.pdf`: optional boundary-enrichment plot when `make_plots=True`.
- `outputs/figure_6f_gene_region_response_heatmap.pdf`: optional gene summary heatmap when `make_plots=True`.
- `results/pretreatment_malignant_sample_region_summary.csv`: sample-region means and medians for tumor inner and tumor margin only.
- `results/pretreatment_malignant_sample_region_summary_long.csv`: tidy long-format sample-region table for plotting.
- `results/pretreatment_malignant_signature_region_summary_long.csv`: tidy signature-only sample-region table.
- `results/pretreatment_malignant_gene_region_summary_long.csv`: tidy gene-only sample-region table.
- `results/pretreatment_malignant_sample_region_summary_all_regions.csv`: all sample-region means and medians, including isolated malignant cells for reference.
- `results/pretreatment_malignant_boundary_delta_by_sample.csv`: paired margin minus inner values per sample.
- `results/pretreatment_malignant_signature_boundary_delta_long.csv`: tidy signature-only boundary delta table.
- `results/pretreatment_malignant_gene_boundary_delta_long.csv`: tidy gene-only boundary delta table.
- `results/pretreatment_malignant_spatial_signature_statistics.csv`: sample-level statistical tests, raw p-values, BH-adjusted p-values within test family, and mean-difference effect sizes.

Statistical note: the script computes tests from sample-level summaries, not individual cells, to avoid pseudoreplication. Isolated malignant cells are excluded from the main comparisons and retained only in the all-region reference table. The large cell-level table is not written by default; set `write_cell_scores=True` only if needed.
