# Figure 7 MP Spatial Proximity

Final candidate right-side panels for Figure 7, designed to sit next to the pretreatment co-occurrence heatmap.

The shared malignant metaprogram palette is saved in `config/palettes.yaml` under `malignant_metaprograms`, and exported as `outputs/figure7_malignant_metaprogram_palette.csv`.

## Main Panels

- `outputs/figure7_main_B_mp2_cd4cxcl13_boxplot.pdf`
- `outputs/figure7_main_C_mp7_myeloid_boxplot.pdf`
- `outputs/figure7_main_D_mp2_cd4cxcl13_cell_cdf.pdf`
- `outputs/figure7_main_E_mp7_myeloid_cell_cdf.pdf`

## Upstream Server Workflow

`code/run_pretreatment_mp_spatial_neighborhood_scan.py` is the upstream Python script for scanning pretreatment spatial neighborhoods after adding malignant-cell MP labels.

The script:

- loads `combined_final_all_cell_types.h5ad`;
- filters to pretreatment samples;
- uses `final_cell_type2_malignant_MP` as the label column;
- computes a spatial graph with a 20 um radius;
- builds a square-root transformed local composition matrix;
- runs NMF stability scans over k = 5 to 12 and grid search over k = 3 to 12;
- writes outputs with the prefix `pre_malignant_added_20um_treatment_neighborhoods`.

This workflow is meant to be run on the server environment where the full AnnData object and `xenium_spatial_neighborhood` package are available.
