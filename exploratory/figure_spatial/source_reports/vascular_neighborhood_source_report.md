# Endothelial and Vascular CAF Neighborhood Source Report

Source analysis folder:

`exploratory/figure_spatial/vascular_neighborhood_abundance/`

Source data:

`data/geo_ready/spatial_neighborhoods/postC_k7_cell_niche_assignments.csv`

Logic:

Each post-chemoradiotherapy cell was assigned to one unique spatial neighborhood by selecting the highest score among `neighborhood_1_k7` through `neighborhood_7_k7`. For each patient and neighborhood, vascular abundance was calculated as the number of `endo_PLVAP` or `vCAF` cells divided by the total number of cells assigned to that neighborhood.

Main source plots:

- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vascular_neighborhood_cell_fraction_boxplots.pdf`
- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_endo_plvap_neighborhood_cell_fraction_boxplot.pdf`
- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vcaf_neighborhood_cell_fraction_boxplot.pdf`

Source tables:

- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vascular_neighborhood_abundance_by_patient.csv`
- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vascular_neighborhood_abundance_summary.csv`
- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vascular_neighborhood_global_friedman_stats.csv`
- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vascular_neighborhood_pairwise_cell_percent_pvalues.csv`
- `exploratory/figure_spatial/vascular_neighborhood_abundance/outputs/postc_vascular_neighborhood_pairwise_area_percent_pvalues.csv`

Summary:

`endo_PLVAP` and `vCAF` abundance differ strongly across spatial neighborhoods. The current canonical niche order is `Niche_iCAF_CXCL5`, `Niche_Mono`, `Niche_Neutro`, `Niche_Epi`, `Niche_Macro_CXCL9`, `Niche_iCAF_CXCL6`, `Niche_T_cell`.

Global statistics:

- endo_PLVAP cell fraction: Friedman P = `1.15e-09`, BH-adjusted P = `1.15e-09`.
- endo_PLVAP area fraction: Friedman P = `2.04e-09`, BH-adjusted P = `2.04e-09`.
- vCAF cell fraction: Friedman P = `1.14e-09`, BH-adjusted P = `1.15e-09`.
- vCAF area fraction: Friedman P = `1.46e-09`, BH-adjusted P = `2.04e-09`.

Both `endo_PLVAP` and `vCAF` are most enriched in `Niche_Macro_CXCL9`, followed by `Niche_iCAF_CXCL6` and `Niche_T_cell`.

Median endo_PLVAP cell fraction by neighborhood:

- Niche_iCAF_CXCL5: 0.19%
- Niche_Mono: 0%
- Niche_Neutro: 0%
- Niche_Epi: 0%
- Niche_Macro_CXCL9: 7.01%
- Niche_iCAF_CXCL6: 3.53%
- Niche_T_cell: 2.05%

Median vCAF cell fraction by neighborhood:

- Niche_iCAF_CXCL5: 0.06%
- Niche_Mono: 0%
- Niche_Neutro: 0%
- Niche_Epi: 0%
- Niche_Macro_CXCL9: 4.05%
- Niche_iCAF_CXCL6: 2.21%
- Niche_T_cell: 1.17%
