# Macrophage Neighborhood Source Report

Source analysis folder:

`exploratory/figure_spatial/macrophage_neighborhood_abundance/`

Source data:

`data/geo_ready/spatial_neighborhoods/postC_k7_cell_niche_assignments.csv`

Logic:

Each post-chemoradiotherapy cell was assigned to one unique spatial neighborhood by selecting the highest score among `neighborhood_1_k7` through `neighborhood_7_k7`. For each patient and neighborhood, macrophage abundance was calculated as the number of `Macro_CXCL9` or `Macro_CXCL5` cells divided by the total number of cells assigned to that neighborhood.

Main source plots:

- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macro_cxcl9_neighborhood_cell_fraction_boxplot.pdf`
- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macro_cxcl5_neighborhood_cell_fraction_boxplot.pdf`
- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macrophage_neighborhood_cell_fraction_boxplots.pdf`

Source tables:

- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macrophage_neighborhood_abundance_by_patient.csv`
- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macrophage_neighborhood_abundance_summary.csv`
- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macrophage_neighborhood_global_friedman_stats.csv`
- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macrophage_neighborhood_pairwise_cell_percent_pvalues.csv`
- `exploratory/figure_spatial/macrophage_neighborhood_abundance/outputs/postc_macrophage_neighborhood_pairwise_area_percent_pvalues.csv`

Summary:

Both macrophage subsets differ across niches.

- `Macro_CXCL9` cell fraction: Friedman P = `5.41e-06`, BH-adjusted P = `5.41e-06`; highest median abundance in `Niche_Macro_CXCL9`.
- `Macro_CXCL5` cell fraction: Friedman P = `5.36e-08`, BH-adjusted P = `1.07e-07`; highest median abundance in `Niche_iCAF_CXCL5`.

Single-feature plots order niches from least to most abundant by median cell fraction.
