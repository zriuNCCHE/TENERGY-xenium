# T-Cell Neighborhood Source Report

Source analysis folder:

`exploratory/figure_spatial/t_cell_neighborhood_abundance/`

Source data:

`data/geo_ready/spatial_neighborhoods/postC_k7_cell_niche_assignments.csv`

Logic:

Each post-chemoradiotherapy cell was assigned to one unique spatial neighborhood by selecting the highest score among `neighborhood_1_k7` through `neighborhood_7_k7`. For each patient and neighborhood, T-cell abundance was calculated as the number of cells in the T-cell subset divided by the total number of cells assigned to that neighborhood.

Main source plots:

- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_individual_t_subset_neighborhood_boxplots.pdf`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_combined_t_groups_neighborhood_cell_fraction_boxplots.pdf`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_total_t_cells_neighborhood_cell_fraction_boxplot.pdf`

The boxplot panels display the global paired Friedman P value and BH-adjusted P value for each T-cell feature across the seven niches.

Source tables:

- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_t_cell_neighborhood_abundance_by_patient.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_t_cell_neighborhood_abundance_median_summary.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_t_cell_neighborhood_abundance_global_friedman_stats.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_t_cell_neighborhood_abundance_pairwise_pvalues.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_total_t_cells_neighborhood_abundance_by_patient.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_total_t_cells_neighborhood_abundance_summary.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_total_t_cells_neighborhood_global_friedman_stats.csv`
- `exploratory/figure_spatial/t_cell_neighborhood_abundance/outputs/postc_total_t_cells_neighborhood_pairwise_pvalues.csv`

Summary:

All individual T-cell subsets and the combined T-cell groups vary significantly across spatial neighborhoods in the exploratory paired Friedman tests. The current canonical niche order is `Niche_iCAF_CXCL5`, `Niche_Mono`, `Niche_Neutro`, `Niche_Epi`, `Niche_Macro_CXCL9`, `Niche_iCAF_CXCL6`, `Niche_T_cell`.

Key global statistics:

- T/NK: Friedman P = `6.53e-07`, BH-adjusted P = `5.22e-06`.
- non-Treg T/NK: Friedman P = `4.49e-08`, BH-adjusted P = `8.98e-08`.
- Treg-like: Friedman P = `2.29e-05`, BH-adjusted P = `2.29e-05`.
- Total T cells: Friedman P = `2.43e-08`, BH-adjusted P = `2.43e-08`.

This source report is exploratory because the neighborhoods were derived from cell-type neighborhood information.
