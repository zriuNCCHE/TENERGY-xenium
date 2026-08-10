# Figure 7 MP2-CD4_CXCL13 Distance Candidates

This folder contains candidate plots generated from server-exported nearest-neighbor distance metrics.

## Source Data

The source CSV files are preserved in:

`data/geo_ready/figure_7_mp_cd4cxcl13_distance/`

## Plot Logic

The preferred main-figure logic is sample-level:

1. For each sample and malignant-cell metaprogram, calculate the median nearest distance from cells in that metaprogram to the closest CD4_CXCL13 T cell.
2. Plot one value per sample per metaprogram.
3. Test global differences across metaprograms using a Friedman test.
4. Test MP2 versus each other MP using paired Wilcoxon signed-rank tests across samples, followed by Benjamini-Hochberg correction.

This avoids treating thousands of cells from the same patient as independent biological replicates.

## Candidate Recommendation

Recommended main panel:

`outputs/figure7_main_candidate_A_sample_median_distance_boxplot.pdf`

Alternative concise main panel:

`outputs/figure7_main_candidate_C_mp2_closeness_delta.pdf`

Recommended supplementary panels:

`outputs/figure7_supp_candidate_A_cell_level_distance_boxplot_zoomed.pdf`

`outputs/figure7_supp_candidate_B_sample_level_distance_heatmap.pdf`

`outputs/figure7_supp_candidate_C_cell_level_cumulative_distribution.pdf`
