# Immune Subset Proportion

This final Figure 3 part plots selected immune-cell subsets as proportions of fibro-immune cells from `data/raw/combined_final_all_cell_types_obs.csv`.

Denominator:

- Fibro-immune cells in each sample, defined as all cells except `A2ML1+ epi` and `malignant`.

Plotted subsets:

- `Macro_CXCL5`
- `Mono_CDC27`
- `Mono_SLC2A3`
- `Neutrophil`

Outputs:

- `outputs/immune_subset_proportion_summary.csv`: sample-level subset percentages among fibro-immune cells.
- `outputs/immune_subset_proportion_unadjusted_stats.csv`: unadjusted paired Wilcoxon tests.
- `outputs/immune_subset_proportion_by_response.pdf`: vector figure.
