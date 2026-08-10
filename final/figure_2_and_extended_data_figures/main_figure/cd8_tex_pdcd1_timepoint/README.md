# CD8 Tex PDCD1 Proportion Across Treatment

This final Figure 3 part plots `CD8_Tex_PDCD1` proportions among refined T cells across treatment timepoints in cCR and non-cCR samples.

Denominator:

- Refined T cells: `CD8_Teff`, `CD8_Tex_PDCD1`, `CD4_Treg_FOXP3`, `CD4_Treg_CCR8`, `CD8_prolif`, `CD4_CXCL13`, and `CD4_prolif`.
- The generic `T/NK` bucket is excluded.

Outputs:

- `outputs/cd8_tex_pdcd1_timepoint_summary.csv`: sample-level percentages.
- `outputs/cd8_tex_pdcd1_timepoint_unadjusted_stats.csv`: unadjusted paired Wilcoxon tests.
- `outputs/cd8_tex_pdcd1_timepoint.pdf`: vector figure.

