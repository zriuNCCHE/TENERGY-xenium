# PostC CCR8+ Treg Proportion By Response

This final Figure 3 part compares postC `CD4_Treg_CCR8` proportions among refined T cells between cCR and non-cCR samples.

Denominator:

- Refined T cells: `CD8_Teff`, `CD8_Tex_PDCD1`, `CD4_Treg_FOXP3`, `CD4_Treg_CCR8`, `CD8_prolif`, `CD4_CXCL13`, and `CD4_prolif`.
- The generic `T/NK` bucket is excluded.

Outputs:

- `outputs/treg_ccr8_postc_by_response_summary.csv`: sample-level percentages.
- `outputs/treg_ccr8_postc_by_response_unadjusted_stats.csv`: unadjusted Wilcoxon rank-sum test.
- `outputs/treg_ccr8_postc_by_response.pdf`: vector figure.

