# Pretreatment cCR versus non-cCR DEG

Candidate Figure 6B.

This part uses the saved DEG table:

`data/raw/figure_6_deg/pre-treatment-cCR_vs_non-cCR.csv`

The input table contains mirrored Scanpy-style rankings for `cCR` and `non-cCR`. The plotting script uses the `non-cCR` rows as the canonical direction, so positive log2 fold change means higher expression in pretreatment non-cCR malignant cells relative to pretreatment cCR malignant cells.

## Files

- `code/plot_pretreatment_ccr_vs_nonccr_volcano.R`: volcano plot and summary table script.
- `outputs/figure_6b_pretreatment_nonccr_vs_ccr_volcano.pdf`: candidate Figure 6B panel.
- `results/pretreatment_nonccr_vs_ccr_volcano_table.csv`: processed table used for plotting.
- `results/pretreatment_nonccr_vs_ccr_summary.csv`: DEG count summary.
- `results/pretreatment_nonccr_vs_ccr_candidate_genes.csv`: candidate ECM/remodeling genes found in this contrast.

## Current Candidate Gene Readout

Detected candidate genes in this DEG table: `MMP14`, `COL17A1`, `LAMB3`, and `LAMA3`.

Not detected in this exported table: `LAMB2` and `LAMC3`.
