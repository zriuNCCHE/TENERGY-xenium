# non-cCR postC versus pre DEG

Candidate Figure 6C.

This part uses the saved DEG table:

`data/raw/figure_6_deg/non-cCR_postC_vs_pre.csv`

The input table contains mirrored Scanpy-style rankings for `postC` and `pre`. The plotting script uses the `postC` rows as the canonical direction, so positive log2 fold change means higher expression in post-chemo-radiotherapy non-cCR malignant cells relative to pretreatment non-cCR malignant cells.

## Files

- `code/plot_nonccr_postc_vs_pre_volcano.R`: volcano plot and summary table script.
- `outputs/figure_6c_nonccr_postc_vs_pre_volcano.pdf`: candidate Figure 6C panel.
- `results/nonccr_postc_vs_pre_volcano_table.csv`: processed table used for plotting.
- `results/nonccr_postc_vs_pre_summary.csv`: DEG count summary.
- `results/nonccr_postc_vs_pre_candidate_genes.csv`: candidate ECM/remodeling genes found in this contrast.

## Current Candidate Gene Readout

Detected candidate genes in this DEG table: `LAMB3`, `COL17A1`, `LAMA3`, and `MMP14`.

Not detected in this exported table: `LAMB2` and `LAMC3`.
