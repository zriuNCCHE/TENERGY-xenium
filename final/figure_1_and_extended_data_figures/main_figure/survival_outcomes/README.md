# Figure 1b-c: survival outcomes

## Current source

The active survival release is 2026-10-08. The patient-level workbook is
`data/_updated_survival/20261008TENERGY初発症例データ一覧.xlsx`, sheet
`TENERGY長期追跡`. The accompanying presentation and figure workbook are
statistical-department reference outputs, not the input for fitting these curves.
All 40 initial-cohort patients are included (16 cCR and 24 non-cCR).

## Reproduce

Requires R, readxl, survival, ggplot2 and gridExtra; the PDF device uses macOS
Quartz and installed Arial. From the repository root:

```sh
Rscript final/figure_1_and_extended_data_figures/main_figure/survival_outcomes/code/plot_figure1_survival_outcomes.R --endpoint=pfs
```

Use `--endpoint=os` for OS only, or omit the flag to regenerate both endpoints.
The October 8 update regenerated PFS only; existing OS artwork was untouched.
Outputs include the combined Figure 1c PDF, separate pooled and response-group
PFS PDFs, and CSVs of curves, risk counts, KM summaries and exact statistics.
Generated outputs are excluded from GitHub under the repository policy.

## Statistical definitions

PFS uses `leng_pfs_obsm` and `event_pfs_obs`, with trial registration as time zero
and 1 indicating an event, 0 censoring. Censor-only times are retained and curves
start at time zero with survival 1. Curves stop at observed follow-up, without
artificial extension. Median-survival confidence intervals use log-log transformed
KM confidence limits, matching the statistical-department report's rounding.

The response-group log-rank test and Cox HR displayed here treat eventual response
as a fixed group. They are descriptive, potentially affected by guarantee-time
bias, and must not be represented as a baseline prognostic or causal effect.
The HR direction in these plots is non-cCR versus cCR.

The statistical department instead reports a time-dependent Cox model: PFS HR
0.39 (95% CI 0.15-1.05), P = 0.063, for cCR versus non-cCR. This is not the
plotted fixed-group estimate. Response-onset dates are not available in the
patient workbook, so the time-dependent model is not reconstructed here.

## October 8 changes

Eight PFS event flags changed from 1 to 0 (seven cCR, one non-cCR), under both
registration-origin and CRT-origin definitions. No PFS durations changed.
There are now 30 PFS events (9 cCR, 21 non-cCR), rather than 38.
Earlier source files were moved outside this Git repository to
`../_archive/survival_before_20261008/`; they were not destroyed.
