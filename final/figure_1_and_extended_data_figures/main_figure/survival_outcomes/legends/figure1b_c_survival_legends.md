# Figure 1b-c survival outcomes

**Figure 1b | Overall survival in the updated TENERGY initial cohort.** Kaplan-Meier estimates of overall survival from clinical-trial registration for all patients (left) and by clinical response (right; cCR, n = 16; non-cCR, n = 24).

**Figure 1c | Progression-free survival in the updated TENERGY initial cohort.** Kaplan-Meier estimates of progression-free survival from clinical-trial registration for all patients (left) and by clinical response (right; cCR, n = 16; non-cCR, n = 24).

Crosses denote censored observations. Numbers below each curve are patients at risk. Between-response comparisons use two-sided log-rank tests; hazard ratios (HRs; non-cCR versus cCR) and 95% confidence intervals are from univariable Cox proportional-hazards models with eventual response treated as a fixed group. These are descriptive post-treatment group comparisons, not time-dependent Cox estimates, and are susceptible to guarantee-time bias. OS uses `leng_os_obsm` and `event_os_obs`; PFS uses `leng_pfs_obsm` and `event_pfs_obs`. Median survival confidence intervals use the log-log transformation.

PFS source updated to `20261008TENERGY初発症例データ一覧.xlsx`. OS artwork was not regenerated in the October 8 PFS-only run. The statistical department's time-dependent Cox estimate for PFS is a different analysis: cCR versus non-cCR HR 0.39 (95% CI 0.15-1.05), P = 0.063, as reported in the October 8 presentation; this script does not reconstruct that model without response-onset dates.
