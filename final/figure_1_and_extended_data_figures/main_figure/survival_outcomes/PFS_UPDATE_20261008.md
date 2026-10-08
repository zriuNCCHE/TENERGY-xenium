# PFS update: 2026-10-08

## Verified changes

Input: `data/_updated_survival/20261008TENERGY初発症例データ一覧.xlsx`.
Compared with the archived September 4 workbook, both PFS event columns changed
from 1 to 0 for eight patients. Both PFS duration columns are unchanged for all
40 patients. There are 30 PFS events: 9/16 cCR and 21/24 non-cCR.

Affected trial IDs: TENERGY-01-004, TENERGY-01-011, TENERGY-01-017,
TENERGY-01-018, TENERGY-04-006, TENERGY-04-011, TENERGY-06-008,
TENERGY-06-009. TENERGY-01-017 is non-cCR; the others are cCR.

Registration-origin KM estimates with log-log confidence limits:

| Group | Median PFS, months | 95% CI | PFS events |
|---|---:|---|---:|
| All patients | 3.2 | 2.4-8.9 | 30/40 |
| cCR | 36.0 | 8.9-not reached | 9/16 |
| non-cCR | 2.3 | 1.8-3.0 | 21/24 |

Two-sided log-rank P = 3.306111980667155e-07.
Fixed-group Cox HR (non-cCR versus cCR) = 8.2347; 95% CI 3.4307-19.7658
(use the unrounded validation CSV for exact machine-readable values).
The statistical department reports a DIFFERENT, time-dependent Cox model:
cCR versus non-cCR HR 0.39 (95% CI 0.15-1.05), P = 0.063.
Do not substitute this into a sentence describing a fixed-group Cox model.

## Manuscript passages requiring changes

Manuscript inspected: `_manuscript/_nature_cancer/Spatially resolved cellular dynamics during chemoradiotherapy and immunotherapy in esophageal squamous cell carcinoma_HB.docx`.
Paragraph numbers below count top-level Word body paragraphs; they are not
rendered Word line numbers. Find the quoted text to locate each passage.
No Word document or assembled Figure 1 TIFF was edited in this update.

### Results: survival paragraph (body paragraph 48)

Find: "3.2 months (95% CI, 2.5-16.6 months)" (the document uses an en dash).
Replace with: "3.2 months (95% CI, 2.4-8.9 months)".
The median is unchanged; the CI is updated using log-log confidence limits to
match the statistical-department report.

Replace the sentence beginning "The estimated median progression-free survival
was 34.7 months" with:

> The estimated median progression-free survival was 36.0 months (95% CI, 8.9 months to not reached) in the cCR group and 2.3 months (95% CI, 1.8-3.0 months) in the non-cCR group (Fig. 1c).

The old fixed-group HR of 4.73 (95% CI, 2.27-9.84) is obsolete. If keeping the
current descriptive fixed-group analysis, the replacement is:

> In the descriptive fixed-group Cox analysis, the hazard ratio for non-cCR relative to cCR was 8.23 (95% CI, 3.43-19.77).

For a time-dependent interpretation, use the department's separate analysis
instead, after confirming its full specification with the statisticians:

> In the statistical department's time-dependent Cox analysis, the hazard ratio for cCR relative to non-cCR was 0.39 (95% CI, 0.15-1.05).

That time-dependent result does not reach the conventional P < 0.05 threshold.
The significant log-rank result must not be described as the time-dependent
model's P value. The fixed-group result remains vulnerable to guarantee-time
bias, because cCR is observed after treatment.

### Methods: survival analysis (body paragraph 111)

If retaining the plotted fixed-group HR, clarify the existing Cox sentence:

> Univariable Cox proportional-hazards models treated eventual cCR status as a fixed group to estimate descriptive hazard ratios and 95% confidence intervals for non-cCR relative to cCR. Median survival confidence intervals were derived from log-log transformed Kaplan-Meier confidence limits. These post-treatment group comparisons are susceptible to guarantee-time bias.

Add the source date (October 8, 2026). If reporting the department's
time-dependent Cox result, separately describe its time-varying response
indicator and response-onset date definition after confirmation; the current
workbook does not contain the onset dates needed to reproduce that analysis.

### Figure 1 legend (body paragraph 199)

Group sizes remain 16 and 24. Keep the trial-registration time origin.
If retaining the fixed-group HR, clarify that "univariable Cox
proportional-hazards models" use eventual cCR as a fixed group, and that the
log-rank P value and HR are descriptive comparisons. Do not call this a
time-dependent HR. Replace the two PFS panels in the assembled Figure 1 artwork.

### Claims not to carry forward

Any statement that all 16 cCR patients experienced a PFS event is incorrect:
the current count is 9 events and 7 censored observations.
No such statement was found in the inspected manuscript's current body text.

## Checks performed

All 40 records retained; 16 cCR and 24 non-cCR; PFS events total 30.
Curves start at survival 1 and time 0, are non-increasing, and stop at each
group's actual last observation. Censor-only times are retained.
Time-zero risk counts are 40 pooled, 16 cCR, 24 non-cCR, without clipping.
36- and 48-month group estimates agree with the October 8 statistical report.
Both PFS panels were visually inspected. Arial PDF rendering is warning-free.
OS artwork was not regenerated. Generated PDFs/CSVs remain local under the
existing GitHub policy excluding figure outputs; source files and code are
versioned.
