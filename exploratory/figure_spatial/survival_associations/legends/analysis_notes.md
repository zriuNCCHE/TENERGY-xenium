# Exploratory postC survival association analysis

This exploratory analysis tests whether postC spatial neighborhood scores and T/fibroblast/endothelial cell-type proportions are associated with progression-related outcomes.

Survival metadata found in the old OBS/spatial metadata: time_to_PD_month and time_to_CR_month. No explicit censoring/follow-up time column was found.

Spatial neighborhood features were calculated per patient as area-weighted continuous NMF scores: sum(neighborhood weight x cell area) divided by total annotated cell area.

Cell-type features were calculated as percent of all non-epithelial, non-malignant cells at postC.

Primary exploratory outputs:
- PD-event association: Wilcoxon rank-sum test comparing patients with versus without recorded PD.
- Time-to-PD association: Spearman correlation among patients with recorded PD only.
- Proxy PD-free curves: median split by feature value; patients without recorded PD were censored at the maximum observed PD time. This is only a rough visual screen because true censoring times are not available.
- Cox PH continuous models: univariable Cox proportional hazards models using each postC feature scaled to mean 0 and SD 1; hazard ratios therefore correspond to a 1 SD increase in the feature.
- Cox PH median-split models: univariable Cox proportional hazards models comparing high versus low feature groups split at the median; this is the same grouping used for the proxy PD-free curves.
- Proportional hazards diagnostics: Schoenfeld residual global p values were calculated for each univariable Cox model when possible.

These results should not be treated as final survival statistics until an explicit event indicator and censor/follow-up time are provided.
