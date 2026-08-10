Extended Figure 6 all malignant metaprogram score legend draft

Cell-level relative malignant metaprogram scores by clinical outcome. Violin plots show all pretreatment malignant cells; embedded boxes show median and interquartile range. Scores below zero are expected for relative module scores and indicate expression below the matched background/reference level, not negative expression.

P values were calculated by two-sided Wilcoxon rank-sum tests comparing cCR and non-cCR cells for each MP, with log-scale p-value storage to avoid numerical underflow. The figure displays threshold-style significance labels only; exact medians, mean scores, log10 p values and BH-adjusted p values across the seven MPs are saved in `results/all_mp_score_cell_level_response_stats.csv`.

Discordance outputs assign each cell to its highest-scoring MP and summarize the top-score margin, MP7 dominance, MP7 rank and pairwise Spearman correlations among MP scores. These are saved in `results/all_mp_score_dominance_discordance_cell_level.csv`, `results/all_mp_score_dominance_discordance_summary.csv` and `results/all_mp_score_spearman_correlation.csv`.
