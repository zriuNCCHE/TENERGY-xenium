# PostC T-cell abundance across highest-score spatial neighborhoods

Exploratory analysis based on the simple highest-score spatial neighborhood assignment.
Each cell was assigned to one unique neighborhood by selecting the largest value among neighborhood_1_k7 through neighborhood_7_k7.
For each patient and neighborhood, T-cell subset abundance was calculated as the number of cells in that subset divided by the total number of cells assigned to that neighborhood.

Individual T subsets: T/NK, CD8_Teff, CD8_Tex_PDCD1, CD8_prolif, CD4_Treg_CCR8, CD4_Treg_FOXP3, CD4_CXCL13 and CD4_prolif.
Combined T groups: Treg-like = CD4_Treg_CCR8 + CD4_Treg_FOXP3; non-Treg T/NK = all other listed T/NK, CD8 and CD4 T subsets.

Global P values use a paired Friedman test across the seven neighborhoods.
Pairwise P values use paired Wilcoxon signed-rank tests between neighborhood pairs, with BH adjustment within each feature.
This is intentionally circular/exploratory because neighborhoods were derived from cell-type neighborhood information.
