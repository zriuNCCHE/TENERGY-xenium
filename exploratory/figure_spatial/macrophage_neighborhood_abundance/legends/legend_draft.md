# PostC macrophage abundance across highest-score spatial neighborhoods

Exploratory analysis based on the simple highest-score spatial neighborhood assignment.
Each cell was assigned to one unique neighborhood by selecting the largest value among neighborhood_1_k7 through neighborhood_7_k7.
For each patient and neighborhood, macrophage abundance was calculated as the number of Macro_CXCL9 or Macro_CXCL5 cells divided by the total number of cells assigned to that neighborhood.
Single-feature plots order niches from least to most abundant by median cell fraction.
Global P values use a paired Friedman test across the seven neighborhoods. Pairwise P values use paired Wilcoxon signed-rank tests between neighborhood pairs with BH adjustment.
