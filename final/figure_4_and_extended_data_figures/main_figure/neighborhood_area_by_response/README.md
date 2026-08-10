# Figure 4 Spatial Neighborhood Area by Response

This folder contains the main Figure 4 spatial-neighborhood area/score comparison between cCR and non-cCR samples.

## Upstream Python Workflow

`code/run_postc_k7_spatial_neighborhood_nmf.py` is the server-side script used to generate the postC k=7 spatial-neighborhood assignments from the combined all-cell AnnData object.

The script:

- loads `combined_final_all_cell_types.h5ad`;
- filters to `sample_timepoint == "postC"`;
- builds a spatial neighbor graph using `final_cell_type2`;
- scans NMF stability over k values;
- runs the selected k=7 NMF model;
- exports `postC_k=7_cell_niche_assignments.csv`.

This output table is treated as source data for the Figure 4 spatial niche analysis and the matched Extended Data Figure 4 niche-composition panels.

## Main R Plotting Workflow

`code/plot_neighborhood_area_by_response.R` uses the exported k=7 neighborhood-assignment table to generate sample-level neighborhood area/score boxplots by clinical response.

The preferred statistical logic is sample-level area-weighted neighborhood score:

`sum(cell neighborhood weight * cell area) / sum(cell area)`

Each plotted point is one patient/sample, not one cell.
