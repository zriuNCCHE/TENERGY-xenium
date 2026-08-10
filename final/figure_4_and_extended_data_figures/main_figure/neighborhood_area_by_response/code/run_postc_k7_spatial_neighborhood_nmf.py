"""Run postC spatial-neighborhood NMF and export k=7 cell assignments.

This server-side workflow produces the postC k=7 neighborhood-assignment table
used by Figure 4 and its Extended Data Figure 4 spatial niche panels.
"""

import os
from pathlib import Path

from collections import Counter

import numpy as np
from scipy import sparse
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad
import squidpy as sq
import spatialdata as sd
import spatialdata_plot
from spatialdata_io import xenium
from spatialdata import bounding_box_query

from xenium_spatial_neighborhood import SpatialNeighborhoodAnalyzer, SpatialNeighborhoodConfig


combined_adata = sc.read_h5ad("../../data/combined_adata/combined_final_all_cell_types.h5ad")
combined_adata_postC = combined_adata[combined_adata.obs["sample_timepoint"] == "postC"].copy()

config = SpatialNeighborhoodConfig(
    cell_type_col="final_cell_type2",
    sample_col="sampleID",
    radius=100.0,
    n_components=6,
    transformation="sqrt",
    output_prefix="postC_treatment_neighborhoods",
)

sn = SpatialNeighborhoodAnalyzer(combined_adata_postC, config=config)

sn.compute_spatial_graph()
sn.build_proportion_matrix(transformation="sqrt")

summary, pair_table = sn.stability_scan_nmf(
    k_range=range(5, 13),
    n_repeats=30,
    sample_fraction=0.9,
)

sn.plot_neighbor_counts()
sn.plot_neighbor_distances()

sn.grid_search_nmf(range(3, 13))

sn.run_nmf(k=7)
sn.plot_sample_level_boxplots()
sn.plot_sample_level_boxplots(agg_mode="area")
sn.export_cell_assignments("postC_k=7_cell_niche_assignments.csv")
