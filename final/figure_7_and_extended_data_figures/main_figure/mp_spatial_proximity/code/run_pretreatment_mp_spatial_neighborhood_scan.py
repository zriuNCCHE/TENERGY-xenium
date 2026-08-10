"""Run pretreatment MP-aware spatial-neighborhood scan.

This server-side workflow uses `final_cell_type2_malignant_MP`, where malignant
cells are represented by metaprogram labels and non-malignant cells retain their
cell-type labels. It was used as an upstream exploratory/selection workflow for
Figure 7 spatial proximity analyses.
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
combined_adata_pre = combined_adata[combined_adata.obs["sample_timepoint"] == "pre"].copy()

config = SpatialNeighborhoodConfig(
    cell_type_col="final_cell_type2_malignant_MP",
    sample_col="sampleID",
    radius=20.0,
    n_components=6,
    transformation="sqrt",
    output_prefix="pre_malignant_added_20um_treatment_neighborhoods",
)

sn = SpatialNeighborhoodAnalyzer(combined_adata_pre, config=config)

sn.compute_spatial_graph()
sn.build_proportion_matrix(transformation="sqrt")

summary, pair_table = sn.stability_scan_nmf(
    k_range=range(5, 13),
    n_repeats=30,
    sample_fraction=0.9,
)

sn.grid_search_nmf(range(3, 13))
