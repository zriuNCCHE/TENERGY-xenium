"""Run pretreatment malignant-distance region annotation.

This server-side workflow annotates pretreatment cells by distance from
malignant regions and exports `malignant_distance_info.csv`, which is used by
Figure 6 spatial boundary / tumor margin analyses.
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

from xenium_malignant_distance import MalignantDistanceAnalyzer, MalignantDistanceConfig


adata = sc.read_h5ad("../../data/combined_adata/combined_final_all_cell_types.h5ad")
pre_treat_adata = adata[adata.obs["sample_timepoint"] == "pre",]

config = MalignantDistanceConfig(
    cell_type_col="final_cell_type",
    malignant_values=("malignant",),
    immune_stromal_values=("A2ML1+ epi", "fibroblast & endothelial cells", "immune cells"),
    spatial_key="spatial",
    sample_col="sampleID",
    neighbor_radius=50.0,
    close_distance=30.0,
    intermediate_distance=100.0,
    core_malignant_fraction=0.9,
    min_neighbors_for_region=5,
)

analyzer = MalignantDistanceAnalyzer(
    pre_treat_adata,
    config=config,
)

results = analyzer.run_pipeline(
    save=True,
    output_dir="results",
    output_prefix="my_malignant_distance",
)

pre_treat_adata.obs.to_csv("malignant_distance_info.csv")

analyzer.plot_spatial_annotations(
    color="malignant_region",
    spot_size=4,
    save_path="results/distance_region_spatial.pdf",
)
