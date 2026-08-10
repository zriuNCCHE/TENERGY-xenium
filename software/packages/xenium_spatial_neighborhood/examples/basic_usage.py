"""Minimal example for Xenium spatial neighborhood analysis."""

import anndata as ad

from xenium_spatial_neighborhood import SpatialNeighborhoodAnalyzer, SpatialNeighborhoodConfig


adata = ad.read_h5ad("xenium_cells.h5ad")

config = SpatialNeighborhoodConfig(
    cell_type_col="cell_type_lvl1",
    sample_col="sampleID",
    radius=50,
    n_components=6,
    transformation="sqrt",
)

analyzer = SpatialNeighborhoodAnalyzer(adata, config=config)
analyzer.compute_spatial_graph()
analyzer.build_proportion_matrix()
analyzer.find_k_elbow(range(2, 12))
weights, signatures = analyzer.run_nmf()
analyzer.plot_signatures()
analyzer.plot_sample_level_boxplots(clinical_key="cCR", agg_mode="area")
