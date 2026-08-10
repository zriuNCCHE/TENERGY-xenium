"""Run pretreatment malignant-cell metaprogram discovery.

This server-side workflow identifies malignant-cell metaprograms from
pretreatment malignant cells. The resulting MP definitions are used by Figure 6
metaprogram analyses and downstream Figure 7 spatial proximity analyses.
"""

import scanpy as sc

from xenium_metaprogram import MalignantMetaProgrammer, MetaProgramConfig


adata = sc.read_h5ad("../../data/combined_adata/combined_all_patients_epi_analyzed_final.h5ad")
adata = adata[adata.obs["sample_timepoint"] == "pre"].copy()
adata = adata[adata.obs["final_cell_type"] == "malignant"].copy()

config = MetaProgramConfig(
    sample_col="sampleID",
    layer="counts",
    k_range=(7, 10),
    top_n_genes=20,
    n_mps=10,
    min_patients=5,
    consensus_fraction=0.20,
)

mp_analyzer = MalignantMetaProgrammer(adata, config=config)
mp_analyzer.run_pipeline(plot_results=True)
