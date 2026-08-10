import scanpy as sc

from xenium_metaprogram import MalignantMetaProgrammer, MetaProgramConfig


adata = sc.read_h5ad("../../data/combined_adata/combined_all_patients_epi_analyzed_final.h5ad")
adata = adata[adata.obs["sample_timepoint"] == "pre"].copy()
adata = adata[adata.obs["final_cell_type"] == "malignant"].copy()

config = MetaProgramConfig(
    sample_col="sampleID",
    layer="counts",
    k_range=(7, 15),
    top_n_genes=50,
    n_mps=15,
    min_patients=3,
    consensus_fraction=0.20,
)

mp_analyzer = MalignantMetaProgrammer(adata, config=config)
mp_analyzer.run_pipeline(
    plot_results=True,
    save=True,
    output_prefix="tenergy_meta_programs",
)
