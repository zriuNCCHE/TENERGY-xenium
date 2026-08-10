import scanpy as sc

from xenium_malignant_distance import MalignantDistanceAnalyzer, MalignantDistanceConfig


adata = sc.read_h5ad("/path/to/xenium_cells.h5ad")

config = MalignantDistanceConfig(
    cell_type_col="final_cell_type",
    malignant_values=("malignant",),
    immune_stromal_values=(
        "T cell",
        "B cell",
        "myeloid",
        "fibroblast",
        "endothelial",
        "stromal",
        "immune",
    ),
    spatial_key="spatial",
    sample_col="sampleID",
    neighbor_radius=50.0,
    close_distance=30.0,
    intermediate_distance=100.0,
    core_malignant_fraction=0.75,
)

distance_analyzer = MalignantDistanceAnalyzer(adata, config=config)
distance_analyzer.run_pipeline(
    save=True,
    output_dir="results",
    output_prefix="example_malignant_distance",
)
