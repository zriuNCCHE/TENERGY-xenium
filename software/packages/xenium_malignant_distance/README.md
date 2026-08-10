# Xenium Malignant Cell Distance Analysis

Reusable Python tools for annotating tumor margin/core structure and
immune-stromal proximity to malignant cells from Xenium/AnnData spatial data.

The package follows the same compact analyzer pattern as the Xenium
Metaprogram analyzer: configure a dataclass, instantiate an analyzer with an
AnnData object, run the pipeline, and receive new per-cell columns in
`adata.obs`.

## What This Package Does

`xenium-malignant-cell-distance-analysis` takes a mixed-cell `AnnData` object and:

1. identifies malignant cells from a configured `adata.obs` cell-type column,
2. computes each cell's distance to the nearest malignant cell,
3. computes local malignant-neighbor fractions within a spatial radius,
4. classifies malignant cells as `tumor_inner`, `tumor_margin`, or
   `isolated_malignant`,
5. classifies immune/stromal cells as `close_to_malignant`,
   `intermediate_to_malignant`, `distant_to_malignant`, or
   `no_malignant_reference`,
6. writes annotations back to `adata.obs`,
7. optionally exports per-cell annotations and summary counts.

The main entry point is:

```python
from xenium_malignant_distance import MalignantDistanceAnalyzer, MalignantDistanceConfig
```

## Quick Start

```python
import scanpy as sc
from xenium_malignant_distance import MalignantDistanceAnalyzer, MalignantDistanceConfig

adata = sc.read_h5ad("/path/to/xenium_cells.h5ad")

config = MalignantDistanceConfig(
    cell_type_col="final_cell_type",
    malignant_values=("malignant",),
    immune_stromal_values=("T cell", "B cell", "myeloid", "fibroblast", "stromal"),
    spatial_key="spatial",
    sample_col="sampleID",
    neighbor_radius=50.0,
    close_distance=30.0,
    intermediate_distance=100.0,
    core_malignant_fraction=0.75,
)

analyzer = MalignantDistanceAnalyzer(adata, config=config)
results = analyzer.run_pipeline(save=True, output_dir="results")

analyzer.plot_spatial_annotations(color="distance_region")
```

## Main Outputs

The analyzer adds these columns to `adata.obs`:

- `distance_to_malignant`: Euclidean distance to the nearest malignant cell.
- `nearest_malignant_cell`: observation name of the nearest malignant cell.
- `neighbor_count`: number of cells within `neighbor_radius`.
- `malignant_neighbor_count`: number of local neighbors labeled malignant.
- `malignant_neighbor_fraction`: malignant fraction among local neighbors.
- `malignant_region`: malignant-cell label: `tumor_inner`, `tumor_margin`, or
  `isolated_malignant`.
- `malignant_proximity`: non-malignant proximity label. Samples with no
  malignant reference cells are labeled `no_malignant_reference`.
- `distance_region`: unified region label across all cells.

The analyzer also stores parameter provenance in:

```python
adata.uns["malignant_distance_analysis"]
```

## Plotting

After `run_pipeline()`, plot the spatial annotations across all samples:

```python
analyzer.plot_spatial_annotations(color="distance_region")
```

Useful color columns include:

```python
analyzer.plot_spatial_annotations(color="malignant_region")
analyzer.plot_spatial_annotations(color="malignant_proximity")
analyzer.plot_spatial_annotations(color="distance_to_malignant", cmap="magma")
analyzer.plot_spatial_annotations(color="malignant_neighbor_fraction", cmap="viridis")
```

If `sample_col` is configured, the plot uses one panel per sample. To plot only
selected samples:

```python
analyzer.plot_spatial_annotations(
    color="distance_region",
    samples=["sample_1", "sample_2"],
)
```

## Algorithm Notes

Nearest malignant-cell distance is most useful for non-malignant cells. For
malignant cells, the analyzer uses local composition: malignant cells surrounded
mostly by malignant neighbors are called tumor-inner, while malignant cells with
more mixed local neighborhoods are called tumor-margin. This makes margin/core
classification depend on tumor-nest context rather than the trivial
self-distance of malignant cells.

## Installation

```bash
cd /path/to/xenium-malignant-cell-distance-analysis
python3 -m pip install -e . --no-build-isolation
```

## Project Layout

```text
src/xenium_malignant_distance/
  malignant_distance_analyzer.py
examples/
  run_example_distance.py
docs/
```
