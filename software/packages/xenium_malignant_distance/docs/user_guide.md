# Xenium Malignant Cell Distance Analysis User Guide

This guide describes the first-pass malignant distance analyzer. The goal is to
produce transparent, editable annotations rather than a hidden black-box tumor
boundary model.

## Required AnnData Fields

The analyzer needs:

- `adata.obs[cell_type_col]`: cell-type labels containing malignant and
  immune/stromal labels.
- `adata.obsm[spatial_key]`: spatial coordinates, usually `adata.obsm["spatial"]`.

If coordinates are stored in observation columns instead, set
`coordinate_cols=("x", "y")` in `MalignantDistanceConfig`.

## Configuration

```python
from xenium_malignant_distance import MalignantDistanceConfig

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
    min_neighbors_for_region=5,
)
```

Use `sample_col` when the AnnData contains multiple tissues, samples, patients,
or fields of view that should not share spatial neighborhoods.

## Classification Logic

For every cell, the analyzer computes distance to the nearest malignant cell
within the same sample group.

For malignant cells, the analyzer counts nearby cells within `neighbor_radius`
and computes:

```text
malignant_neighbor_fraction = malignant_neighbor_count / neighbor_count
```

Malignant cells are assigned:

- `isolated_malignant` when there are fewer than `min_neighbors_for_region`
  local neighbors.
- `tumor_inner` when the local malignant fraction is at least
  `core_malignant_fraction`.
- `tumor_margin` otherwise.

For immune/stromal cells, the analyzer uses nearest malignant-cell distance:

- `close_to_malignant` when distance is at most `close_distance`.
- `intermediate_to_malignant` when distance is above `close_distance` and at
  most `intermediate_distance`.
- `distant_to_malignant` when distance is above `intermediate_distance`.
- `no_malignant_reference` when the sample group has no malignant cells.

Other non-malignant cells are labeled `other_non_malignant` by default.

## Plotting Spatial Annotations

After running the pipeline, use:

```python
analyzer.plot_spatial_annotations(color="distance_region")
```

This plots all samples by default. If `sample_col` is set, each sample is shown
in its own panel. The `color` argument can be any column in `adata.obs`,
including:

```python
"distance_region"
"malignant_region"
"malignant_proximity"
"distance_to_malignant"
"malignant_neighbor_fraction"
```

For continuous columns:

```python
analyzer.plot_spatial_annotations(
    color="distance_to_malignant",
    cmap="magma",
    spot_size=4,
)
```

For selected samples:

```python
analyzer.plot_spatial_annotations(
    color="distance_region",
    samples=["sample_1", "sample_2"],
)
```

To save the plot:

```python
analyzer.plot_spatial_annotations(
    color="distance_region",
    save_path="results/distance_region_spatial.png",
)
```
