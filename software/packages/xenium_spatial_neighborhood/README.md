# Xenium Spatial Neighborhood Analysis

Reusable Python tools for identifying cellular neighborhoods from Xenium/AnnData
spatial data using local cell-type composition and NMF.

The package follows the same compact analyzer pattern as the Xenium malignant
cell distance package: configure a dataclass, instantiate an analyzer with an
AnnData object, run the pipeline, and receive new per-cell columns in
`adata.obs`.

## What This Package Does

`xenium-spatial-neighborhood-analysis` takes a mixed-cell `AnnData` object and:

1. builds a Squidpy spatial neighbor graph,
2. counts each cell's neighboring cell types,
3. converts those counts into a cell-by-cell-type neighborhood matrix,
4. optionally applies `sqrt` or `log1p` transformation,
5. runs NMF to identify recurrent cellular-neighborhood signatures,
6. writes dominant neighborhood labels and continuous NMF weights to `adata.obs`,
7. optionally summarizes neighborhood abundance by sample and clinical group.
8. scans NMF ranks using repeated cell-type co-occurrence stability.

The main entry point is:

```python
from xenium_spatial_neighborhood import SpatialNeighborhoodAnalyzer, SpatialNeighborhoodConfig
```

## Quick Start

```python
import scanpy as sc
from xenium_spatial_neighborhood import SpatialNeighborhoodAnalyzer, SpatialNeighborhoodConfig

adata = sc.read_h5ad("/path/to/xenium_cells.h5ad")

config = SpatialNeighborhoodConfig(
    cell_type_col="cell_type_lvl1",
    sample_col="sampleID",
    radius=50.0,
    n_components=6,
    transformation="sqrt",
    figure_dir="figures",
)

analyzer = SpatialNeighborhoodAnalyzer(adata, config=config)
weights, signatures = analyzer.run_pipeline(save=True, output_dir="results")

analyzer.plot_signatures()
analyzer.plot_neighborhood_pie_charts()
analyzer.plot_sample_level_boxplots(clinical_key="cCR", agg_mode="area")

pvalues = analyzer.compare_neighborhood_proportions(
    clinical_key="cCR",
    save=True,
    output_dir="results",
)

analyzer.export_cell_assignments("results/cell_niche_assignments.csv")
```

The exported assignment table is intended to be self-contained for downstream
sample/group boxplots. By default it includes `cell_id`, all `adata.obs`
metadata columns, the dominant niche label, and all continuous neighborhood
weight columns. If `adata.obs["cell_id"]` exists, that exact column is used as
the exported `cell_id`; otherwise the AnnData index is used.

## Main Outputs

The analyzer adds these columns to `adata.obs`:

- `nmf_dominant`: dominant NMF neighborhood for each cell.
- `neighborhood_1_k6`, `neighborhood_2_k6`, ...: continuous per-cell
  neighborhood weights.

The analyzer also stores parameter provenance in:

```python
adata.uns["spatial_neighborhood_analysis"]
```

If `save=True`, the pipeline writes:

- `*_proportion_matrix.csv`: local cell-type neighborhood matrix.
- `*_nmf_weights.csv`: per-cell NMF neighborhood weights.
- `*_nmf_signatures.csv`: NMF component-by-cell-type signatures.
- `*_cell_niche_assignments.csv`: per-cell dominant niche and NMF weights.
- `*_neighborhood_proportion_pvalues.csv`: optional two-group p-value table
  from sample-level neighborhood proportions.

Plotting methods save PDFs to `config.figure_dir` by default, using
`config.output_prefix` in the filename. For example:

```python
analyzer.plot_neighbor_counts()
analyzer.plot_neighbor_distances()
```

These two QC violin plots are sorted from lowest to highest median value by
default. Use `sort=False` if you want the original cell-type order.

writes:

```text
figures/spatial_neighborhood_neighbor_counts.pdf
figures/spatial_neighborhood_neighbor_distances.pdf
```

To keep every plot in one analysis folder, set it once:

```python
config = SpatialNeighborhoodConfig(
    output_prefix="my_analysis",
    figure_dir="results/figures",
)
```

## Interactive Rank Selection

```python
analyzer.compute_spatial_graph()
analyzer.build_proportion_matrix(transformation="sqrt")
errors = analyzer.find_k_elbow(range(2, 12))
signature_grid = analyzer.grid_search_nmf(range(3, 9))
weights, signatures = analyzer.run_nmf(k=6)
```

`grid_search_nmf()` now draws vertical warm heatmaps by default, with cell types
as rows and NMF components as columns. Cell types are sorted into
component-specific blocks so the heatmaps are easier to read than the raw input
cell-type order.

Heatmap-style plots use the same Warm amber palette by default, including
grid-search signatures, K-stability co-occurrence heatmaps, final NMF
signatures, and patient heatmaps.

For a stability-based K scan inspired by repeated NMF co-occurrence analysis:

```python
summary, pair_table = analyzer.stability_scan_nmf(
    k_range=range(2, 11),
    n_repeats=50,
    sample_fraction=0.8,
)

analyzer.plot_k_cooccurrence_heatmap(k=6)
```

The scan automatically plots and saves:

- `figures/spatial_neighborhood_combined_cooccurrence_heatmap.pdf`
- `figures/spatial_neighborhood_k_scale_cooccurrence_heatmap.pdf`

The combined heatmap is a square cell-type-by-cell-type matrix using relative
co-occurrence frequency from 0 to 1, where 1 means 100% of scanned K/repeat runs
placed the pair in the same NMF component. The K-scale heatmap shows each
cell-type pair across candidate ranks.

## Installation

On an HPC or server, copy this source folder and install it directly:

```bash
cd /path/to/xenium-spatialneighborhood-analysis
python3 -m pip install -e . --no-build-isolation
```

If dependencies are provided by a local pip mirror/image, install them from that
mirror as usual. No wheel build is required for using this package.

## Project Layout

```text
src/xenium_spatial_neighborhood/
  spatial_neighborhood_analyzer.py
examples/
  basic_usage.py
docs/
  user_guide.md
  figures/spatial_neighborhood_workflow_schematic.pdf
```
