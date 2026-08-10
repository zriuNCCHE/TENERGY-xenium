# User Guide

## Purpose

`xenium-spatial-neighborhood-analysis` identifies recurrent cellular
neighborhoods from Xenium-style spatial transcriptomics data. Instead of using
gene-expression programs, it uses the local cell-type composition around every
cell. NMF decomposes these local compositions into interpretable neighborhood
signatures.

The package mirrors `xenium-malignant-cell-distance-analysis`: one config
dataclass, one analyzer class, notebook-friendly methods, and an optional
`run_pipeline()` helper.

## Input Requirements

The input should be an `AnnData` object with:

- cell metadata in `adata.obs`,
- cell-type labels in a configured column such as `cell_type_lvl1`,
- spatial coordinates compatible with `squidpy.gr.spatial_neighbors`,
- sample IDs when using multi-sample graphs or patient-level summaries.

Typical columns:

```text
adata.obs["cell_type_lvl1"]
adata.obs["sampleID"]
adata.obs["cCR"]
adata.obs["cell_area"]
```

## Standard Workflow

A schematic overview of the analysis design is available at:

```text
docs/figures/spatial_neighborhood_workflow_schematic.pdf
```

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
    output_prefix="pre_treatment_neighborhoods",
    figure_dir="results/figures",
)

sn = SpatialNeighborhoodAnalyzer(adata, config=config)
weights, signatures = sn.run_pipeline(save=True, output_dir="results")
```

This writes:

```text
results/pre_treatment_neighborhoods_proportion_matrix.csv
results/pre_treatment_neighborhoods_nmf_weights.csv
results/pre_treatment_neighborhoods_nmf_signatures.csv
```

## Interactive Rank Selection

For exploratory analysis, choose the NMF rank manually:

```python
sn.compute_spatial_graph()
sn.build_proportion_matrix(transformation="sqrt")
errors = sn.find_k_elbow(range(2, 12))
signature_grid = sn.grid_search_nmf(range(3, 9))
weights, signatures = sn.run_nmf(k=6)
```

Use the elbow plot together with component heatmaps. Prefer a rank that gives
interpretable cell-type signatures without splitting one obvious biological
pattern into many small components.

`grid_search_nmf()` draws vertical warm heatmaps by default, with cell types as
rows and NMF components as columns. For each K, every cell type is
assigned to the component where it has the highest loading, then cell types are
ordered by component and loading strength. This makes the heatmap easier to
interpret than the original arbitrary cell-type order.

## Stability-Based K Scan

Elbow plots are useful but often insufficient for NMF. A more biological
diagnostic is to ask whether the same cell types repeatedly group together in
NMF signatures across candidate K values. This follows the same logic as the
immune-cell co-occurrence analysis used in the anti-PD-1 NSCLC atlas: repeated
NMF runs are used to identify stable co-occurrence modules rather than relying
on reconstruction error alone.

In this package, each NMF run produces an `H` matrix with neighborhoods as rows
and neighboring cell types as columns. For each run, every cell type is assigned
to the NMF component where it has the strongest loading. The stability scan then
counts how often pairs of cell types are assigned to the same component.

```python
sn.compute_spatial_graph()
sn.build_proportion_matrix(transformation="sqrt")

summary, pair_table = sn.stability_scan_nmf(
    k_range=range(2, 11),
    n_repeats=50,
    sample_fraction=0.8,
)
```

`sample_fraction` defaults to `0.8`, so each repeat uses a random 80% subset of
cells. The scan also uses random NMF initialization by default. Together, these
settings test whether the same cell-type modules are recovered despite modest
data perturbation and different random starts.

Because this scan can take a while, the function prints progress for each K and
every 10 repeats by default. To change the reporting frequency:

```python
summary, pair_table = sn.stability_scan_nmf(log_every=5)
```

The scan automatically plots and saves two figures under `config.figure_dir`:

```text
results/figures/pre_treatment_neighborhoods_combined_cooccurrence_heatmap.pdf
results/figures/pre_treatment_neighborhoods_k_scale_cooccurrence_heatmap.pdf
```

The combined co-occurrence heatmap is the square, diagonal-block style plot. It
combines all scanned K values and repeats into one cell-type-by-cell-type
matrix, with the colorbar showing relative co-occurrence frequency from 0 to 1.
A value of 1 means that the pair was assigned to the same NMF component in 100%
of scanned K/repeat runs.

To plot it again:

```python
sn.plot_combined_cooccurrence_heatmap()
```

To use raw absolute counts instead:

```python
sn.plot_combined_cooccurrence_heatmap(value="count")
```

The K-scale co-occurrence heatmap has cell-type pairs as rows and K values as
columns. To plot it again or adjust filtering:

```python
sn.plot_k_stability_heatmap()
```

Values near 1 mean that a pair repeatedly appears in the same NMF signature.
Stable vertical patterns across adjacent K values indicate robust modules;
sudden splitting or merging suggests over- or under-classification.

To inspect the full cell-type co-occurrence matrix at one K:

```python
sn.plot_k_cooccurrence_heatmap(k=6)
```

The returned `summary` table contains reconstruction error and average
co-occurrence statistics by K. The returned `pair_table` contains one row per
cell-type pair per K, so you can filter or export it:

```python
pair_table.sort_values("cooccurrence_frequency", ascending=False).head()
```

If you want a different output folder just for this scan:

```python
summary, pair_table = sn.stability_scan_nmf(
    figure_dir="results/figures",
    output_prefix="my_neighborhood_scan",
)
```

After choosing K from the stability heatmap and biology, run the final model:

```python
weights, signatures = sn.run_nmf(k=6)
sn.plot_signatures()
```

## Graph Construction

Use either a fixed radius or a fixed number of nearest neighbors:

```python
SpatialNeighborhoodConfig(radius=50.0)
SpatialNeighborhoodConfig(radius=None, n_neighbors=20)
```

For multi-sample AnnData objects, set `sample_col`. It is passed to Squidpy as
`library_key` by default so neighborhoods are not connected across samples.

## Main Outputs

- `analyzer.prop_df`: per-cell cell-type neighborhood percentages.
- `weights`: per-cell NMF weights, one column per neighborhood.
- `signatures`: NMF component-by-cell-type signature matrix.
- `adata.obs["nmf_dominant"]`: dominant neighborhood assignment per cell.
- `adata.obs["neighborhood_*"]`: continuous NMF weights per cell.
- `adata.uns["spatial_neighborhood_analysis"]`: parameter provenance.

## Function Reference

All plotting functions save a PDF by default under `figures/` using
`config.output_prefix`. The recommended way to keep all plots in one place is
to set `figure_dir` once:

```python
config = SpatialNeighborhoodConfig(
    output_prefix="pre_treatment_neighborhoods",
    figure_dir="results/figures",
)
sn = SpatialNeighborhoodAnalyzer(adata, config=config)
```

Then all plot calls save into that same folder, including neighbor counts,
neighbor distances, K scans, signatures, neighborhood pie charts, and
sample-level boxplots.

To disable saving for any plot:

```python
sn.plot_signatures(save=False)
```

The package uses a consistent visual style by default. Heatmaps use the Warm
amber palette, where low values are pale and high values are darker. Violin and
box plots use matching muted colors.

### Setup

```python
config = SpatialNeighborhoodConfig(
    cell_type_col="cell_type_lvl1",
    cell_id_col="cell_id",
    sample_col="sampleID",
    radius=50.0,
    n_components=6,
)
sn = SpatialNeighborhoodAnalyzer(adata, config=config)
```

`SpatialNeighborhoodConfig` stores reusable parameters: cell-type column, cell
ID column, sample column, graph radius or nearest-neighbor count, NMF rank,
transformation, random seed, output prefix, and plotting/export defaults.

### Graph And Matrix

`compute_spatial_graph()`

Builds a Squidpy spatial graph and writes it to:

```text
adata.obsp["spatial_connectivities"]
adata.obsp["spatial_distances"]
```

Use this when the AnnData object does not already contain a spatial graph.

`build_proportion_matrix(transformation=None)`

Builds the main NMF input matrix. Rows are cells, columns are cell types, and
values are the percentage of each neighboring cell type around each cell. The
result is stored as:

```python
sn.prop_df
```

`transformation="sqrt"` or `"log1p"` can be used to reduce dominance by very
abundant neighbor cell types. Leave it as `None` for raw percentages.

`plot_neighbor_counts()`

Plots the number of spatial neighbors around each cell, grouped by cell type.
By default, cell types are sorted from lowest to highest median neighbor count.

`plot_neighbor_distances()`

Plots the mean distance from each cell to its spatial neighbors, grouped by
cell type. By default, cell types are sorted from lowest to highest median mean
neighbor distance.

Use `sort=False` for either plot if you want to keep the original cell-type
order:

```python
sn.plot_neighbor_counts(sort=False)
sn.plot_neighbor_distances(sort=False)
```

### Choosing K

`find_k_elbow(k_range=range(2, 12))`

Runs NMF for each candidate K and plots reconstruction error. It returns a table
with `k` and `reconstruction_error`. The plot is saved as
`<output_prefix>_nmf_elbow.pdf`.

`grid_search_nmf(k_range=range(3, 8))`

Runs NMF for each candidate K and plots normalized cell-type signature heatmaps.
Use this to inspect whether each K creates interpretable neighborhoods. By
default, the plot is vertical, with cell types on rows and components on
columns. To use the older horizontal orientation:

```python
signature_grid = sn.grid_search_nmf(
    k_range=range(3, 9),
    orientation="horizontal",
)
```

To keep the original cell-type order:

```python
signature_grid = sn.grid_search_nmf(
    k_range=range(3, 9),
    order_cell_types=False,
)
```

The returned dictionary contains the plotted signature table for each K:

```python
signature_grid[6]
```

Each grid-search heatmap is saved as
`<output_prefix>_grid_search_nmf_k<K>.pdf`.

`stability_scan_nmf(k_range=range(2, 11), n_repeats=50, sample_fraction=0.8)`

Runs repeated NMF across K values. For each repeat, cell types are assigned to
their strongest NMF component, and the package records how often cell-type pairs
co-occur in the same component.

This stores:

```python
sn.k_stability_summary
sn.k_stability_pair_table
sn.k_stability_cooccurrence
sn.k_stability_combined_cooccurrence
sn.k_stability_combined_frequency
```

It automatically saves:

```text
figures/<output_prefix>_combined_cooccurrence_heatmap.pdf
figures/<output_prefix>_k_scale_cooccurrence_heatmap.pdf
```

`plot_combined_cooccurrence_heatmap(value="frequency")`

Plots the square cell-type-by-cell-type co-occurrence heatmap. The default
`value="frequency"` uses relative frequency from 0 to 1. Use `value="count"` for
absolute co-occurrence counts.

`plot_k_stability_heatmap()`

Plots a cell-type-pair-by-K heatmap. This is useful for seeing whether a pair is
stable across a range of K values or only appears at one rank.

`plot_k_cooccurrence_heatmap(k=6)`

Plots the square co-occurrence matrix for a single K.

Saved co-occurrence figure names are:

```text
<output_prefix>_combined_cooccurrence_heatmap.pdf
<output_prefix>_k_scale_cooccurrence_heatmap.pdf
<output_prefix>_k6_cooccurrence_heatmap.pdf
```

### Final NMF

`run_nmf(k=None)`

Runs the final NMF model. If `k` is not provided, it uses
`config.n_components`. It returns:

```python
weights, signatures = sn.run_nmf(k=6)
```

It also writes to `adata.obs`:

```text
nmf_dominant
neighborhood_1_k6
neighborhood_2_k6
...
```

`plot_signatures(normalize=False)`

Plots the final NMF signature matrix. Rows are neighborhoods and columns are
neighboring cell types. The plot is saved as `<output_prefix>_nmf_signatures.pdf`
or `<output_prefix>_nmf_signatures_normalized.pdf`.

`plot_neighborhood_pie_charts(top_n=8, min_fraction=0.02)`

Plots pie charts showing the cell-type composition of each NMF neighborhood
signature. The function normalizes each row of the final NMF `H` matrix to sum
to 1, optionally keeps the top contributors, and groups smaller contributors as
`Other` for readability.

```python
pie_table = sn.plot_neighborhood_pie_charts(
    top_n=8,
    min_fraction=0.02,
)
```

The plot is saved as:

```text
<output_prefix>_neighborhood_pie_charts.pdf
```

The returned `pie_table` contains the full ungrouped contribution table:

```text
neighborhood
cell_type
fraction
percent
```

### Sample-Level Summaries

`summarize_samples(clinical_key="cCR", agg_mode="count")`

Aggregates cell-level NMF weights to sample level. Use `agg_mode="count"` for a
simple mean across cells, or `agg_mode="area"` to weight each cell by area.

```python
sample_df = sn.summarize_samples(
    clinical_key="cCR",
    agg_mode="area",
    area_col="cell_area",
)
```

`compare_clinical_groups(sample_summary=None, clinical_key="cCR")`

Backward-compatible wrapper around `compare_neighborhood_proportions()`.

`compare_neighborhood_proportions(sample_summary=None, clinical_key="cCR")`

If the clinical column has exactly two groups, this runs Mann-Whitney U tests
for each neighborhood and applies Benjamini-Hochberg FDR correction. The output
table includes group sizes, means, medians, mean/median differences, U
statistics, raw p-values, and FDR-adjusted p-values.

```python
sample_df = sn.summarize_samples(clinical_key="cCR", agg_mode="area")
pvalues = sn.compare_neighborhood_proportions(
    sample_summary=sample_df,
    clinical_key="cCR",
    group_order=["cCR", "non-cCR"],
    save=True,
    output_dir="results",
)
```

This writes:

```text
results/<output_prefix>_neighborhood_proportion_pvalues.csv
```

`group_order` controls the direction of the reported differences:

```text
mean_diff_group2_minus_group1
median_diff_group2_minus_group1
```

For example, `group_order=["cCR", "non-cCR"]` means positive differences are
higher in `non-cCR`.

`print_neighborhood_proportion_pvalues()`

Convenience method that computes, prints, optionally saves, and returns the same
p-value table.

`plot_sample_level_boxplots(clinical_key="cCR", agg_mode="count")`

Plots sample-level neighborhood weights by clinical group and prints the same
two-group statistical comparison table when applicable. The plot is saved as
`<output_prefix>_sample_level_boxplots_count.pdf` or
`<output_prefix>_sample_level_boxplots_area.pdf`.

### Export Helpers

`save_pipeline_outputs(output_dir="results")`

Saves the proportion matrix, NMF weights, and NMF signatures as CSV files.

Individual export helpers:

```python
sn.export_proportion_matrix("proportion_matrix.csv")
sn.export_nmf_weights("nmf_weights.csv")
sn.export_nmf_signatures("nmf_signatures.csv")
sn.export_cell_assignments("cell_niche_assignments.csv")
```

`cell_assignment_table()`

Returns a per-cell table designed to be useful outside the package. By default,
it contains:

- `cell_id`
- all `adata.obs` metadata columns except duplicated NMF output columns
- the dominant niche assignment
- continuous NMF neighborhood weights

```python
cell_df = sn.cell_assignment_table()
```

If `adata.obs["cell_id"]` exists, that exact column is used for the exported
`cell_id`, so it should match your original Xenium/cell metadata table. If no
`cell_id` column exists, the AnnData index is used.

Because the exported table includes sample/group metadata and the continuous
neighborhood weights, you can make sample-level boxplots from the CSV alone by
grouping the neighborhood weight columns by `sampleID`, `cCR`, or another
clinical column.

`export_cell_assignments(path)`

Writes the same table to CSV:

```python
sn.export_cell_assignments("results/cell_niche_assignments.csv")
```

To export a compact table with only the configured cell type, configured sample
ID, requested metadata, dominant niche, and neighborhood weights:

```python
sn.export_cell_assignments(
    "results/cell_niche_assignments.csv",
    metadata_cols=["cCR", "sample_timepoint"],
    include_obs_metadata=False,
)
```

If your cell ID column has a different name, set it once in the config:

```python
config = SpatialNeighborhoodConfig(
    cell_type_col="cell_type_lvl1",
    cell_id_col="cell_ID",
    sample_col="sampleID",
)
```

### One-Call Pipeline

`run_pipeline(save=True, output_dir="results")`

Runs graph construction, proportion matrix creation, final NMF, metadata
writing, and optional CSV export. Use this after you have already chosen a final
K, or pass `k` directly:

```python
weights, signatures = sn.run_pipeline(k=6, save=True, output_dir="results")
```

## Sample-Level Summaries

After `run_nmf`, aggregate neighborhoods by sample:

```python
sample_df = sn.plot_sample_level_boxplots(
    clinical_key="cCR",
    agg_mode="area",
    area_col="cell_area",
)
```

`agg_mode="count"` computes a simple mean of cell-level weights.
`agg_mode="area"` weights each cell by area, using `cell_area` or falling back
to `area_microns2`.

When the clinical key has exactly two groups, the package prints Mann-Whitney U
tests with Benjamini-Hochberg FDR correction.

## Method Text

Spatial cellular neighborhoods were inferred from cell-type composition around
each cell. A spatial neighbor graph was constructed using Squidpy with either a
fixed radius or fixed nearest-neighbor count, stratified by sample for
multi-sample objects. For each cell, neighboring cells were counted by annotated
cell type and normalized to percentages, producing a cell-by-cell-type
neighborhood composition matrix. Optional square-root or log1p transformation
was applied before min-max scaling. Non-negative matrix factorization was then
used to decompose the matrix into per-cell neighborhood weights and cell-type
signatures. Each cell was assigned to the neighborhood with maximal NMF weight,
while continuous component weights were retained for downstream visualization
and sample-level aggregation. Sample-level neighborhood abundance was computed
as either mean cell-level weight or area-weighted mean and compared across
clinical groups using two-sided Mann-Whitney U tests with FDR correction.

For rank selection, repeated NMF was optionally performed across candidate K
values using random 80% cell subsampling by default. In each run, cell types
were assigned to the component with maximal loading in the NMF signature matrix,
and cell-type-pair co-occurrence frequencies were calculated across repeated
runs. K was selected by jointly considering reconstruction error, the K-scale
co-occurrence heatmap, and biological interpretability of the resulting
neighborhood signatures.
