# Xenium Metaprogram Analysis

This project refactors the previous `MalignantMetaProgrammer` script into a small reusable Python package for ESCC Xenium 5k malignant-cell meta-program analysis.

## Context

The method is inspired by the multi-K NMF meta-program workflow from the Nature paper referenced in this project. The key Xenium-specific adaptation is per-sample positive centering:

```text
V_centered = max(0, V - mean(V_sample))
```

This is intended to suppress patient-specific CNV/clonal background so NMF focuses on intra-tumour functional variation.

## Project Layout

```text
src/xenium_metaprogram/
  malignant_metaprogrammer.py   # main analysis class
examples/
  run_tenergy_mp.py             # example runner for the current ESCC dataset
references/
  original_MP.py                # untouched copy of previous script
  gemini_handoff.txt            # Gemini development notes
  example_similarity_heatmap.png
```

## Recommended TENERGY Defaults

- `k_range=(7, 15)`
- `top_n_genes=50`
- `n_mps=10` to `15`
- `consensus_fraction=0.20`
- `center_first=True`

The expected biological anchor checks are:

- Proliferation: `MKI67`, `TOP2A`
- Keratinization: `KRT1`, `KRT10`, `IVL`
- Invasive/resistance front: `MMP14`, `COL17A1`, `LAMA3`, `LAMC2`, `CD44`

## Example

```python
import scanpy as sc
from xenium_metaprogram import MalignantMetaProgrammer, MetaProgramConfig

adata = sc.read_h5ad("../../data/combined_adata/combined_all_patients_epi_analyzed_final.h5ad")
adata = adata[adata.obs["sample_timepoint"] == "pre"].copy()
adata = adata[adata.obs["final_cell_type"] == "malignant"].copy()

config = MetaProgramConfig(sample_col="sampleID", layer="counts")
mp = MalignantMetaProgrammer(adata, config=config)
mp.run_pipeline(k_range=(7, 15), top_n_genes=50)
```

By default, `run_pipeline(save=True)` writes these files to the current working
directory:

```text
xenium_metaprogram_genes.csv
xenium_metaprogram_annotations.csv
xenium_metaprogram_similarity_heatmap.pdf
```

You can customize the output location and prefix:

```python
mp.run_pipeline(
    k_range=(7, 15),
    top_n_genes=50,
    save=True,
    output_dir="results",
    output_prefix="tenergy_pre_malignant",
)
```

## Similarity Heatmap

`run_pipeline(plot_results=True)` now plots a labeled similarity heatmap by default
using the paper-like `"nature"` colormap. The available named colormap options are:

```python
mp.available_heatmap_colormaps()
```

You can replot the heatmap without rerunning NMF:

```python
mp.plot_similarity_heatmap(cmap="nature")
mp.plot_similarity_heatmap(cmap="mako")
mp.plot_similarity_heatmap(cmap="magma", vmax=0.4)
```

The colored side/top bars mark retained consensus MPs. Gray blocks are raw NMF
clusters that were present in the similarity matrix but did not pass the
meta-program filters.

## Naming Meta-Programs

After running the pipeline, suggested biological names can be generated from
overlap with built-in marker panels:

```python
mp_annotations = mp.annotate_meta_programs()
mp_annotations
```

The table includes the suggested label, marker overlap genes, and the first
consensus genes for each MP. Export it with:

```python
mp.export_meta_program_annotations("tenergy_mp_annotations.csv")
```

When `run_pipeline(save=True)`, this annotation table is exported automatically
as `{output_prefix}_annotations.csv`.

These labels are intended as starting points. For manuscript figures, confirm
each name by inspecting consensus genes, spatial distribution, and clinical or
histological associations.

## Algorithmic Methodology

This package identifies recurrent malignant-cell expression programs from
Xenium targeted-transcriptomics data using a multi-rank, per-sample NMF
meta-program workflow. The approach is designed for cohorts in which malignant
cells contain both sample-specific clonal/CNV-like expression structure and
within-tumor functional variation.

### Input and malignant-cell subset

The method expects an `AnnData` object containing malignant cells, gene names in
`adata.var_names`, sample identifiers in `adata.obs`, and a non-negative
expression matrix in either `adata.layers[layer]` or `adata.X`. In the TENERGY
analysis, the input is restricted to pre-treatment malignant epithelial cells
before running the pipeline.

### Expression scale and layer choice

For meta-program discovery, the input AnnData should already have been
normalized and log-transformed with `log1p`. In the current TENERGY objects, the
configured layer is named `"counts"`, but this layer contains log1p-transformed
values rather than raw molecule counts. This naming convention is easy to forget,
so check the object provenance before running the pipeline.

The recommended configuration for the current project is:

```python
config = MetaProgramConfig(sample_col="sampleID", layer="counts")
```

With `center_first=True`, the package creates `adata.layers["centered_pos"]`
from this log1p layer and runs NMF on the positive-centered matrix. The resulting
programs therefore reflect positive within-sample deviations in log-transformed
expression.

Raw counts are not the intended input for this implementation unless the method
is deliberately being re-benchmarked. In the current implementation, final
meta-program gene-set scores are computed with Scanpy's `score_genes` default
expression source, which is typically `adata.X` unless `adata.raw` is used by
Scanpy. For consistency, `adata.X` should also represent normalized/log1p
expression when interpreting the final MP scores.

### Per-sample positive centering

Before NMF, the package optionally constructs a sample-centered positive layer:

```text
V_centered[s, g] = max(0, V[s, g] - mean_sample(V[:, g]))
```

where `V[s, g]` is the expression of gene `g` in cell `s`, and
`mean_sample(V[:, g])` is the mean expression of gene `g` across malignant cells
from the same sample. Negative centered values are truncated to zero so that the
matrix remains compatible with NMF. This step is intended to reduce
patient-specific baseline programs, including broad clonal or CNV-associated
expression shifts, while preserving positive within-sample deviations that may
represent functional malignant-cell states.

### Multi-rank NMF within each sample

For each sample independently, NMF is run across a user-specified range of
factorization ranks, typically `k=7..15` for the current Xenium ESCC analysis.
For each rank `k`, the expression matrix is decomposed as:

```text
V_sample ~= W_sample,k H_sample,k
```

where rows of `H_sample,k` represent gene-loading vectors for NMF components.
For every component, the top `top_n_genes` genes by loading are retained as an
NMF program. Genes with zero expression within a sample are removed before NMF,
and samples with too few cells or too few nonzero genes are skipped.

### Program similarity

All NMF programs from all samples and all tested ranks are compared by Jaccard
similarity over their top-gene sets:

```text
similarity(A, B) = |genes(A) intersect genes(B)| / |genes(A) union genes(B)|
```

This produces a program-by-program similarity matrix. The matrix is used both
for clustering recurrent programs and for the diagnostic similarity heatmap.

### Meta-program clustering and filtering

The package converts similarity to distance as `1 - similarity` and performs
average-linkage hierarchical clustering. Clusters are cut into `n_mps` candidate
groups using a maximum-cluster criterion. Each candidate group is then filtered
for recurrence across samples: it must contain NMF programs from at least
`min_patients` distinct samples.

For retained clusters, consensus genes are defined by recurrence within the
cluster. A gene is included in the meta-program if it appears in at least:

```text
ceil(consensus_fraction * number_of_programs_in_cluster)
```

programs from that cluster. Candidate clusters without consensus genes after
this filter are discarded. Final retained clusters are renamed compactly as
`MP_1`, `MP_2`, etc.; the raw hierarchical cluster IDs are preserved in
`program_table["raw_cluster"]`.

### Cell-level scoring

Each retained meta-program is scored across all cells with Scanpy's
`score_genes`, using the consensus gene set present in `adata.var_names`. Scores
are written to `adata.obs` under the corresponding meta-program names, such as
`MP_1`, `MP_2`, and so on. These scores can be used for downstream spatial
visualization, clinical-response comparisons, differential abundance-style
summaries, or association with histological regions.

### Recommended reporting text

For manuscript drafting, the method can be summarized as follows:

> We identified recurrent malignant-cell meta-programs using a sample-aware
> multi-rank NMF workflow. For each sample, malignant-cell expression values were
> first taken from a normalized log1p expression layer and centered gene-wise by
> subtracting the sample-specific malignant-cell mean, followed by truncation of
> negative values to zero. This positive-centered matrix was used to reduce
> sample-specific baseline expression programs while retaining within-sample
> positive deviations in log-transformed expression. NMF was then performed
> independently within each sample over multiple ranks. For each NMF component,
> the top-loading genes were retained as a candidate program. Candidate programs
> from all samples and ranks were compared using Jaccard similarity of their
> top-gene sets, clustered by average-linkage hierarchical clustering, and
> filtered for recurrence across samples. Consensus genes for each retained
> cluster were defined by recurrent membership among programs in that cluster.
> The resulting consensus gene sets were used as malignant-cell meta-programs and
> scored in each cell using gene-set scoring.

## First Improvements Over The Original Script

- Removed executable usage code from the package module.
- Added `MetaProgramConfig` for reproducible parameters.
- Added minimum sample/gene checks before NMF.
- Handles sparse matrices explicitly.
- Filters zero-count genes per sample before NMF.
- Adds a reusable zero-variance filter before correlation on sparse Xenium subsets.
- Adds optional RAPIDS/cuML NMF support when available.
- Keeps program metadata in `program_table` for later debugging and summaries.

## Next Development Targets

1. Add a notebook that compares raw counts vs. `centered_pos` on one test cohort.
2. Add stability summaries across `k`, sample, and cluster.
3. Add export of MP score plots, violin plots by clinical response, and spatial panels.
4. Add tests using a small synthetic `AnnData` object so refactors cannot silently break the pipeline.

## Troubleshooting

On restricted machines, Scanpy/Matplotlib/Numba may need writable cache directories:

```bash
export MPLCONFIGDIR=/tmp/matplotlib-cache
export NUMBA_CACHE_DIR=/tmp/numba-cache
```
