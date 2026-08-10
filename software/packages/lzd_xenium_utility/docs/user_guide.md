# User Guide

## Purpose

`lzd-xenium-utility` provides reusable helper classes for Xenium processing.
The first utility aligns an H&E image to an existing Xenium `SpatialData`
object using the alignment matrix produced by the Xenium workflow. The second
utility transfers annotations from a filtered AnnData object into the cell
table of a raw Xenium `SpatialData` object.

## Input Requirements

For H&E alignment, provide:

- an existing `SpatialData` object with an `.images` mapping,
- an H&E OME-TIFF image,
- either an alignment CSV such as `matrix.csv` or a ZIP containing that file.

For AnnData annotation transfer, provide:

- an existing `SpatialData` object with a `.tables` mapping,
- a table such as `sdata.tables["table"]` with `.obs`,
- a filtered, annotated `AnnData` object with `.obs`,
- matching cell IDs, defaulting to a `cell_id` column in both objects.

## Standard Workflow

```python
import spatialdata as sd
from lzd_xenium_utility import HEAlignmentConfig, HEImageAligner

sdata = sd.read_zarr("/path/to/xenium_spatialdata.zarr")

aligner = HEImageAligner(
    sdata,
    config=HEAlignmentConfig(image_key="pre_treatment_he"),
)

aligned_he = aligner.align_from_zip(
    he_image_path="/path/to/pre_treatment_he.ome.tif",
    alignment_zip="/path/to/pre_treatment_he_alignment_files.zip",
)
```

This modifies `sdata` in place and adds:

```python
sdata.images["pre_treatment_he"]
```

The returned `aligned_he` is the image element that was added to the existing
`SpatialData` object:

```python
aligned_he is sdata.images["pre_treatment_he"]  # True
```

It is not a new full `SpatialData` object. Keep using the original `sdata`
object for downstream SpatialData plotting, saving, or analysis.

## Return Value And Data Structure

`HEImageAligner` has two effects:

1. it returns the aligned image element created by
   `spatialdata_io.xenium_aligned_image()`,
2. it stores that same image element in `sdata.images[image_key]`.

After alignment, the relevant `SpatialData` structure is:

```text
sdata
  images
    pre_treatment_he  # aligned H&E image element
  shapes
  points
  tables
  labels
```

The exact image element type is controlled by `spatialdata_io` and may vary
with the installed SpatialData version and the input image. In normal use, treat
it as a SpatialData-compatible image element with the transformation needed to
place the H&E image into the Xenium coordinate system.

## Direct Alignment Matrix

If the matrix has already been extracted:

```python
aligned_he = aligner.align_from_files(
    he_image_path="/path/to/pre_treatment_he.ome.tif",
    alignment_file="/path/to/matrix.csv",
    image_key="pre_treatment_he",
)
```

## Overwriting Existing Images

By default, the aligner refuses to overwrite an existing image key:

```python
aligner.align_from_zip(..., image_key="aligned_he")
```

If replacement is intentional:

```python
aligner.align_from_zip(..., image_key="aligned_he", overwrite=True)
```

## Provenance

The aligner stores lightweight records in:

```python
aligner.provenance
```

If `sdata.attrs` exists, the same records are appended under:

```python
sdata.attrs["lzd_xenium_utility"]["he_alignment"]
```

## Troubleshooting

`FileNotFoundError` means the image, ZIP, or alignment file path is wrong, or
the ZIP does not contain the configured alignment filename.

`KeyError` for an image key means `sdata.images` already contains that name.
Use a new `image_key` or pass `overwrite=True`.

`ImportError` for `spatialdata-io` means the server environment needs the
Xenium SpatialData dependencies installed.

## AnnData Annotation Transfer

Use `SpatialDataAnnotationMatcher` when the raw Xenium Zarr contains all cells
from the instrument output, but your annotated AnnData contains only filtered
cells.

```python
import anndata as ad
import spatialdata as sd
from lzd_xenium_utility import AnnDataAnnotationConfig, SpatialDataAnnotationMatcher

sdata = sd.read_zarr("/path/to/raw_xenium_spatialdata.zarr")
adata = ad.read_h5ad("/path/to/filtered_annotated_cells.h5ad")

matcher = SpatialDataAnnotationMatcher(
    sdata=sdata,
    adata=adata,
    config=AnnDataAnnotationConfig(
        table_key="table",
        sdata_cell_id_key="cell_id",
        adata_cell_id_key="cell_id",
        annotation_columns=None,
    ),
)

print(matcher.inspect_overlap())
summary = matcher.transfer_annotations()
```

The default behavior is:

```python
unmatched_policy="drop"
```

This means all `adata.obs` columns except the AnnData cell ID column are
transferred into `sdata.tables["table"].obs`, and raw SpatialData cells that do
not appear in the filtered AnnData are removed from `sdata.tables["table"]`.
The resulting table follows the filtered, annotated AnnData cell set while
preserving the raw SpatialData table order for matched cells.

`sdata` is modified in place. In default `drop` mode, the utility replaces
`sdata.tables["table"]` with a filtered table copy. Write to a new Zarr path if
you want to keep the original raw Zarr unchanged.

To keep raw cells and mark cells absent from the filtered AnnData:

```python
matcher = SpatialDataAnnotationMatcher(
    sdata=sdata,
    adata=adata,
    config=AnnDataAnnotationConfig(
        annotation_columns=("cell_type", "cell_subtype"),
        unmatched_policy="keep",
    ),
)

matcher.transfer_annotations()
```

In keep mode, unmatched raw cells remain in the table and receive:

```python
sdata.tables["table"].obs["annotation_match_status"] == "filtered_out"
```

Matched cells receive:

```python
sdata.tables["table"].obs["annotation_match_status"] == "annotated"
```

## Combined Multi-Sample AnnData

If the AnnData object contains multiple samples, specify the sample column and
sample ID. The utility filters AnnData first, then matches cell IDs.

```python
matcher = SpatialDataAnnotationMatcher(
    sdata=sdata,
    adata=combined_adata,
    config=AnnDataAnnotationConfig(
        sample_key="sample_id",
        sample_id="sample_01",
        annotation_columns=None,
    ),
)

summary = matcher.transfer_annotations()
```

If `sample_id=None`, the AnnData is treated as a single-sample object and is
not filtered by sample. If duplicated AnnData cell IDs remain after sample
filtering, the utility raises an error rather than guessing which annotation to
use.

If the filtered AnnData contains cell IDs that are not present in the SpatialData
table, the utility emits a warning. Those AnnData-only cells cannot be
transferred because there is no destination SpatialData cell for them.

## Server Copy-Paste Examples

The `examples/` folder contains one runnable server-style script per utility:

```text
examples/align_he_image.py
examples/transfer_anndata_annotations.py
```

Each script keeps editable paths and keys at the top, modifies `sdata` in
memory, then writes a new output Zarr.
