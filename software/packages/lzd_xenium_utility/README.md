# LZD Xenium Utility

Reusable Python utilities for Xenium `SpatialData` processing.

This package is designed to grow as a collection of focused Xenium helper
modules. It follows the same compact server-friendly style as the Xenium
spatial-neighborhood package: configure a dataclass, instantiate a utility
class, run notebook-friendly methods, and keep installation simple.

## What This Package Does

The package currently provides two focused utilities.

### Add an aligned H&E image

The first utility adds an aligned H&E image to an existing Xenium `SpatialData`
object:

1. accepts an H&E OME-TIFF,
2. accepts either an alignment CSV or an alignment ZIP containing `matrix.csv`,
3. calls `spatialdata_io.xenium_aligned_image`,
4. writes the aligned image to `sdata.images`,
5. records lightweight provenance when possible.

Main entry point:

```python
from lzd_xenium_utility import HEAlignmentConfig, HEImageAligner
```

### Transfer filtered AnnData annotations

The second utility transfers cell annotations from a filtered, annotated
`AnnData` object into the table of a raw Xenium `SpatialData` object. By
default, it uses all `adata.obs` columns except the cell ID column, overwrites
matching columns in `sdata.tables["table"].obs`, and drops raw SpatialData
cells that are absent from the annotated AnnData. Set `unmatched_policy="keep"`
to retain raw cells and mark them as `filtered_out`.

Main entry point:

```python
from lzd_xenium_utility import AnnDataAnnotationConfig, SpatialDataAnnotationMatcher
```

## Quick Start

```python
import spatialdata as sd
from lzd_xenium_utility import HEAlignmentConfig, HEImageAligner

sdata = sd.read_zarr("/path/to/xenium_spatialdata.zarr")

config = HEAlignmentConfig(image_key="aligned_he")
aligner = HEImageAligner(sdata, config=config)

aligned_he = aligner.align_from_zip(
    he_image_path="/path/to/pre_he.ome.tif",
    alignment_zip="/path/to/pre_he_alignment_files.zip",
)

print(sdata.images["aligned_he"])
print(aligned_he is sdata.images["aligned_he"])
```

If the alignment matrix is already extracted:

```python
aligned_he = aligner.align_from_files(
    he_image_path="/path/to/pre_he.ome.tif",
    alignment_file="/path/to/matrix.csv",
    image_key="aligned_he",
)
```

Transfer annotations from a filtered AnnData object:

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
        annotation_columns=None,  # default: transfer all AnnData obs columns except cell_id
    ),
)

print(matcher.inspect_overlap())
summary = matcher.transfer_annotations()
print(summary)
```

For a combined multi-sample AnnData, specify the sample to pull annotations
from:

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

matcher.transfer_annotations()
```

To keep raw cells that were filtered out of the annotated AnnData:

```python
matcher = SpatialDataAnnotationMatcher(
    sdata=sdata,
    adata=adata,
    config=AnnDataAnnotationConfig(
        annotation_columns=("cell_type",),
        unmatched_policy="keep",
    ),
)

matcher.transfer_annotations()
```

## Main Outputs

The H&E aligner modifies the provided `SpatialData` object in place. It does not
create a second full `SpatialData` object. The main result is:

```python
sdata.images["aligned_he"]
```

The return value from `align_from_zip()` or `align_from_files()` is the aligned
image element that was inserted into `sdata.images`:

```python
aligned_he = aligner.align_from_zip(...)
aligned_he is sdata.images["aligned_he"]  # True
```

So downstream workflows should usually keep using the same `sdata` object, now
with the aligned H&E image available under the configured image key.

It also keeps a copy of the latest records in:

```python
aligner.provenance
```

When `sdata.attrs` is available, provenance is also appended under:

```python
sdata.attrs["lzd_xenium_utility"]["he_alignment"]
```

The annotation matcher also modifies `sdata` in place. In default
`unmatched_policy="drop"` mode, it replaces `sdata.tables[table_key]` with a
filtered table containing only matched cells. In `unmatched_policy="keep"` mode,
it keeps the table rows and updates `.obs` columns in place. Write to a new Zarr
path if you want to keep the raw input unchanged.

## Installation

On an HPC or server, copy this source folder and install it directly:

```bash
cd /path/to/lzd-xenium-utility
python3 -m pip install -e . --no-build-isolation
```

No wheel build is required for normal use.

## Project Layout

```text
src/lzd_xenium_utility/
  he_alignment.py
  anndata_annotation.py
docs/
  user_guide.md
  design_principles.md
examples/
  align_he_image.py
  transfer_anndata_annotations.py
tests/
  test_he_alignment.py
  test_anndata_annotation.py
```
