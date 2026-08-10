"""Server example: transfer filtered AnnData annotations into raw Xenium SpatialData."""

import anndata as ad
import spatialdata as sd

from lzd_xenium_utility import AnnDataAnnotationConfig, SpatialDataAnnotationMatcher


RAW_SPATIALDATA_ZARR = "/path/to/raw_xenium_spatialdata.zarr"
FILTERED_ANNOTATED_H5AD = "/path/to/filtered_annotated_combined.h5ad"
OUTPUT_SPATIALDATA_ZARR = "/path/to/filtered_annotated_spatialdata.zarr"

TABLE_KEY = "table"
CELL_ID_KEY = "cell_id"
SAMPLE_KEY = "sample_id"
SAMPLE_ID = "sample_01"  # Use None if the AnnData contains only one sample.


sdata = sd.read_zarr(RAW_SPATIALDATA_ZARR)
combined_adata = ad.read_h5ad(FILTERED_ANNOTATED_H5AD)

matcher = SpatialDataAnnotationMatcher(
    sdata=sdata,
    adata=combined_adata,
    config=AnnDataAnnotationConfig(
        table_key=TABLE_KEY,
        sdata_cell_id_key=CELL_ID_KEY,
        adata_cell_id_key=CELL_ID_KEY,
        sample_key=SAMPLE_KEY,
        sample_id=SAMPLE_ID,
        annotation_columns=None,  # Default: transfer all AnnData obs columns except cell_id.
        unmatched_policy="drop",  # Default: filter SpatialData table to annotated cells.
        overwrite=True,  # Default: replace existing columns with AnnData obs values.
    ),
)

print(matcher.inspect_overlap())
summary = matcher.transfer_annotations()
print(summary)

# sdata is modified in place. Write it to a new Zarr path to keep the raw input.
sdata.write(OUTPUT_SPATIALDATA_ZARR)
