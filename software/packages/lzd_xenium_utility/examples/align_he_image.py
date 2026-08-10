"""Server example: add an aligned H&E image to Xenium SpatialData."""

import spatialdata as sd

from lzd_xenium_utility import HEAlignmentConfig, HEImageAligner


INPUT_SPATIALDATA_ZARR = "/path/to/xenium_spatialdata.zarr"
HE_IMAGE = "/path/to/pre_he.ome.tif"
ALIGNMENT_ZIP = "/path/to/pre_he_alignment_files.zip"
OUTPUT_SPATIALDATA_ZARR = "/path/to/xenium_spatialdata_with_he.zarr"


sdata = sd.read_zarr(INPUT_SPATIALDATA_ZARR)

config = HEAlignmentConfig(
    image_key="aligned_he",
    alignment_filename="matrix.csv",
    extract_dir="/path/to/alignment_metadata",
    overwrite=False,
)

aligner = HEImageAligner(sdata, config=config)
aligned_he = aligner.align_from_zip(
    he_image_path=HE_IMAGE,
    alignment_zip=ALIGNMENT_ZIP,
)

print(aligned_he)

# sdata is modified in place. Write it to a new Zarr path to keep the raw input.
sdata.write(OUTPUT_SPATIALDATA_ZARR)
