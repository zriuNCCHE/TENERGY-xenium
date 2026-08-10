import sys
import types
import zipfile

import pytest

from lzd_xenium_utility import HEAlignmentConfig, HEImageAligner


class MiniSpatialData:
    def __init__(self):
        self.images = {}
        self.attrs = {}


def install_fake_spatialdata_io(monkeypatch):
    module = types.ModuleType("spatialdata_io")

    def xenium_aligned_image(image_path, alignment_file):
        return {
            "image_path": str(image_path),
            "alignment_file": str(alignment_file),
        }

    module.xenium_aligned_image = xenium_aligned_image
    monkeypatch.setitem(sys.modules, "spatialdata_io", module)


def test_align_from_files_adds_image_and_provenance(tmp_path, monkeypatch):
    install_fake_spatialdata_io(monkeypatch)
    image = tmp_path / "he.ome.tif"
    matrix = tmp_path / "matrix.csv"
    image.write_text("fake image")
    matrix.write_text("1,0,0\n0,1,0\n0,0,1\n")
    sdata = MiniSpatialData()

    aligner = HEImageAligner(sdata, config=HEAlignmentConfig(image_key="aligned_he"))
    aligned = aligner.align_from_files(image, matrix)

    assert sdata.images["aligned_he"] == aligned
    assert aligned["image_path"] == str(image)
    assert aligned["alignment_file"] == str(matrix)
    assert aligner.provenance[-1]["image_key"] == "aligned_he"
    assert sdata.attrs["lzd_xenium_utility"]["he_alignment"]


def test_align_from_zip_extracts_matrix_and_adds_image(tmp_path, monkeypatch):
    install_fake_spatialdata_io(monkeypatch)
    image = tmp_path / "he.ome.tif"
    image.write_text("fake image")
    alignment_zip = tmp_path / "alignment.zip"
    with zipfile.ZipFile(alignment_zip, "w") as zip_ref:
        zip_ref.writestr("nested/matrix.csv", "1,0,0\n0,1,0\n0,0,1\n")

    sdata = MiniSpatialData()
    aligner = HEImageAligner(sdata)
    aligned = aligner.align_from_zip(
        he_image_path=image,
        alignment_zip=alignment_zip,
        extract_dir=tmp_path / "extracted",
    )

    assert (tmp_path / "extracted" / "matrix.csv").exists()
    assert sdata.images["aligned_he"] == aligned
    assert aligner.provenance[-1]["source"] == "zip"
    assert aligner.provenance[-1]["alignment_zip"] == str(alignment_zip)


def test_duplicate_image_key_requires_overwrite(tmp_path, monkeypatch):
    install_fake_spatialdata_io(monkeypatch)
    image = tmp_path / "he.ome.tif"
    matrix = tmp_path / "matrix.csv"
    image.write_text("fake image")
    matrix.write_text("matrix")
    sdata = MiniSpatialData()
    sdata.images["aligned_he"] = "existing"

    aligner = HEImageAligner(sdata)

    with pytest.raises(KeyError, match="overwrite=True"):
        aligner.align_from_files(image, matrix)

    aligned = aligner.align_from_files(image, matrix, overwrite=True)
    assert sdata.images["aligned_he"] == aligned


def test_missing_matrix_inside_zip_raises_clear_error(tmp_path):
    alignment_zip = tmp_path / "alignment.zip"
    with zipfile.ZipFile(alignment_zip, "w") as zip_ref:
        zip_ref.writestr("other.csv", "content")

    aligner = HEImageAligner(MiniSpatialData())

    with pytest.raises(FileNotFoundError, match="matrix.csv"):
        aligner.extract_alignment_file(alignment_zip, extract_dir=tmp_path / "out")


def test_requires_spatialdata_like_images_mapping():
    with pytest.raises(TypeError, match="images"):
        HEImageAligner(object())

