import pandas as pd
import pytest

from lzd_xenium_utility import AnnDataAnnotationConfig, SpatialDataAnnotationMatcher


class MiniAnnData:
    def __init__(self, obs):
        self.obs = obs

    def __getitem__(self, item):
        if isinstance(item, tuple):
            row_selector = item[0]
        else:
            row_selector = item
        return MiniAnnData(self.obs.iloc[row_selector].copy())

    def copy(self):
        return MiniAnnData(self.obs.copy())


class MiniSpatialData:
    def __init__(self, table):
        self.tables = {"table": table}
        self.attrs = {}


def make_table():
    obs = pd.DataFrame(
        {
            "cell_id": ["cell_1", "cell_2", "cell_3", "cell_4"],
            "area": [10, 20, 30, 40],
        },
        index=["row_1", "row_2", "row_3", "row_4"],
    )
    return MiniAnnData(obs)


def make_adata():
    obs = pd.DataFrame(
        {
            "cell_id": ["cell_1", "cell_3"],
            "cell_type": ["T cell", "Epithelial"],
            "leiden": ["0", "1"],
        },
        index=["a", "b"],
    )
    return MiniAnnData(obs)


def test_transfer_annotations_drops_unmatched_cells_by_default():
    sdata = MiniSpatialData(make_table())
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        make_adata(),
        config=AnnDataAnnotationConfig(annotation_columns=("cell_type",)),
    )

    summary = matcher.transfer_annotations()

    assert list(sdata.tables["table"].obs["cell_id"]) == ["cell_1", "cell_3"]
    assert list(sdata.tables["table"].obs["cell_type"]) == ["T cell", "Epithelial"]
    assert list(sdata.tables["table"].obs["annotation_match_status"]) == ["annotated", "annotated"]
    assert summary["sdata_cells"] == 4
    assert summary["matched_cells"] == 2
    assert summary["sdata_only_cells"] == 2
    assert summary["unmatched_policy"] == "drop"
    assert sdata.attrs["lzd_xenium_utility"]["anndata_annotation"]


def test_default_transfers_all_source_obs_columns_except_cell_id():
    sdata = MiniSpatialData(make_table())
    matcher = SpatialDataAnnotationMatcher(sdata, make_adata())

    matcher.transfer_annotations()

    assert list(sdata.tables["table"].obs["cell_type"]) == ["T cell", "Epithelial"]
    assert list(sdata.tables["table"].obs["leiden"]) == ["0", "1"]


def test_transfer_annotations_can_keep_and_mark_filtered_out_cells():
    sdata = MiniSpatialData(make_table())
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        make_adata(),
        config=AnnDataAnnotationConfig(
            annotation_columns=("cell_type",),
            unmatched_policy="keep",
        ),
    )

    matcher.transfer_annotations()

    assert list(sdata.tables["table"].obs["cell_id"]) == ["cell_1", "cell_2", "cell_3", "cell_4"]
    assert pd.isna(sdata.tables["table"].obs.loc["row_2", "cell_type"])
    assert sdata.tables["table"].obs.loc["row_2", "annotation_match_status"] == "filtered_out"
    assert sdata.tables["table"].obs.loc["row_4", "annotation_match_status"] == "filtered_out"


def test_sample_id_filters_combined_adata_before_matching():
    sdata = MiniSpatialData(make_table())
    combined_obs = pd.DataFrame(
        {
            "cell_id": ["cell_1", "cell_1", "cell_3"],
            "sample_id": ["sample_a", "sample_b", "sample_a"],
            "cell_type": ["T cell", "B cell", "Epithelial"],
        }
    )
    combined_adata = MiniAnnData(combined_obs)
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        combined_adata,
        config=AnnDataAnnotationConfig(
            sample_id="sample_a",
            annotation_columns=("cell_type",),
        ),
    )

    matcher.transfer_annotations()

    assert list(sdata.tables["table"].obs["cell_id"]) == ["cell_1", "cell_3"]
    assert list(sdata.tables["table"].obs["cell_type"]) == ["T cell", "Epithelial"]


def test_duplicate_source_cell_ids_raise_clear_error():
    sdata = MiniSpatialData(make_table())
    adata = MiniAnnData(
        pd.DataFrame(
            {
                "cell_id": ["cell_1", "cell_1"],
                "cell_type": ["T cell", "B cell"],
            }
        )
    )
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        adata,
        config=AnnDataAnnotationConfig(annotation_columns=("cell_type",)),
    )

    with pytest.raises(ValueError, match="Duplicate examples"):
        matcher.transfer_annotations()


def test_existing_destination_column_is_overwritten_by_default():
    table = make_table()
    table.obs["cell_type"] = "old"
    sdata = MiniSpatialData(table)
    matcher = SpatialDataAnnotationMatcher(sdata, make_adata())

    matcher.transfer_annotations()

    assert list(sdata.tables["table"].obs["cell_type"]) == ["T cell", "Epithelial"]


def test_existing_destination_column_can_require_explicit_overwrite():
    table = make_table()
    table.obs["cell_type"] = "old"
    sdata = MiniSpatialData(table)
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        make_adata(),
        config=AnnDataAnnotationConfig(annotation_columns=("cell_type",), overwrite=False),
    )

    with pytest.raises(KeyError, match="overwrite=True"):
        matcher.transfer_annotations()

    matcher.transfer_annotations(overwrite=True)
    assert list(sdata.tables["table"].obs["cell_type"]) == ["T cell", "Epithelial"]


def test_can_match_by_obs_names_when_cell_id_key_is_none():
    table = MiniAnnData(pd.DataFrame({"area": [10, 20, 30]}, index=["cell_1", "cell_2", "cell_3"]))
    adata = MiniAnnData(pd.DataFrame({"cell_type": ["T cell", "B cell"]}, index=["cell_1", "cell_3"]))
    sdata = MiniSpatialData(table)
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        adata,
        config=AnnDataAnnotationConfig(
            sdata_cell_id_key=None,
            adata_cell_id_key=None,
            annotation_columns=("cell_type",),
        ),
    )

    matcher.transfer_annotations()

    assert list(sdata.tables["table"].obs.index) == ["cell_1", "cell_3"]
    assert list(sdata.tables["table"].obs["cell_type"]) == ["T cell", "B cell"]


def test_warns_when_adata_has_cells_missing_from_sdata():
    sdata = MiniSpatialData(make_table())
    adata = MiniAnnData(
        pd.DataFrame(
            {
                "cell_id": ["cell_1", "cell_3", "cell_999"],
                "cell_type": ["T cell", "Epithelial", "Unknown"],
            }
        )
    )
    matcher = SpatialDataAnnotationMatcher(
        sdata,
        adata,
        config=AnnDataAnnotationConfig(annotation_columns=("cell_type",)),
    )

    with pytest.warns(UserWarning, match="not present in the SpatialData table"):
        summary = matcher.transfer_annotations()

    assert summary["adata_only_cells"] == 1
    assert list(sdata.tables["table"].obs["cell_id"]) == ["cell_1", "cell_3"]
