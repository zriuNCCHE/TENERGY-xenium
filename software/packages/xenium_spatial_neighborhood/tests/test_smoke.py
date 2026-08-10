import numpy as np
import pandas as pd
from scipy import sparse

from xenium_spatial_neighborhood import SpatialNeighborhoodAnalyzer, SpatialNeighborhoodConfig


class MiniAnnData:
    def __init__(self):
        self.obs = pd.DataFrame(
            {
                "cell_id": [f"xenium_cell_{i}" for i in range(6)],
                "cell_type_lvl1": ["A", "A", "B", "B", "C", "C"],
                "sampleID": ["S1", "S1", "S1", "S2", "S2", "S2"],
                "cCR": ["cCR", "cCR", "cCR", "non-cCR", "non-cCR", "non-cCR"],
                "cell_area": [1, 2, 1, 2, 1, 2],
            },
            index=[f"cell_{i}" for i in range(6)],
        )
        self.obs_names = self.obs.index
        self.n_obs = len(self.obs)
        rows = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5])
        cols = np.array([1, 2, 0, 2, 0, 1, 4, 5, 3, 5, 3, 4])
        self.obsp = {"spatial_connectivities": sparse.csr_matrix((np.ones_like(rows), (rows, cols)), shape=(6, 6))}
        distances = np.array([1, 2, 1, 2, 2, 2, 4, 5, 4, 5, 5, 5])
        self.obsp["spatial_distances"] = sparse.csr_matrix((distances, (rows, cols)), shape=(6, 6))


def test_build_proportion_matrix_and_run_nmf():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)

    prop = calc.build_proportion_matrix()
    assert prop.shape == (6, 3)
    assert list(prop.columns) == ["A", "B", "C"]

    weights, signatures = calc.run_nmf()
    assert weights.shape == (6, 2)
    assert signatures.shape == (2, 3)
    assert "nmf_dominant" in adata.obs


def test_stability_scan_nmf():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)
    calc.build_proportion_matrix()

    summary, pair_table = calc.stability_scan_nmf(k_range=range(2, 4), n_repeats=2, plot=False, save_figures=False)

    assert summary["k"].tolist() == [2, 3]
    assert {"k", "pair", "cooccurrence_frequency"}.issubset(pair_table.columns)
    assert set(calc.k_stability_cooccurrence) == {2, 3}
    assert calc.k_stability_combined_cooccurrence.shape == (3, 3)


def test_grid_search_orders_cell_types():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)
    calc.build_proportion_matrix()

    signature_grid = calc.grid_search_nmf(k_range=range(2, 3), save=False, show=False)

    assert list(signature_grid) == [2]
    assert sorted(signature_grid[2].columns) == ["A", "B", "C"]


def test_compare_neighborhood_proportions_table():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)
    calc.build_proportion_matrix()
    calc.run_nmf()
    sample_summary = calc.summarize_samples(clinical_key="cCR")

    pvalues = calc.compare_neighborhood_proportions(sample_summary, clinical_key="cCR")

    expected_columns = {
        "neighborhood",
        "clinical_key",
        "group_1",
        "group_2",
        "n_group_1",
        "n_group_2",
        "mean_group_1",
        "mean_group_2",
        "median_group_1",
        "median_group_2",
        "mean_diff_group2_minus_group1",
        "median_diff_group2_minus_group1",
        "mannwhitney_u",
        "p_value",
        "p_adj_fdr_bh",
    }
    assert expected_columns.issubset(pvalues.columns)
    assert len(pvalues) == 2
    assert set(pvalues["clinical_key"]) == {"cCR"}


def test_plot_neighborhood_pie_charts_returns_table():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)
    calc.build_proportion_matrix()
    calc.run_nmf()

    pie_table = calc.plot_neighborhood_pie_charts(save=False, show=False)

    assert {"neighborhood", "cell_type", "fraction", "percent"}.issubset(pie_table.columns)
    assert set(pie_table["neighborhood"]) == set(calc.comp_names)


def test_cell_assignment_table_contains_niche_columns():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)
    calc.build_proportion_matrix()
    calc.run_nmf()

    table = calc.cell_assignment_table()

    assert "cell_id" in table.columns
    assert table["cell_id"].tolist() == adata.obs["cell_id"].tolist()
    assert "cell_type_lvl1" in table.columns
    assert "sampleID" in table.columns
    assert "cCR" in table.columns
    assert "cell_area" in table.columns
    assert "nmf_dominant" in table.columns
    assert set(calc.comp_names).issubset(table.columns)
    assert len(table) == adata.n_obs


def test_cell_assignment_table_can_keep_minimal_metadata():
    adata = MiniAnnData()
    config = SpatialNeighborhoodConfig(sample_col="sampleID", n_components=2, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)
    calc.build_proportion_matrix()
    calc.run_nmf()

    table = calc.cell_assignment_table(metadata_cols=["cCR"], include_obs_metadata=False)

    assert table["cell_id"].tolist() == adata.obs["cell_id"].tolist()
    assert "cell_type_lvl1" in table.columns
    assert "sampleID" in table.columns
    assert "cCR" in table.columns
    assert "cell_area" not in table.columns


def test_neighbor_quality_plots_sort_lowest_left():
    adata = MiniAnnData()
    rows = np.array([0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5])
    cols = np.array([1, 0, 2, 0, 1, 2, 4, 5, 3, 2, 5, 3, 4, 2])
    adata.obsp["spatial_connectivities"] = sparse.csr_matrix((np.ones_like(rows), (rows, cols)), shape=(6, 6))
    config = SpatialNeighborhoodConfig(sample_col="sampleID", min_cells_for_plot=1, random_state=0)
    calc = SpatialNeighborhoodAnalyzer(adata, config=config)

    count_ax = calc.plot_neighbor_counts(save=False, show=False)
    distance_ax = calc.plot_neighbor_distances(save=False, show=False)

    assert [tick.get_text() for tick in count_ax.get_xticklabels()] == ["A", "B", "C"]
    assert [tick.get_text() for tick in distance_ax.get_xticklabels()] == ["A", "B", "C"]
