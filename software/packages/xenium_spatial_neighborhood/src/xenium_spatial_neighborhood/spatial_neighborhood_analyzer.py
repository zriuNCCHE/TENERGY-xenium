from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy import sparse, stats
from sklearn.decomposition import NMF
from sklearn.preprocessing import MinMaxScaler
from statsmodels.stats.multitest import multipletests


Transformation = Literal["sqrt", "log1p"] | None
AggregationMode = Literal["count", "area"]

DEFAULT_HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "warm_amber",
    ["#fffdf0", "#fee8a8", "#fdbb6d", "#e34a33", "#7f0000"],
)
DEFAULT_LINE_COLOR = "teal"
DEFAULT_POINT_COLOR = "#2a9d8f"
DEFAULT_GROUP_PALETTE = "Set2"


@dataclass
class SpatialNeighborhoodConfig:
    """Default parameters for Xenium spatial cellular-neighborhood analysis."""

    cell_type_col: str = "cell_type_lvl1"
    cell_id_col: str | None = "cell_id"
    sample_col: str | None = None
    radius: float | None = 50.0
    n_neighbors: int | None = None
    coord_type: str = "generic"
    library_key: str | None = None
    transformation: Transformation = None
    n_components: int = 6
    random_state: int = 42
    nmf_init: str = "nndsvda"
    nmf_max_iter: int = 2000
    elbow_max_iter: int = 1000
    dominant_col: str = "nmf_dominant"
    min_cells_for_plot: int = 20
    max_plot_cells_per_type: int = 2000
    area_col: str = "cell_area"
    alternative_area_col: str = "area_microns2"
    output_prefix: str = "spatial_neighborhood"
    figure_dir: str = "figures"


class SpatialNeighborhoodAnalyzer:
    """
    Identify spatial cellular neighborhoods from local cell-type composition.

    The analyzer expects an AnnData object with cell-type labels in ``adata.obs``
    and spatial coordinates compatible with ``squidpy.gr.spatial_neighbors``.
    NMF weights and dominant neighborhood labels are written back to
    ``adata.obs``.
    """

    def __init__(
        self,
        adata,
        cell_type_col: str | None = None,
        sample_col: str | None = None,
        config: SpatialNeighborhoodConfig | None = None,
    ):
        self.config = config or SpatialNeighborhoodConfig()
        if cell_type_col is not None:
            self.config.cell_type_col = cell_type_col
        if sample_col is not None:
            self.config.sample_col = sample_col

        self.adata = adata
        self.cell_type_col = self.config.cell_type_col
        self.sample_col = self.config.sample_col

        self.prop_df: pd.DataFrame | None = None
        self.nmf_model: NMF | None = None
        self.W: np.ndarray | None = None
        self.H: np.ndarray | None = None
        self.comp_names: list[str] = []
        self.signature_table: pd.DataFrame | None = None
        self.weight_table: pd.DataFrame | None = None
        self.sample_summary_table: pd.DataFrame | None = None
        self.k_stability_summary: pd.DataFrame | None = None
        self.k_stability_pair_table: pd.DataFrame | None = None
        self.k_stability_cooccurrence: dict[int, pd.DataFrame] = {}
        self.k_stability_cooccurrence_counts: dict[int, pd.DataFrame] = {}
        self.k_stability_combined_cooccurrence: pd.DataFrame | None = None
        self.k_stability_combined_frequency: pd.DataFrame | None = None
        self.figure_paths: dict[str, Path] = {}
        self.k_stability_figure_paths: dict[str, Path] = {}

        self._validate_inputs(require_graph=False)

    def _validate_inputs(self, require_graph: bool = False) -> None:
        if self.cell_type_col not in self.adata.obs:
            raise KeyError(f"Cell-type column '{self.cell_type_col}' is not present in adata.obs.")
        if self.sample_col is not None and self.sample_col not in self.adata.obs:
            raise KeyError(f"Sample column '{self.sample_col}' is not present in adata.obs.")
        if self.config.radius is not None and self.config.n_neighbors is not None:
            raise ValueError("Set only one of radius or n_neighbors, not both.")
        if self.config.radius is None and self.config.n_neighbors is None:
            raise ValueError("Set either radius or n_neighbors.")
        if self.config.n_components < 1:
            raise ValueError("n_components must be at least 1.")
        if require_graph and "spatial_connectivities" not in self.adata.obsp:
            raise KeyError("spatial_connectivities not found. Run compute_spatial_graph first.")

    def run_pipeline(
        self,
        k: int | None = None,
        transformation: Transformation = None,
        build_graph: bool = True,
        save: bool = False,
        output_dir: str | Path = ".",
        output_prefix: str | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Run graph construction, neighborhood matrix creation, and NMF."""

        if build_graph:
            print("\n[1/4] Building spatial neighbor graph")
            self.compute_spatial_graph()
        else:
            print("\n[1/4] Using existing spatial neighbor graph")
            self._validate_inputs(require_graph=True)

        print("[2/4] Building local cell-type proportion matrix")
        self.build_proportion_matrix(transformation=transformation)

        print("[3/4] Running NMF")
        weights, signatures = self.run_nmf(k=k)

        print("[4/4] Writing metadata")
        self._write_uns_metadata()

        if save:
            self.save_pipeline_outputs(output_dir=output_dir, output_prefix=output_prefix)

        print(f"Pipeline complete: {len(self.comp_names)} spatial neighborhoods added to adata.obs.")
        return weights, signatures

    def compute_spatial_graph(
        self,
        radius: float | None = None,
        n_neighbors: int | None = None,
        coord_type: str | None = None,
        library_key: str | None = None,
    ) -> Mapping[str, object]:
        """Build a Squidpy spatial-neighbor graph in ``adata.obsp``."""

        try:
            import squidpy as sq
        except ImportError as exc:
            raise ImportError("squidpy is required for compute_spatial_graph().") from exc

        graph_radius = self.config.radius if radius is None else radius
        graph_neighbors = self.config.n_neighbors if n_neighbors is None else n_neighbors
        if graph_radius is None and graph_neighbors is None:
            raise ValueError("Provide radius or n_neighbors.")
        if graph_radius is not None and graph_neighbors is not None:
            raise ValueError("Provide only one of radius or n_neighbors.")

        params: dict[str, object] = {
            "coord_type": coord_type or self.config.coord_type,
            "library_key": library_key or self.config.library_key or self.sample_col,
        }
        if graph_neighbors is not None:
            params["n_neighs"] = graph_neighbors
        else:
            params["radius"] = graph_radius

        sq.gr.spatial_neighbors(self.adata, **params)
        return params

    def build_proportion_matrix(self, transformation: Transformation = None) -> pd.DataFrame:
        """Return a cell-by-cell-type matrix of local neighborhood percentages."""

        conn = self._connectivities()
        cat = pd.Categorical(self.adata.obs[self.cell_type_col])
        if np.any(cat.codes < 0):
            raise ValueError(f"Column '{self.cell_type_col}' contains missing values.")

        cell_types = list(cat.categories)
        n_obs = self.adata.n_obs
        indicator = sparse.csr_matrix(
            (np.ones(n_obs), (np.arange(n_obs), cat.codes)),
            shape=(n_obs, len(cell_types)),
        )
        counts = conn.dot(indicator).toarray().astype(float)
        totals = counts.sum(axis=1)

        with np.errstate(divide="ignore", invalid="ignore"):
            props = np.nan_to_num((counts / totals[:, None]) * 100.0)

        transform = self.config.transformation if transformation is None else transformation
        if transform == "sqrt":
            props = np.sqrt(props)
        elif transform == "log1p":
            props = np.log1p(props)
        elif transform is not None:
            raise ValueError("transformation must be None, 'sqrt', or 'log1p'.")

        self.prop_df = pd.DataFrame(props, index=self.adata.obs_names, columns=cell_types)
        return self.prop_df

    def find_k_elbow(
        self,
        k_range: Iterable[int] = range(2, 12),
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Plot NMF reconstruction error over candidate ranks."""

        prop_df = self._proportion_matrix()
        x_scaled = MinMaxScaler().fit_transform(prop_df.values)
        rows = []
        for k in k_range:
            model = NMF(
                n_components=k,
                init=self.config.nmf_init,
                random_state=self.config.random_state,
                max_iter=self.config.elbow_max_iter,
            )
            model.fit(x_scaled)
            rows.append({"k": k, "reconstruction_error": model.reconstruction_err_})

        result = pd.DataFrame(rows)
        plt.figure(figsize=(7, 4))
        plt.plot(result["k"], result["reconstruction_error"], marker="o", color=DEFAULT_LINE_COLOR)
        plt.xlabel("NMF components")
        plt.ylabel("Reconstruction error")
        plt.title("NMF Elbow Plot")
        plt.tight_layout()
        path = self._resolve_figure_path(save, figure_dir, output_prefix, "nmf_elbow", save_path)
        self._save_current_figure(path, "nmf_elbow")
        if show:
            plt.show()
        else:
            plt.close()
        return result

    def grid_search_nmf(
        self,
        k_range: Iterable[int] = range(3, 8),
        order_cell_types: bool = True,
        show_component_boundaries: bool = True,
        orientation: Literal["vertical", "horizontal"] = "vertical",
        cmap: str = DEFAULT_HEATMAP_CMAP,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        show: bool = True,
    ) -> dict[int, pd.DataFrame]:
        """Plot normalized NMF signatures for each candidate rank.

        By default, cell types are sorted into component-specific blocks:
        each cell type is assigned to the NMF component where it has the highest
        loading, then ordered by assignment and loading strength. The default
        orientation places cell types on rows and NMF components on columns,
        similar to compact manuscript heatmaps.
        """

        if orientation not in {"vertical", "horizontal"}:
            raise ValueError("orientation must be 'vertical' or 'horizontal'.")

        prop_df = self._proportion_matrix()
        x_scaled = MinMaxScaler().fit_transform(prop_df.values)
        signature_tables = {}
        for k in k_range:
            model = NMF(
                n_components=k,
                init=self.config.nmf_init,
                random_state=self.config.random_state,
                max_iter=self.config.elbow_max_iter,
            )
            model.fit_transform(x_scaled)
            h_table = pd.DataFrame(model.components_, columns=prop_df.columns)
            h_norm = h_table.div(h_table.max(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
            if order_cell_types:
                ordered_columns, boundaries = self._ordered_signature_columns(h_norm)
                h_plot = h_norm.loc[:, ordered_columns]
            else:
                h_plot = h_norm
                boundaries = []
            signature_tables[k] = h_plot

            component_labels = [f"N{i + 1}" for i in range(k)]
            if orientation == "vertical":
                display_table = h_plot.T
                figsize = (max(2.4, 0.42 * k), max(5.5, 0.18 * len(display_table)))
                plt.figure(figsize=figsize)
                ax = sns.heatmap(
                    display_table,
                    cmap=cmap,
                    vmin=0,
                    vmax=1,
                    xticklabels=component_labels,
                    yticklabels=True,
                    cbar=False,
                )
                ax.xaxis.tick_top()
                ax.tick_params(axis="x", rotation=90, labelsize=8, length=0)
                ax.tick_params(axis="y", rotation=0, labelsize=7, length=0)
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                ax.set_xlabel("")
                ax.set_ylabel("")
                if show_component_boundaries:
                    for boundary in boundaries:
                        ax.axhline(boundary, color="white", linewidth=0.8)
            else:
                plt.figure(figsize=(max(8, 0.3 * len(h_plot.columns)), 1.5 + k * 0.4))
                ax = sns.heatmap(
                    h_plot,
                    cmap=cmap,
                    vmin=0,
                    vmax=1,
                    yticklabels=component_labels,
                    cbar=False,
                )
                if show_component_boundaries:
                    for boundary in boundaries:
                        ax.axvline(boundary, color="white", linewidth=0.8)
            plt.title(f"Normalized Neighborhood Signatures (k={k})")
            plt.tight_layout()
            path = self._resolve_figure_path(save, figure_dir, output_prefix, f"grid_search_nmf_k{k}", None)
            self._save_current_figure(path, f"grid_search_nmf_k{k}")
            if show:
                plt.show()
            else:
                plt.close(ax.figure)
        return signature_tables

    def stability_scan_nmf(
        self,
        k_range: Iterable[int] = range(2, 11),
        n_repeats: int = 50,
        sample_fraction: float = 0.8,
        nmf_init: str = "random",
        plot: bool = True,
        save_figures: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        log_every: int = 10,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Evaluate NMF ranks by repeated cell-type co-occurrence in signatures.

        For each rank, NMF is run ``n_repeats`` times on a random subset of
        cells. By default this follows the paper-style 80% subsampling strategy
        and uses random initialization. In each run, every cell type is assigned
        to the component where it has the strongest loading in the ``H`` matrix.
        The returned pair table reports how often each pair of cell types is
        assigned to the same component.
        """

        if n_repeats < 1:
            raise ValueError("n_repeats must be at least 1.")
        if not 0 < sample_fraction <= 1:
            raise ValueError("sample_fraction must be greater than 0 and less than or equal to 1.")
        if log_every < 1:
            raise ValueError("log_every must be at least 1.")

        prop_df = self._proportion_matrix()
        x_full = prop_df.values
        cell_types = prop_df.columns.to_list()
        n_cell_types = len(cell_types)
        if n_cell_types < 2:
            raise ValueError("At least two cell types are required for co-occurrence analysis.")

        rng = np.random.default_rng(self.config.random_state)
        summary_rows = []
        pair_rows = []
        cooccurrence_by_k = {}
        cooccurrence_counts_by_k = {}
        k_values = list(k_range)
        combined_counts = np.zeros((n_cell_types, n_cell_types), dtype=float)

        print(
            "\nStarting NMF stability scan: "
            f"{len(k_values)} K values, {n_repeats} repeats per K, "
            f"sample_fraction={sample_fraction}."
        )

        for k_index, k in enumerate(k_values, start=1):
            if k < 1:
                raise ValueError("All k values must be at least 1.")
            if k > n_cell_types:
                raise ValueError(
                    f"k={k} is larger than the number of cell types ({n_cell_types}); "
                    "co-occurrence rank scans should use k <= number of cell types."
                )

            co_counts = np.zeros((n_cell_types, n_cell_types), dtype=float)
            errors = []
            print(f"[K {k_index}/{len(k_values)}] k={k}: starting {n_repeats} repeats")
            for repeat in range(n_repeats):
                if sample_fraction < 1:
                    n_rows = max(2, int(np.ceil(x_full.shape[0] * sample_fraction)))
                    row_idx = rng.choice(x_full.shape[0], size=n_rows, replace=False)
                    x_input = x_full[row_idx]
                else:
                    x_input = x_full

                x_scaled = MinMaxScaler().fit_transform(x_input)
                model = NMF(
                    n_components=k,
                    init=nmf_init,
                    random_state=self.config.random_state + repeat,
                    max_iter=self.config.elbow_max_iter,
                )
                model.fit(x_scaled)
                errors.append(model.reconstruction_err_)

                component_assignment = model.components_.argmax(axis=0)
                same_component = component_assignment[:, None] == component_assignment[None, :]
                co_counts += same_component.astype(float)
                if (repeat + 1) % log_every == 0 or repeat + 1 == n_repeats:
                    print(f"  k={k}: finished repeat {repeat + 1}/{n_repeats}")

            co_matrix = co_counts / n_repeats
            co_table = pd.DataFrame(co_matrix, index=cell_types, columns=cell_types)
            co_count_table = pd.DataFrame(co_counts, index=cell_types, columns=cell_types)
            cooccurrence_by_k[k] = co_table
            cooccurrence_counts_by_k[k] = co_count_table
            combined_counts += co_counts

            off_diag = co_matrix[~np.eye(n_cell_types, dtype=bool)]
            summary_rows.append(
                {
                    "k": k,
                    "n_repeats": n_repeats,
                    "sample_fraction": sample_fraction,
                    "mean_reconstruction_error": float(np.mean(errors)),
                    "sd_reconstruction_error": float(np.std(errors, ddof=1)) if n_repeats > 1 else 0.0,
                    "mean_pair_cooccurrence": float(np.mean(off_diag)),
                    "sd_pair_cooccurrence": float(np.std(off_diag, ddof=1)) if off_diag.size > 1 else 0.0,
                }
            )

            for i, cell_type_1 in enumerate(cell_types):
                for j in range(i + 1, n_cell_types):
                    cell_type_2 = cell_types[j]
                    pair_rows.append(
                        {
                            "k": k,
                            "cell_type_1": cell_type_1,
                            "cell_type_2": cell_type_2,
                            "pair": f"{cell_type_1} | {cell_type_2}",
                            "cooccurrence_frequency": co_matrix[i, j],
                        }
                    )

        self.k_stability_summary = pd.DataFrame(summary_rows)
        self.k_stability_pair_table = pd.DataFrame(pair_rows)
        self.k_stability_cooccurrence = cooccurrence_by_k
        self.k_stability_cooccurrence_counts = cooccurrence_counts_by_k
        self.k_stability_combined_cooccurrence = pd.DataFrame(combined_counts, index=cell_types, columns=cell_types)
        self.k_stability_combined_frequency = self.k_stability_combined_cooccurrence / (n_repeats * len(k_values))
        print("NMF stability scan complete.")
        self.k_stability_figure_paths = {}
        figure_dir = Path(figure_dir or self.config.figure_dir)
        prefix = output_prefix or self.config.output_prefix
        combined_path = figure_dir / f"{prefix}_combined_cooccurrence_heatmap.pdf" if save_figures else None
        k_scale_path = figure_dir / f"{prefix}_k_scale_cooccurrence_heatmap.pdf" if save_figures else None
        if plot:
            self.plot_combined_cooccurrence_heatmap(save=save_figures, save_path=combined_path)
            self.plot_k_stability_heatmap(self.k_stability_pair_table, save=save_figures, save_path=k_scale_path)
        elif save_figures:
            self.plot_combined_cooccurrence_heatmap(save=True, save_path=combined_path, show=False)
            self.plot_k_stability_heatmap(self.k_stability_pair_table, save=True, save_path=k_scale_path, show=False)
        if self.k_stability_figure_paths:
            print("\nSaved stability figures:")
            for label, path in self.k_stability_figure_paths.items():
                print(f"  {label}: {path}")
        return self.k_stability_summary, self.k_stability_pair_table

    def plot_k_stability_heatmap(
        self,
        pair_table: pd.DataFrame | None = None,
        min_max_frequency: float = 0.0,
        figsize: tuple[float, float] | None = None,
        cmap: str = DEFAULT_HEATMAP_CMAP,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ):
        """Plot a K-scale heatmap of cell-type-pair co-occurrence frequencies."""

        table = pair_table if pair_table is not None else self.k_stability_pair_table
        if table is None:
            raise RuntimeError("Run stability_scan_nmf() before plotting the K stability heatmap.")

        matrix = table.pivot(index="pair", columns="k", values="cooccurrence_frequency").fillna(0.0)
        if min_max_frequency > 0:
            keep = matrix.max(axis=1) >= min_max_frequency
            matrix = matrix.loc[keep]
        if matrix.empty:
            raise ValueError("No cell-type pairs remain after applying min_max_frequency.")

        if figsize is None:
            figsize = (max(6, 0.6 * matrix.shape[1]), max(5, 0.22 * matrix.shape[0]))
        plt.figure(figsize=figsize)
        ax = sns.heatmap(
            matrix,
            cmap=cmap,
            vmin=0,
            vmax=1,
            cbar_kws={"label": "Co-occurrence frequency"},
        )
        ax.set_xlabel("NMF rank k")
        ax.set_ylabel("Cell-type pair")
        ax.set_title("NMF K Stability: Cell-Type Co-occurrence Across Ranks")
        plt.tight_layout()
        path = self._resolve_figure_path(save, figure_dir, output_prefix, "k_scale_cooccurrence_heatmap", save_path)
        self._save_current_figure(path, "k_scale_cooccurrence")
        if show:
            plt.show()
        else:
            plt.close(ax.figure)
        return ax

    def plot_combined_cooccurrence_heatmap(
        self,
        matrix: pd.DataFrame | None = None,
        value: Literal["count", "frequency"] = "frequency",
        cluster: bool = True,
        figsize: tuple[float, float] | None = None,
        cmap: str = DEFAULT_HEATMAP_CMAP,
        label_fontsize: float | None = None,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ):
        """
        Plot one square co-occurrence heatmap combining all scanned K values.

        ``value="frequency"`` shows the fraction of all NMF runs where each
        cell-type pair was assigned to the same component, from 0 to 1.
        ``value="count"`` shows the corresponding absolute number of
        co-occurrences.
        """

        if matrix is None:
            if value == "count":
                matrix = self.k_stability_combined_cooccurrence
            elif value == "frequency":
                matrix = self.k_stability_combined_frequency
            else:
                raise ValueError("value must be 'count' or 'frequency'.")
        if matrix is None:
            raise RuntimeError("Run stability_scan_nmf() before plotting the combined co-occurrence heatmap.")

        plot_matrix = self._cluster_square_matrix(matrix) if cluster else matrix
        n_labels = plot_matrix.shape[0]
        if figsize is None:
            heatmap_size = max(6.5, 0.24 * n_labels)
            label_space = 2.5 if n_labels <= 35 else 3.4
            figsize = (heatmap_size + label_space + 0.9, heatmap_size)
        if label_fontsize is None:
            label_fontsize = self._heatmap_label_fontsize(n_labels)

        cbar_label = "Number of co-occurrences" if value == "count" else "Co-occurrence frequency"
        vmax = float(np.nanmax(plot_matrix.to_numpy())) if value == "count" else 1
        fig = plt.figure(figsize=figsize)
        grid = fig.add_gridspec(1, 2, width_ratios=[1, 0.045], wspace=0.55)
        ax = fig.add_subplot(grid[0, 0])
        cax = fig.add_subplot(grid[0, 1])
        ax = sns.heatmap(
            plot_matrix,
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            square=True,
            ax=ax,
            cbar_ax=cax,
            cbar_kws={"label": cbar_label},
            xticklabels=False,
            yticklabels=True,
        )
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.tick_params(axis="y", labelsize=label_fontsize, labelrotation=0, length=0, pad=4)
        ax.tick_params(axis="x", bottom=False, top=False, labelbottom=False)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title("Combined Cell-Type Co-occurrence Across NMF K Scan")
        path = self._resolve_figure_path(save, figure_dir, output_prefix, "combined_cooccurrence_heatmap", save_path)
        self._save_current_figure(path, "combined_cooccurrence")
        if show:
            plt.show()
        else:
            plt.close(ax.figure)
        return ax

    def plot_k_cooccurrence_heatmap(
        self,
        k: int,
        cooccurrence: dict[int, pd.DataFrame] | None = None,
        figsize: tuple[float, float] | None = None,
        cmap: str = DEFAULT_HEATMAP_CMAP,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ):
        """Plot the cell-type co-occurrence matrix for one NMF rank."""

        matrices = cooccurrence if cooccurrence is not None else self.k_stability_cooccurrence
        if k not in matrices:
            raise KeyError(f"No co-occurrence matrix found for k={k}. Run stability_scan_nmf() first.")
        matrix = matrices[k]
        if figsize is None:
            figsize = (max(6, 0.35 * matrix.shape[1]), max(5, 0.35 * matrix.shape[0]))
        plt.figure(figsize=figsize)
        ax = sns.heatmap(
            matrix,
            cmap=cmap,
            vmin=0,
            vmax=1,
            square=True,
            cbar_kws={"label": "Co-occurrence frequency"},
        )
        ax.set_title(f"Cell-Type Co-occurrence in NMF Signatures (k={k})")
        plt.tight_layout()
        path = self._resolve_figure_path(save, figure_dir, output_prefix, f"k{k}_cooccurrence_heatmap", save_path)
        self._save_current_figure(path, f"k_{k}_cooccurrence")
        if show:
            plt.show()
        else:
            plt.close(ax.figure)
        return ax

    def run_nmf(self, k: int | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Run final NMF and write component weights to ``adata.obs``."""

        rank = k or self.config.n_components
        if rank < 1:
            raise ValueError("k must be at least 1.")

        prop_df = self._proportion_matrix()
        x_scaled = MinMaxScaler().fit_transform(prop_df.values)
        self.nmf_model = NMF(
            n_components=rank,
            init=self.config.nmf_init,
            random_state=self.config.random_state,
            max_iter=self.config.nmf_max_iter,
        )
        self.W = self.nmf_model.fit_transform(x_scaled)
        self.H = self.nmf_model.components_
        self.comp_names = [f"neighborhood_{i + 1}_k{rank}" for i in range(rank)]

        self.weight_table = pd.DataFrame(self.W, index=prop_df.index, columns=self.comp_names)
        self.signature_table = pd.DataFrame(self.H, columns=prop_df.columns, index=self.comp_names)

        self.adata.obs[self.config.dominant_col] = self.weight_table.idxmax(axis=1).values
        for column in self.comp_names:
            self.adata.obs[column] = self.weight_table[column].values
        return self.weight_table, self.signature_table

    def summarize_samples(
        self,
        clinical_key: str = "cCR",
        agg_mode: AggregationMode = "count",
        area_col: str | None = None,
    ) -> pd.DataFrame:
        """Aggregate NMF weights to sample level."""

        self._require_nmf_results()
        if self.sample_col is None:
            raise ValueError("sample_col must be configured for sample-level summaries.")

        group_cols = [self.sample_col, clinical_key]
        if agg_mode == "count":
            summary = self.adata.obs.groupby(group_cols, observed=False)[self.comp_names].mean().reset_index()
        elif agg_mode == "area":
            resolved_area_col = self._area_col(area_col)

            def area_weighted_average(group: pd.DataFrame) -> pd.Series:
                weighted = group[self.comp_names].multiply(group[resolved_area_col], axis=0)
                return weighted.sum() / group[resolved_area_col].sum()

            summary = self.adata.obs.groupby(group_cols, observed=False).apply(area_weighted_average).reset_index()
        else:
            raise ValueError("agg_mode must be 'count' or 'area'.")

        self.sample_summary_table = summary
        return summary

    def compare_clinical_groups(
        self,
        sample_summary: pd.DataFrame | None = None,
        clinical_key: str = "cCR",
        group_order: Sequence[object] | None = None,
        save: bool = False,
        output_dir: str | Path = "results",
        output_prefix: str | None = None,
    ) -> pd.DataFrame | None:
        """Compare two clinical groups with Mann-Whitney U tests and FDR correction."""

        return self.compare_neighborhood_proportions(
            sample_summary=sample_summary,
            clinical_key=clinical_key,
            group_order=group_order,
            save=save,
            output_dir=output_dir,
            output_prefix=output_prefix,
        )

    def compare_neighborhood_proportions(
        self,
        sample_summary: pd.DataFrame | None = None,
        clinical_key: str = "cCR",
        group_order: Sequence[object] | None = None,
        save: bool = False,
        output_dir: str | Path = "results",
        output_prefix: str | None = None,
    ) -> pd.DataFrame | None:
        """Return a clear two-group p-value table for sample-level neighborhoods.

        The comparison is performed on sample-level neighborhood weights from
        ``summarize_samples``. For each neighborhood, the table reports group
        sizes, means, medians, the median difference, Mann-Whitney U statistic,
        raw p-value, and Benjamini-Hochberg FDR-adjusted p-value.
        """

        self._require_nmf_results()
        summary = sample_summary if sample_summary is not None else self.sample_summary_table
        if summary is None:
            summary = self.summarize_samples(clinical_key=clinical_key)

        if clinical_key not in summary:
            raise KeyError(f"Clinical/group column '{clinical_key}' is not present in sample_summary.")

        groups = list(pd.Series(summary[clinical_key].dropna().unique()))
        if group_order is not None:
            missing = [group for group in group_order if group not in groups]
            if missing:
                raise ValueError(f"Groups from group_order are not present in sample_summary: {missing}")
            groups = list(group_order)

        if len(groups) != 2:
            return None

        group_1_label, group_2_label = groups[0], groups[1]
        rows = []
        for component in self.comp_names:
            group_1 = summary.loc[summary[clinical_key] == group_1_label, component].dropna().astype(float).values
            group_2 = summary.loc[summary[clinical_key] == group_2_label, component].dropna().astype(float).values
            if len(group_1) == 0 or len(group_2) == 0:
                u_stat, p_value = np.nan, np.nan
            else:
                u_stat, p_value = stats.mannwhitneyu(group_1, group_2, alternative="two-sided")

            group_1_mean = float(np.mean(group_1)) if len(group_1) else np.nan
            group_2_mean = float(np.mean(group_2)) if len(group_2) else np.nan
            group_1_median = float(np.median(group_1)) if len(group_1) else np.nan
            group_2_median = float(np.median(group_2)) if len(group_2) else np.nan
            rows.append(
                {
                    "neighborhood": component,
                    "clinical_key": clinical_key,
                    "group_1": group_1_label,
                    "group_2": group_2_label,
                    "n_group_1": len(group_1),
                    "n_group_2": len(group_2),
                    "mean_group_1": group_1_mean,
                    "mean_group_2": group_2_mean,
                    "median_group_1": group_1_median,
                    "median_group_2": group_2_median,
                    "mean_diff_group2_minus_group1": group_2_mean - group_1_mean,
                    "median_diff_group2_minus_group1": group_2_median - group_1_median,
                    "mannwhitney_u": float(u_stat) if np.isfinite(u_stat) else np.nan,
                    "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                }
            )

        result = pd.DataFrame(rows)
        valid_p = result["p_value"].notna()
        result["p_adj_fdr_bh"] = np.nan
        if valid_p.any():
            _, p_adj, _, _ = multipletests(result.loc[valid_p, "p_value"], method="fdr_bh")
            result.loc[valid_p, "p_adj_fdr_bh"] = p_adj
        result = result.sort_values(["p_adj_fdr_bh", "p_value"], na_position="last").reset_index(drop=True)

        if save:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            prefix = output_prefix or self.config.output_prefix
            path = output_dir / f"{prefix}_neighborhood_proportion_pvalues.csv"
            result.to_csv(path, index=False)
            print(f"Saved p-value table: {path}")

        return result

    def print_neighborhood_proportion_pvalues(
        self,
        sample_summary: pd.DataFrame | None = None,
        clinical_key: str = "cCR",
        group_order: Sequence[object] | None = None,
        save: bool = False,
        output_dir: str | Path = "results",
        output_prefix: str | None = None,
    ) -> pd.DataFrame | None:
        """Print and return the neighborhood proportion p-value table."""

        result = self.compare_neighborhood_proportions(
            sample_summary=sample_summary,
            clinical_key=clinical_key,
            group_order=group_order,
            save=save,
            output_dir=output_dir,
            output_prefix=output_prefix,
        )
        if result is None:
            print(f"Need exactly two groups in '{clinical_key}' to compute p-values.")
            return None
        print("\n--- Neighborhood Proportion Comparison ---")
        print(result.to_string(index=False))
        return result

    def plot_neighbor_counts(
        self,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
        sort: bool = True,
    ):
        """Plot spatial-neighbor counts by cell type.

        By default, cell types are ordered from lowest to highest median
        neighbor count to make density differences easier to scan.
        """

        conn = self._connectivities()
        neighbor_counts = np.diff(conn.indptr)
        diagonal = conn.diagonal()
        if np.any(diagonal != 0):
            neighbor_counts = neighbor_counts - (diagonal != 0).astype(int)

        plot_df = pd.DataFrame(
            {"cell_type": self.adata.obs[self.cell_type_col].values, "n_neighbors": neighbor_counts},
            index=self.adata.obs_names,
        )
        sampled, valid_types = self._sample_by_cell_type(plot_df, "cell_type")
        if sort:
            valid_types = self._order_categories_by_median(plot_df, "cell_type", "n_neighbors", valid_types)

        plt.figure(figsize=(max(8, 0.3 * len(valid_types)), 5))
        ax = sns.violinplot(
            data=sampled,
            x="cell_type",
            y="n_neighbors",
            order=valid_types,
            cut=0,
            inner=None,
            color=DEFAULT_POINT_COLOR,
        )
        ax.tick_params(axis="x", rotation=90)
        ax.set_title("Neighbor Density (Counts)")
        plt.tight_layout()
        path = self._resolve_figure_path(save, figure_dir, output_prefix, "neighbor_counts", save_path)
        self._save_current_figure(path, "neighbor_counts")
        if show:
            plt.show()
        else:
            plt.close(ax.figure)
        return ax

    def plot_neighbor_distances(
        self,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
        sort: bool = True,
    ):
        """Plot mean spatial-neighbor distance by cell type.

        By default, cell types are ordered from lowest to highest median mean
        distance to neighbors.
        """

        dist_matrix = self.adata.obsp.get("spatial_distances")
        if dist_matrix is None:
            raise KeyError("spatial_distances not found. Run compute_spatial_graph first.")
        dist_matrix = sparse.csr_matrix(dist_matrix)

        nnz_per_row = np.diff(dist_matrix.indptr)
        row_sums = np.asarray(dist_matrix.sum(axis=1)).ravel()
        with np.errstate(divide="ignore", invalid="ignore"):
            mean_distance = np.nan_to_num(row_sums / nnz_per_row)

        plot_df = pd.DataFrame(
            {"cell_type": self.adata.obs[self.cell_type_col].values, "avg_distance": mean_distance},
            index=self.adata.obs_names,
        )
        sampled, valid_types = self._sample_by_cell_type(plot_df, "cell_type")
        if sort:
            valid_types = self._order_categories_by_median(plot_df, "cell_type", "avg_distance", valid_types)

        plt.figure(figsize=(max(8, 0.3 * len(valid_types)), 5))
        ax = sns.violinplot(
            data=sampled,
            x="cell_type",
            y="avg_distance",
            order=valid_types,
            cut=0,
            inner=None,
            color=DEFAULT_POINT_COLOR,
        )
        ax.tick_params(axis="x", rotation=90)
        ax.set_ylabel("Mean distance to neighbors")
        ax.set_title("Spatial Dispersion by Cell Type")
        plt.tight_layout()
        path = self._resolve_figure_path(save, figure_dir, output_prefix, "neighbor_distances", save_path)
        self._save_current_figure(path, "neighbor_distances")
        if show:
            plt.show()
        else:
            plt.close(ax.figure)
        return ax

    def plot_signatures(
        self,
        normalize: bool = False,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ):
        """Plot final NMF cell-type signatures."""

        self._require_nmf_results()
        signature_table = self.signature_table.copy()
        if normalize:
            signature_table = signature_table.div(signature_table.max(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

        plt.figure(figsize=(max(8, 0.3 * len(signature_table.columns)), 2 + len(signature_table) * 0.5))
        ax = sns.heatmap(
            signature_table,
            cmap=DEFAULT_HEATMAP_CMAP,
            vmin=0,
            vmax=float(np.nanmax(signature_table.to_numpy())) if signature_table.size else None,
            yticklabels=self.comp_names,
        )
        ax.set_title("Final NMF Neighborhood Signatures")
        plt.tight_layout()
        suffix = "nmf_signatures_normalized" if normalize else "nmf_signatures"
        path = self._resolve_figure_path(save, figure_dir, output_prefix, suffix, save_path)
        self._save_current_figure(path, suffix)
        if show:
            plt.show()
        else:
            plt.close(ax.figure)
        return ax

    def plot_neighborhood_pie_charts(
        self,
        top_n: int | None = 8,
        min_fraction: float = 0.02,
        include_other: bool = True,
        ncols: int = 3,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Plot cell-type contribution pies for each final NMF neighborhood.

        Contributions are calculated from the final NMF ``H`` signature matrix
        by normalizing each neighborhood signature to sum to 1. Small
        contributors can be grouped into ``Other`` for readability.
        """

        self._require_nmf_results()
        if ncols < 1:
            raise ValueError("ncols must be at least 1.")
        if top_n is not None and top_n < 1:
            raise ValueError("top_n must be at least 1 or None.")
        if min_fraction < 0:
            raise ValueError("min_fraction must be non-negative.")

        proportions = self.signature_table.div(self.signature_table.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
        rows = []
        plot_rows = []
        for neighborhood, values in proportions.iterrows():
            ordered = values.sort_values(ascending=False)
            for cell_type, fraction in ordered.items():
                rows.append(
                    {
                        "neighborhood": neighborhood,
                        "cell_type": cell_type,
                        "fraction": float(fraction),
                        "percent": float(fraction * 100),
                    }
                )

            keep = ordered[ordered >= min_fraction]
            if top_n is not None:
                keep = keep.iloc[:top_n]
            remainder = float(ordered.drop(index=keep.index).sum())
            for cell_type, fraction in keep.items():
                plot_rows.append(
                    {
                        "neighborhood": neighborhood,
                        "cell_type": cell_type,
                        "plot_cell_type": cell_type,
                        "fraction": float(fraction),
                    }
                )
            if include_other and remainder > 0:
                plot_rows.append(
                    {
                        "neighborhood": neighborhood,
                        "cell_type": "Other",
                        "plot_cell_type": "Other",
                        "fraction": remainder,
                    }
                )

        contribution_table = pd.DataFrame(rows)
        plot_table = pd.DataFrame(plot_rows)
        plot_cell_types = plot_table["plot_cell_type"].drop_duplicates().tolist()
        palette = self._categorical_palette(plot_cell_types)

        n_neighborhoods = len(self.comp_names)
        ncols = min(ncols, n_neighborhoods)
        nrows = int(np.ceil(n_neighborhoods / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.5 * nrows), squeeze=False)
        axes_flat = axes.ravel()

        for ax, neighborhood in zip(axes_flat, self.comp_names):
            subset = plot_table.loc[plot_table["neighborhood"] == neighborhood]
            colors = [palette[cell_type] for cell_type in subset["plot_cell_type"]]
            ax.pie(
                subset["fraction"],
                colors=colors,
                startangle=90,
                counterclock=False,
                wedgeprops={"linewidth": 0.6, "edgecolor": "white"},
            )
            ax.set_title(neighborhood, fontsize=10)
            ax.set_aspect("equal")

        for ax in axes_flat[n_neighborhoods:]:
            ax.axis("off")

        handles = [
            plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[cell_type], markersize=8)
            for cell_type in plot_cell_types
        ]
        fig.legend(
            handles,
            plot_cell_types,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            title="Cell type",
            fontsize=8,
            title_fontsize=9,
        )
        fig.suptitle("Cell-Type Composition of NMF Neighborhood Signatures", y=0.995, fontsize=12)
        fig.tight_layout(rect=(0, 0, 0.86, 0.96))

        path = self._resolve_figure_path(save, figure_dir, output_prefix, "neighborhood_pie_charts", save_path)
        self._save_current_figure(path, "neighborhood_pie_charts")
        if show:
            plt.show()
        else:
            plt.close(fig)
        return contribution_table

    def plot_sample_level_boxplots(
        self,
        clinical_key: str = "cCR",
        agg_mode: AggregationMode = "count",
        area_col: str | None = None,
        save: bool = True,
        figure_dir: str | Path | None = None,
        output_prefix: str | None = None,
        save_path: str | Path | None = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Plot sample-level neighborhood weights by clinical group."""

        summary = self.summarize_samples(clinical_key=clinical_key, agg_mode=agg_mode, area_col=area_col)
        ylabel = "Mean NMF Weight (by Cell Count)"
        if agg_mode == "area":
            ylabel = f"Mean NMF Weight (Weighted by {self._area_col(area_col)})"

        melt_df = summary.melt(
            id_vars=[clinical_key],
            value_vars=self.comp_names,
            var_name="component",
            value_name="weight",
        )

        plt.figure(figsize=(len(self.comp_names) * 2.5, 6))
        ax = sns.boxplot(
            data=melt_df,
            x="component",
            y="weight",
            hue=clinical_key,
            showfliers=False,
            palette=DEFAULT_GROUP_PALETTE,
            boxprops={"alpha": 0.4},
        )
        sns.stripplot(
            data=melt_df,
            x="component",
            y="weight",
            hue=clinical_key,
            dodge=True,
            alpha=0.8,
            size=6,
            palette=DEFAULT_GROUP_PALETTE,
            legend=False,
            linewidth=1,
            ax=ax,
        )
        ax.tick_params(axis="x", rotation=45)
        ax.set_ylabel(ylabel)
        ax.set_title(f"Patient-Level Neighborhood Comparison ({agg_mode})")
        plt.tight_layout()
        suffix = f"sample_level_boxplots_{agg_mode}"
        path = self._resolve_figure_path(save, figure_dir, output_prefix, suffix, save_path)
        self._save_current_figure(path, suffix)
        if show:
            plt.show()
        else:
            plt.close(ax.figure)

        stats_df = self.compare_neighborhood_proportions(summary, clinical_key=clinical_key)
        if stats_df is not None:
            print("\n--- Statistical Summary ---")
            print(stats_df.to_string(index=False))
        return summary

    def export_proportion_matrix(self, path: str | Path) -> None:
        self._proportion_matrix().to_csv(path)

    def export_nmf_weights(self, path: str | Path) -> None:
        self._require_nmf_results()
        self.weight_table.to_csv(path)

    def export_nmf_signatures(self, path: str | Path) -> None:
        self._require_nmf_results()
        self.signature_table.to_csv(path)

    def cell_assignment_table(
        self,
        metadata_cols: Sequence[str] | None = None,
        include_weights: bool = True,
        dominant_col: str | None = None,
        cell_id_col: str | None = None,
        include_obs_metadata: bool = True,
    ) -> pd.DataFrame:
        """Return per-cell niche assignments, metadata, and optional NMF weights.

        The table always contains ``cell_id`` and the dominant niche column. If
        ``adata.obs`` contains the configured cell ID column, that column is used
        for ``cell_id``; otherwise ``adata.obs_names`` are used. By default, all
        ``adata.obs`` metadata columns are included so the exported CSV can be
        used directly for downstream sample/group boxplots.
        """

        self._require_nmf_results()
        dominant_col = dominant_col or self.config.dominant_col
        if dominant_col not in self.adata.obs:
            raise KeyError(f"Dominant niche column '{dominant_col}' is not present in adata.obs.")

        resolved_cell_id_col = self.config.cell_id_col if cell_id_col is None else cell_id_col
        if resolved_cell_id_col is not None and resolved_cell_id_col in self.adata.obs:
            cell_ids = self.adata.obs[resolved_cell_id_col].to_numpy()
        else:
            cell_ids = self.adata.obs_names.to_numpy()

        cols = []
        if include_obs_metadata:
            default_metadata = list(self.adata.obs.columns)
        else:
            default_metadata = [self.cell_type_col]
            if self.sample_col is not None:
                default_metadata.append(self.sample_col)
        if metadata_cols is not None:
            default_metadata.extend(metadata_cols)

        reserved_cols = {dominant_col, *self.comp_names}
        if resolved_cell_id_col is not None:
            reserved_cols.add(resolved_cell_id_col)
        for col in default_metadata:
            if col is not None and col in self.adata.obs and col not in cols and col not in reserved_cols:
                cols.append(col)
        if dominant_col not in cols:
            cols.append(dominant_col)
        if include_weights:
            cols.extend([col for col in self.comp_names if col not in cols])

        table = self.adata.obs.loc[:, cols].copy()
        table.insert(0, "cell_id", cell_ids)
        return table

    def export_cell_assignments(
        self,
        path: str | Path,
        metadata_cols: Sequence[str] | None = None,
        include_weights: bool = True,
        dominant_col: str | None = None,
        cell_id_col: str | None = None,
        include_obs_metadata: bool = True,
    ) -> pd.DataFrame:
        """Export per-cell niche assignments to CSV and return the table."""

        table = self.cell_assignment_table(
            metadata_cols=metadata_cols,
            include_weights=include_weights,
            dominant_col=dominant_col,
            cell_id_col=cell_id_col,
            include_obs_metadata=include_obs_metadata,
        )
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path, index=False)
        return table

    def save_pipeline_outputs(
        self,
        output_dir: str | Path = ".",
        output_prefix: str | None = None,
    ) -> dict[str, Path]:
        """Save the neighborhood matrix, NMF outputs, and cell assignments."""

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_prefix = output_prefix or self.config.output_prefix

        paths = {
            "proportion_matrix": output_dir / f"{output_prefix}_proportion_matrix.csv",
            "weights": output_dir / f"{output_prefix}_nmf_weights.csv",
            "signatures": output_dir / f"{output_prefix}_nmf_signatures.csv",
            "cell_assignments": output_dir / f"{output_prefix}_cell_niche_assignments.csv",
        }
        self.export_proportion_matrix(paths["proportion_matrix"])
        self.export_nmf_weights(paths["weights"])
        self.export_nmf_signatures(paths["signatures"])
        self.export_cell_assignments(paths["cell_assignments"])

        print("\nSaved pipeline outputs:")
        for label, path in paths.items():
            print(f"  {label}: {path}")
        return paths

    def _connectivities(self) -> sparse.csr_matrix:
        self._validate_inputs(require_graph=True)
        return sparse.csr_matrix(self.adata.obsp["spatial_connectivities"])

    def _proportion_matrix(self) -> pd.DataFrame:
        if self.prop_df is None:
            return self.build_proportion_matrix()
        return self.prop_df

    def _require_nmf_results(self) -> None:
        if self.weight_table is None or self.signature_table is None or not self.comp_names:
            raise RuntimeError("Run run_nmf() before using this method.")

    def _area_col(self, area_col: str | None) -> str:
        requested = area_col or self.config.area_col
        if requested in self.adata.obs:
            return requested
        if self.config.alternative_area_col in self.adata.obs:
            return self.config.alternative_area_col
        raise KeyError(f"Area column '{requested}' is not present in adata.obs.")

    def _sample_by_cell_type(self, df: pd.DataFrame, cell_type_col: str) -> tuple[pd.DataFrame, Sequence[str]]:
        valid_types = (
            df[cell_type_col]
            .value_counts()[lambda counts: counts >= self.config.min_cells_for_plot]
            .index.tolist()
        )
        plot_df = df[df[cell_type_col].isin(valid_types)]

        rng = np.random.default_rng(self.config.random_state)
        sampled_indices = []
        for cell_type in valid_types:
            idx = plot_df.index[plot_df[cell_type_col] == cell_type].to_numpy()
            sampled_indices.extend(
                rng.choice(idx, size=min(len(idx), self.config.max_plot_cells_per_type), replace=False)
            )
        return plot_df.loc[sampled_indices], valid_types

    @staticmethod
    def _order_categories_by_median(
        df: pd.DataFrame,
        category_col: str,
        value_col: str,
        valid_categories: Sequence[str],
    ) -> list[str]:
        """Order categories from lowest to highest median plotted value."""

        order_table = (
            df[df[category_col].isin(valid_categories)]
            .groupby(category_col, observed=True)[value_col]
            .median()
            .sort_values(kind="mergesort")
        )
        return order_table.index.tolist()

    def _cluster_square_matrix(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """Order a symmetric co-occurrence matrix to reveal diagonal blocks."""

        if matrix.shape[0] <= 2:
            return matrix
        values = matrix.to_numpy(dtype=float)
        max_value = np.nanmax(values)
        if not np.isfinite(max_value) or max_value <= 0:
            return matrix
        similarity = values / max_value
        similarity = np.clip(similarity, 0, 1)
        distance = 1 - similarity
        np.fill_diagonal(distance, 0)
        try:
            order = leaves_list(linkage(squareform(distance, checks=False), method="average"))
        except ValueError:
            return matrix
        return matrix.iloc[order, order]

    @staticmethod
    def _ordered_signature_columns(signature_table: pd.DataFrame) -> tuple[list[str], list[int]]:
        """Order cell types by strongest NMF component and loading strength."""

        assignments = signature_table.idxmax(axis=0)
        component_rank = assignments.map({label: i for i, label in enumerate(signature_table.index)})
        max_loading = signature_table.max(axis=0)
        order_df = pd.DataFrame(
            {
                "cell_type": signature_table.columns,
                "component_rank": component_rank.to_numpy(),
                "max_loading": max_loading.to_numpy(),
            }
        )
        order_df = order_df.sort_values(
            ["component_rank", "max_loading", "cell_type"],
            ascending=[True, False, True],
            kind="mergesort",
        )
        ordered_columns = order_df["cell_type"].tolist()
        block_sizes = order_df.groupby("component_rank", sort=False).size().to_numpy()
        boundaries = np.cumsum(block_sizes)[:-1].tolist()
        return ordered_columns, boundaries

    @staticmethod
    def _heatmap_label_fontsize(n_labels: int) -> float:
        if n_labels <= 25:
            return 8
        if n_labels <= 40:
            return 6.5
        if n_labels <= 60:
            return 5.5
        return 4.5

    @staticmethod
    def _categorical_palette(labels: Sequence[object]) -> dict[object, tuple[float, float, float]]:
        colors = sns.color_palette("tab20", n_colors=max(1, len(labels)))
        palette = {label: colors[i] for i, label in enumerate(labels)}
        if "Other" in palette:
            palette["Other"] = (0.72, 0.72, 0.72)
        return palette

    def _save_current_figure(self, save_path: str | Path | None, label: str) -> None:
        if save_path is None:
            return
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        plt.gcf().savefig(path, bbox_inches="tight")
        self.figure_paths[label] = path
        self.k_stability_figure_paths[label] = path

    def _resolve_figure_path(
        self,
        save: bool,
        figure_dir: str | Path | None,
        output_prefix: str | None,
        suffix: str,
        save_path: str | Path | None,
    ) -> Path | None:
        if save_path is not None:
            return Path(save_path)
        if not save:
            return None
        prefix = output_prefix or self.config.output_prefix
        resolved_figure_dir = figure_dir or self.config.figure_dir
        return Path(resolved_figure_dir) / f"{prefix}_{suffix}.pdf"

    def _write_uns_metadata(self) -> None:
        self.adata.uns["spatial_neighborhood_analysis"] = {
            "cell_type_col": self.cell_type_col,
            "cell_id_col": self.config.cell_id_col,
            "sample_col": self.sample_col,
            "radius": self.config.radius,
            "n_neighbors": self.config.n_neighbors,
            "coord_type": self.config.coord_type,
            "library_key": self.config.library_key or self.sample_col,
            "transformation": self.config.transformation,
            "n_components": self.config.n_components,
            "dominant_col": self.config.dominant_col,
            "component_columns": list(self.comp_names),
        }

    @staticmethod
    def default_output_columns() -> Mapping[str, str]:
        """Describe columns added to ``adata.obs`` by the analyzer."""

        return {
            "nmf_dominant": "Dominant NMF neighborhood for each cell.",
            "neighborhood_*": "Continuous per-cell NMF neighborhood weights.",
        }


SpatialNeighborhoodCalculator = SpatialNeighborhoodAnalyzer
