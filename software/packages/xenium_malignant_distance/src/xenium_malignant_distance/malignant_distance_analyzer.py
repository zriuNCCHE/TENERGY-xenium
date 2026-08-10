from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


@dataclass
class MalignantDistanceConfig:
    """Default parameters for Xenium malignant-cell distance analysis."""

    cell_type_col: str = "final_cell_type"
    malignant_values: tuple[str, ...] = ("malignant",)
    immune_stromal_values: tuple[str, ...] = (
        "immune",
        "stromal",
        "T cell",
        "B cell",
        "myeloid",
        "fibroblast",
        "endothelial",
    )
    spatial_key: str = "spatial"
    sample_col: str | None = None
    coordinate_cols: tuple[str, str] | None = None
    neighbor_radius: float = 50.0
    close_distance: float = 30.0
    intermediate_distance: float | None = 100.0
    core_malignant_fraction: float = 0.75
    min_neighbors_for_region: int = 5
    include_self_in_malignant_neighbors: bool = False
    output_prefix: str = "malignant_distance"


class MalignantDistanceAnalyzer:
    """
    Annotate malignant tumor margin/core regions and non-malignant proximity.

    The analyzer expects an AnnData object with cell-type labels in ``adata.obs``
    and spatial coordinates either in ``adata.obsm[spatial_key]`` or in two
    explicit observation columns. Results are written back to ``adata.obs``.
    """

    def __init__(
        self,
        adata,
        cell_type_col: str | None = None,
        spatial_key: str | None = None,
        sample_col: str | None = None,
        config: MalignantDistanceConfig | None = None,
    ):
        self.config = config or MalignantDistanceConfig()
        if cell_type_col is not None:
            self.config.cell_type_col = cell_type_col
        if spatial_key is not None:
            self.config.spatial_key = spatial_key
        if sample_col is not None:
            self.config.sample_col = sample_col

        self.adata = adata
        self.cell_type_col = self.config.cell_type_col
        self.spatial_key = self.config.spatial_key
        self.sample_col = self.config.sample_col

        self.result_table: pd.DataFrame | None = None
        self.summary_table: pd.DataFrame | None = None

        self._validate_inputs()

    def _validate_inputs(self) -> None:
        if self.cell_type_col not in self.adata.obs:
            raise KeyError(f"Cell-type column '{self.cell_type_col}' is not present in adata.obs.")
        if self.sample_col is not None and self.sample_col not in self.adata.obs:
            raise KeyError(f"Sample column '{self.sample_col}' is not present in adata.obs.")
        if self.config.coordinate_cols is not None:
            missing = [col for col in self.config.coordinate_cols if col not in self.adata.obs]
            if missing:
                raise KeyError(f"Coordinate columns are missing from adata.obs: {missing}")
        elif self.spatial_key not in self.adata.obsm:
            raise KeyError(f"Spatial key '{self.spatial_key}' is not present in adata.obsm.")

        coords = self._coordinates()
        if coords.shape[1] < 2:
            raise ValueError("Spatial coordinates must have at least two columns.")
        if not np.isfinite(coords[:, :2]).all():
            raise ValueError("Spatial coordinates contain NaN or infinite values.")
        if self.config.neighbor_radius <= 0:
            raise ValueError("neighbor_radius must be positive.")
        if self.config.close_distance < 0:
            raise ValueError("close_distance must be non-negative.")
        if self.config.intermediate_distance is not None and (
            self.config.intermediate_distance < self.config.close_distance
        ):
            raise ValueError("intermediate_distance must be greater than or equal to close_distance.")
        if not 0 <= self.config.core_malignant_fraction <= 1:
            raise ValueError("core_malignant_fraction must be between 0 and 1.")

        malignant_mask = self._malignant_mask()
        if malignant_mask.sum() == 0:
            raise ValueError("No malignant cells were found with the configured malignant_values.")

    def run_pipeline(
        self,
        neighbor_radius: float | None = None,
        close_distance: float | None = None,
        intermediate_distance: float | None = None,
        core_malignant_fraction: float | None = None,
        save: bool = False,
        output_dir: str | Path = ".",
        output_prefix: str | None = None,
    ) -> pd.DataFrame:
        """Run distance calculation, neighborhood scoring, and region annotation."""

        if neighbor_radius is not None:
            self.config.neighbor_radius = neighbor_radius
        if close_distance is not None:
            self.config.close_distance = close_distance
        if intermediate_distance is not None:
            self.config.intermediate_distance = intermediate_distance
        if core_malignant_fraction is not None:
            self.config.core_malignant_fraction = core_malignant_fraction
        self._validate_inputs()

        print("\n[1/4] Computing nearest malignant-cell distances")
        result = self._compute_nearest_malignant_distances()

        print("[2/4] Computing local malignant-neighbor fractions")
        neighborhood = self._compute_neighborhood_composition()
        result = result.join(neighborhood)

        print("[3/4] Classifying malignant regions and non-malignant proximity")
        result["malignant_region"] = self._classify_malignant_regions(result)
        result["malignant_proximity"] = self._classify_non_malignant_proximity(result)
        result["distance_region"] = self._combine_region_labels(result)

        print("[4/4] Writing annotations to adata.obs")
        self.result_table = result
        self._write_obs_columns(result)
        self.summary_table = self.summarize()
        self._write_uns_metadata()

        if save:
            self.save_pipeline_outputs(output_dir=output_dir, output_prefix=output_prefix)

        print("Pipeline complete: malignant distance annotations added to adata.obs.")
        return result

    def _coordinates(self) -> np.ndarray:
        if self.config.coordinate_cols is not None:
            coords = self.adata.obs.loc[:, list(self.config.coordinate_cols)].to_numpy()
        else:
            coords = self.adata.obsm[self.spatial_key]
        if sparse.issparse(coords):
            coords = coords.toarray()
        return np.asarray(coords, dtype=np.float64)[:, :2]

    def _malignant_mask(self) -> np.ndarray:
        values = self.adata.obs[self.cell_type_col].astype(str)
        return values.isin(self.config.malignant_values).to_numpy()

    def _immune_stromal_mask(self) -> np.ndarray:
        values = self.adata.obs[self.cell_type_col].astype(str)
        return values.isin(self.config.immune_stromal_values).to_numpy()

    def _group_indices(self) -> list[np.ndarray]:
        if self.sample_col is None:
            return [np.arange(self.adata.n_obs)]
        groups = []
        samples = self.adata.obs[self.sample_col]
        for sample_id in samples.drop_duplicates():
            groups.append(np.flatnonzero(samples.to_numpy() == sample_id))
        return groups

    def _compute_nearest_malignant_distances(self) -> pd.DataFrame:
        coords = self._coordinates()
        malignant_mask = self._malignant_mask()
        nearest_distance = np.full(self.adata.n_obs, np.nan, dtype=np.float64)
        nearest_index = np.full(self.adata.n_obs, -1, dtype=np.int64)

        for group_idx in self._group_indices():
            group_malignant_idx = group_idx[malignant_mask[group_idx]]
            if group_malignant_idx.size == 0:
                continue
            model = NearestNeighbors(n_neighbors=1)
            model.fit(coords[group_malignant_idx])
            distances, local_indices = model.kneighbors(coords[group_idx])
            nearest_distance[group_idx] = distances[:, 0]
            nearest_index[group_idx] = group_malignant_idx[local_indices[:, 0]]

        nearest_cell = pd.Series(pd.NA, index=self.adata.obs_names, dtype="object")
        valid = nearest_index >= 0
        nearest_cell.iloc[np.flatnonzero(valid)] = self.adata.obs_names.to_numpy()[nearest_index[valid]]

        return pd.DataFrame(
            {
                "distance_to_malignant": nearest_distance,
                "nearest_malignant_cell": nearest_cell,
            },
            index=self.adata.obs_names,
        )

    def _compute_neighborhood_composition(self) -> pd.DataFrame:
        coords = self._coordinates()
        malignant_mask = self._malignant_mask()
        radius = self.config.neighbor_radius

        neighbor_count = np.zeros(self.adata.n_obs, dtype=np.int64)
        malignant_neighbor_count = np.zeros(self.adata.n_obs, dtype=np.int64)

        for group_idx in self._group_indices():
            if group_idx.size == 0:
                continue
            model = NearestNeighbors(radius=radius)
            model.fit(coords[group_idx])
            local_neighbors = model.radius_neighbors(coords[group_idx], return_distance=False)
            for local_pos, local_neighbor_idx in enumerate(local_neighbors):
                global_cell = group_idx[local_pos]
                global_neighbors = group_idx[local_neighbor_idx]
                if not self.config.include_self_in_malignant_neighbors:
                    global_neighbors = global_neighbors[global_neighbors != global_cell]
                neighbor_count[global_cell] = global_neighbors.size
                malignant_neighbor_count[global_cell] = int(malignant_mask[global_neighbors].sum())

        malignant_fraction = np.full(self.adata.n_obs, np.nan, dtype=np.float64)
        valid = neighbor_count > 0
        malignant_fraction[valid] = malignant_neighbor_count[valid] / neighbor_count[valid]

        return pd.DataFrame(
            {
                "neighbor_count": neighbor_count,
                "malignant_neighbor_count": malignant_neighbor_count,
                "malignant_neighbor_fraction": malignant_fraction,
            },
            index=self.adata.obs_names,
        )

    def _classify_malignant_regions(self, result: pd.DataFrame) -> pd.Series:
        malignant_mask = self._malignant_mask()
        labels = pd.Series(pd.NA, index=self.adata.obs_names, dtype="object")
        enough_neighbors = result["neighbor_count"] >= self.config.min_neighbors_for_region
        core_like = result["malignant_neighbor_fraction"] >= self.config.core_malignant_fraction

        labels.loc[malignant_mask & ~enough_neighbors.to_numpy()] = "isolated_malignant"
        labels.loc[malignant_mask & enough_neighbors.to_numpy() & core_like.to_numpy()] = "tumor_inner"
        labels.loc[malignant_mask & enough_neighbors.to_numpy() & ~core_like.to_numpy()] = "tumor_margin"
        return labels

    def _classify_non_malignant_proximity(self, result: pd.DataFrame) -> pd.Series:
        malignant_mask = self._malignant_mask()
        immune_stromal_mask = self._immune_stromal_mask()
        labels = pd.Series(pd.NA, index=self.adata.obs_names, dtype="object")
        candidate_mask = ~malignant_mask

        labels.loc[candidate_mask] = "other_non_malignant"
        labels.loc[candidate_mask & immune_stromal_mask] = "distant_to_malignant"

        distance = result["distance_to_malignant"]
        no_reference = candidate_mask & immune_stromal_mask & distance.isna().to_numpy()
        labels.loc[no_reference] = "no_malignant_reference"

        close = candidate_mask & immune_stromal_mask & (distance <= self.config.close_distance).to_numpy()
        labels.loc[close] = "close_to_malignant"

        if self.config.intermediate_distance is not None:
            intermediate = (
                candidate_mask
                & immune_stromal_mask
                & (distance > self.config.close_distance).to_numpy()
                & (distance <= self.config.intermediate_distance).to_numpy()
            )
            labels.loc[intermediate] = "intermediate_to_malignant"

        return labels

    def _combine_region_labels(self, result: pd.DataFrame) -> pd.Series:
        labels = pd.Series(pd.NA, index=self.adata.obs_names, dtype="object")
        malignant_mask = self._malignant_mask()
        labels.loc[malignant_mask] = result.loc[malignant_mask, "malignant_region"]
        labels.loc[~malignant_mask] = result.loc[~malignant_mask, "malignant_proximity"]
        return labels

    def _write_obs_columns(self, result: pd.DataFrame) -> None:
        for column in result.columns:
            self.adata.obs[column] = result[column]

    def _write_uns_metadata(self) -> None:
        self.adata.uns["malignant_distance_analysis"] = {
            "cell_type_col": self.cell_type_col,
            "malignant_values": list(self.config.malignant_values),
            "immune_stromal_values": list(self.config.immune_stromal_values),
            "spatial_key": self.spatial_key,
            "sample_col": self.sample_col,
            "coordinate_cols": list(self.config.coordinate_cols) if self.config.coordinate_cols else None,
            "neighbor_radius": self.config.neighbor_radius,
            "close_distance": self.config.close_distance,
            "intermediate_distance": self.config.intermediate_distance,
            "core_malignant_fraction": self.config.core_malignant_fraction,
            "min_neighbors_for_region": self.config.min_neighbors_for_region,
            "include_self_in_malignant_neighbors": self.config.include_self_in_malignant_neighbors,
        }

    def summarize(self) -> pd.DataFrame:
        """Return counts by combined distance-region labels."""

        if self.result_table is None:
            raise RuntimeError("Run the pipeline before summarizing results.")

        group_cols = ["distance_region"]
        if self.sample_col is not None:
            summary_source = self.adata.obs[[self.sample_col]].join(self.result_table[["distance_region"]])
            group_cols = [self.sample_col, "distance_region"]
        else:
            summary_source = self.result_table[["distance_region"]]

        summary = (
            summary_source.groupby(group_cols, dropna=False)
            .size()
            .rename("n_cells")
            .reset_index()
            .sort_values(group_cols)
            .reset_index(drop=True)
        )
        return summary

    def export_annotations(self, path: str | Path) -> None:
        if self.result_table is None:
            raise RuntimeError("Run the pipeline before exporting annotations.")
        export_df = self.result_table.copy()
        export_df.insert(0, "cell_id", export_df.index)
        export_df.to_csv(path, index=False)

    def export_summary(self, path: str | Path) -> None:
        if self.summary_table is None:
            raise RuntimeError("Run the pipeline before exporting a summary.")
        self.summary_table.to_csv(path, index=False)

    def save_pipeline_outputs(
        self,
        output_dir: str | Path = ".",
        output_prefix: str | None = None,
    ) -> dict[str, Path]:
        """Save per-cell annotations and region-count summary CSVs."""

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_prefix = output_prefix or self.config.output_prefix

        paths = {
            "annotations": output_dir / f"{output_prefix}_annotations.csv",
            "summary": output_dir / f"{output_prefix}_summary.csv",
        }
        self.export_annotations(paths["annotations"])
        self.export_summary(paths["summary"])

        print("\nSaved pipeline outputs:")
        for label, path in paths.items():
            print(f"  {label}: {path}")
        return paths

    def region_counts(self) -> pd.Series:
        """Convenience accessor for the combined annotation counts."""

        if self.result_table is None:
            raise RuntimeError("Run the pipeline before requesting region counts.")
        return self.result_table["distance_region"].value_counts(dropna=False)

    def plot_spatial_annotations(
        self,
        color: str = "distance_region",
        samples: Sequence[object] | None = None,
        spot_size: float = 6,
        alpha: float = 0.85,
        ncols: int = 3,
        figsize: tuple[float, float] | None = None,
        palette: Mapping[object, str] | None = None,
        cmap: str = "viridis",
        invert_y: bool = True,
        show: bool = True,
        save_path: str | Path | None = None,
    ):
        """
        Plot spatial annotation results for all samples or a selected subset.

        ``color`` can be any column in ``adata.obs`` after running the pipeline,
        such as ``distance_region``, ``malignant_region``,
        ``malignant_proximity``, ``distance_to_malignant``, or
        ``malignant_neighbor_fraction``.
        """

        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        if color not in self.adata.obs:
            raise KeyError(f"Column '{color}' is not present in adata.obs.")
        if ncols < 1:
            raise ValueError("ncols must be at least 1.")

        coords = self._coordinates()
        plot_df = pd.DataFrame(
            {
                "x": coords[:, 0],
                "y": coords[:, 1],
                color: self.adata.obs[color].to_numpy(),
            },
            index=self.adata.obs_names,
        )
        if self.sample_col is not None:
            plot_df[self.sample_col] = self.adata.obs[self.sample_col].to_numpy()
            sample_values = list(pd.unique(plot_df[self.sample_col])) if samples is None else list(samples)
        else:
            plot_df["_sample"] = "all_cells"
            sample_values = ["all_cells"]

        n_samples = len(sample_values)
        ncols = min(ncols, n_samples)
        nrows = int(np.ceil(n_samples / ncols))
        if figsize is None:
            figsize = (5 * ncols, 5 * nrows)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        axes_flat = axes.ravel()

        color_values = plot_df[color]
        numeric_color = pd.api.types.is_numeric_dtype(color_values)
        category_color_map = None
        if not numeric_color:
            categories = pd.Series(color_values, dtype="object").dropna().unique().tolist()
            default_colors = plt.get_cmap("tab20").colors
            category_color_map = {
                category: default_colors[i % len(default_colors)] for i, category in enumerate(categories)
            }
            if palette is not None:
                category_color_map.update(palette)

        scatter_for_colorbar = None
        for ax, sample_id in zip(axes_flat, sample_values):
            if self.sample_col is not None:
                sample_df = plot_df.loc[plot_df[self.sample_col] == sample_id]
            else:
                sample_df = plot_df

            if numeric_color:
                scatter_for_colorbar = ax.scatter(
                    sample_df["x"],
                    sample_df["y"],
                    c=sample_df[color],
                    s=spot_size,
                    alpha=alpha,
                    cmap=cmap,
                    linewidths=0,
                )
            else:
                point_colors = sample_df[color].map(category_color_map).fillna("#d0d0d0")
                ax.scatter(
                    sample_df["x"],
                    sample_df["y"],
                    c=point_colors,
                    s=spot_size,
                    alpha=alpha,
                    linewidths=0,
                )

            ax.set_title(str(sample_id))
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_aspect("equal", adjustable="box")
            if invert_y:
                ax.invert_yaxis()

        for ax in axes_flat[n_samples:]:
            ax.axis("off")

        if numeric_color and scatter_for_colorbar is not None:
            fig.colorbar(scatter_for_colorbar, ax=axes_flat[:n_samples], label=color, shrink=0.75)
        elif category_color_map is not None:
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="none",
                    markerfacecolor=category_color_map[category],
                    markeredgewidth=0,
                    markersize=7,
                    label=str(category),
                )
                for category in category_color_map
            ]
            fig.legend(handles=handles, title=color, loc="center right", frameon=False)
            fig.subplots_adjust(right=0.82)

        fig.suptitle(color, y=0.995)
        fig.tight_layout()
        if save_path is not None:
            fig.savefig(save_path, bbox_inches="tight", dpi=300)
        if show:
            plt.show()
        return fig

    @staticmethod
    def default_output_columns() -> Mapping[str, str]:
        """Describe columns added to ``adata.obs`` by the analyzer."""

        return {
            "distance_to_malignant": "Euclidean distance to the nearest malignant cell.",
            "nearest_malignant_cell": "Observation name of the nearest malignant cell.",
            "neighbor_count": "Number of cells within the configured radius.",
            "malignant_neighbor_count": "Number of local neighbors labeled malignant.",
            "malignant_neighbor_fraction": "Fraction of local neighbors labeled malignant.",
            "malignant_region": "Tumor-inner, tumor-margin, or isolated label for malignant cells.",
            "malignant_proximity": "Close/intermediate/distant label for immune-stromal cells.",
            "distance_region": "Unified region label for malignant and non-malignant cells.",
        }

    @staticmethod
    def recommended_starting_parameters() -> Mapping[str, float | int]:
        """Starting point for Xenium coordinates measured in microns."""

        return {
            "neighbor_radius": 50.0,
            "close_distance": 30.0,
            "intermediate_distance": 100.0,
            "core_malignant_fraction": 0.75,
            "min_neighbors_for_region": 5,
        }
