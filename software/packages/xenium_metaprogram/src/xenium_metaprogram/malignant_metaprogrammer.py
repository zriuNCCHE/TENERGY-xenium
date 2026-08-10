from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy import sparse
from sklearn.decomposition import NMF as SklearnNMF


DEFAULT_MP_MARKER_SETS: dict[str, tuple[str, ...]] = {
    "Cell cycle / proliferation": (
        "MKI67",
        "TOP2A",
        "UBE2C",
        "CENPF",
        "PCNA",
        "MCM2",
        "MCM5",
        "MCM6",
        "CDK1",
        "CCNB1",
        "CCNB2",
    ),
    "Hypoxia / stress": (
        "VEGFA",
        "CA9",
        "ENO1",
        "LDHA",
        "SLC2A1",
        "DDIT4",
        "HILPDA",
        "ATF3",
        "JUN",
        "FOS",
        "HSPA1A",
        "HSPA1B",
    ),
    "Squamous differentiation / keratinization": (
        "KRT1",
        "KRT4",
        "KRT5",
        "KRT10",
        "KRT13",
        "KRT14",
        "KRTDAP",
        "IVL",
        "SPRR1A",
        "SPRR1B",
        "SPRR2A",
        "TGM1",
    ),
    "Basal / invasive matrix remodeling": (
        "MMP14",
        "MMP10",
        "MMP1",
        "ITGA6",
        "ITGB4",
        "COL17A1",
        "LAMA3",
        "LAMC2",
        "CD44",
        "VIM",
        "FN1",
        "SERPINE1",
    ),
    "Interferon / antigen presentation": (
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "ISG15",
        "MX1",
        "OAS1",
        "STAT1",
        "IRF1",
        "B2M",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-DRA",
    ),
    "MHC-II / immune-interacting": (
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-DQA1",
        "HLA-DQB1",
        "CD74",
        "CIITA",
    ),
    "Respiration / oxidative phosphorylation": (
        "MT-CO1",
        "MT-CO2",
        "MT-CO3",
        "MT-CYB",
        "MT-ND1",
        "MT-ND2",
        "MT-ND3",
        "MT-ND4",
        "MT-ATP6",
        "NDUFA4",
        "COX6C",
    ),
    "Ribosome / translation": (
        "RPL3",
        "RPL5",
        "RPL7",
        "RPL10",
        "RPL13",
        "RPLP0",
        "RPS3",
        "RPS6",
        "RPS8",
        "RPS18",
        "EEF1A1",
        "EFTUD2",
    ),
    "Secretory / mucinous": (
        "MUC1",
        "MUC4",
        "MUC5B",
        "AGR2",
        "AGR3",
        "TFF1",
        "TFF3",
        "SPINK5",
        "SLPI",
        "LCN2",
    ),
    "Senescence / inflammatory": (
        "CDKN1A",
        "CDKN2A",
        "GADD45A",
        "GADD45B",
        "IL6",
        "CXCL1",
        "CXCL2",
        "CXCL8",
        "SERPINE1",
        "CEBPB",
    ),
    "Cilia / motility": (
        "FOXJ1",
        "PIFO",
        "DNAH5",
        "DNAH11",
        "DNAI1",
        "RSPH1",
        "RSPH4A",
        "TPPP3",
        "CAPS",
    ),
    "Neural / lineage-related": (
        "CHGA",
        "CHGB",
        "SYP",
        "ENO2",
        "SOX2",
        "SOX9",
        "NGFR",
        "TP63",
    ),
}


@dataclass
class MetaProgramConfig:
    """Default parameters tuned for targeted Xenium malignant-cell data."""

    sample_col: str = "patient_id"
    layer: str = "counts"
    centered_layer: str = "centered_pos"
    k_range: tuple[int, int] = (7, 15)
    top_n_genes: int = 50
    n_mps: int = 15
    min_patients: int = 3
    consensus_fraction: float = 0.20
    min_cells_per_sample: int = 30
    min_nonzero_genes_per_sample: int = 25
    nmf_max_iter: int = 1000
    random_state: int = 42
    nmf_backend: str = "auto"
    heatmap_cmap: str = "nature"
    heatmap_vmax: float = 0.3
    output_prefix: str = "xenium_metaprogram"


class MalignantMetaProgrammer:
    """
    Discover conserved malignant-cell meta-programs across samples.

    This is a Xenium-adapted version of the multi-K NMF meta-program workflow
    inspired by Gavish et al. The main deviation is per-sample positive centering,
    which reduces patient-specific CNV/clonal background before NMF.
    """

    def __init__(
        self,
        adata,
        sample_col: str | None = None,
        layer: str | None = None,
        config: MetaProgramConfig | None = None,
    ):
        self.config = config or MetaProgramConfig()
        if sample_col is not None:
            self.config.sample_col = sample_col
        if layer is not None:
            self.config.layer = layer

        self.adata = adata
        self.sample_col = self.config.sample_col
        self.layer = self.config.layer

        self.all_programs: list[dict] = []
        self.similarity_matrix: pd.DataFrame | None = None
        self.meta_programs: dict[str, list[str]] = {}
        self.program_assignments: np.ndarray | None = None
        self.program_table: pd.DataFrame | None = None
        self.meta_program_cluster_map: dict[int, str] = {}

        self._validate_inputs()

    def _validate_inputs(self) -> None:
        if self.sample_col not in self.adata.obs:
            raise KeyError(f"Sample column '{self.sample_col}' is not present in adata.obs.")
        if self.layer not in self.adata.layers and self.layer != "X":
            raise KeyError(f"Layer '{self.layer}' is not present in adata.layers.")

    def preprocess_centered_pos(self, new_layer: str | None = None) -> None:
        """
        Create V_centered = max(0, V - mean(V_sample)).

        For malignant ESCC Xenium data, this suppresses patient-specific clonal/CNV
        baselines so NMF sees intra-tumour functional variation.
        """

        new_layer = new_layer or self.config.centered_layer
        print(f"Neutralizing sample clonal signals from layer '{self.layer}'...")

        v_matrix = self._matrix_for_layer(self.layer, dtype=np.float32)
        for sample in self.adata.obs[self.sample_col].unique():
            idx = np.asarray(self.adata.obs[self.sample_col] == sample)
            if idx.sum() == 0:
                continue
            v_matrix[idx] -= v_matrix[idx].mean(axis=0)

        v_matrix[v_matrix < 0] = 0
        self.adata.layers[new_layer] = v_matrix
        self.layer = new_layer
        self.config.layer = new_layer
        print(f"Created '{new_layer}'. NMF will run on sample-relative positive signal.")

    def run_pipeline(
        self,
        k_range: tuple[int, int] | None = None,
        top_n_genes: int | None = None,
        n_mps: int | None = None,
        min_patients: int | None = None,
        consensus_fraction: float | None = None,
        center_first: bool = True,
        plot_results: bool = True,
        heatmap_cmap: str | None = None,
        save: bool = True,
        output_dir: str | Path = ".",
        output_prefix: str | None = None,
    ) -> Mapping[str, list[str]]:
        """Run centering, multi-K NMF, program clustering, and cell scoring."""

        self._reset_results()
        k_range = k_range or self.config.k_range
        top_n_genes = top_n_genes or self.config.top_n_genes
        n_mps = n_mps or self.config.n_mps
        min_patients = min_patients or self.config.min_patients
        consensus_fraction = consensus_fraction or self.config.consensus_fraction

        if center_first:
            self.preprocess_centered_pos()

        print(f"\n[1/4] Running multi-K NMF for k={k_range[0]}..{k_range[1]}")
        for k in range(k_range[0], k_range[1] + 1):
            print(f"  > rank k={k}")
            self._run_sample_nmf(k=k, top_n_genes=top_n_genes)

        if len(self.all_programs) < 2:
            raise RuntimeError("Fewer than two NMF programs were produced; check filters and input data.")

        print(f"[2/4] Computing similarity across {len(self.all_programs)} programs")
        self._compute_similarity()

        print("[3/4] Clustering meta-programs")
        self._identify_meta_programs(
            n_clusters=n_mps,
            min_patient_threshold=min_patients,
            consensus_fraction=consensus_fraction,
        )

        print("[4/4] Scoring cells")
        self.score_cells()
        self.print_hallmark_report()

        print(f"Pipeline complete: {len(self.meta_programs)} MPs identified.")
        grid = None
        if plot_results or save:
            grid = self.plot_similarity_heatmap(
                cmap=heatmap_cmap or self.config.heatmap_cmap,
                show=plot_results,
            )
        if save:
            self.save_pipeline_outputs(
                output_dir=output_dir,
                output_prefix=output_prefix,
                heatmap_grid=grid,
                heatmap_cmap=heatmap_cmap or self.config.heatmap_cmap,
            )
        return self.meta_programs

    def _reset_results(self) -> None:
        self.all_programs = []
        self.similarity_matrix = None
        self.meta_programs = {}
        self.program_assignments = None
        self.program_table = None
        self.meta_program_cluster_map = {}

    def _matrix_for_layer(self, layer: str, dtype=np.float32) -> np.ndarray:
        data = self.adata.X if layer == "X" else self.adata.layers[layer]
        if sparse.issparse(data):
            data = data.toarray()
        return np.asarray(data, dtype=dtype).copy()

    def _sample_matrix(self, sample_id) -> tuple[np.ndarray, pd.Index]:
        sample_adata = self.adata[self.adata.obs[self.sample_col] == sample_id].copy()
        data = sample_adata.X if self.layer == "X" else sample_adata.layers[self.layer]
        if sparse.issparse(data):
            data = data.toarray()
        data = np.asarray(data, dtype=np.float32)

        nonzero_gene_mask = np.asarray((data > 0).sum(axis=0)).ravel() > 0
        data = data[:, nonzero_gene_mask]
        feature_names = sample_adata.var_names[nonzero_gene_mask]
        return data, feature_names

    def _get_nmf_loadings(self, sample_id, k: int) -> tuple[np.ndarray, pd.Index]:
        data, feature_names = self._sample_matrix(sample_id)
        if data.shape[0] < self.config.min_cells_per_sample:
            raise ValueError(f"Only {data.shape[0]} cells available.")
        if data.shape[1] < max(k, self.config.min_nonzero_genes_per_sample):
            raise ValueError(f"Only {data.shape[1]} nonzero genes available.")

        backend = self._resolve_nmf_backend()
        if isinstance(backend, dict):
            model = backend["class"](
                n_components=k,
                init="nndsvd",
                random_state=self.config.random_state,
                max_iter=self.config.nmf_max_iter,
            )
            model.fit(backend["cupy"].asarray(data, dtype=backend["cupy"].float32))
            components = backend["cupy"].asnumpy(model.components_)
        else:
            model = SklearnNMF(
                n_components=k,
                init="nndsvda",
                random_state=self.config.random_state,
                max_iter=self.config.nmf_max_iter,
            )
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                model.fit(data)
            components = model.components_

        return components, feature_names

    def _resolve_nmf_backend(self):
        if self.config.nmf_backend == "sklearn":
            return "sklearn"
        if self.config.nmf_backend in {"auto", "cuml"}:
            try:
                import cupy as cp
                from cuml.decomposition import NMF as CuNMF

                return {"class": CuNMF, "cupy": cp}
            except Exception:
                if self.config.nmf_backend == "cuml":
                    raise
        return "sklearn"

    def _run_sample_nmf(self, k: int, top_n_genes: int) -> None:
        for sample_id in self.adata.obs[self.sample_col].unique():
            try:
                h_matrix, feature_names = self._get_nmf_loadings(sample_id, k)
                for i in range(k):
                    top_indices = h_matrix[i].argsort()[-top_n_genes:][::-1]
                    genes = feature_names[top_indices].tolist()
                    self.all_programs.append(
                        {
                            "sample": sample_id,
                            "k": k,
                            "component": i,
                            "genes": set(genes),
                            "program_id": f"{sample_id}_k{k}_P{i}",
                        }
                    )
            except Exception as exc:
                print(f"  ! Skipping {sample_id} at k={k}: {exc}")

    def _compute_similarity(self) -> None:
        n_programs = len(self.all_programs)
        sim_mat = np.zeros((n_programs, n_programs), dtype=np.float32)
        for i in range(n_programs):
            sim_mat[i, i] = 1.0
            for j in range(i + 1, n_programs):
                genes_i = self.all_programs[i]["genes"]
                genes_j = self.all_programs[j]["genes"]
                union = len(genes_i | genes_j)
                score = len(genes_i & genes_j) / union if union else 0
                sim_mat[i, j] = sim_mat[j, i] = score

        ids = [p["program_id"] for p in self.all_programs]
        self.similarity_matrix = pd.DataFrame(sim_mat, index=ids, columns=ids)

    def _identify_meta_programs(
        self,
        n_clusters: int,
        min_patient_threshold: int,
        consensus_fraction: float,
    ) -> None:
        if self.similarity_matrix is None:
            raise RuntimeError("Similarity matrix has not been computed.")

        distance = 1 - self.similarity_matrix.to_numpy()
        np.fill_diagonal(distance, 0)
        condensed = squareform(distance, checks=False)
        z_matrix = linkage(condensed, method="average")
        self.program_assignments = fcluster(z_matrix, t=n_clusters, criterion="maxclust")

        prog_df = pd.DataFrame(self.all_programs)
        prog_df["cluster"] = self.program_assignments
        retained_rows = []
        for cluster_id, subset in prog_df.groupby("cluster", sort=True):
            if subset["sample"].nunique() < min_patient_threshold:
                continue
            all_genes = [gene for gene_set in subset["genes"] for gene in gene_set]
            counts = pd.Series(all_genes).value_counts()
            min_count = max(1, int(np.ceil(len(subset) * consensus_fraction)))
            consensus = counts[counts >= min_count].index.tolist()
            if consensus:
                retained_rows.append((cluster_id, consensus))

        for mp_idx, (cluster_id, consensus) in enumerate(retained_rows, start=1):
            mp_name = f"MP_{mp_idx}"
            self.meta_programs[mp_name] = consensus
            self.meta_program_cluster_map[cluster_id] = mp_name

        prog_df["raw_cluster"] = prog_df["cluster"]
        prog_df["meta_program"] = prog_df["cluster"].map(self.meta_program_cluster_map)
        self.program_table = prog_df

    def score_cells(self) -> None:
        for mp_name, genes in self.meta_programs.items():
            genes_present = [gene for gene in genes if gene in self.adata.var_names]
            if not genes_present:
                print(f"  ! Skipping {mp_name}: no genes present in adata.var_names")
                continue
            sc.tl.score_genes(self.adata, gene_list=genes_present, score_name=mp_name)

    def print_hallmark_report(
        self,
        markers: Sequence[str] = ("MKI67", "TOP2A", "VEGFA", "MMP14", "COL17A1", "LAMA3"),
    ) -> None:
        print("\n--- Hallmark Gene Report ---")
        for marker in markers:
            found_in = [mp for mp, genes in self.meta_programs.items() if marker in genes]
            status = f"found in {found_in}" if found_in else "not found"
            print(f"{marker:8}: {status}")

    def export_meta_programs(self, path: str) -> None:
        rows = [
            {"meta_program": mp, "rank": rank + 1, "gene": gene}
            for mp, genes in self.meta_programs.items()
            for rank, gene in enumerate(genes)
        ]
        pd.DataFrame(rows).to_csv(path, index=False)

    def annotate_meta_programs(
        self,
        marker_sets: Mapping[str, Sequence[str]] | None = None,
        min_overlap: int = 1,
    ) -> pd.DataFrame:
        """
        Suggest biological names for MPs from overlaps with marker gene sets.

        The output is intentionally a suggestion table, not a hard relabeling step.
        Inspect the top genes and spatial/clinical patterns before using names in a manuscript.
        """

        marker_sets = marker_sets or DEFAULT_MP_MARKER_SETS
        rows = []
        for mp_name, genes in self.meta_programs.items():
            gene_set = set(genes)
            best_label = "Unassigned"
            best_overlap = 0
            best_fraction = 0.0
            best_genes: list[str] = []
            for label, markers in marker_sets.items():
                markers_present = [marker for marker in markers if marker in gene_set]
                overlap = len(markers_present)
                fraction = overlap / len(markers) if markers else 0.0
                if (overlap, fraction) > (best_overlap, best_fraction):
                    best_label = label
                    best_overlap = overlap
                    best_fraction = fraction
                    best_genes = markers_present

            if best_overlap < min_overlap:
                best_label = "Unassigned"
                best_genes = []
                best_fraction = 0.0

            rows.append(
                {
                    "meta_program": mp_name,
                    "suggested_label": best_label,
                    "marker_overlap": best_overlap,
                    "marker_fraction": best_fraction,
                    "overlap_genes": ", ".join(best_genes),
                    "n_consensus_genes": len(genes),
                    "top_consensus_genes": ", ".join(genes[:15]),
                }
            )

        return pd.DataFrame(rows)

    def export_meta_program_annotations(
        self,
        path: str,
        marker_sets: Mapping[str, Sequence[str]] | None = None,
        min_overlap: int = 1,
    ) -> pd.DataFrame:
        annotations = self.annotate_meta_programs(marker_sets=marker_sets, min_overlap=min_overlap)
        annotations.to_csv(path, index=False)
        return annotations

    def save_pipeline_outputs(
        self,
        output_dir: str | Path = ".",
        output_prefix: str | None = None,
        heatmap_grid=None,
        heatmap_cmap: str | None = None,
    ) -> dict[str, Path]:
        """Save MP genes, suggested annotations, and the similarity heatmap PDF."""

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_prefix = output_prefix or self.config.output_prefix

        paths = {
            "genes": output_dir / f"{output_prefix}_genes.csv",
            "annotations": output_dir / f"{output_prefix}_annotations.csv",
            "heatmap": output_dir / f"{output_prefix}_similarity_heatmap.pdf",
        }

        self.export_meta_programs(str(paths["genes"]))
        self.export_meta_program_annotations(str(paths["annotations"]))

        if heatmap_grid is None:
            heatmap_grid = self.plot_similarity_heatmap(
                cmap=heatmap_cmap or self.config.heatmap_cmap,
                show=False,
            )
        heatmap_grid.fig.savefig(paths["heatmap"], bbox_inches="tight")
        print("\nSaved pipeline outputs:")
        for label, path in paths.items():
            print(f"  {label}: {path}")
        return paths

    @staticmethod
    def available_heatmap_colormaps() -> dict[str, str]:
        """Named colormap options for similarity heatmaps."""

        return {
            "nature": "rocket_r",
            "mako": "mako",
            "magma": "magma",
            "rocket": "rocket",
            "viridis": "viridis",
            "icefire": "icefire",
            "flare": "flare",
        }

    def plot_similarity_heatmap(
        self,
        cmap: str = "nature",
        vmax: float | None = None,
        figsize: tuple[float, float] = (10, 10),
        label_mps: bool = True,
        show_filtered_clusters: bool = False,
        show: bool = True,
    ):
        if self.similarity_matrix is None or self.program_assignments is None:
            raise RuntimeError("Run the pipeline before plotting the similarity heatmap.")

        if self.program_table is not None and "meta_program" in self.program_table:
            program_meta = self.program_table.set_index("program_id")["meta_program"]
        else:
            program_meta = pd.Series(index=self.similarity_matrix.index, dtype=object)

        order_df = (
            pd.DataFrame(
                {
                    "program_id": self.similarity_matrix.index,
                    "cluster": self.program_assignments,
                }
            )
            .assign(meta_program=lambda df: df["program_id"].map(program_meta))
            .sort_values("cluster")
            .reset_index(drop=True)
        )
        reordered_mat = self.similarity_matrix.loc[order_df["program_id"], order_df["program_id"]]

        mp_names = list(self.meta_programs.keys())
        mp_palette = sns.color_palette("tab20", max(len(mp_names), 1))
        mp_color_map = dict(zip(mp_names, mp_palette))
        filtered_color = (0.82, 0.82, 0.82)
        if show_filtered_clusters:
            filtered_clusters = sorted(
                set(order_df["cluster"].unique()) - set(self.meta_program_cluster_map)
            )
            filtered_palette = sns.color_palette("Greys", max(len(filtered_clusters) + 2, 3))[2:]
            filtered_color_map = dict(zip(filtered_clusters, filtered_palette))
            network_colors = [
                mp_color_map.get(row.meta_program, filtered_color_map.get(row.cluster, filtered_color))
                for row in order_df.itertuples()
            ]
        else:
            network_colors = [
                mp_color_map.get(mp_name, filtered_color) for mp_name in order_df["meta_program"]
            ]

        cmap_name = self.available_heatmap_colormaps().get(cmap, cmap)
        vmax = self.config.heatmap_vmax if vmax is None else vmax

        grid = sns.clustermap(
            reordered_mat,
            row_cluster=False,
            col_cluster=False,
            row_colors=network_colors,
            col_colors=network_colors,
            cmap=cmap_name,
            vmin=0,
            vmax=vmax,
            figsize=figsize,
            xticklabels=False,
            yticklabels=False,
            cbar_kws={"label": "Similarity (Jaccard index)"},
        )
        self._decorate_similarity_heatmap(grid, order_df, mp_color_map, label_mps)
        if show:
            plt.show()
        return grid

    def plot_nature_style_heatmap(self, vmax: float | None = None, cmap: str = "nature"):
        """Backward-compatible alias for the standalone similarity heatmap plotter."""

        return self.plot_similarity_heatmap(cmap=cmap, vmax=vmax)

    def _decorate_similarity_heatmap(
        self,
        grid,
        order_df: pd.DataFrame,
        mp_color_map: Mapping[str, tuple[float, float, float]],
        label_mps: bool,
    ) -> None:
        ax = grid.ax_heatmap
        n_programs = len(order_df)
        if n_programs == 0:
            return

        boundaries = np.flatnonzero(order_df["cluster"].to_numpy()[1:] != order_df["cluster"].to_numpy()[:-1]) + 1
        for boundary in boundaries:
            ax.axhline(boundary, color="white", linewidth=0.35, alpha=0.8)
            ax.axvline(boundary, color="white", linewidth=0.35, alpha=0.8)

        if not label_mps:
            return

        mp_blocks = order_df.dropna(subset=["meta_program"]).groupby("meta_program", sort=False)
        for mp_name, subset in mp_blocks:
            start = subset.index.min()
            end = subset.index.max() + 1
            center = start + (end - start) / 2
            color = mp_color_map.get(mp_name, "black")
            ax.text(
                n_programs + max(n_programs * 0.012, 2),
                center,
                mp_name,
                va="center",
                ha="left",
                color=color,
                fontsize=9,
                fontweight="bold",
                clip_on=False,
            )
            ax.text(
                center,
                -max(n_programs * 0.012, 2),
                mp_name,
                va="bottom",
                ha="center",
                rotation=90,
                color=color,
                fontsize=9,
                fontweight="bold",
                clip_on=False,
            )

    def plot_spatial_samples(
        self,
        mp_ids: Iterable[str] | None = None,
        spot_size: float = 30,
        cmap: str = "magma",
    ):
        mps_to_plot = list(mp_ids) if mp_ids else list(self.meta_programs.keys())
        samples = self.adata.obs[self.sample_col].unique()
        n_samples = len(samples)
        n_mps = len(mps_to_plot)
        fig, axes = plt.subplots(n_samples, n_mps, figsize=(n_mps * 5, n_samples * 5))
        axes = np.asarray(axes)
        if axes.ndim == 0:
            axes = axes.reshape(1, 1)
        elif n_samples == 1:
            axes = axes.reshape(1, -1)
        elif n_mps == 1:
            axes = axes.reshape(-1, 1)

        for i, sample_id in enumerate(samples):
            sample_adata = self.adata[self.adata.obs[self.sample_col] == sample_id]
            for j, mp_id in enumerate(mps_to_plot):
                sc.pl.spatial(
                    sample_adata,
                    color=mp_id,
                    cmap=cmap,
                    spot_size=spot_size,
                    ax=axes[i, j],
                    show=False,
                    title=f"{sample_id} | {mp_id}",
                    frameon=False,
                )
                axes[i, j].set_xlabel("")
                axes[i, j].set_ylabel("")

        plt.tight_layout()
        plt.show()
        return fig

    def test_sample_run(self, sample_id=None, k: int = 8, top_n_genes: int = 30):
        if sample_id is None:
            sample_id = self.adata.obs[self.sample_col].unique()[0]
        h_matrix, names = self._get_nmf_loadings(sample_id, k)
        fig = plt.figure(figsize=(12, 6))
        for i in range(min(k, 4)):
            top = h_matrix[i].argsort()[-top_n_genes:][::-1]
            plt.subplot(2, 2, i + 1)
            plt.bar(names[top], h_matrix[i][top], color="gray")
            plt.xticks(rotation=90, fontsize=7)
            plt.title(f"Program {i}")
        plt.tight_layout()
        plt.show()
        return fig

    @staticmethod
    def drop_zero_variance_columns(df: pd.DataFrame) -> pd.DataFrame:
        """Guard before correlation/clustering on sparse Xenium subsets."""

        return df.loc[:, df.std(axis=0) > 0]
