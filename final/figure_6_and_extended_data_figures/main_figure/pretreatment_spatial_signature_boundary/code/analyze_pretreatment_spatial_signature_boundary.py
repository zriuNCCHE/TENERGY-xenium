#!/usr/bin/env python3
"""
Pretreatment malignant-cell spatial signature analysis for Figure 6.

This script tests whether MP7-related malignant-cell programs are enriched at
the tumor margin versus tumor inner compartment, and whether that spatial
enrichment differs between cCR and non-cCR patients.

Use from the command line on the server:

    python analyze_pretreatment_spatial_signature_boundary.py \
        --input-h5ad /path/to/pre_treat_adata.h5ad \
        --output-root /path/to/final/figure_6/pretreatment_spatial_signature_boundary

Use inside a notebook where pre_treat_epi_adata already exists:

    from pathlib import Path
    from analyze_pretreatment_spatial_signature_boundary import run_analysis

    run_analysis(
        pre_treat_epi_adata,
        output_root=Path("final/figure_6/pretreatment_spatial_signature_boundary"),
        layer=None,
        make_plots=False,
    )

Important statistical choice:
    P values are computed on sample-level summaries, not cell-level values.
    This avoids treating thousands of cells from the same sample as independent
    biological replicates.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Iterable
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.stats import mannwhitneyu, wilcoxon


MP7_FALLBACK = [
    "LAMB3",
    "LAMC2",
    "PLEC",
    "COL17A1",
    "LAMA3",
    "ITGA6",
    "ITGB4",
    "APP",
    "TNC",
    "PLAU",
    "INHBA",
    "SDC4",
    "ITGA3",
    "CDH3",
    "PTHLH",
    "NDRG1",
    "MYH9",
    "DST",
    "SLC2A1",
    "PRNP",
    "SERPINE1",
    "ANO1",
    "PTGFRN",
    "SLC7A5",
    "MMP14",
    "CDKN1A",
    "ANXA2",
]

SIGNATURES = {
    "MP7 all": MP7_FALLBACK,
    "Laminin-332": ["LAMA3", "LAMB3", "LAMC2"],
    "Adhesion / anchoring": ["COL17A1", "ITGA6", "ITGB4", "ITGA3", "PLEC"],
    "Remodeling / invasion": ["MMP14", "PLAU", "SERPINE1", "TNC", "INHBA"],
    "Stress / adaptation": ["NDRG1", "MYH9", "SLC2A1", "SLC7A5", "CDKN1A"],
}

PLOT_SIGNATURES = [
    "MP7 all",
    "Laminin-332",
    "Adhesion / anchoring",
    "Remodeling / invasion",
    "Stress / adaptation",
]

GENE_PANEL = [
    "LAMA3",
    "LAMB3",
    "LAMC2",
    "COL17A1",
    "ITGB4",
    "ITGA3",
    "ITGA6",
    "PLEC",
    "MMP14",
    "PLAU",
    "SERPINE1",
    "TNC",
    "INHBA",
    "NDRG1",
    "MYH9",
    "SLC2A1",
    "SLC7A5",
    "CDKN1A",
]

COMPARISON_REGION_ORDER = ["tumor_inner", "tumor_margin"]
REGION_ORDER = ["tumor_inner", "tumor_margin", "isolated_malignant"]
REGION_LABELS = {
    "tumor_inner": "Inner",
    "tumor_margin": "Margin",
    "isolated_malignant": "Isolated",
}
RESPONSE_ORDER = ["cCR", "non-cCR"]
# Matched to config/palettes.yaml clinical_outcome.
RESPONSE_COLORS = {"cCR": "#7BB6A4", "non-cCR": "#F2A38A"}
SIGNATURE_COLORS = {
    "MP7 all": "#6F6F6F",
    "Laminin-332": "#5B8DB8",
    "Adhesion / anchoring": "#7E6AAD",
    "Remodeling / invasion": "#C46A4A",
    "Stress / adaptation": "#6AA378",
}


@dataclass
class AnalysisConfig:
    layer: str | None = None
    log1p_transform: bool = False
    major_cell_type_col: str = "major_cell_type"
    malignant_label: str = "malignant"
    timepoint_col: str = "sample_timepoint"
    pre_label: str = "pre"
    region_col: str = "malignant_region"
    response_col: str = "cCR"
    sample_col: str = "sampleID"
    patient_col: str = "patientID"
    min_cells_per_region: int = 20


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def load_mp7_genes(project_root: Path | None = None) -> list[str]:
    """Load MP7 genes from the project result table if available."""
    if project_root is None:
        project_root = project_root_from_script()
    mp_file = (
        project_root
        / "final"
        / "figure_6"
        / "metaprogram_similarity_heatmap"
        / "results"
        / "xenium_metaprogram_genes.csv"
    )
    if not mp_file.exists():
        return MP7_FALLBACK
    table = pd.read_csv(mp_file)
    required = {"meta_program", "gene"}
    if not required.issubset(table.columns):
        return MP7_FALLBACK
    genes = (
        table.loc[table["meta_program"].astype(str).eq("MP_7"), "gene"]
        .dropna()
        .astype(str)
        .tolist()
    )
    return genes if genes else MP7_FALLBACK


def sanitize_response(value) -> str:
    text = str(value).strip()
    low = text.lower()
    if low in {"ccr", "true", "1", "yes", "y", "complete_response"}:
        return "cCR"
    if low in {"non-ccr", "nonccr", "false", "0", "no", "n", "incomplete_response"}:
        return "non-cCR"
    return text


def get_matrix(adata, layer: str | None = None):
    if layer is None:
        mat = adata.X
    else:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' was not found in adata.layers.")
        mat = adata.layers[layer]
    if sparse.issparse(mat):
        return mat.tocsr()
    return np.asarray(mat)


def vector_from_gene(adata, gene: str, matrix) -> np.ndarray:
    gene_idx = adata.var_names.get_loc(gene)
    column = matrix[:, gene_idx]
    if sparse.issparse(column):
        return np.asarray(column.toarray()).ravel()
    return np.asarray(column).ravel()


def zscore(values: np.ndarray) -> np.ndarray:
    values = values.astype(float)
    mean = np.nanmean(values)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros(values.shape[0], dtype=float)
    return (values - mean) / sd


def subset_pretreatment_malignant(adata, config: AnalysisConfig):
    obs = adata.obs
    keep = pd.Series(True, index=obs.index)
    if config.major_cell_type_col in obs.columns:
        keep &= obs[config.major_cell_type_col].astype(str).eq(config.malignant_label)
    if config.timepoint_col in obs.columns:
        keep &= obs[config.timepoint_col].astype(str).eq(config.pre_label)
    return adata[keep.to_numpy()].copy()


def compute_cell_scores(adata, config: AnalysisConfig, project_root: Path | None = None):
    mp7_genes = load_mp7_genes(project_root)
    signatures = dict(SIGNATURES)
    signatures["MP7 all"] = mp7_genes

    available_genes = set(map(str, adata.var_names))
    matrix = get_matrix(adata, config.layer)
    obs = adata.obs.copy()

    score_table = obs[
        [
            col
            for col in [
                config.sample_col,
                config.patient_col,
                config.response_col,
                config.region_col,
            ]
            if col in obs.columns
        ]
    ].copy()
    score_table[config.response_col] = score_table[config.response_col].map(sanitize_response)

    missing_records = []
    gene_expression = {}

    genes_to_extract = []
    for genes in signatures.values():
        genes_to_extract.extend(genes)
    genes_to_extract.extend(GENE_PANEL)
    genes_to_extract = list(dict.fromkeys(genes_to_extract))

    for gene in genes_to_extract:
        if gene not in available_genes:
            missing_records.append({"feature": gene, "type": "gene", "status": "missing"})
            continue
        values = vector_from_gene(adata, gene, matrix)
        if config.log1p_transform:
            values = np.log1p(values)
        gene_expression[gene] = values.astype(float)
        score_table[f"gene__{gene}"] = values.astype(float)

    for signature_name, genes in signatures.items():
        genes_present = [gene for gene in genes if gene in gene_expression]
        genes_missing = [gene for gene in genes if gene not in gene_expression]
        for gene in genes_missing:
            missing_records.append(
                {"feature": signature_name, "type": "signature_gene", "status": f"missing:{gene}"}
            )
        if len(genes_present) == 0:
            score_table[f"score__{signature_name}"] = np.nan
            continue
        z_vectors = np.vstack([zscore(gene_expression[gene]) for gene in genes_present])
        score_table[f"score__{signature_name}"] = np.nanmean(z_vectors, axis=0)

    missing_table = pd.DataFrame(missing_records)
    return score_table, missing_table


def summarize_by_sample_region(cell_scores: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    feature_cols = [col for col in cell_scores.columns if col.startswith(("score__", "gene__"))]
    group_cols = [config.sample_col, config.patient_col, config.response_col, config.region_col]
    group_cols = [col for col in group_cols if col in cell_scores.columns]

    summary = (
        cell_scores.groupby(group_cols, observed=True)[feature_cols]
        .agg(["mean", "median"])
        .reset_index()
    )
    summary.columns = [
        "__".join([part for part in col if part]) if isinstance(col, tuple) else col
        for col in summary.columns
    ]
    counts = (
        cell_scores.groupby(group_cols, observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    summary = summary.merge(counts, on=group_cols, how="left")
    return summary


def filter_min_cells(sample_region: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    return sample_region.loc[sample_region["n_cells"] >= config.min_cells_per_region].copy()


def comparison_regions_only(sample_region: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    """Keep only inner and margin regions for publication comparisons."""
    return sample_region.loc[
        sample_region[config.region_col].isin(COMPARISON_REGION_ORDER)
    ].copy()


def boundary_delta_table(sample_region: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    feature_cols = [
        col
        for col in sample_region.columns
        if col.startswith(("score__", "gene__")) and col.endswith("__mean")
    ]
    metadata = [
        col
        for col in [config.sample_col, config.patient_col, config.response_col]
        if col in sample_region.columns
    ]
    rows = []
    for feature_col in feature_cols:
        wide = sample_region.pivot_table(
            index=metadata,
            columns=config.region_col,
            values=feature_col,
            aggfunc="mean",
        ).reset_index()
        if "tumor_margin" not in wide.columns or "tumor_inner" not in wide.columns:
            continue
        tmp = wide[metadata].copy()
        tmp["feature"] = feature_col.replace("__mean", "")
        tmp["margin_mean"] = wide["tumor_margin"]
        tmp["inner_mean"] = wide["tumor_inner"]
        tmp["margin_minus_inner"] = wide["tumor_margin"] - wide["tumor_inner"]
        tmp["margin_over_inner"] = (wide["tumor_margin"] + 1e-9) / (wide["tumor_inner"] + 1e-9)
        rows.append(tmp)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def safe_mannwhitney(x: Iterable[float], y: Iterable[float]) -> float:
    x = pd.Series(x).dropna().astype(float)
    y = pd.Series(y).dropna().astype(float)
    if len(x) < 2 or len(y) < 2:
        return np.nan
    try:
        return mannwhitneyu(x, y, alternative="two-sided").pvalue
    except ValueError:
        return np.nan


def safe_wilcoxon(x: Iterable[float]) -> float:
    x = pd.Series(x).dropna().astype(float)
    if len(x) < 2:
        return np.nan
    if np.allclose(x, 0):
        return np.nan
    try:
        return wilcoxon(x, alternative="two-sided").pvalue
    except ValueError:
        return np.nan


def compute_statistics(
    sample_region: pd.DataFrame,
    boundary_delta: pd.DataFrame,
    config: AnalysisConfig,
) -> pd.DataFrame:
    rows = []
    feature_cols = [
        col
        for col in sample_region.columns
        if col.startswith(("score__", "gene__")) and col.endswith("__mean")
    ]

    for feature_col in feature_cols:
        feature = feature_col.replace("__mean", "")

        margin = sample_region.loc[
            sample_region[config.region_col].eq("tumor_margin"), feature_col
        ]
        inner = sample_region.loc[
            sample_region[config.region_col].eq("tumor_inner"), feature_col
        ]
        rows.append(
            {
                "feature": feature,
                "test": "margin_vs_inner_unpaired_samples",
                "comparison": "tumor_margin vs tumor_inner",
                "group_1": "tumor_margin",
                "group_2": "tumor_inner",
                "n_1": int(margin.dropna().shape[0]),
                "n_2": int(inner.dropna().shape[0]),
                "mean_1": float(margin.dropna().mean()) if margin.dropna().shape[0] else np.nan,
                "mean_2": float(inner.dropna().mean()) if inner.dropna().shape[0] else np.nan,
                "effect_mean_difference": float(margin.dropna().mean() - inner.dropna().mean())
                if margin.dropna().shape[0] and inner.dropna().shape[0]
                else np.nan,
                "p_value": safe_mannwhitney(margin, inner),
            }
        )

        margin_table = sample_region.loc[
            sample_region[config.region_col].eq("tumor_margin")
        ]
        ccr_margin = margin_table.loc[margin_table[config.response_col].eq("cCR"), feature_col]
        nonccr_margin = margin_table.loc[
            margin_table[config.response_col].eq("non-cCR"), feature_col
        ]
        rows.append(
            {
                "feature": feature,
                "test": "margin_nonccr_vs_ccr",
                "comparison": "non-cCR margin vs cCR margin",
                "group_1": "non-cCR",
                "group_2": "cCR",
                "n_1": int(nonccr_margin.dropna().shape[0]),
                "n_2": int(ccr_margin.dropna().shape[0]),
                "mean_1": float(nonccr_margin.dropna().mean()) if nonccr_margin.dropna().shape[0] else np.nan,
                "mean_2": float(ccr_margin.dropna().mean()) if ccr_margin.dropna().shape[0] else np.nan,
                "effect_mean_difference": float(nonccr_margin.dropna().mean() - ccr_margin.dropna().mean())
                if nonccr_margin.dropna().shape[0] and ccr_margin.dropna().shape[0]
                else np.nan,
                "p_value": safe_mannwhitney(nonccr_margin, ccr_margin),
            }
        )

        for region in COMPARISON_REGION_ORDER:
            region_table = sample_region.loc[sample_region[config.region_col].eq(region)]
            ccr_values = region_table.loc[
                region_table[config.response_col].eq("cCR"), feature_col
            ]
            nonccr_values = region_table.loc[
                region_table[config.response_col].eq("non-cCR"), feature_col
            ]
            rows.append(
                {
                    "feature": feature,
                    "test": "region_nonccr_vs_ccr",
                    "comparison": f"non-cCR vs cCR in {region}",
                    "region": region,
                    "group_1": "non-cCR",
                    "group_2": "cCR",
                    "n_1": int(nonccr_values.dropna().shape[0]),
                    "n_2": int(ccr_values.dropna().shape[0]),
                    "mean_1": float(nonccr_values.dropna().mean()) if nonccr_values.dropna().shape[0] else np.nan,
                    "mean_2": float(ccr_values.dropna().mean()) if ccr_values.dropna().shape[0] else np.nan,
                    "effect_mean_difference": float(nonccr_values.dropna().mean() - ccr_values.dropna().mean())
                    if nonccr_values.dropna().shape[0] and ccr_values.dropna().shape[0]
                    else np.nan,
                    "p_value": safe_mannwhitney(nonccr_values, ccr_values),
                }
            )

    if not boundary_delta.empty:
        for feature, feature_df in boundary_delta.groupby("feature", observed=True):
            rows.append(
                {
                    "feature": feature,
                    "test": "paired_boundary_enrichment",
                    "comparison": "margin_minus_inner != 0",
                    "group_1": "margin_minus_inner",
                    "group_2": "0",
                    "n_1": int(feature_df["margin_minus_inner"].dropna().shape[0]),
                    "n_2": np.nan,
                    "mean_1": float(feature_df["margin_minus_inner"].dropna().mean())
                    if feature_df["margin_minus_inner"].dropna().shape[0]
                    else np.nan,
                    "mean_2": 0.0,
                    "effect_mean_difference": float(feature_df["margin_minus_inner"].dropna().mean())
                    if feature_df["margin_minus_inner"].dropna().shape[0]
                    else np.nan,
                    "p_value": safe_wilcoxon(feature_df["margin_minus_inner"]),
                }
            )
            ccr_delta = feature_df.loc[
                feature_df[config.response_col].eq("cCR"), "margin_minus_inner"
            ]
            nonccr_delta = feature_df.loc[
                feature_df[config.response_col].eq("non-cCR"), "margin_minus_inner"
            ]
            rows.append(
                {
                    "feature": feature,
                    "test": "boundary_enrichment_by_response",
                    "comparison": "non-cCR delta vs cCR delta",
                    "group_1": "non-cCR",
                    "group_2": "cCR",
                    "n_1": int(nonccr_delta.dropna().shape[0]),
                    "n_2": int(ccr_delta.dropna().shape[0]),
                    "mean_1": float(nonccr_delta.dropna().mean()) if nonccr_delta.dropna().shape[0] else np.nan,
                    "mean_2": float(ccr_delta.dropna().mean()) if ccr_delta.dropna().shape[0] else np.nan,
                    "effect_mean_difference": float(nonccr_delta.dropna().mean() - ccr_delta.dropna().mean())
                    if nonccr_delta.dropna().shape[0] and ccr_delta.dropna().shape[0]
                    else np.nan,
                    "p_value": safe_mannwhitney(nonccr_delta, ccr_delta),
                }
            )
    stats = pd.DataFrame(rows)
    if not stats.empty:
        stats["p_value"] = stats["p_value"].astype(float)
        stats["p_adj_bh_within_test"] = np.nan
        for test_name, idx in stats.groupby("test", observed=True).groups.items():
            pvals = stats.loc[idx, "p_value"]
            stats.loc[idx, "p_adj_bh_within_test"] = benjamini_hochberg(pvals)
    return stats


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    p = pd.Series(p_values, dtype=float)
    adjusted = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return adjusted.to_numpy()
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    n = len(ranked)
    raw_adj = ranked * n / np.arange(1, n + 1)
    monotone = np.minimum.accumulate(raw_adj[::-1])[::-1]
    adjusted.loc[order] = np.minimum(monotone, 1.0)
    return adjusted.to_numpy()


def clean_feature_label(feature: str) -> str:
    return feature.replace("score__", "").replace("gene__", "")


def p_value_label(p_value: float | None) -> str:
    if p_value is None or not np.isfinite(p_value):
        return "p = NA"
    if p_value < 0.001:
        return "p < 0.001"
    return f"p = {p_value:.3f}"


def lookup_p_value(
    stats: pd.DataFrame,
    feature: str,
    test: str,
    region: str | None = None,
) -> float:
    if stats.empty:
        return np.nan
    mask = stats["feature"].eq(feature) & stats["test"].eq(test)
    if region is not None and "region" in stats.columns:
        mask &= stats["region"].eq(region)
    values = stats.loc[mask, "p_value"].dropna()
    if values.empty:
        return np.nan
    return float(values.iloc[0])


def setup_plotting():
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "Helvetica",
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.5,
        }
    )


def plot_signature_region_response(
    sample_region: pd.DataFrame,
    stats: pd.DataFrame,
    outputs_dir: Path,
    config: AnalysisConfig,
):
    setup_plotting()
    plot_df = sample_region.copy()
    plot_df = plot_df.loc[plot_df[config.region_col].isin(COMPARISON_REGION_ORDER)]
    plot_df["region_label"] = plot_df[config.region_col].map(REGION_LABELS)
    region_order = [
        REGION_LABELS[r]
        for r in COMPARISON_REGION_ORDER
        if r in plot_df[config.region_col].unique()
    ]

    n_cols = 3
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.1, 4.6), constrained_layout=True)
    axes = axes.ravel()

    for ax, signature in zip(axes, PLOT_SIGNATURES):
        col = f"score__{signature}__mean"
        if col not in plot_df.columns:
            ax.axis("off")
            continue
        sns.boxplot(
            data=plot_df,
            x="region_label",
            y=col,
            hue=config.response_col,
            order=region_order,
            hue_order=RESPONSE_ORDER,
            palette=RESPONSE_COLORS,
            width=0.7,
            fliersize=0,
            linewidth=0.6,
            ax=ax,
        )
        sns.stripplot(
            data=plot_df,
            x="region_label",
            y=col,
            hue=config.response_col,
            order=region_order,
            hue_order=RESPONSE_ORDER,
            palette=RESPONSE_COLORS,
            dodge=True,
            size=2.5,
            alpha=0.85,
            linewidth=0.2,
            edgecolor="white",
            ax=ax,
        )
        ax.set_title(signature)
        ax.set_xlabel("")
        ax.set_ylabel("Mean z-scored signature")
        ax.tick_params(axis="x", rotation=25)
        y_values = plot_df[col].dropna()
        if not y_values.empty:
            y_min, y_max = float(y_values.min()), float(y_values.max())
            y_span = max(y_max - y_min, 0.1)
            ax.set_ylim(y_min - 0.12 * y_span, y_max + 0.26 * y_span)
            for i, region in enumerate(COMPARISON_REGION_ORDER):
                if REGION_LABELS[region] not in region_order:
                    continue
                p = lookup_p_value(stats, f"score__{signature}", "region_nonccr_vs_ccr", region)
                ax.text(
                    i,
                    y_max + (0.08 + 0.08 * i) * y_span,
                    p_value_label(p),
                    ha="center",
                    va="bottom",
                    fontsize=5.5,
                )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend_.remove()
    for ax in axes[len(PLOT_SIGNATURES) :]:
        ax.axis("off")

    fig.legend(
        handles[: len(RESPONSE_ORDER)],
        labels[: len(RESPONSE_ORDER)],
        loc="lower center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.savefig(outputs_dir / "figure_6d_signature_region_by_response_boxplot.pdf")
    plt.close(fig)


def plot_boundary_delta(
    boundary_delta: pd.DataFrame,
    stats: pd.DataFrame,
    outputs_dir: Path,
    config: AnalysisConfig,
):
    if boundary_delta.empty:
        return
    setup_plotting()
    features = [f"score__{signature}" for signature in PLOT_SIGNATURES]
    plot_df = boundary_delta.loc[boundary_delta["feature"].isin(features)].copy()
    if plot_df.empty:
        return
    plot_df["signature"] = plot_df["feature"].map(clean_feature_label)

    fig, ax = plt.subplots(figsize=(6.8, 2.4), constrained_layout=True)
    ax.axhline(0, color="#999999", linewidth=0.5, linestyle="--")
    sns.boxplot(
        data=plot_df,
        x="signature",
        y="margin_minus_inner",
        hue=config.response_col,
        order=PLOT_SIGNATURES,
        hue_order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        width=0.7,
        fliersize=0,
        linewidth=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="signature",
        y="margin_minus_inner",
        hue=config.response_col,
        order=PLOT_SIGNATURES,
        hue_order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        dodge=True,
        size=2.8,
        alpha=0.9,
        linewidth=0.2,
        edgecolor="white",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Margin - inner\nsignature score")
    ax.tick_params(axis="x", rotation=25)
    y_values = plot_df["margin_minus_inner"].dropna()
    if not y_values.empty:
        y_min, y_max = float(y_values.min()), float(y_values.max())
        y_span = max(y_max - y_min, 0.1)
        ax.set_ylim(y_min - 0.12 * y_span, y_max + 0.30 * y_span)
        for i, signature in enumerate(PLOT_SIGNATURES):
            feature = f"score__{signature}"
            p = lookup_p_value(stats, feature, "boundary_enrichment_by_response")
            ax.text(
                i,
                y_max + 0.08 * y_span,
                p_value_label(p),
                ha="center",
                va="bottom",
                fontsize=5.5,
            )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[: len(RESPONSE_ORDER)],
        labels[: len(RESPONSE_ORDER)],
        frameon=False,
        title="",
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
    )
    fig.savefig(outputs_dir / "figure_6e_signature_boundary_enrichment_by_response.pdf")
    plt.close(fig)


def plot_gene_region_heatmap(
    sample_region: pd.DataFrame,
    outputs_dir: Path,
    config: AnalysisConfig,
):
    setup_plotting()
    plot_rows = []
    for response in RESPONSE_ORDER:
        for region in ["tumor_inner", "tumor_margin"]:
            subset = sample_region.loc[
                sample_region[config.response_col].eq(response)
                & sample_region[config.region_col].eq(region)
            ]
            label = f"{response}\n{REGION_LABELS[region]}"
            row = {"group": label}
            for gene in GENE_PANEL:
                col = f"gene__{gene}__mean"
                if col in subset.columns:
                    row[gene] = subset[col].mean()
                else:
                    row[gene] = np.nan
            plot_rows.append(row)
    heat = pd.DataFrame(plot_rows).set_index("group")
    heat = heat.dropna(axis=1, how="all")
    if heat.empty:
        return

    centered = (heat - heat.mean(axis=0)) / heat.std(axis=0, ddof=0).replace(0, np.nan)
    centered = centered.fillna(0)

    fig, ax = plt.subplots(figsize=(7.1, 2.0), constrained_layout=True)
    sns.heatmap(
        centered,
        cmap="vlag",
        center=0,
        linewidths=0.25,
        linecolor="#D8D8D8",
        cbar_kws={"label": "Column z-score"},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)
    fig.savefig(outputs_dir / "figure_6f_gene_region_response_heatmap.pdf")
    plt.close(fig)


def feature_kind(feature: str) -> str:
    if feature.startswith("score__"):
        return "signature"
    if feature.startswith("gene__"):
        return "gene"
    return "other"


def sample_region_to_long(sample_region: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    id_cols = [
        col
        for col in [config.sample_col, config.patient_col, config.response_col, config.region_col, "n_cells"]
        if col in sample_region.columns
    ]
    feature_cols = [
        col
        for col in sample_region.columns
        if col.startswith(("score__", "gene__")) and col.endswith(("__mean", "__median"))
    ]
    if not feature_cols:
        return pd.DataFrame()
    long = sample_region.melt(
        id_vars=id_cols,
        value_vars=feature_cols,
        var_name="feature_statistic",
        value_name="value",
    )
    parsed = long["feature_statistic"].str.extract(r"^(score|gene)__(.+)__(mean|median)$")
    long["feature_type"] = parsed[0].map({"score": "signature", "gene": "gene"})
    long["feature"] = parsed[0] + "__" + parsed[1]
    long["feature_name"] = parsed[1]
    long["statistic"] = parsed[2]
    long = long.drop(columns=["feature_statistic"])
    return long[
        id_cols + ["feature_type", "feature", "feature_name", "statistic", "value"]
    ]


def boundary_delta_to_long(boundary_delta: pd.DataFrame, config: AnalysisConfig) -> pd.DataFrame:
    if boundary_delta.empty:
        return boundary_delta.copy()
    long = boundary_delta.copy()
    long["feature_type"] = long["feature"].map(feature_kind)
    long["feature_name"] = long["feature"].map(clean_feature_label)
    id_cols = [
        col
        for col in [config.sample_col, config.patient_col, config.response_col]
        if col in long.columns
    ]
    value_cols = ["margin_mean", "inner_mean", "margin_minus_inner", "margin_over_inner"]
    return long[id_cols + ["feature_type", "feature", "feature_name"] + value_cols]


def write_outputs(
    cell_scores: pd.DataFrame,
    sample_region_all: pd.DataFrame,
    sample_region: pd.DataFrame,
    boundary_delta: pd.DataFrame,
    stats: pd.DataFrame,
    missing_table: pd.DataFrame,
    results_dir: Path,
    config: AnalysisConfig,
    write_cell_scores: bool = False,
):
    if write_cell_scores:
        cell_scores.to_csv(results_dir / "pretreatment_malignant_cell_signature_scores.csv", index=False)
    sample_region_all.to_csv(
        results_dir / "pretreatment_malignant_sample_region_summary_all_regions.csv",
        index=False,
    )
    sample_region.to_csv(results_dir / "pretreatment_malignant_sample_region_summary.csv", index=False)
    boundary_delta.to_csv(results_dir / "pretreatment_malignant_boundary_delta_by_sample.csv", index=False)
    stats.to_csv(results_dir / "pretreatment_malignant_spatial_signature_statistics.csv", index=False)
    missing_table.to_csv(results_dir / "pretreatment_malignant_signature_missing_genes.csv", index=False)

    sample_region_long = sample_region_to_long(sample_region, config)
    sample_region_long.to_csv(
        results_dir / "pretreatment_malignant_sample_region_summary_long.csv",
        index=False,
    )
    sample_region_long.loc[sample_region_long["feature_type"].eq("signature")].to_csv(
        results_dir / "pretreatment_malignant_signature_region_summary_long.csv",
        index=False,
    )
    sample_region_long.loc[sample_region_long["feature_type"].eq("gene")].to_csv(
        results_dir / "pretreatment_malignant_gene_region_summary_long.csv",
        index=False,
    )

    boundary_delta_long = boundary_delta_to_long(boundary_delta, config)
    boundary_delta_long.to_csv(
        results_dir / "pretreatment_malignant_boundary_delta_by_sample_long.csv",
        index=False,
    )
    if not boundary_delta_long.empty:
        boundary_delta_long.loc[boundary_delta_long["feature_type"].eq("signature")].to_csv(
            results_dir / "pretreatment_malignant_signature_boundary_delta_long.csv",
            index=False,
        )
        boundary_delta_long.loc[boundary_delta_long["feature_type"].eq("gene")].to_csv(
            results_dir / "pretreatment_malignant_gene_boundary_delta_long.csv",
            index=False,
        )


def run_analysis(
    adata,
    output_root: str | Path,
    layer: str | None = None,
    log1p_transform: bool = False,
    min_cells_per_region: int = 20,
    project_root: str | Path | None = None,
    make_plots: bool = False,
    write_cell_scores: bool = False,
):
    output_root = Path(output_root)
    outputs_dir = output_root / "outputs"
    results_dir = output_root / "results"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    if project_root is not None:
        project_root = Path(project_root)

    config = AnalysisConfig(
        layer=layer,
        log1p_transform=log1p_transform,
        min_cells_per_region=min_cells_per_region,
    )

    required_cols = [
        config.region_col,
        config.response_col,
        config.sample_col,
    ]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise KeyError(f"Missing required adata.obs columns: {missing_cols}")

    adata_epi = subset_pretreatment_malignant(adata, config)
    if adata_epi.n_obs == 0:
        raise ValueError("No pretreatment malignant cells remained after subsetting.")

    cell_scores, missing_table = compute_cell_scores(adata_epi, config, project_root)
    sample_region_all = summarize_by_sample_region(cell_scores, config)
    sample_region_all = filter_min_cells(sample_region_all, config)
    sample_region = comparison_regions_only(sample_region_all, config)
    boundary_delta = boundary_delta_table(sample_region, config)
    stats = compute_statistics(sample_region, boundary_delta, config)

    write_outputs(
        cell_scores=cell_scores,
        sample_region_all=sample_region_all,
        sample_region=sample_region,
        boundary_delta=boundary_delta,
        stats=stats,
        missing_table=missing_table,
        results_dir=results_dir,
        config=config,
        write_cell_scores=write_cell_scores,
    )
    if make_plots:
        plot_signature_region_response(sample_region, stats, outputs_dir, config)
        plot_boundary_delta(boundary_delta, stats, outputs_dir, config)
        plot_gene_region_heatmap(sample_region, outputs_dir, config)

    return {
        "sample_region": sample_region,
        "sample_region_all": sample_region_all,
        "sample_region_long": sample_region_to_long(sample_region, config),
        "boundary_delta": boundary_delta,
        "boundary_delta_long": boundary_delta_to_long(boundary_delta, config),
        "statistics": stats,
        "missing_genes": missing_table,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze pretreatment malignant-cell MP7 signatures across tumor margin and inner regions."
    )
    parser.add_argument("--input-h5ad", required=True, help="Input AnnData .h5ad file.")
    parser.add_argument(
        "--output-root",
        default=str(Path(__file__).resolve().parents[1]),
        help="Output folder with results/ and outputs/ subfolders.",
    )
    parser.add_argument(
        "--project-root",
        default=str(project_root_from_script()),
        help="Project root used to find existing MP7 gene tables.",
    )
    parser.add_argument("--layer", default=None, help="AnnData layer to use. Default uses adata.X.")
    parser.add_argument(
        "--log1p-transform",
        action="store_true",
        help="Apply log1p to extracted expression values before scoring.",
    )
    parser.add_argument(
        "--min-cells-per-region",
        type=int,
        default=20,
        help="Minimum malignant cells required per sample-region summary.",
    )
    parser.add_argument(
        "--make-plots",
        action="store_true",
        help="Also render exploratory PDFs on the server. Default only writes CSV summaries.",
    )
    parser.add_argument(
        "--write-cell-scores",
        action="store_true",
        help="Write the large cell-level score table. Default writes compact sample-level CSVs only.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("Please install anndata or run this inside the existing Scanpy environment.") from exc

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = ad.read_h5ad(args.input_h5ad)

    run_analysis(
        adata,
        output_root=args.output_root,
        layer=args.layer,
        log1p_transform=args.log1p_transform,
        min_cells_per_region=args.min_cells_per_region,
        project_root=args.project_root,
        make_plots=args.make_plots,
        write_cell_scores=args.write_cell_scores,
    )


if __name__ == "__main__":
    main()
