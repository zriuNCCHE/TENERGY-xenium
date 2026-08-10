#!/usr/bin/env python3
"""
Create Figure 6 pretreatment spatial-signature plots from server-exported CSVs.

This local plotting script expects compact sample-level outputs from
analyze_pretreatment_spatial_signature_boundary.py. It recomputes p-values from
sample-level summaries so plot annotations match the rendered data.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import os
import re
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
warnings.filterwarnings("ignore", category=FutureWarning)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon


RESPONSE_ORDER = ["cCR", "non-cCR"]
RESPONSE_COLORS = {"cCR": "#7BB6A4", "non-cCR": "#F2A38A"}
REGION_ORDER = ["tumor_inner", "tumor_margin"]
REGION_LABELS = {"tumor_inner": "Inner", "tumor_margin": "Margin"}

SIGNATURE_ORDER = [
    "MP7 all",
    "Laminin-332",
    "Adhesion / anchoring",
    "Remodeling / invasion",
    "Stress / adaptation",
]

CORE_GENE_ORDER = [
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


def p_value(x, y) -> float:
    x = pd.Series(x).dropna().astype(float)
    y = pd.Series(y).dropna().astype(float)
    if len(x) < 2 or len(y) < 2:
        return np.nan
    try:
        return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def paired_p_value(delta) -> float:
    delta = pd.Series(delta).dropna().astype(float)
    if len(delta) < 2 or np.allclose(delta, 0):
        return np.nan
    try:
        return float(wilcoxon(delta, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def p_label(p: float) -> str:
    if not np.isfinite(p):
        return "p = NA"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def feature_label(feature: str) -> str:
    return feature.replace("score__", "").replace("gene__", "")


def safe_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return text.lower()


def sample_region_to_long(sample_region: pd.DataFrame) -> pd.DataFrame:
    id_cols = ["sampleID", "patientID", "cCR", "malignant_region", "n_cells"]
    id_cols = [col for col in id_cols if col in sample_region.columns]
    value_cols = [
        col
        for col in sample_region.columns
        if col.startswith(("score__", "gene__")) and col.endswith(("__mean", "__median"))
    ]
    long = sample_region.melt(
        id_vars=id_cols,
        value_vars=value_cols,
        var_name="feature_statistic",
        value_name="value",
    )
    parsed = long["feature_statistic"].str.extract(r"^(score|gene)__(.+)__(mean|median)$")
    long["feature_type"] = parsed[0].map({"score": "signature", "gene": "gene"})
    long["feature"] = parsed[0] + "__" + parsed[1]
    long["feature_name"] = parsed[1]
    long["statistic"] = parsed[2]
    long = long.drop(columns=["feature_statistic"])
    long["region_label"] = long["malignant_region"].map(REGION_LABELS)
    return long


def boundary_to_long(boundary_delta: pd.DataFrame) -> pd.DataFrame:
    out = boundary_delta.copy()
    out["feature_type"] = np.where(out["feature"].str.startswith("score__"), "signature", "gene")
    out["feature_name"] = out["feature"].map(feature_label)
    return out


def order_existing(items, available):
    available = list(dict.fromkeys(available))
    ordered = [item for item in items if item in available]
    ordered.extend([item for item in available if item not in ordered])
    return ordered


def region_response_p_table(long_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    df = long_df.loc[long_df["statistic"].eq("mean")].copy()
    for (feature_type, feature, feature_name, region), sub in df.groupby(
        ["feature_type", "feature", "feature_name", "malignant_region"], observed=True
    ):
        ccr = sub.loc[sub["cCR"].eq("cCR"), "value"]
        nonccr = sub.loc[sub["cCR"].eq("non-cCR"), "value"]
        rows.append(
            {
                "feature_type": feature_type,
                "feature": feature,
                "feature_name": feature_name,
                "test": "region_nonccr_vs_ccr",
                "region": region,
                "comparison": f"non-cCR vs cCR in {region}",
                "n_nonccr": int(nonccr.dropna().shape[0]),
                "n_ccr": int(ccr.dropna().shape[0]),
                "mean_nonccr": float(nonccr.mean()) if len(nonccr.dropna()) else np.nan,
                "mean_ccr": float(ccr.mean()) if len(ccr.dropna()) else np.nan,
                "effect_nonccr_minus_ccr": float(nonccr.mean() - ccr.mean())
                if len(nonccr.dropna()) and len(ccr.dropna())
                else np.nan,
                "p_value": p_value(nonccr, ccr),
            }
        )
    return pd.DataFrame(rows)


def boundary_p_table(boundary_long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (feature_type, feature, feature_name), sub in boundary_long.groupby(
        ["feature_type", "feature", "feature_name"], observed=True
    ):
        ccr = sub.loc[sub["cCR"].eq("cCR"), "margin_minus_inner"]
        nonccr = sub.loc[sub["cCR"].eq("non-cCR"), "margin_minus_inner"]
        all_delta = sub["margin_minus_inner"]
        rows.append(
            {
                "feature_type": feature_type,
                "feature": feature,
                "feature_name": feature_name,
                "test": "boundary_enrichment_by_response",
                "comparison": "non-cCR delta vs cCR delta",
                "n_nonccr": int(nonccr.dropna().shape[0]),
                "n_ccr": int(ccr.dropna().shape[0]),
                "mean_nonccr": float(nonccr.mean()) if len(nonccr.dropna()) else np.nan,
                "mean_ccr": float(ccr.mean()) if len(ccr.dropna()) else np.nan,
                "effect_nonccr_minus_ccr": float(nonccr.mean() - ccr.mean())
                if len(nonccr.dropna()) and len(ccr.dropna())
                else np.nan,
                "p_value": p_value(nonccr, ccr),
                "paired_margin_minus_inner_p_value": paired_p_value(all_delta),
            }
        )
    return pd.DataFrame(rows)


def add_region_p_labels(ax, data, feature, region_pvals, y_col="value"):
    y = data[y_col].dropna()
    if y.empty:
        return
    y_min, y_max = float(y.min()), float(y.max())
    y_span = max(y_max - y_min, 0.1)
    ax.set_ylim(y_min - 0.12 * y_span, y_max + 0.30 * y_span)
    for i, region in enumerate(REGION_ORDER):
        if region not in data["malignant_region"].unique():
            continue
        p = region_pvals.loc[
            region_pvals["feature"].eq(feature) & region_pvals["region"].eq(region), "p_value"
        ]
        p = float(p.iloc[0]) if not p.empty else np.nan
        ax.text(i, y_max + 0.08 * y_span, p_label(p), ha="center", va="bottom", fontsize=5.5)


def add_boundary_p_label(ax, data, feature, boundary_pvals):
    y = data["margin_minus_inner"].dropna()
    if y.empty:
        return
    y_min, y_max = float(y.min()), float(y.max())
    y_span = max(y_max - y_min, 0.1)
    ax.set_ylim(y_min - 0.12 * y_span, y_max + 0.28 * y_span)
    p = boundary_pvals.loc[boundary_pvals["feature"].eq(feature), "p_value"]
    p = float(p.iloc[0]) if not p.empty else np.nan
    ax.text(0.5, y_max + 0.08 * y_span, p_label(p), ha="center", va="bottom", fontsize=5.5)


def draw_region_box(ax, data, feature, region_pvals, title=None, ylabel="Mean expression"):
    sub = data.loc[data["feature"].eq(feature) & data["statistic"].eq("mean")].copy()
    sub = sub.loc[sub["malignant_region"].isin(REGION_ORDER)]
    sns.boxplot(
        data=sub,
        x="malignant_region",
        y="value",
        hue="cCR",
        order=REGION_ORDER,
        hue_order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        width=0.7,
        fliersize=0,
        linewidth=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=sub,
        x="malignant_region",
        y="value",
        hue="cCR",
        order=REGION_ORDER,
        hue_order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        dodge=True,
        size=2.4,
        alpha=0.9,
        linewidth=0.2,
        edgecolor="white",
        ax=ax,
    )
    ax.set_xticks(range(len(REGION_ORDER)))
    ax.set_xticklabels([REGION_LABELS[x] for x in REGION_ORDER], rotation=25)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title or feature_label(feature))
    add_region_p_labels(ax, sub, feature, region_pvals)
    if ax.legend_:
        ax.legend_.remove()


def draw_boundary_box(ax, data, feature, boundary_pvals, title=None, ylabel="Margin - inner"):
    sub = data.loc[data["feature"].eq(feature)].copy()
    ax.axhline(0, color="#999999", linewidth=0.5, linestyle="--")
    sns.boxplot(
        data=sub,
        x="cCR",
        y="margin_minus_inner",
        order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        width=0.55,
        fliersize=0,
        linewidth=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=sub,
        x="cCR",
        y="margin_minus_inner",
        order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        size=2.6,
        alpha=0.9,
        linewidth=0.2,
        edgecolor="white",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_title(title or feature_label(feature))
    add_boundary_p_label(ax, sub, feature, boundary_pvals)


def plot_signature_region(long_df, region_pvals, output_dir):
    setup_plotting()
    features = [f"score__{x}" for x in SIGNATURE_ORDER]
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.6), constrained_layout=True)
    axes = axes.ravel()
    for ax, feature in zip(axes, features):
        draw_region_box(
            ax,
            long_df,
            feature,
            region_pvals,
            title=feature_label(feature),
            ylabel="Mean z-scored signature",
        )
    handles, labels = axes[0].get_legend_handles_labels()
    for ax in axes[len(features) :]:
        ax.axis("off")
        ax.legend(handles[:2], labels[:2], loc="center", ncol=1, frameon=False, title="")
    fig.savefig(output_dir / "pretreatment_signature_region_response_boxplots.pdf")
    plt.close(fig)


def plot_signature_boundary(boundary_long, boundary_pvals, output_dir):
    setup_plotting()
    features = [f"score__{x}" for x in SIGNATURE_ORDER]
    plot_df = boundary_long.loc[boundary_long["feature"].isin(features)].copy()
    plot_df["feature_name"] = pd.Categorical(
        plot_df["feature_name"], categories=SIGNATURE_ORDER, ordered=True
    )
    fig, ax = plt.subplots(figsize=(6.9, 2.6), constrained_layout=True)
    ax.axhline(0, color="#999999", linewidth=0.5, linestyle="--")
    sns.boxplot(
        data=plot_df,
        x="feature_name",
        y="margin_minus_inner",
        hue="cCR",
        order=SIGNATURE_ORDER,
        hue_order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        width=0.7,
        fliersize=0,
        linewidth=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="feature_name",
        y="margin_minus_inner",
        hue="cCR",
        order=SIGNATURE_ORDER,
        hue_order=RESPONSE_ORDER,
        palette=RESPONSE_COLORS,
        dodge=True,
        size=2.6,
        alpha=0.9,
        linewidth=0.2,
        edgecolor="white",
        ax=ax,
    )
    y = plot_df["margin_minus_inner"].dropna()
    y_min, y_max = float(y.min()), float(y.max())
    y_span = max(y_max - y_min, 0.1)
    ax.set_ylim(y_min - 0.12 * y_span, y_max + 0.30 * y_span)
    for i, feature in enumerate(features):
        p = boundary_pvals.loc[boundary_pvals["feature"].eq(feature), "p_value"]
        p = float(p.iloc[0]) if not p.empty else np.nan
        ax.text(i, y_max + 0.08 * y_span, p_label(p), ha="center", va="bottom", fontsize=5.5)
    ax.set_xlabel("")
    ax.set_ylabel("Margin - inner\nsignature score")
    ax.tick_params(axis="x", rotation=25)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, title="", loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(output_dir / "pretreatment_signature_boundary_delta_boxplots.pdf")
    plt.close(fig)


def plot_gene_region_grid(long_df, region_pvals, genes, output_path, n_cols=4):
    setup_plotting()
    features = [f"gene__{gene}" for gene in genes]
    n_rows = int(np.ceil(len(features) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.1, max(2.0, 1.9 * n_rows)), constrained_layout=True)
    axes = np.ravel(axes)
    for ax, feature in zip(axes, features):
        draw_region_box(ax, long_df, feature, region_pvals, title=feature_label(feature), ylabel="Mean expression")
    for ax in axes[len(features) :]:
        ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles[:2],
        labels[:2],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_gene_boundary_grid(boundary_long, boundary_pvals, genes, output_path, n_cols=4):
    setup_plotting()
    features = [f"gene__{gene}" for gene in genes]
    n_rows = int(np.ceil(len(features) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.1, max(2.0, 1.8 * n_rows)), constrained_layout=True)
    axes = np.ravel(axes)
    for ax, feature in zip(axes, features):
        draw_boundary_box(ax, boundary_long, feature, boundary_pvals, title=feature_label(feature))
    for ax in axes[len(features) :]:
        ax.axis("off")
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=RESPONSE_COLORS["cCR"], label="cCR"),
        plt.Line2D([0], [0], marker="s", linestyle="", color=RESPONSE_COLORS["non-cCR"], label="non-cCR"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.01),
        ncol=2,
        frameon=False,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_individual_gene_files(long_df, boundary_long, region_pvals, boundary_pvals, genes, output_dir):
    region_dir = output_dir / "individual_gene_region_boxplots"
    boundary_dir = output_dir / "individual_gene_boundary_delta_boxplots"
    region_dir.mkdir(parents=True, exist_ok=True)
    boundary_dir.mkdir(parents=True, exist_ok=True)
    setup_plotting()
    for gene in genes:
        feature = f"gene__{gene}"
        fig, ax = plt.subplots(figsize=(2.2, 2.4), constrained_layout=True)
        draw_region_box(ax, long_df, feature, region_pvals, title=gene, ylabel="Mean expression")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], frameon=False, title="", loc="upper left", bbox_to_anchor=(1.02, 1))
        fig.savefig(region_dir / f"{gene}_region_response_boxplot.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(2.0, 2.4), constrained_layout=True)
        draw_boundary_box(ax, boundary_long, feature, boundary_pvals, title=gene)
        fig.savefig(boundary_dir / f"{gene}_boundary_delta_boxplot.pdf")
        plt.close(fig)


def plot_gene_matrix(sample_region, genes, output_dir):
    setup_plotting()
    rows = []
    for response in RESPONSE_ORDER:
        for region in REGION_ORDER:
            sub = sample_region.loc[
                sample_region["cCR"].eq(response) & sample_region["malignant_region"].eq(region)
            ]
            row = {"group": f"{response}\n{REGION_LABELS[region]}"}
            for gene in genes:
                col = f"gene__{gene}__mean"
                row[gene] = sub[col].mean() if col in sub.columns else np.nan
            rows.append(row)
    heat = pd.DataFrame(rows).set_index("group")
    heat = heat.dropna(axis=1, how="all")
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
    fig.savefig(output_dir / "pretreatment_gene_region_response_matrix.pdf")
    plt.close(fig)


def parse_args():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Plot pretreatment spatial signature/gene boxplots from sample-level CSVs."
    )
    parser.add_argument(
        "--results-dir",
        default=str(root / "results"),
        help="Folder containing pretreatment malignant CSV summaries.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "outputs" / "csv_derived_plots"),
        help="Folder where PDF plots should be written.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_region = pd.read_csv(results_dir / "pretreatment_malignant_sample_region_summary.csv")
    boundary_delta = pd.read_csv(results_dir / "pretreatment_malignant_boundary_delta_by_sample.csv")

    sample_region = sample_region.loc[sample_region["malignant_region"].isin(REGION_ORDER)].copy()
    sample_region["cCR"] = pd.Categorical(sample_region["cCR"], categories=RESPONSE_ORDER, ordered=True)
    boundary_delta["cCR"] = pd.Categorical(boundary_delta["cCR"], categories=RESPONSE_ORDER, ordered=True)

    long_df = sample_region_to_long(sample_region)
    boundary_long = boundary_to_long(boundary_delta)

    all_genes = sorted(long_df.loc[long_df["feature_type"].eq("gene"), "feature_name"].dropna().unique())
    gene_order = order_existing(CORE_GENE_ORDER, all_genes)
    core_genes = [gene for gene in CORE_GENE_ORDER if gene in all_genes]

    region_pvals = region_response_p_table(long_df)
    boundary_pvals = boundary_p_table(boundary_long)
    region_pvals.to_csv(results_dir / "pretreatment_spatial_region_response_pvalues_from_csv.csv", index=False)
    boundary_pvals.to_csv(results_dir / "pretreatment_spatial_boundary_delta_pvalues_from_csv.csv", index=False)
    long_df.to_csv(results_dir / "pretreatment_malignant_sample_region_summary_long_from_csv.csv", index=False)
    boundary_long.to_csv(results_dir / "pretreatment_malignant_boundary_delta_by_sample_long_from_csv.csv", index=False)

    plot_signature_region(long_df, region_pvals, output_dir)
    plot_signature_boundary(boundary_long, boundary_pvals, output_dir)
    plot_gene_region_grid(
        long_df,
        region_pvals,
        gene_order,
        output_dir / "pretreatment_gene_region_response_boxplots_all_genes.pdf",
        n_cols=4,
    )
    plot_gene_region_grid(
        long_df,
        region_pvals,
        core_genes,
        output_dir / "pretreatment_gene_region_response_boxplots_core_genes.pdf",
        n_cols=3,
    )
    plot_gene_boundary_grid(
        boundary_long,
        boundary_pvals,
        gene_order,
        output_dir / "pretreatment_gene_boundary_delta_boxplots_all_genes.pdf",
        n_cols=4,
    )
    plot_gene_boundary_grid(
        boundary_long,
        boundary_pvals,
        core_genes,
        output_dir / "pretreatment_gene_boundary_delta_boxplots_core_genes.pdf",
        n_cols=3,
    )
    plot_individual_gene_files(long_df, boundary_long, region_pvals, boundary_pvals, gene_order, output_dir)
    plot_gene_matrix(sample_region, core_genes, output_dir)

    print(f"Wrote plots to: {output_dir}")
    print(f"Genes plotted: {len(gene_order)}")


if __name__ == "__main__":
    main()
