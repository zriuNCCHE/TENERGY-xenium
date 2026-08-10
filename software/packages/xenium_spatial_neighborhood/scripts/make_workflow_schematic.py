from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle


OUTPUT_DIR = Path(__file__).resolve().parents[1] / "docs" / "figures"


CELL_COLORS = {
    "Tumor": "#ef9a9a",
    "T cell": "#7aa6d9",
    "CAF": "#f4c46a",
    "Myeloid": "#91c9a4",
    "Endothelial": "#b79adf",
}


BAR_COLORS = ["#355c9a", "#e91e63", "#ffd11a", "#7b4a12", "#00796b", "#f2c879"]


def arrow(ax, start, end, color="#f2c300"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.4,
            color=color,
            shrinkA=3,
            shrinkB=3,
        )
    )


def draw_cell(ax, x, y, label, r=0.012, alpha=0.95, edge="#b67b7b"):
    color = CELL_COLORS[label]
    ax.add_patch(Circle((x, y), r, facecolor=color, edgecolor=edge, linewidth=0.55, alpha=alpha))
    ax.add_patch(Circle((x, y), r * 0.33, facecolor=edge, edgecolor="none", alpha=0.55))


def draw_spatial_map(ax):
    ax.text(0.09, 0.86, "Spatial omics data", ha="center", fontsize=10)
    tissue = Polygon(
        [(0.015, 0.18), (0.17, 0.12), (0.205, 0.74), (0.04, 0.78)],
        closed=True,
        facecolor="#fbf7ef",
        edgecolor="#c9c2b8",
        linewidth=1.0,
    )
    ax.add_patch(tissue)
    rng = np.random.default_rng(4)
    for _ in range(115):
        x = rng.uniform(0.035, 0.175)
        y = rng.uniform(0.19, 0.73)
        if y < 0.12 + 3.0 * (x - 0.015) or y > 0.81 - 0.4 * x:
            continue
        label = rng.choice(list(CELL_COLORS), p=[0.45, 0.18, 0.18, 0.14, 0.05])
        draw_cell(ax, x, y, label, r=0.0058, alpha=0.75)
    ax.add_patch(Rectangle((0.09, 0.43), 0.055, 0.115, facecolor="#f8d84a", edgecolor="none", alpha=0.55))
    ax.add_patch(Rectangle((0.018, 0.16), 0.005, 0.63, facecolor="#9ca3af", edgecolor="none"))


def draw_zoomed_neighbors(ax):
    ax.text(0.31, 0.86, "Cells and their neighbors", ha="center", fontsize=10)
    ax.add_patch(Circle((0.31, 0.47), 0.145, facecolor="#fbf7ef", edgecolor="#222222", linewidth=0.65))
    rng = np.random.default_rng(8)
    for _ in range(62):
        angle = rng.uniform(0, 2 * np.pi)
        rad = 0.138 * np.sqrt(rng.uniform())
        x = 0.31 + rad * np.cos(angle)
        y = 0.47 + rad * np.sin(angle)
        label = rng.choice(list(CELL_COLORS), p=[0.35, 0.2, 0.2, 0.18, 0.07])
        draw_cell(ax, x, y, label, r=rng.uniform(0.006, 0.012), alpha=0.58)
    for center, width, height, angle in [
        ((0.26, 0.55), 0.105, 0.135, -15),
        ((0.35, 0.55), 0.12, 0.13, 18),
        ((0.31, 0.37), 0.16, 0.11, 8),
    ]:
        ax.add_patch(Ellipse(center, width, height, angle=angle, facecolor="none", edgecolor="#376da3", linewidth=0.9))


def stacked_bar(ax, x, y, width, height, fractions):
    left = x
    for frac, color in zip(fractions, BAR_COLORS):
        ax.add_patch(Rectangle((left, y), width * frac, height, facecolor=color, edgecolor="white", linewidth=0.25))
        left += width * frac


def draw_composition_matrix(ax):
    ax.text(0.52, 0.86, "Neighboring cell composition", ha="center", fontsize=10)
    rng = np.random.default_rng(12)
    y0 = 0.66
    for i in range(9):
        fractions = rng.dirichlet([1.5, 1.2, 1.0, 0.9, 0.8, 0.5])
        stacked_bar(ax, 0.455, y0 - i * 0.038, 0.125, 0.017, fractions)
    ax.text(0.515, 0.315, "......", ha="center", fontsize=8, color="#555555")
    for i in range(3):
        fractions = rng.dirichlet([1.5, 1.2, 1.0, 0.9, 0.8, 0.5])
        stacked_bar(ax, 0.455, 0.255 - i * 0.038, 0.125, 0.017, fractions)


def draw_pattern_extraction(ax):
    ax.text(0.70, 0.69, "Pattern extraction", ha="center", fontsize=10)
    ax.text(0.66, 0.37, "NMF", ha="center", fontsize=7.5, rotation=15)
    patterns = [
        [0.40, 0.25, 0.20, 0.10, 0.05, 0.00],
        [0.05, 0.45, 0.25, 0.05, 0.15, 0.05],
        [0.10, 0.05, 0.10, 0.20, 0.45, 0.10],
    ]
    for i, fractions in enumerate(patterns):
        stacked_bar(ax, 0.67, 0.58 - i * 0.085, 0.135, 0.027, fractions)
    for y1, y2 in [(0.66, 0.59), (0.55, 0.505), (0.43, 0.42)]:
        ax.plot([0.58, 0.67], [y1, y2], color="#cfcfcf", linewidth=0.6, linestyle="--")


def draw_neighborhood_outputs(ax):
    ax.text(0.91, 0.86, "Spatial neighborhood\nprograms", ha="center", fontsize=10, linespacing=0.9)
    centers = [(0.86, 0.60), (0.95, 0.60), (0.86, 0.36), (0.95, 0.36)]
    rng = np.random.default_rng(22)
    dominant = ["Tumor", "Myeloid", "CAF", "T cell"]
    for center, dom in zip(centers, dominant):
        ax.add_patch(Circle(center, 0.053, facecolor="#fbf7ef", edgecolor="#222222", linewidth=0.5))
        for _ in range(24):
            angle = rng.uniform(0, 2 * np.pi)
            rad = 0.045 * np.sqrt(rng.uniform())
            x = center[0] + rad * np.cos(angle)
            y = center[1] + rad * np.sin(angle)
            label = dom if rng.uniform() < 0.52 else rng.choice(list(CELL_COLORS))
            draw_cell(ax, x, y, label, r=rng.uniform(0.004, 0.008), alpha=0.72)


def build_schematic():
    fig, ax = plt.subplots(figsize=(12, 3.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    draw_spatial_map(ax)
    draw_zoomed_neighbors(ax)
    draw_composition_matrix(ax)
    draw_pattern_extraction(ax)
    draw_neighborhood_outputs(ax)

    arrow(ax, (0.19, 0.48), (0.235, 0.48))
    arrow(ax, (0.41, 0.48), (0.445, 0.48))
    arrow(ax, (0.585, 0.48), (0.66, 0.48))
    arrow(ax, (0.81, 0.48), (0.845, 0.48))

    return fig


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build_schematic()
    pdf_path = OUTPUT_DIR / "spatial_neighborhood_workflow_schematic.pdf"
    png_path = OUTPUT_DIR / "spatial_neighborhood_workflow_schematic.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(pdf_path)
    print(png_path)


if __name__ == "__main__":
    main()
