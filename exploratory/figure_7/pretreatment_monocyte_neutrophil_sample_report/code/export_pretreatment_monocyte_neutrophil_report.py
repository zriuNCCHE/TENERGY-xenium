from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[4]
DATA_PATH = PROJECT_ROOT / "data" / "raw" / "combined_final_all_cell_types_obs.csv"
OUT_DIR = (
    PROJECT_ROOT
    / "exploratory"
    / "figure_7"
    / "pretreatment_monocyte_neutrophil_sample_report"
    / "outputs"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_CELL_TYPES = ["Mono_CDC27", "Mono_SLC2A3", "Neutrophil"]


def pct(x: float) -> str:
    return f"{x:.2f}%"


df = pd.read_csv(
    DATA_PATH,
    usecols=[
        "sampleID",
        "sample_timepoint",
        "patientID",
        "cCR",
        "final_cell_type2",
        "cell_area",
    ],
)

pre = df[df["sample_timepoint"].astype(str).str.lower().eq("pre")].copy()

total_by_sample = (
    pre.groupby(["sampleID", "patientID", "cCR"], dropna=False)
    .agg(
        total_cells=("final_cell_type2", "size"),
        total_cell_area=("cell_area", "sum"),
    )
    .reset_index()
)

target = pre[pre["final_cell_type2"].isin(TARGET_CELL_TYPES)].copy()

counts = (
    target.groupby(["sampleID", "final_cell_type2"], dropna=False)
    .size()
    .unstack(fill_value=0)
    .reindex(columns=TARGET_CELL_TYPES, fill_value=0)
    .reset_index()
)

areas = (
    target.groupby(["sampleID", "final_cell_type2"], dropna=False)["cell_area"]
    .sum()
    .unstack(fill_value=0)
    .reindex(columns=TARGET_CELL_TYPES, fill_value=0)
    .add_suffix("_area")
    .reset_index()
)

report = total_by_sample.merge(counts, on="sampleID", how="left").merge(
    areas, on="sampleID", how="left"
)

for col in TARGET_CELL_TYPES + [f"{c}_area" for c in TARGET_CELL_TYPES]:
    report[col] = report[col].fillna(0)

report["combined_monocyte_neutrophil_cells"] = report[TARGET_CELL_TYPES].sum(axis=1)
report["combined_monocyte_cells"] = report[["Mono_CDC27", "Mono_SLC2A3"]].sum(axis=1)
report["combined_monocyte_neutrophil_area"] = report[
    [f"{c}_area" for c in TARGET_CELL_TYPES]
].sum(axis=1)
report["combined_monocyte_area"] = report[["Mono_CDC27_area", "Mono_SLC2A3_area"]].sum(axis=1)

report["combined_monocyte_neutrophil_fraction_percent"] = (
    report["combined_monocyte_neutrophil_cells"] / report["total_cells"] * 100
)
report["combined_monocyte_fraction_percent"] = (
    report["combined_monocyte_cells"] / report["total_cells"] * 100
)
report["neutrophil_fraction_percent"] = report["Neutrophil"] / report["total_cells"] * 100
report["mono_cdc27_fraction_percent"] = report["Mono_CDC27"] / report["total_cells"] * 100
report["mono_slc2a3_fraction_percent"] = report["Mono_SLC2A3"] / report["total_cells"] * 100
report["combined_monocyte_neutrophil_area_fraction_percent"] = (
    report["combined_monocyte_neutrophil_area"] / report["total_cell_area"] * 100
)

report = report.sort_values(
    ["combined_monocyte_neutrophil_cells", "combined_monocyte_neutrophil_fraction_percent"],
    ascending=[False, False],
)
report.insert(0, "rank_by_combined_count", range(1, len(report) + 1))

count_path = OUT_DIR / "pretreatment_monocyte_neutrophil_sample_report.csv"
report.to_csv(count_path, index=False)

rank_tables = {}
rank_specs = {
    "combined_count": "combined_monocyte_neutrophil_cells",
    "combined_fraction": "combined_monocyte_neutrophil_fraction_percent",
    "neutrophil_count": "Neutrophil",
    "neutrophil_fraction": "neutrophil_fraction_percent",
    "mono_cdc27_count": "Mono_CDC27",
    "mono_slc2a3_count": "Mono_SLC2A3",
    "combined_monocyte_count": "combined_monocyte_cells",
}

for name, col in rank_specs.items():
    ranked = report.drop(columns=[c for c in report.columns if c.startswith("rank_by_")]).sort_values(
        col, ascending=False
    ).copy()
    ranked.insert(0, f"rank_by_{name}", range(1, len(ranked) + 1))
    ranked.to_csv(OUT_DIR / f"pretreatment_monocyte_neutrophil_rank_by_{name}.csv", index=False)
    rank_tables[name] = ranked

top_combined = rank_tables["combined_count"].iloc[0]
top_fraction = rank_tables["combined_fraction"].iloc[0]
top_neut = rank_tables["neutrophil_count"].iloc[0]
top_mono = rank_tables["combined_monocyte_count"].iloc[0]

summary_lines = [
    "# Pretreatment Monocyte/Neutrophil Sample Report",
    "",
    "Cell types included: `Mono_CDC27`, `Mono_SLC2A3`, and `Neutrophil`.",
    "Macrophage subsets were not included in this report.",
    "",
    "## Main Results",
    "",
    (
        f"- Highest absolute combined monocyte/neutrophil count: "
        f"`{top_combined.sampleID}` ({int(top_combined.combined_monocyte_neutrophil_cells)} cells; "
        f"{pct(top_combined.combined_monocyte_neutrophil_fraction_percent)} of all pretreatment cells in that sample)."
    ),
    (
        f"- Highest combined monocyte/neutrophil fraction: "
        f"`{top_fraction.sampleID}` ({pct(top_fraction.combined_monocyte_neutrophil_fraction_percent)}; "
        f"{int(top_fraction.combined_monocyte_neutrophil_cells)} cells)."
    ),
    (
        f"- Highest neutrophil count: "
        f"`{top_neut.sampleID}` ({int(top_neut.Neutrophil)} neutrophils; "
        f"{pct(top_neut.neutrophil_fraction_percent)} of all cells)."
    ),
    (
        f"- Highest monocyte count: "
        f"`{top_mono.sampleID}` ({int(top_mono.combined_monocyte_cells)} monocytes; "
        f"{pct(top_mono.combined_monocyte_fraction_percent)} of all cells)."
    ),
    "",
    "## Exported Files",
    "",
    f"- `{count_path.name}`: full sample-level report.",
    "- `pretreatment_monocyte_neutrophil_rank_by_combined_count.csv`",
    "- `pretreatment_monocyte_neutrophil_rank_by_combined_fraction.csv`",
    "- `pretreatment_monocyte_neutrophil_rank_by_neutrophil_count.csv`",
    "- `pretreatment_monocyte_neutrophil_rank_by_combined_monocyte_count.csv`",
]

(OUT_DIR / "pretreatment_monocyte_neutrophil_sample_report.md").write_text(
    "\n".join(summary_lines) + "\n"
)

print("\n".join(summary_lines[:18]))
print(f"\nSaved report to: {count_path}")
