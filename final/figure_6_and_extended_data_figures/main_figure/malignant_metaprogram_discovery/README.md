# Figure 6 Malignant Metaprogram Discovery

This folder stores the upstream server-side workflow used to discover malignant-cell metaprograms from pretreatment malignant cells.

## Upstream Python Workflow

`code/run_pretreatment_malignant_metaprogram_discovery.py`:

- loads `combined_all_patients_epi_analyzed_final.h5ad`;
- filters to `sample_timepoint == "pre"`;
- filters to malignant cells using `final_cell_type == "malignant"`;
- runs `xenium_metaprogram.MalignantMetaProgrammer`;
- scans `k_range = (7, 10)`;
- uses `top_n_genes = 20`, `n_mps = 10`, `min_patients = 5`, and `consensus_fraction = 0.20`;
- generates metaprogram plots and output tables in the server working directory.

The resulting malignant metaprogram definitions are used by downstream Figure 6 metaprogram DEG-overlap analyses and Figure 7 malignant-MP spatial proximity analyses.

## Notes

This script is intended to be run on the server environment where the full AnnData object and `xenium_metaprogram` package are available. It is preserved here for reproducibility and future manuscript/GEO reference.
