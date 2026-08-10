# QC Overview

This final Supplementary Figure 1 part summarizes Xenium sample scale and QC metrics from `data/raw/combined_final_all_cell_types_obs.csv`.

Panels:

- Total cells per sample.
- Sample-level median transcript counts, detected genes, cell area and nucleus ratio.
- Summed segmented cell area and cell density per segmented cell area.

Note: true tissue/section mask area is not available in the current metadata table. Area and density panels therefore use summed segmented cell area, not whole-section area.

Outputs:

- `outputs/qc_overview_sample_summary.csv`: sample-level QC and area summaries.
- `outputs/qc_overview.pdf`: vector figure.

