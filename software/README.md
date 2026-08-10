# ESCC Xenium Supplementary Software

This directory collects reusable research utilities that were previously prepared as separate supplementary software archives. They are now organized as source packages within this project so the analysis code can be shared through GitHub while large data files remain outside the code release.

## Packages

- `packages/xenium_malignant_distance`: annotates malignant cells by distance to tumour boundary or adjacent immune/stromal regions. Used for Figure 6 tumour-margin and boundary analyses.
- `packages/xenium_spatial_neighborhood`: builds spatial neighborhood / niche programs from Xenium cell-type coordinates. Used for Figure 4 niche analysis.
- `packages/xenium_metaprogram`: discovers malignant-cell metaprograms from single-cell expression data. Used for Figure 6 malignant metaprogram analysis.
- `packages/lzd_xenium_utility`: lightweight utilities for Xenium H&E alignment and transferring AnnData annotations.

## Installation For Local Reuse

Install only the package needed for the analysis you are running, for example:

```bash
python -m pip install -e software/packages/xenium_spatial_neighborhood
python -m pip install -e software/packages/xenium_malignant_distance
python -m pip install -e software/packages/xenium_metaprogram
python -m pip install -e software/packages/lzd_xenium_utility
```

These packages are research utilities for reproducing the ESCC Xenium analyses. They do not include the manuscript data; exported CSV tables and submission-ready data are kept in the project `data/` area and should be deposited separately, for example through GEO or another journal-approved repository.

## Repository Hygiene

The integrated copies intentionally omit macOS metadata, Python bytecode caches, and generated output folders from the original supplementary archives. Figure-specific wrapper scripts remain under `final/figure_*_and_extended_data_figures/.../code/`, while reusable package logic lives here.
