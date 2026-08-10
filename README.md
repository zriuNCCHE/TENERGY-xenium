# ESCC Xenium Figure Analysis

This repository contains code, configuration files, reusable research software, and figure-specific scripts for the ESCC Xenium manuscript analysis.

Large exported data tables and raw analysis inputs are intentionally not tracked in Git. Submission-ready data should be deposited separately, for example through GEO or another journal-approved repository.

## Structure

- `config/`: shared plotting parameters, color palettes, and spatial niche names.
- `final/`: figure-specific code and documentation for final manuscript figures and extended data figures.
- `exploratory/`: exploratory analysis code retained for internal scientific review.
- `software/`: reusable Python packages that were previously prepared as supplementary software archives.
- `design_principles.md`: project-level plotting and organization principles.

## Reusable Software

The reusable packages are integrated under `software/packages/` and can be installed independently:

```bash
python -m pip install -e software/packages/xenium_spatial_neighborhood
python -m pip install -e software/packages/xenium_malignant_distance
python -m pip install -e software/packages/xenium_metaprogram
python -m pip install -e software/packages/lzd_xenium_utility
```

Figure wrapper scripts remain in the corresponding `final/figure_*_and_extended_data_figures/.../code/` folders.

## GitHub Scope

This repository is intended for sharing code and lightweight documentation. Generated outputs, result tables, and data files are excluded from version control to keep the repository suitable for GitHub.
