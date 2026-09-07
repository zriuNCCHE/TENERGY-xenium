# Data Upload Plan

## Current inventory

The local `data/` directory is approximately 1.9 GB and is currently excluded from Git by `.gitignore`. It contains large cell-level exports, analysis matrices, clinical follow-up files, and smaller reference tables. None of these files are currently tracked by GitHub.

## Recommended GitHub scope

Keep only lightweight, non-sensitive reproducibility metadata in the code repository:

- `data/geo_ready/metaprograms/xenium_malignant_metaprogram_genes.json`
- `data/geo_ready/metaprograms/xenium_malignant_metaprogram_genes.txt`
- `data/geo_ready/tenergy_rnaseq/TENERGY_TPM_sample_manifest_inferred.csv`
- `data/geo_ready/tenergy_rnaseq/xenium_tenergy_patient_id_crosswalk.csv`
- small configuration and reference tables in `config/`
- `data/_updated_survival/20260904TENERGY初発症例データ一覧.xlsx`
- `data/_updated_survival/20260904EPOC1802長期予後解析_draft1.pptx`

These files should be added deliberately with `git add -f`, because the repository currently ignores `data/` as a whole.

## Recommended GEO or journal data deposit

Deposit the cell-level and expression data separately rather than placing them in the Git repository:

- `combined_adata_obs_with_all_obsm.csv` and the immune, fibroblast, CD4 and CD8 exports
- `combined_final_all_cell_types_obs.csv`
- `postC_k7_cell_niche_assignments.csv`
- malignant-cell MP score tables, including region annotations
- MP2/CD4-CXCL13 and MP7/myeloid distance cell-level tables
- `TENERGY_TPM.tsv`
- marker and differential-expression matrices used as source data
- the updated clinical survival workbook and any supplementary clinical tables, according to the journal's clinical-data policy

The current `_updated_survival` folder is small and has been allowlisted for Git tracking. Temporary Excel lock files are excluded.

The largest files are not appropriate for ordinary GitHub storage: the combined observation export is about 883 MB, the immune export about 297 MB, the fibroblast export about 171 MB, and the spatial neighborhood table about 87 MB. GitHub also imposes practical repository-size limits even when Git LFS is used.

## Reproducibility requirement

The final plotting scripts currently expect the deposited files at paths beginning with `data/geo_ready/` or `data/raw/`. Before public release, add a small manifest containing the data-deposit accession, filename, checksum, and expected local path. The repository should then contain code, configuration, package documentation, and the manifest; the large data should be downloaded from the external repository into the documented paths.

## Files intentionally not recommended for GitHub

- raw or full cell-level CSV exports
- duplicate derived tables that can be regenerated from a primary source
- PPTX or spreadsheet working files used only during analysis
- generated PDFs, PNGs, TIFFs and intermediate result tables
- private notebook outputs and operating-system metadata
