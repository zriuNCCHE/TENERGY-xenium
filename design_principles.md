# eSCC Xenium Figure Project Design Principles

This document summarizes the intended organization and working principles for the eSCC Xenium figure project. It is meant to be a shared reference for future plotting, exploratory analysis, final figure assembly, and discussion of biological findings.

## Project Goals

- Build publication-quality figures for an eSCC Xenium study, initially targeting Nature Cancer or another Nature-family journal.
- Keep exploratory plotting separate from final manuscript figures.
- Allow final figures to be reorganized later, for example moving panels between Figure 1, Figure 2, supplementary figures, or renamed figure folders.
- Store reusable plotting parameters such as colors, labels, ordering, font sizes, figure widths, and sample metadata in shared configuration files rather than hard-coding them inside individual plotting scripts.
- Keep data tables exported from the server in a central data folder so exploratory and final analyses can reference the same inputs.
- Use table-based local inputs rather than storing AnnData `.h5ad` files in this project.
- Support a workflow where some figures or code are produced on the no-internet Qashiwab server and then pasted back into this project.
- Generate figure legends together with figure outputs so that plots, code, statistics, and written interpretation remain connected.
- Keep a discussion area for notes and chats about biological findings.
- Support both Python and R plotting/analysis code, with R preferred for publication plots when practical.
- Keep version control lightweight. The repository is mainly a scientific figure workspace, not a software package.
- Assume the project can eventually be public; privacy restrictions are not expected to constrain file organization.

## Proposed Folder Structure

The initial structure should be flexible and easy to rename:

```text
figures/
  design_principles.md
  data/
    raw/
    processed/
    external_from_qashiwab/
  config/
    cell_types.yaml
    palettes.yaml
    plotting.yaml
    samples.yaml
  exploratory/
    README.md
    2026-xx-xx_topic_name/
      code/
      outputs/
      legends/
      notes.md
  final/
    figure_1/
      part_name/
        code/
        outputs/
        legends/
        stats/
        qashiwab_imports/
        notes.md
    figure_2/
      malignant_cell_proportion/
        code/
        outputs/
        legends/
    figure_3/
    figure_4/
    extended_figure_1/
    extended_figure_2/
    supplementary_figure_1/
    supplementary_figure_2/
  discussion/
    figure_1.md
    figure_2.md
    xenium_biology_notes.md
  scripts/
    io/
    plotting/
    stats/
    python/
    r/
  templates/
    figure_legend_template.md
    figure_folder_template/
```

## Folder Responsibilities

- `data/raw/`: CSV, TXT, TSV, or other exported tables copied from the server or generated upstream. These files should be treated as read-only once added.
- `data/processed/`: cleaned or merged tables generated from raw data. Scripts should document how these are made.
- `data/external_from_qashiwab/`: files produced on Qashiwab, especially when plots require AnnData or other server-only resources.
- `config/`: shared metadata and parameters used by all plotting scripts. Cell-type colors, group ordering, sample labels, cohort definitions, and journal plotting defaults should live here.
- `exploratory/`: temporary or iterative analysis. Exploratory outputs may be rough, but should still use the same data and config system as final figures.
- `final/`: manuscript-facing figures. Each major figure can contain multiple part folders. Each part folder should keep the same internal hierarchy: `code/`, `outputs/`, `legends/`, and optional `stats/`, `qashiwab_imports/`, and `notes.md`.
- `discussion/`: biological interpretation, hypotheses, unresolved questions, and summary notes. This should be separated from plotting code so ideas remain easy to read.
- `scripts/`: reusable parsing, plotting, export, and statistics helpers, including both Python and R code when useful.
- `templates/`: reusable structure for new figures, legends, and notes.

## Figure Naming and Renaming

Figure folders can begin with provisional names such as `figure_1`, `figure_2`, or `supplementary_figure_1`. To make future renaming safe:

- Each final figure folder should include a short `notes.md` file describing the current figure concept and panel list.
- Related parts can be grouped under a larger figure folder. For example, malignant-cell proportion and malignant-cell state plots can both live under `final/figure_2/`.
- Extended figure folders use the same structure as main figure folders, for example `final/extended_figure_2/a2ml1_positive_cell_proportion/`.
- Panel files should use panel-aware names, for example `fig1a_cell_type_map.pdf`, `fig1b_region_composition.pdf`, or `panel_a_cell_type_map.pdf`.
- When figure numbering is uncertain, prefer semantic names inside the files, such as `panel_a_tumor_microenvironment_map.pdf`, and let the folder name carry the current figure number.
- Avoid hard-coded paths inside scripts. Scripts should resolve project paths from a small shared path utility or configuration file.

## Shared Configuration Principles

Plotting choices should be centralized so all figures remain consistent.

Recommended shared configuration:

- `config/cell_types.yaml`: canonical cell-type names, display labels, order, parent categories, and colors.
- `config/palettes.yaml`: color palettes for cell types, tissue regions, samples, patients, treatment groups, and continuous gradients.
- `config/plotting.yaml`: journal style defaults such as font family, font size, DPI, figure widths, file formats, and panel label style.
- `config/samples.yaml`: sample IDs, display names, patient IDs, experimental groups, tissue regions, and any batch or cohort metadata.

YAML is a good default for human-readable configuration. JSON is also acceptable when it is easier for a script or tool to consume. The important rule is that shared parameters should live in `config/`, not inside individual figure scripts.

Suggested rule: individual figure scripts may choose which cell types or groups to display, but should not define new colors or labels locally unless there is a documented reason.

Fixed shared encodings:

- Timepoints must use the colors in `config/palettes.yaml`: `pre`, `postC`, and `postA` should not be recolored locally.
- Clinical outcome groups must use the colors in `config/palettes.yaml`: `cCR` and `non-cCR` should not be recolored locally.
- Figure 2 major-cell-type UMAP and H&E/spatial overlays must use the colors in `config/palettes.yaml` under `figure2_major_cell_types`.
- If a new repeated biological category appears, add its colors or ordering to `config/` before using it in multiple figures.

Reusable figure sizing:

- Store reusable part dimensions in `config/plotting.yaml`.
- Size final parts at the physical dimensions they should occupy in the assembled figure instead of making oversized canvases that must be scaled down later.
- Keep font sizes fixed at journal size when changing the canvas. Do not compensate for layout by shrinking text below the 5-7 pt standard.

## Nature-Family Figure Requirements To Track

Based on the supplied requirements:

- Figures should be cited in manuscript order as `Fig. 1`, `Fig. 2`, and so on.
- Figure panels should be prepared at minimum 300 dpi and maximum width 180 mm.
- Standard text labels should use 5-7 pt sans-serif font at final size.
- Standard plot lines, box outlines, bracket lines, axes and panel borders should use 0.5 pt strokes unless a specific figure type requires a different width.
- Greek characters should use Symbol font where required.
- Use scale bars rather than magnification factors.
- Include error bars where appropriate.
- Keep labels, scale bars, and error bars editable when possible. Avoid flattening them into low-resolution raster images.
- Legends should include a brief title and a short description of panels in sequence.
- Legends should avoid excess methodological detail.
- Legends should describe keys verbally, for example "open red triangles" rather than relying only on symbols.
- Legends should describe center values, error bars, sample size, statistical test, and P values where applicable.

Practical export preference:

- Save final vector plots as PDF only by default, to keep the project small and manuscript-ready.
- Save microscopy or spatial raster images as high-resolution TIFF/PNG.
- Save dense UMAP or spatial point overlays as high-resolution raster PNG/TIFF rather than vector PDF, so individual cells are not stored as editable objects.
- Create PNG/SVG preview files only when specifically needed for review, debugging, or downstream assembly.

## Figure Legend Workflow

Each final figure folder should contain a `legends/` folder with:

- `legend_draft.md`: editable working legend.
- `legend_final.md`: polished version intended for manuscript transfer.
- `panel_descriptions.md`: optional panel-by-panel notes, useful before the final legend is polished.
- `stats_summary.md` or `stats_summary.csv`: sample sizes, statistical tests, P values, and error-bar definitions used in the legend.

Each exploratory folder may also contain `legends/` or `notes.md`, but exploratory legends can be lighter and focused on interpretation.

Legends should be tight. For each figure part, keep the legend draft in the same part folder hierarchy under `legends/`.

Recommended wording style:

- Brief title.
- Panel-wise sentences, following the style of prior manuscript legends: `(A) ...`, `(B and C) ...`.
- Marker/line/box definitions.
- Center value, error definition, sample size, test, and P value when applicable.
- Avoid repeating sample sizes in the legend if they are already stated in the manuscript text or statistics table, unless required for submission.

## Qashiwab Server Handoff

Some analysis involving AnnData or server-only resources may need to run on Qashiwab, which has no internet access. This project should support pasting those outputs back cleanly:

- Put imported Qashiwab code, logs, exported tables, and figures into the relevant `qashiwab_imports/` folder.
- Also copy shared config snapshots used by Qashiwab when relevant, especially cell-type colors and sample metadata.
- Prefer server outputs as CSV/TSV plus figure files, so plots can be inspected and regenerated locally when possible.
- Do not store `.h5ad` files in this project. Code pasted from Qashiwab may reference `.h5ad` paths on the server, but local reproducible inputs should be exported tables or figures.
- Add a short `README.md` or note with each Qashiwab import describing the source path, date, input data, and command or notebook used.

## Exploratory Analysis Principles

- Before creating a new part, decide whether it is exploratory or final. If the user has not specified, ask or make a conservative folder choice based on context.
- Exploratory analyses should be fast to create and easy to discard.
- Exploratory code should still read from `data/` and `config/`, not from hidden local paths.
- Each exploratory topic should live in a dated folder such as `exploratory/2026-06-12_cell_type_composition/`.
- Exploratory outputs do not need to satisfy all final figure requirements, but should keep enough provenance to reproduce promising plots.

## Final Figure Principles

- Final figure scripts should be reproducible from project data and config.
- Python and R are both acceptable, but prefer R for final publication plots when practical because it is often simpler for statistical graphics and journal-style PDF output.
- Each final figure should have clearly separated code, panel outputs, assembled figure, legend, stats, and notes.
- Final panels should use consistent fonts, colors, naming, and sample ordering.
- Each panel should have enough associated statistics and metadata to support the figure legend.
- Final part outputs should be PDF only by default. Add other formats only when there is a clear need.
- For Nature-family submissions, use 5-7 pt sans-serif text at final figure size unless a specific journal requirement says otherwise.
- Use 0.5 pt strokes for standard plot lines, statistical brackets, axes, box outlines and panel borders. In R/ggplot2, convert point units to the package's linewidth units rather than assuming `linewidth = 0.5` equals 0.5 pt.
- Match each part's physical PDF size to its intended size in the assembled multi-panel figure. Adjust canvas width and height before changing fonts.
- For screening plots with many possible comparisons, save unadjusted P values in the part's statistics table and plot unadjusted P values unless the user requests multiplicity adjustment. The user may select biologically relevant comparisons for later adjustment/reporting.
- For Figure 3 and Extended Figure 3 immune-composition plots, use fibro-immune proportions and unadjusted P values by default. Fibro-immune cells are defined as all cells except malignant and A2ML1+ epithelial cells.
- For Figure 3 T-cell subtype plots, use refined T-cell subtype proportions by default. Refined T cells are `CD8_Teff`, `CD8_Tex_PDCD1`, `CD4_Treg_FOXP3`, `CD4_Treg_CCR8`, `CD8_prolif`, `CD4_CXCL13`, and `CD4_prolif`; exclude the generic `T/NK` bucket unless explicitly requested.
- For Figure 4 and Extended Figure 4 stromal/fibroblast plots, use fibro-immune proportions and unadjusted P values by default. Fibro-immune cells are defined as all cells except malignant and A2ML1+ epithelial cells.
- For Figure 2 UMAP parts, use equal x/y coordinates, no axis boundary line, no ticks, and no grid. Use `config/palettes.yaml` for major cell type colors and raster PNG/TIFF output at publication resolution.

## Version Control Principles

Use light version control mainly to keep track of important code, configuration, legends, and manuscript-facing changes.

- Track plotting scripts, reusable helpers, config files, figure legends, notes, and small example tables.
- Do not over-engineer the repository as a software package unless the project later needs it.
- Large generated figures, large raw exports, and temporary exploratory outputs can be ignored if they make the repository noisy.
- Because the work is intended to become public, default to clear names and reproducible paths rather than private abbreviations.

## Discussion Principles

The `discussion/` folder should capture biological thinking separately from code:

- Main observation from each figure.
- Possible biological interpretation.
- Alternative explanations.
- Follow-up analysis ideas.
- Reviewer-risk notes.
- Manuscript wording candidates.

## Suggested Improvements To The Original Plan

1. Add a `config/` folder from the beginning. This will prevent color and label drift across figures.
2. Add `stats/` or `stats_summary` files beside legends. Nature Cancer legends require sample size, center values, error bars, statistical tests, and P values, so the statistics should be stored near the figure.
3. Use semantic panel names when figure numbering is uncertain. This makes it easier to move panels between figures later.
4. Treat raw data as read-only. If we transform data, write the transformed version to `data/processed/` with a script that explains how it was made.
5. Keep Qashiwab imports explicitly separated from locally generated outputs. This will make provenance much clearer when mixing server-side and local work.
6. Use templates for new final figure folders so every figure starts with the same reproducible structure.

## Open Questions Before We Start Data And Code

1. Final assembly: do you want assembled multi-panel figures created by code, by Illustrator/Inkscape/PowerPoint, or both?
2. Data scale: will most local data tables be small enough for ordinary CSV parsing, or should we expect very large spatial tables that need Parquet?
3. Git ignore policy: which generated outputs should be kept in Git, and which should be treated as local/publication artifacts only?

## Decisions Already Made

- Main plotting workflow: both Python and R, with R preferred for final publication graphics when practical.
- Config format: YAML or JSON are both acceptable; YAML is the default for human-readable metadata.
- AnnData boundary: do not store `.h5ad` files in this project. Export tables or figures from Qashiwab instead.
- Version control: use lightweight Git, focused on code, config, legends, and important notes.
- Privacy: no special privacy restriction is expected because the project is intended to be public.
- Final plot output: PDF only by default.
- Repeated colors: keep shared encodings, such as timepoint colors, in `config/`.
- Current malignant-cell proportion plot: final part under `final/figure_2/malignant_cell_proportion/`.
