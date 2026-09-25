# Human Pancreatic Islet scRNA-seq: Beta Cell Analysis in Type 2 Diabetes

Analysis of public human pancreatic islet single-cell RNA-seq data
([GEO GSE153855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153855),
control vs. type II diabetes donors) with Scanpy. The pipeline first
annotates all islet cell types, then focuses specifically on **Beta cells**
to characterize what changes within the insulin-producing population in
type 2 diabetes.

## Pipeline overview

1. **Load & build AnnData** — read raw count matrix + per-cell annotation, align to build an `AnnData` object.
2. **QC & filtering** — mitochondrial gene flagging, QC metrics, filter cells with `< 200` genes and genes detected in `< 3` cells.
3. **Normalization & clustering** — `normalize_total` → `log1p` → HVG selection → scale → PCA → neighbors → UMAP → Leiden clustering.
4. **Cell-type annotation** — automated annotation via [CellTypist](https://www.celltypist.org/) (`Adult_Human_PancreaticIslet` model, majority voting), cross-checked against Leiden clusters, then collapsed into broad islet cell types (Alpha / Beta / Delta / PP / Other).
5. **Islet-level composition** — cell-type proportions compared between control and T2D donors (cohort-level and per-donor).
6. **Beta cell deep-dive**:
   - DEG between T2D and normal Beta cells (Wilcoxon)
   - Beta cell sub-clustering (Leiden)
   - Removal of non-Beta contamination sub-clusters (identified via marker genes)
   - Diffusion pseudotime trajectory from a healthy root sub-cluster toward the diseased state
   - Expression of Beta-cell identity TFs (`PDX1`, `MAFA`, `NKX6-1`, `NEUROD1`) vs. ER-stress/UPR TFs (`ATF4`, `XBP1`, `DDIT3`) across sub-clusters and along the trajectory

This last step is how the analysis addresses "what's mainly disrupted in
Beta cells in diabetes": identity-TF expression drops and stress-TF
expression rises along the healthy → diseased trajectory, alongside the
P4HB (stress) vs. INS (function) marker comparison in step 6.

## Repository structure

```
diabetes-beta-cell-scrnaseq/
├── README.md
├── requirements.txt
├── .gitignore
└── beta_cell_analysis.py
```

`beta_cell_analysis.py` uses `# %%` cell markers (Jupytext-style) — it opens
as a notebook, cell-by-cell, in Jupyter or VS Code, or runs top-to-bottom as
a plain script.

## Data

Download the two GSE153855 supplementary files from GEO and place them
under `data/` (or update the paths at the top of the script):
- `GSE153855_Expression_counts_HQ_allsamples.txt.gz`
- `GSE153855_Cell_annotation.txt.gz`

## Installation

```bash
git clone https://github.com/<your-username>/diabetes-beta-cell-scrnaseq.git
cd diabetes-beta-cell-scrnaseq
pip install -r requirements.txt
```

CellTypist downloads its reference model (`Adult_Human_PancreaticIslet.pkl`)
on first run via `models.download_models()`.

## Usage

```bash
python beta_cell_analysis.py
```

Or open it in Jupyter/VS Code and run cell-by-cell.

## Notes

- This dataset's gene symbols didn't follow the standard `MT-` mitochondrial
  naming convention, so no explicit %MT-based filtering was applied — only
  the minimum genes/cells thresholds.
- The contamination sub-clusters excluded before trajectory analysis
  (`leiden_beta` clusters 4, 6, 9) were identified by cross-referencing
  sub-cluster marker genes against non-Beta islet cell markers; re-run the
  marker gene step on your own clustering if the resolution or sub-cluster
  numbering changes.

## License

Add a license (e.g. MIT) if you plan to share this publicly.
