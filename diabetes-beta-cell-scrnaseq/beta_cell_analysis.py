#!/usr/bin/env python
# coding: utf-8

# %% [markdown]
# # Human Pancreatic Islet scRNA-seq: Beta Cell Analysis in Type 2 Diabetes
#
# Data: GEO accession GSE153855 (human pancreatic islet scRNA-seq,
# control vs. type II diabetes donors).
#
# Goal: identify which pancreatic islet cell types shift in proportion
# between control and T2D donors, then zoom into Beta cells specifically
# to characterize *what changes within Beta cells* in T2D — marker/DEG
# analysis, sub-clustering, a pseudotime trajectory from a healthy root
# state toward the diseased state, and identity vs. stress transcription
# factor expression along that trajectory.

# %% [code cell 0] Setup
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Update these to point at your local copies of the GSE153855 files.
COUNTS_PATH = "./data/GSE153855_Expression_counts_HQ_allsamples.txt.gz"
ANNO_PATH = "./data/GSE153855_Cell_annotation.txt.gz"

# %% [code cell 1] Load counts + cell annotation and build AnnData
counts_df = pd.read_csv(COUNTS_PATH, sep="\t", index_col=0)
anno_df = pd.read_csv(ANNO_PATH, sep="\t")

# Counts columns (cells) and annotation rows are expected to be in the same
# order for this dataset; align annotation to the counts matrix explicitly.
n_cells = counts_df.shape[1]
anno_df = anno_df.iloc[:n_cells]
anno_df.index = counts_df.columns

adata = sc.AnnData(X=counts_df.T, obs=anno_df)
print(adata)

# %% [code cell 2] QC metrics
# Note: this dataset's gene symbols did not match the standard "MT-" mito
# naming convention, so no cells/genes were filtered on %MT here — only on
# minimum genes/cells detected (see below).
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-", "Mt-"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

# %% [code cell 3] Filter low-quality cells and rarely-detected genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata)

# %% [code cell 4] Normalize, log-transform, find highly variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.012, max_mean=3, min_disp=0.5)
print(f"Highly variable genes: {sum(adata.var['highly_variable'])}")

# %% [code cell 5] Subset to HVGs, scale, PCA
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# %% [code cell 6] Neighbors, UMAP, Leiden clustering
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=0.5)
sc.pl.umap(adata, color=["leiden"])

# %% [code cell 7] Cluster marker genes (Wilcoxon)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# %% [code cell 8] Automated cell-type annotation with CellTypist
import celltypist
from celltypist import models

models.download_models()
model = models.Model.load("Adult_Human_PancreaticIslet.pkl")

predictions = celltypist.annotate(adata, model=model, majority_voting=True)
adata.obs["predicted_cell_type"] = predictions.predicted_labels["predicted_labels"]

sc.pl.umap(adata, color=["predicted_cell_type"])

# %% [code cell 9] Cross-check: Leiden clusters vs. predicted cell types
cluster_cell_mapping = pd.crosstab(adata.obs["leiden"], adata.obs["predicted_cell_type"])
print(cluster_cell_mapping)

# %% [code cell 10] Collapse fine-grained labels into broad islet cell types
def simplify_cell_type(label: str) -> str:
    if "alpha" in label:
        return "Alpha cell"
    elif "beta" in label:
        return "Beta cell"
    elif "delta" in label:
        return "Delta cell"
    elif "PP" in label:
        return "PP cell"
    else:
        return "Other"


adata.obs["broad_cell_type"] = adata.obs["predicted_cell_type"].apply(simplify_cell_type)
sc.pl.umap(adata, color=["broad_cell_type"])

# %% [code cell 11] Cell-type proportions: control vs. type II diabetes
condition_col = "Disease"
print(adata.obs[condition_col].value_counts())

cell_counts = pd.crosstab(adata.obs[condition_col], adata.obs["broad_cell_type"], normalize="index") * 100
print("\n--- Cell Type Proportions (%) ---")
print(cell_counts)

ax = cell_counts.plot(kind="bar", stacked=True, figsize=(10, 6), colormap="Set2")
plt.title("Cell Type Proportions: Control vs Type 2 Diabetes", fontsize=14)
plt.xlabel("Condition / Group", fontsize=12)
plt.ylabel("Percentage (%)", fontsize=12)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title="Cell Type")
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

# %% [code cell 12] Per-donor cell-type proportions (sanity check across individuals)
donor_proportions = pd.crosstab(adata.obs["Donor"], adata.obs["broad_cell_type"], normalize="index") * 100
donor_info = adata.obs[["Donor", "Disease"]].drop_duplicates().set_index("Donor")
donor_summary = donor_proportions.join(donor_info)
print(donor_summary)

# %% [markdown]
# ## Beta cell-focused analysis
#
# Beta cells are the insulin-producing population destroyed/dysfunctional
# in diabetes, so from here on the analysis subsets to Beta cells only and
# asks: which genes/pathways differ between T2D and normal Beta cells, and
# is there a trajectory from a "healthy" to a "stressed/diseased" state?

# %% [code cell 13] Subset Beta cells, DEG: T2D vs. normal
adata_beta = adata[adata.obs["broad_cell_type"] == "Beta cell"].copy()

sc.tl.rank_genes_groups(
    adata_beta,
    groupby="Disease",
    reference="normal",
    method="wilcoxon",
)
sc.pl.rank_genes_groups(adata_beta, n_genes=20, sharey=False)

deg_df = sc.get.rank_genes_groups_df(adata_beta, group="type II diabetes")
print(deg_df.head(20))

# %% [code cell 14] Beta cell sub-clustering
sc.pp.pca(adata_beta, svd_solver="arpack")
sc.pp.neighbors(adata_beta, n_neighbors=15, n_pcs=10)
sc.tl.leiden(adata_beta, resolution=0.5, key_added="leiden_beta")

sc.pl.umap(
    adata_beta,
    color=["leiden_beta", "Disease"],
    palette="viridis",
    size=50,
    legend_loc="right margin",
)
plt.suptitle("Beta Cell Sub-clustering: Leiden Clusters & Disease Condition", fontsize=14)
plt.tight_layout()
plt.show()

# Stress marker (P4HB, ER stress/UPR) vs. functional marker (INS)
sc.pl.umap(
    adata_beta,
    color=["leiden_beta", "P4HB", "INS"],
    color_map="RdYlBu_r",
    size=50,
    ncols=3,
)
plt.suptitle("Stress Marker (P4HB) and Insulin (INS) Expression in Beta Cells", fontsize=14)
plt.tight_layout()
plt.show()

# %% [code cell 15] Sub-cluster marker genes
sc.tl.dendrogram(adata_beta, groupby="leiden_beta")
sc.pl.rank_genes_groups_dotplot(adata_beta, n_genes=3, dendrogram=True)

marker_genes_df = sc.get.rank_genes_groups_df(adata_beta, group=None)
print(marker_genes_df.head(15))

# %% [code cell 16] Disease composition of each Beta sub-cluster
cluster_disease_prop = pd.crosstab(adata_beta.obs["leiden_beta"], adata_beta.obs["Disease"], normalize="index") * 100
print(cluster_disease_prop)

# %% [code cell 17] Remove non-Beta contamination sub-clusters
# Clusters 4, 6, 9 showed marker profiles inconsistent with Beta cell identity
# (cross-checked against cluster_cell_mapping / marker_genes_df above) and
# were excluded before trajectory analysis.
contamination_clusters = ["4", "6", "9"]
adata_beta_clean = adata_beta[~adata_beta.obs["leiden_beta"].isin(contamination_clusters)].copy()

sc.pp.pca(adata_beta_clean, svd_solver="arpack")
sc.pp.neighbors(adata_beta_clean, n_neighbors=15, n_pcs=10)
sc.tl.umap(adata_beta_clean)

# %% [code cell 18] Pseudotime trajectory: healthy -> diseased Beta cell state
sc.tl.diffmap(adata_beta_clean)

# Root the trajectory in a predominantly healthy/normal sub-cluster.
healthy_cluster = "0"
root_cells = adata_beta_clean.obs["leiden_beta"] == healthy_cluster
adata_beta_clean.uns["iroot"] = np.flatnonzero(root_cells)[0]

sc.tl.dpt(adata_beta_clean)

sc.pl.umap(
    adata_beta_clean,
    color=["dpt_pseudotime", "leiden_beta", "Disease"],
    color_map="viridis",
    size=50,
    ncols=3,
)
plt.suptitle("Pseudotime Trajectory Analysis of Beta Cells", fontsize=14)
plt.tight_layout()
plt.show()

# %% [code cell 19] Identity vs. stress transcription factors along the trajectory
# PDX1 / MAFA / NKX6-1 / NEUROD1: core Beta cell identity TFs.
# ATF4 / XBP1 / DDIT3 (CHOP): unfolded protein response / ER stress TFs.
# Comparing these across sub-clusters is how "what's being lost/disrupted
# in Beta cells" was characterized — identity TF expression dropping while
# stress TF expression rises marks the diseased trajectory.
target_tfs = ["PDX1", "MAFA", "NKX6-1", "NEUROD1", "ATF4", "XBP1", "DDIT3"]
valid_tfs = [tf for tf in target_tfs if tf in adata_beta_clean.var_names]
print(f"Target TFs found in dataset: {valid_tfs}")

sc.tl.dendrogram(adata_beta_clean, groupby="leiden_beta")
sc.pl.dotplot(
    adata_beta_clean,
    valid_tfs,
    groupby="leiden_beta",
    dendrogram=True,
    cmap="Reds",
    standard_scale="var",
)
plt.suptitle("Expression of Beta-Cell Identity vs. Stress TFs across Sub-clusters", fontsize=14)
plt.show()
