import json
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for scripts
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path("data/raw/tcga")
PROC_DIR = Path("data/processed")
RES_DIR  = Path("data/results")
PROC_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

def load_star_counts(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", comment="#",
                     names=["gene_id","gene_name","gene_type",
                            "unstranded","stranded_first","stranded_second"])
    df = df[~df["gene_id"].str.startswith("N_")]
    df = df.set_index("gene_id")
    df = df[~df.index.duplicated(keep="first")]
    return df["unstranded"]

def load_cancer_type(ct):
    ct_dir = DATA_DIR / ct
    manifest = json.load(open(ct_dir / "manifest.json"))
    frames = {}
    for file_id, meta in manifest.items():
        file_dir = ct_dir / file_id
        tsvs = list(file_dir.glob("*.tsv"))
        if not tsvs:
            continue
        sample_id = meta["cases"][0]["submitter_id"]
        if sample_id in frames:
            continue
        frames[sample_id] = load_star_counts(tsvs[0])
    df = pd.DataFrame(frames)
    return df, manifest

# Load
counts, manifests = {}, {}
for ct in ["BRCA", "LUAD", "GBM"]:
    print(f"Loading {ct}...")
    counts[ct], manifests[ct] = load_cancer_type(ct)
    df = counts[ct]
    print(f"  {df.shape[0]:,} genes × {df.shape[1]} samples")

# Merge
dfs = []
for ct, df in counts.items():
    df = df.copy()
    df.columns = [f"{ct}_{c}" for c in df.columns]
    dfs.append(df)

merged = pd.concat(dfs, axis=1).fillna(0)
merged = merged.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)

# Obs metadata
obs_rows = []
for ct, manifest in manifests.items():
    for fid, meta in manifest.items():
        case = meta["cases"][0]
        sid  = f"{ct}_{case['submitter_id']}"
        if sid not in merged.columns:
            continue
        obs_rows.append({
            "sample_id":   sid,
            "cancer_type": ct,
            "case_id":     case["case_id"],
            "gender":      case.get("demographic", {}).get("gender", "unknown"),
            "sample_type": case["samples"][0]["sample_type"] if case.get("samples") else "unknown",
        })

obs = pd.DataFrame(obs_rows).set_index("sample_id")
obs = obs.loc[obs.index.isin(merged.columns)]
merged = merged[obs.index]

adata = ad.AnnData(X=merged.T.values, obs=obs,
                   var=pd.DataFrame(index=merged.index))
adata.var_names = merged.index.tolist()
adata.obs_names_make_unique()

# QC plot
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, ct in zip(axes, ["BRCA", "LUAD", "GBM"]):
    sub = adata[adata.obs.cancer_type == ct]
    lib_sizes = sub.X.sum(axis=1)
    ax.hist(lib_sizes / 1e6, bins=15, color="#1a5c8a", alpha=0.8, edgecolor="white")
    ax.set_title(ct)
    ax.set_xlabel("Library size (M reads)")
    ax.set_ylabel("Samples")
plt.suptitle("Library size distribution by cancer type", y=1.02)
plt.tight_layout()
plt.savefig(RES_DIR / "qc_library_sizes.png", dpi=150, bbox_inches="tight")
plt.close()

# Filter + normalize
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)

# UMAP plot
sc.pl.umap(adata, color="cancer_type",
           palette={"BRCA": "#1a5c8a", "LUAD": "#c45c2e", "GBM": "#2e7d4f"},
           title="UMAP — TCGA pan-cancer (BRCA / LUAD / GBM)",
           show=False, save=False)
plt.savefig(RES_DIR / "umap_cancer_type.png", dpi=150, bbox_inches="tight")
plt.close()

adata.write_h5ad(PROC_DIR / "tcga_pancancer_raw.h5ad")
print(f"Saved → {PROC_DIR / 'tcga_pancancer_raw.h5ad'}")
print(f"Shape: {adata.shape[0]} samples × {adata.shape[1]} genes")