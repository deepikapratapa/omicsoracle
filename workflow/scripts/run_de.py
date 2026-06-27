import json
import requests
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from tqdm import tqdm

PROC_DIR = Path("data/processed")
RES_DIR  = Path("data/results")

adata = ad.read_h5ad(PROC_DIR / "tcga_pancancer_raw.h5ad")
adata.obs_names_make_unique()

counts_matrix = pd.DataFrame(
    adata.layers["counts"],
    index=adata.obs_names,
    columns=adata.var_names
).astype(int)

def run_deseq2(counts, obs, target_cancer):
    tumor_mask = obs["sample_type"].str.contains("Primary", case=False, na=False)
    obs_filtered = obs[tumor_mask]
    counts_filtered = counts.loc[obs_filtered.index]
    meta = obs_filtered[["cancer_type"]].copy()
    meta["condition"] = (meta["cancer_type"] == target_cancer).map(
        {True: target_cancer, False: "other"}
    )
    dds = DeseqDataSet(counts=counts_filtered, metadata=meta,
                       design="~condition", refit_cooks=True, quiet=True)
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", target_cancer, "other"], quiet=True)
    stat_res.summary()
    results = stat_res.results_df.copy()
    results["cancer_type"] = target_cancer
    results["gene_id"] = results.index
    results = results.dropna(subset=["padj"]).sort_values("padj")
    return results

de_results = {}
for ct in ["BRCA", "LUAD", "GBM"]:
    print(f"Running DESeq2: {ct} vs rest...")
    de_results[ct] = run_deseq2(counts_matrix, adata.obs, ct)
    sig = de_results[ct][
        (de_results[ct]["padj"] < 0.05) &
        (de_results[ct]["log2FoldChange"].abs() > 1)
    ]
    print(f"  Significant DEGs: {len(sig):,}")
    de_results[ct].to_csv(RES_DIR / f"de_{ct}_vs_rest.csv")

# Entrez → symbol map
all_gene_ids = set()
for df in de_results.values():
    all_gene_ids.update(df["gene_id"].astype(str).tolist())

def entrez_to_symbol(gene_ids, batch_size=200):
    mapping = {}
    ids = [str(g) for g in gene_ids]
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        resp = requests.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            data={"db": "gene", "id": ",".join(batch), "retmode": "json"}
        )
        data = resp.json().get("result", {})
        for uid in batch:
            mapping[uid] = data.get(uid, {}).get("name", uid)
    return mapping

print(f"Mapping {len(all_gene_ids):,} genes to symbols...")
symbol_map = entrez_to_symbol(list(all_gene_ids))
with open(PROC_DIR / "entrez_symbol_map.json", "w") as f:
    json.dump(symbol_map, f)

# Add symbols
for ct, df in de_results.items():
    de_results[ct]["symbol"] = df["gene_id"].astype(str).map(symbol_map)

# Volcano plots
def volcano_plot(results, cancer_type, fdr=0.05, lfc=1.0, top_n=10, ax=None):
    df = results.copy()
    df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["color"] = "gray"
    df.loc[(df["padj"] < fdr) & (df["log2FoldChange"] >  lfc), "color"] = "up"
    df.loc[(df["padj"] < fdr) & (df["log2FoldChange"] < -lfc), "color"] = "down"
    palette = {"up": "#c45c2e", "down": "#1a5c8a", "gray": "#cccccc"}
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    for group, color in palette.items():
        sub = df[df["color"] == group]
        ax.scatter(sub["log2FoldChange"], sub["-log10padj"],
                   c=color, s=8, alpha=0.6, linewidths=0)
    ax.axhline(-np.log10(fdr), color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axvline(lfc,  color="black", lw=0.8, ls="--", alpha=0.5)
    ax.axvline(-lfc, color="black", lw=0.8, ls="--", alpha=0.5)
    up_n   = (df["color"] == "up").sum()
    down_n = (df["color"] == "down").sum()
    ax.set_title(f"{cancer_type} vs rest  |  ↑{up_n:,}  ↓{down_n:,}", fontsize=11)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(adj p-value)")
    return ax

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, ct in zip(axes, ["BRCA", "LUAD", "GBM"]):
    volcano_plot(de_results[ct], ct, ax=ax)
plt.suptitle("Differential expression — TCGA pan-cancer", y=1.02, fontsize=13)
plt.tight_layout()
plt.savefig(RES_DIR / "volcano_plots.png", dpi=150, bbox_inches="tight")
plt.close()
print("Volcano plots saved.")

# Heatmap of top DEGs
import anndata as ad

# Add symbols to de_results
for ct in ["BRCA", "LUAD", "GBM"]:
    de_results[ct]["symbol"] = de_results[ct]["gene_id"].astype(str).map(symbol_map)

top_genes = []
for ct, df in de_results.items():
    sig = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)]
    top = sig.nsmallest(20, "padj").index.tolist()
    top_genes.extend(top)
top_genes = list(dict.fromkeys(top_genes))

adata2 = ad.read_h5ad(PROC_DIR / "tcga_pancancer_raw.h5ad")
adata2.obs_names_make_unique()

expr = pd.DataFrame(adata2.X, index=adata2.obs_names, columns=adata2.var_names)
top_genes = [g for g in top_genes if g in expr.columns]

if top_genes:
    heat_data   = expr[top_genes].T
    sample_order = adata2.obs.sort_values("cancer_type").index
    heat_data    = heat_data[sample_order]
    ct_colors    = {"BRCA": "#1a5c8a", "LUAD": "#c45c2e", "GBM": "#2e7d4f"}
    col_colors   = adata2.obs.loc[sample_order, "cancer_type"].map(ct_colors)

    g = sns.clustermap(
        heat_data, col_colors=col_colors, col_cluster=False,
        row_cluster=True, cmap="RdBu_r", center=0,
        figsize=(14, 8),
        yticklabels=[g.split(".")[0] for g in top_genes],
        xticklabels=False,
        cbar_kws={"label": "log1p(CPM)"},
    )
    patches = [mpatches.Patch(color=c, label=l) for l, c in ct_colors.items()]
    plt.legend(handles=patches, bbox_to_anchor=(1.15, 1), loc="upper left")
    plt.savefig(RES_DIR / "deg_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Heatmap saved.")
else:
    # Create empty placeholder so Snakemake doesn't fail
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(0.5, 0.5, "No significant DEGs for heatmap",
            ha="center", va="center", transform=ax.transAxes)
    plt.savefig(RES_DIR / "deg_heatmap.png", dpi=150)
    plt.close()
    print("Placeholder heatmap saved.")

print("DE analysis complete.")