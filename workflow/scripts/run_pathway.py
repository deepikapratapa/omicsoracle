import json
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

PROC_DIR = Path("data/processed")
RES_DIR  = Path("data/results")
ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

with open(PROC_DIR / "entrez_symbol_map.json") as f:
    symbol_map = json.load(f)

de_results = {}
for ct in ["BRCA", "LUAD", "GBM"]:
    df = pd.read_csv(RES_DIR / f"de_{ct}_vs_rest.csv", index_col=0)
    df["symbol"] = df["gene_id"].astype(str).map(symbol_map)
    de_results[ct] = df

sig_genes = {}
for ct, df in de_results.items():
    sig = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)]["symbol"].dropna().tolist()
    sig_genes[ct] = [str(g).strip().upper() for g in sig if g and str(g) != "nan"]

def enrichr_query(gene_list, databases=None):
    if databases is None:
        databases = ["KEGG_2021_Human", "Reactome_2022",
                     "GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]
    if not gene_list:
        return {}
    resp = requests.post(f"{ENRICHR_BASE}/addList",
                         files={"list": (None, "\n".join(gene_list)),
                                "description": (None, "OmicsOracle")})
    if not resp.ok:
        return {}
    user_list_id = resp.json()["userListId"]
    results = {}
    for db in databases:
        r = requests.get(f"{ENRICHR_BASE}/enrich",
                         params={"userListId": user_list_id, "backgroundType": db})
        if not r.ok:
            continue
        data = r.json().get(db, [])
        if not data:
            continue
        df = pd.DataFrame(data, columns=["rank","term","pval","zscore","combined_score",
                                          "genes","adj_pval","old_pval","old_adj_pval"])
        df["database"] = db
        results[db] = df.sort_values("adj_pval")
    return results

enrichment_results = {}
for ct in ["LUAD", "GBM"]:
    print(f"Running pathway enrichment for {ct}...")
    enrichment_results[ct] = enrichr_query(sig_genes[ct])
    for db, df in enrichment_results[ct].items():
        df.to_csv(RES_DIR / f"pathway_{ct}_{db}.csv", index=False)

# KEGG bubble plot
fig, axes = plt.subplots(1, 2, figsize=(20, 7))
for ax, ct in zip(axes, ["LUAD", "GBM"]):
    res = enrichment_results[ct]
    db  = "KEGG_2021_Human"
    if db not in res:
        continue
    df = res[db].head(15).copy()
    df["-log10padj"] = -np.log10(df["adj_pval"].clip(lower=1e-50))
    df["n_genes"]    = df["genes"].apply(len)
    df["term_short"] = df["term"].str[:50]
    df = df.sort_values("-log10padj")
    scatter = ax.scatter(df["-log10padj"], range(len(df)),
                         s=df["n_genes"]*20, c=df["combined_score"],
                         cmap="RdYlBu_r", alpha=0.85,
                         edgecolors="white", linewidths=0.5)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["term_short"], fontsize=8)
    ax.set_xlabel("-log10(adj p-value)")
    ax.set_title(f"{ct} — KEGG 2021 Human", fontsize=10)
    plt.colorbar(scatter, ax=ax, label="Combined score")
plt.suptitle("KEGG pathway enrichment — TCGA pan-cancer", y=1.02, fontsize=13)
plt.tight_layout()
plt.savefig(RES_DIR / "pathway_bubble_kegg.png", dpi=150, bbox_inches="tight")
plt.close()

# Hallmark heatmap
hallmark_dfs = []
for ct in ["LUAD", "GBM"]:
    db = "MSigDB_Hallmark_2020"
    if db not in enrichment_results[ct]:
        continue
    top = enrichment_results[ct][db].head(20).copy()
    top["cancer_type"]  = ct
    top["term_short"]   = top["term"].str[:45]
    top["-log10padj"]   = -np.log10(top["adj_pval"].clip(lower=1e-50))
    hallmark_dfs.append(top)

if hallmark_dfs:
    hall  = pd.concat(hallmark_dfs)
    pivot = hall.pivot_table(index="term_short", columns="cancer_type",
                              values="-log10padj", fill_value=0)
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(8, 10))
    sns.heatmap(pivot, cmap="YlOrRd", linewidths=0.4, annot=True,
                fmt=".1f", ax=ax, cbar_kws={"label": "-log10(adj p-value)"})
    ax.set_title("MSigDB Hallmark pathways — LUAD vs GBM", fontsize=12)
    plt.tight_layout()
    plt.savefig(RES_DIR / "hallmark_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()

print("Pathway enrichment complete.")