import json
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROC_DIR = Path("data/processed")
RES_DIR  = Path("data/results")

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
    sig_genes[ct] = sig

def query_dgidb(gene_symbols):
    if not gene_symbols:
        return pd.DataFrame()
    query = """
    query GetInteractions($genes: [String!]!) {
      genes(names: $genes) {
        nodes {
          name
          interactions {
            drug { name approved }
            interactionScore
            interactionTypes { type }
          }
        }
      }
    }
    """
    try:
        resp = requests.post("https://dgidb.org/api/graphql",
                             json={"query": query, "variables": {"genes": gene_symbols}},
                             timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"  DGIdb error: {e}")
        return pd.DataFrame()
    rows = []
    for gene_node in data.get("data", {}).get("genes", {}).get("nodes", []):
        gene = gene_node["name"]
        for interaction in gene_node.get("interactions", []):
            drug = interaction.get("drug", {})
            rows.append({
                "gene":        gene,
                "drug":        drug.get("name", ""),
                "approved":    drug.get("approved", False),
                "score":       interaction.get("interactionScore", 0),
                "interaction": ", ".join([t["type"] for t in interaction.get("interactionTypes", [])]),
            })
    return pd.DataFrame(rows)

def build_target_table(ct, dgi):
    rows = []
    if not dgi.empty:
        for _, r in dgi.iterrows():
            rows.append({
                "gene":      r["gene"],
                "drug":      r["drug"],
                "source":    "DGIdb",
                "phase":     4 if r["approved"] else 0,
                "mechanism": r["interaction"],
            })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).drop_duplicates(subset=["gene", "drug"])
    df["actionability_score"] = df["phase"].apply(
        lambda p: 4 if p >= 4 else (3 if p == 3 else (2 if p == 2 else 1))
    )
    return df.sort_values(["actionability_score", "gene"], ascending=[False, True])

dgidb_results, target_tables = {}, {}
for ct in ["LUAD", "GBM"]:
    print(f"Querying DGIdb for {ct}...")
    dgidb_results[ct] = query_dgidb(sig_genes[ct])
    target_tables[ct] = build_target_table(ct, dgidb_results[ct])
    if not target_tables[ct].empty:
        target_tables[ct].to_csv(RES_DIR / f"drug_targets_{ct}.csv", index=False)
        print(f"  Saved {len(target_tables[ct])} drug-gene pairs")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(20, 7))
colors = {4: "#c45c2e", 3: "#e8954a", 2: "#f0c070", 1: "#cccccc"}
for ax, ct in zip(axes, ["LUAD", "GBM"]):
    df = target_tables[ct]
    if df.empty:
        continue
    top = df.head(15).copy()
    top["label"] = top["gene"] + " → " + top["drug"].str[:20]
    bar_colors = top["phase"].map(colors).fillna("#cccccc")
    ax.barh(range(len(top)), top["actionability_score"],
            color=bar_colors, edgecolor="white", height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["label"], fontsize=8)
    ax.set_xlabel("Actionability score")
    ax.set_title(f"{ct} — ranked drug targets", fontsize=11)
    ax.invert_yaxis()
    patches = [
        mpatches.Patch(color="#c45c2e", label="FDA approved"),
        mpatches.Patch(color="#e8954a", label="Phase 3"),
        mpatches.Patch(color="#cccccc", label="Preclinical"),
    ]
    ax.legend(handles=patches, fontsize=7, loc="lower right")
plt.suptitle("Drug target actionability — TCGA pan-cancer", y=1.02, fontsize=13)
plt.tight_layout()
plt.savefig(RES_DIR / "drug_targets.png", dpi=150, bbox_inches="tight")
plt.close()
print("Drug target mining complete.")