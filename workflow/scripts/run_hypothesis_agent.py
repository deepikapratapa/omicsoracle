import json
import requests
import pandas as pd
from pathlib import Path

PROC_DIR   = Path("data/processed")
RES_DIR    = Path("data/results")
OLLAMA_URL = "http://localhost:11434/api/chat"
MODEL      = "mistral"

def ollama_chat(system, user):
    payload = {
        "model": MODEL, "stream": False,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user",   "content": user},
        ],
        "options": {"temperature": 0.2}
    }
    resp = requests.post(OLLAMA_URL, json=payload, timeout=300)
    resp.raise_for_status()
    return resp.json()["message"]["content"]

SYSTEM_PROMPT = """You are OmicsOracle, an expert computational biologist and
translational oncologist. You receive structured multi-omics pipeline outputs
and synthesize them into ranked biological hypotheses.
Generate exactly 3 ranked hypotheses. Return valid JSON only, no markdown."""

with open(PROC_DIR / "entrez_symbol_map.json") as f:
    symbol_map = json.load(f)

de_results, sig_genes, pathway_results, drug_results = {}, {}, {}, {}
for ct in ["LUAD", "GBM"]:
    df = pd.read_csv(RES_DIR / f"de_{ct}_vs_rest.csv", index_col=0)
    df["symbol"] = df["gene_id"].astype(str).map(symbol_map)
    de_results[ct] = df
    sig = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)][
        ["symbol", "log2FoldChange", "padj"]].dropna()
    sig_genes[ct] = sig
    pathway_results[ct] = {}
    for db in ["KEGG_2021_Human", "Reactome_2022", "MSigDB_Hallmark_2020"]:
        path = RES_DIR / f"pathway_{ct}_{db}.csv"
        if path.exists():
            pathway_results[ct][db] = pd.read_csv(path)
    path = RES_DIR / f"drug_targets_{ct}.csv"
    if path.exists():
        drug_results[ct] = pd.read_csv(path)

def build_context(ct):
    degs = sig_genes.get(ct, pd.DataFrame())
    top_degs = []
    if not degs.empty:
        for _, row in degs.head(15).iterrows():
            top_degs.append({
                "gene": row["symbol"],
                "log2fc": round(float(row["log2FoldChange"]), 2),
                "padj": float(row["padj"]),
                "direction": "down" if row["log2FoldChange"] < 0 else "up"
            })
    pathways = {}
    for db, df in pathway_results.get(ct, {}).items():
        if df is not None and not df.empty:
            pathways[db] = df.head(5)[["term","adj_pval","combined_score"]].to_dict("records")
    drugs = drug_results.get(ct, pd.DataFrame())
    top_drugs = []
    if not drugs.empty:
        approved = drugs[drugs["phase"] >= 4].drop_duplicates(subset=["gene","drug"])
        for _, row in approved.head(15).iterrows():
            top_drugs.append({
                "gene": row["gene"], "drug": row["drug"],
                "mechanism": row.get("mechanism", "")
            })
    return {"cancer_type": ct, "n_degs": len(degs),
            "top_degs": top_degs, "pathways": pathways, "drug_targets": top_drugs}

for ct in ["LUAD", "GBM"]:
    print(f"Running hypothesis agent for {ct}...")
    context = build_context(ct)
    user_prompt = f"""
Analyze multi-omics results for {ct} and generate 3 ranked hypotheses.

DEGs: {json.dumps(context['top_degs'], indent=2)}
Pathways: {json.dumps(context['pathways'], indent=2)}
Drug targets: {json.dumps(context['drug_targets'], indent=2)}

Return this exact JSON structure, no markdown:
{{
  "cancer_type": "{ct}",
  "summary": "2-3 sentence narrative",
  "hypotheses": [
    {{
      "rank": 1,
      "title": "title",
      "mechanism": "mechanism",
      "supporting_genes": ["gene1"],
      "supporting_pathways": ["pathway1"],
      "drug_opportunity": "drug and rationale",
      "confidence": "High or Medium or Low",
      "validation_experiment": "experiment"
    }}
  ]
}}"""

    content = ollama_chat(SYSTEM_PROMPT, user_prompt)
    if "```" in content:
        content = content.split("```")[1]
        if content.startswith("json"):
            content = content[4:]
        content = content.strip()
    start = content.find("{")
    end   = content.rfind("}") + 1
    content = content[start:end]

    result = json.loads(content)
    result["reasoning_trace"] = {"model": MODEL, "n_degs": context["n_degs"]}

    with open(RES_DIR / f"hypotheses_{ct}.json", "w") as f:
        json.dump(result, f, indent=2)

    lines = [f"# OmicsOracle — {ct}\n", f"## Summary\n{result['summary']}\n"]
    for h in result["hypotheses"]:
        lines += [
            f"### #{h['rank']} — {h['title']}",
            f"**Confidence:** {h['confidence']}",
            f"**Mechanism:** {h['mechanism']}",
            f"**Genes:** {', '.join(h['supporting_genes'])}",
            f"**Drug:** {h['drug_opportunity']}\n"
        ]
    (RES_DIR / f"hypotheses_{ct}.md").write_text("\n".join(lines))
    print(f"  {ct} done — {len(result['hypotheses'])} hypotheses saved")

print("Hypothesis agent complete.")