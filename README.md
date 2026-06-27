# OmicsOracle

**Agentic multi-omics intelligence pipeline — from raw RNA-seq to ranked biological hypotheses.**

OmicsOracle is a production-grade bioinformatics platform that autonomously processes TCGA RNA-seq data through differential expression, pathway enrichment, and drug target mining, then synthesizes findings into ranked biological hypotheses using a local LLM agent. Built to demonstrate the full stack of modern computational biology: reproducible workflows, external API integration, agentic AI reasoning, and interactive clinical dashboards.

---

## What it does

```
Raw TCGA RNA-seq counts
        ↓
QC + normalization (150 samples, 23K genes, AnnData)
        ↓
DESeq2 differential expression (PyDESeq2, tumor-filtered)
        ↓
Pathway enrichment (Enrichr → KEGG / Reactome / GO / MSigDB Hallmarks)
        ↓
Drug target mining (DGIdb v5 GraphQL → 237 associations, 120 FDA-approved)
        ↓
LLM hypothesis agent (Mistral 7B via Ollama → ranked hypotheses + reasoning trace)
        ↓
Interactive Streamlit dashboard (volcano, pathways, drug targets, hypotheses)
```

---

## Key results

| Cancer type | DEGs | Top pathway | FDA-approved drug hits |
|---|---|---|---|
| GBM | 42 | ABC transporters (score: 127) | 63 |
| LUAD | 17 | p53 signaling / Death receptor | 58 |
| BRCA | — | Insufficient signal at n=44 | — |

**Selected biological findings:**
- XRCC3 downregulation in GBM → DNA repair deficiency → synthetic lethality opportunity with PARP inhibitors
- ABCB11 loss in GBM → ABC transporter dysregulation → Lapatinib and Cyclosporine as repurposing candidates
- PPM1D downregulation in LUAD → p53 pathway disruption → Doxorubicin HCl as top FDA-approved hit
- ECM1 extreme downregulation in LUAD (log2FC = −43) → extracellular matrix remodeling signature

---

## Tech stack

| Layer | Tools |
|---|---|
| Data ingestion | GDC API, TCGA RNA-seq (STAR counts), AnnData |
| Differential expression | PyDESeq2, DESeq2, scanpy |
| Pathway enrichment | Enrichr API, KEGG, Reactome, GO, MSigDB Hallmarks |
| Drug targets | DGIdb v5 GraphQL, OpenTargets |
| Gene annotation | NCBI E-utilities (Entrez → HGNC symbol mapping) |
| LLM agent | Mistral 7B via Ollama (local), LangGraph-ready |
| Workflow | Snakemake DAG (6 rules, parallel execution) |
| Dashboard | Streamlit, Plotly |
| Infrastructure | Docker, AWS S3/EC2, GitHub Actions CI |
| Language | Python 3.11, R (DESeq2 compatible) |

---

## Project structure

```
omicsoracle/
├── configs/
│   └── config.yaml              # Pipeline parameters
├── data/
│   ├── raw/tcga/                # TCGA downloads (gitignored)
│   ├── processed/               # AnnData h5ad, symbol maps
│   └── results/                 # DE tables, pathway CSVs, hypotheses
├── notebooks/
│   ├── 01_eda/                  # TCGA ingestion + UMAP
│   ├── 02_de_analysis/          # PyDESeq2 + volcano plots
│   ├── 03_pathway/              # Enrichr multi-database enrichment
│   ├── 04_drug_targets/         # DGIdb v5 GraphQL mining
│   └── 05_llm_agent/            # Mistral hypothesis generation
├── src/omicsoracle/
│   ├── agents/                  # LangGraph agent definitions
│   ├── pipeline/                # Core pipeline modules
│   ├── api/                     # FastAPI endpoints
│   └── utils/                   # TCGA download, gene mapping
├── workflow/
│   ├── Snakefile                # Full DAG orchestration
│   └── scripts/                 # Standalone pipeline scripts
├── app/
│   └── streamlit_app.py         # Interactive dashboard
└── tests/                       # pytest test suite
```

---

## Quickstart

### 1. Clone and set up environment

```bash
git clone https://github.com/deepikapratapa/omicsoracle.git
cd omicsoracle

conda create -n omicsoracle python=3.11 -y
conda activate omicsoracle

conda install -c conda-forge -c bioconda \
  jupyter anndata scanpy snakemake -y

pip install -r requirements.txt
pip install -e .
```

### 2. Download TCGA data

```bash
python src/omicsoracle/utils/tcga_download.py
```

Downloads 150 TCGA RNA-seq files (BRCA, LUAD, GBM) via GDC API. No account required. ~160MB.

### 3. Run the full pipeline

```bash
# Dry run to preview DAG
snakemake --dry-run

# Execute all 6 rules with 4 cores
snakemake --cores 4
```

Pipeline completes in ~15 minutes. Outputs saved to `data/results/`.

### 4. Start the LLM agent (requires Ollama)

```bash
# Install and start Ollama
brew install ollama
ollama pull mistral
ollama serve
```

### 5. Launch the dashboard

```bash
streamlit run app/streamlit_app.py
```

Opens at `http://localhost:8501`

---

## Pipeline DAG

```
download_tcga
      ↓
    eda
      ↓
differential_expression
    ↓         ↓
pathway    drug_targets
    ↓         ↓
    hypothesis_agent
          ↓
         all
```

Pathway enrichment and drug target mining run in parallel after DE completes.

---

## Dashboard

Five interactive pages:

- **Overview** — key metrics, mini volcano, top DEGs, KEGG pathway ranking
- **Expression** — full interactive volcano with hover tooltips, filterable DEG table
- **Pathways** — bubble plot across KEGG / Reactome / GO / Hallmarks, database switcher
- **Drug targets** — ranked drug-gene pairs, phase filter, FDA-approved highlights
- **Hypotheses** — LLM-generated ranked hypotheses with genes, pathways, drug opportunities, and validation experiments

All views are parameterized by cancer type (BRCA / LUAD / GBM) and adjustable FDR / log2FC thresholds via sidebar.

---

## Reproducing the analysis

All results are fully reproducible from public data:

```bash
# Wipe results and rerun from scratch
rm -rf data/processed data/results
snakemake --cores 4
```

Expected outputs after full run:
- `data/results/de_{BRCA,LUAD,GBM}_vs_rest.csv`
- `data/results/pathway_{LUAD,GBM}_*.csv` (4 databases each)
- `data/results/drug_targets_{LUAD,GBM}.csv`
- `data/results/hypotheses_{LUAD,GBM}.json`
- `data/results/hypotheses_{LUAD,GBM}.md`

---

## Configuration

Edit `configs/config.yaml` to change:

```yaml
data:
  tcga_cancers: ["BRCA", "LUAD", "GBM"]  # cancer types to analyze

de:
  fdr_threshold: 0.05       # adjusted p-value cutoff
  log2fc_threshold: 1.0     # minimum fold change

pathway:
  databases:                # Enrichr databases to query
    - KEGG_2021_Human
    - Reactome_2022
    - GO_Biological_Process_2023
    - MSigDB_Hallmark_2020
  top_n: 20

drug_targets:
  top_n: 10                 # top drug-gene pairs per cancer type
```

---

## Requirements

```
python>=3.11
pydeseq2
anndata
scanpy
snakemake
streamlit
plotly
biopython
langchain
langgraph
ollama
fastapi
requests
pandas
numpy
scipy
seaborn
matplotlib
tqdm
```

Full list: `requirements.txt`

---

## Roadmap

- [ ] scRNA-seq integration (cell type deconvolution)
- [ ] ATAC-seq chromatin accessibility layer
- [ ] Snakemake cloud execution (AWS Batch)
- [ ] FastAPI REST endpoints for programmatic access
- [ ] Docker image for zero-setup deployment
- [ ] Expanded to 10 TCGA cancer types

---

## Author

**Deepika Sarala Pratapa**
M.S. Applied Data Science, University of Florida (2026)
B.Tech Biotechnology and Bioinformatics, KL University

[GitHub](https://github.com/deepikapratapa) · [Portfolio](https://deepikapratapa.github.io) · [LinkedIn](https://linkedin.com/in/deepikapratapa)

---

## License

MIT License — open source, free to use and extend.