import streamlit as st
import pandas as pd
import numpy as np
import json
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path

st.set_page_config(
    page_title="OmicsOracle",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────
st.markdown("""
<style>
/* Global */
[data-testid="stAppViewContainer"] { background: #0f1117; }
[data-testid="stSidebar"] { background: #161b27; border-right: 1px solid #1e2535; }
[data-testid="stSidebar"] * { color: #a0aec0; }
.block-container { padding: 1.5rem 2rem; }

/* Metric cards */
[data-testid="metric-container"] {
    background: #161b27;
    border: 1px solid #1e2535;
    border-radius: 10px;
    padding: 1rem 1.25rem;
}
[data-testid="stMetricValue"] { font-size: 28px; color: #e2e8f0; }
[data-testid="stMetricLabel"] { color: #718096; font-size: 12px; }

/* Section headers */
.section-header {
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 0.1em;
    color: #4a5568;
    text-transform: uppercase;
    margin: 1.5rem 0 0.75rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid #1e2535;
}

/* Hypothesis cards */
.hyp-card {
    background: #161b27;
    border: 1px solid #1e2535;
    border-radius: 10px;
    padding: 1.25rem;
    margin-bottom: 1rem;
}
.hyp-card:hover { border-color: #2d3748; }
.hyp-title { font-size: 15px; font-weight: 600; color: #e2e8f0; margin-bottom: 0.5rem; }
.hyp-body  { font-size: 13px; color: #a0aec0; line-height: 1.6; margin-bottom: 0.75rem; }
.tag {
    display: inline-block;
    font-size: 11px;
    padding: 2px 8px;
    border-radius: 4px;
    margin: 2px;
    background: #1e2d4a;
    color: #63b3ed;
}
.tag-drug  { background: #1a2f2a; color: #68d391; }
.tag-path  { background: #2d2020; color: #fc8181; }
.conf-high   { background: #2d1a1a; color: #fc8181; font-size:11px; padding:2px 8px; border-radius:4px; }
.conf-medium { background: #2d2a1a; color: #f6ad55; font-size:11px; padding:2px 8px; border-radius:4px; }
.conf-low    { background: #1a2d1a; color: #68d391; font-size:11px; padding:2px 8px; border-radius:4px; }

/* Sidebar nav */
.nav-label {
    font-size: 10px;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: #4a5568;
    padding: 1rem 0 0.25rem;
}
</style>
""", unsafe_allow_html=True)

# ── Data loading ──────────────────────────────────────────────
RES_DIR  = Path("data/results")
PROC_DIR = Path("data/processed")

@st.cache_data
def load_de(ct):
    path = RES_DIR / f"de_{ct}_vs_rest.csv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, index_col=0)
    if PROC_DIR.joinpath("entrez_symbol_map.json").exists():
        with open(PROC_DIR / "entrez_symbol_map.json") as f:
            sym = json.load(f)
        df["symbol"] = df["gene_id"].astype(str).map(sym).fillna(df["gene_id"].astype(str))
    return df

@st.cache_data
def load_pathway(ct, db="KEGG_2021_Human"):
    path = RES_DIR / f"pathway_{ct}_{db}.csv"
    return pd.read_csv(path) if path.exists() else pd.DataFrame()

@st.cache_data
def load_drugs(ct):
    path = RES_DIR / f"drug_targets_{ct}.csv"
    return pd.read_csv(path) if path.exists() else pd.DataFrame()

@st.cache_data
def load_hypotheses(ct):
    path = RES_DIR / f"hypotheses_{ct}.json"
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)

# ── Sidebar ───────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🧬 OmicsOracle")
    st.markdown("<div style='font-size:12px;color:#4a5568;margin-bottom:1rem'>pan-cancer intelligence</div>",
                unsafe_allow_html=True)

    page = st.radio("Navigate", [
        "Overview",
        "Expression",
        "Pathways",
        "Drug targets",
        "Hypotheses",
    ], label_visibility="collapsed")

    st.markdown("<div class='nav-label'>Cancer type</div>", unsafe_allow_html=True)
    selected_ct = st.selectbox("Cancer type", ["GBM", "LUAD", "BRCA"],
                               label_visibility="collapsed")

    st.markdown("<div class='nav-label'>Filters</div>", unsafe_allow_html=True)
    fdr_thresh = st.slider("FDR threshold", 0.01, 0.10, 0.05, 0.01)
    lfc_thresh = st.slider("|log2FC| threshold", 0.5, 3.0, 1.0, 0.25)

# ── Load data for selected cancer type ───────────────────────
de_df   = load_de(selected_ct)
path_df = load_pathway(selected_ct)
drug_df = load_drugs(selected_ct)
hyp     = load_hypotheses(selected_ct)

sig_df = pd.DataFrame()
if not de_df.empty:
    sig_df = de_df[
        (de_df["padj"] < fdr_thresh) &
        (de_df["log2FoldChange"].abs() > lfc_thresh)
    ]

# ── OVERVIEW PAGE ─────────────────────────────────────────────
if page == "Overview":
    st.markdown(f"## {selected_ct} — pan-cancer overview")

    # Metrics
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples", "150", "TCGA primary tumor")
    c2.metric("DEGs", len(sig_df), f"FDR<{fdr_thresh}, |log2FC|>{lfc_thresh}")
    c3.metric("Drug associations", len(drug_df) if not drug_df.empty else 0)
    c4.metric("FDA-approved hits",
              int((drug_df["phase"] >= 4).sum()) if not drug_df.empty else 0)

    # Mini volcano + top genes side by side
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("<div class='section-header'>Volcano plot</div>",
                    unsafe_allow_html=True)
        if not de_df.empty:
            vdf = de_df.copy()
            vdf["-log10padj"] = -np.log10(vdf["padj"].clip(lower=1e-300))
            vdf["color"] = "Non-significant"
            vdf.loc[(vdf["padj"] < fdr_thresh) & (vdf["log2FoldChange"] >  lfc_thresh), "color"] = "Up"
            vdf.loc[(vdf["padj"] < fdr_thresh) & (vdf["log2FoldChange"] < -lfc_thresh), "color"] = "Down"
            color_map = {"Up": "#e05c2e", "Down": "#3b82f6", "Non-significant": "#374151"}

            fig = px.scatter(
                vdf, x="log2FoldChange", y="-log10padj",
                color="color", color_discrete_map=color_map,
                hover_data={"symbol": True, "padj": ":.2e",
                            "log2FoldChange": ":.2f", "color": False},
                template="plotly_dark",
            )
            fig.add_hline(y=-np.log10(fdr_thresh), line_dash="dash",
                          line_color="#4a5568", line_width=0.8)
            fig.add_vline(x= lfc_thresh, line_dash="dash",
                          line_color="#4a5568", line_width=0.8)
            fig.add_vline(x=-lfc_thresh, line_dash="dash",
                          line_color="#4a5568", line_width=0.8)
            fig.update_traces(marker_size=4, marker_opacity=0.7)
            fig.update_layout(
                paper_bgcolor="#0f1117", plot_bgcolor="#0f1117",
                margin=dict(l=0, r=0, t=10, b=0), height=300,
                showlegend=False,
                xaxis=dict(gridcolor="#1e2535", title="log2 fold change"),
                yaxis=dict(gridcolor="#1e2535", title="-log10(adj p-value)"),
            )
            st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown("<div class='section-header'>Top DEGs</div>",
                    unsafe_allow_html=True)
        if not sig_df.empty:
            top = sig_df.nsmallest(12, "padj").copy()
            top["gene"] = top["symbol"] if "symbol" in top.columns else top["gene_id"]
            top["abs_lfc"] = top["log2FoldChange"].abs()
            top["direction"] = top["log2FoldChange"].apply(
                lambda x: "Down" if x < 0 else "Up"
            )
            fig2 = px.bar(
                top.sort_values("log2FoldChange"),
                x="log2FoldChange", y="gene",
                color="direction",
                color_discrete_map={"Down": "#3b82f6", "Up": "#e05c2e"},
                template="plotly_dark",
                orientation="h",
            )
            fig2.update_layout(
                paper_bgcolor="#0f1117", plot_bgcolor="#0f1117",
                margin=dict(l=0, r=0, t=10, b=0), height=300,
                showlegend=False,
                xaxis=dict(gridcolor="#1e2535", title="log2 fold change"),
                yaxis=dict(gridcolor="#1e2535", title=""),
            )
            st.plotly_chart(fig2, use_container_width=True)
        else:
            st.info("No significant DEGs at current thresholds.")

    # Top pathway
    st.markdown("<div class='section-header'>Top KEGG pathways</div>",
                unsafe_allow_html=True)
    if not path_df.empty:
        top_p = path_df.head(12).copy()
        top_p["-log10padj"] = -np.log10(top_p["adj_pval"].clip(lower=1e-50))
        top_p["term_short"] = top_p["term"].str[:55]
        fig3 = px.bar(
            top_p.sort_values("-log10padj"),
            x="-log10padj", y="term_short",
            color="combined_score",
            color_continuous_scale=["#1e2535", "#3b82f6", "#e05c2e"],
            template="plotly_dark", orientation="h",
        )
        fig3.update_layout(
            paper_bgcolor="#0f1117", plot_bgcolor="#0f1117",
            margin=dict(l=0, r=0, t=10, b=0), height=350,
            coloraxis_showscale=False,
            xaxis=dict(gridcolor="#1e2535"),
            yaxis=dict(gridcolor="#1e2535", title=""),
        )
        st.plotly_chart(fig3, use_container_width=True)

# ── EXPRESSION PAGE ───────────────────────────────────────────
elif page == "Expression":
    st.markdown(f"## {selected_ct} — differential expression")

    if de_df.empty:
        st.warning("No DE results found. Run the pipeline first.")
    else:
        # Full volcano
        st.markdown("<div class='section-header'>Volcano plot</div>",
                    unsafe_allow_html=True)
        vdf = de_df.copy()
        vdf["-log10padj"] = -np.log10(vdf["padj"].clip(lower=1e-300))
        vdf["color"] = "Non-significant"
        vdf.loc[(vdf["padj"] < fdr_thresh) & (vdf["log2FoldChange"] >  lfc_thresh), "color"] = "Up"
        vdf.loc[(vdf["padj"] < fdr_thresh) & (vdf["log2FoldChange"] < -lfc_thresh), "color"] = "Down"
        color_map = {"Up": "#e05c2e", "Down": "#3b82f6", "Non-significant": "#374151"}

        fig = px.scatter(
            vdf, x="log2FoldChange", y="-log10padj",
            color="color", color_discrete_map=color_map,
            hover_data={"symbol": True, "padj": ":.2e", "log2FoldChange": ":.2f", "color": False},
            template="plotly_dark",
        )
        fig.add_hline(y=-np.log10(fdr_thresh), line_dash="dash",
                      line_color="#4a5568", line_width=0.8)
        fig.add_vline(x= lfc_thresh, line_dash="dash", line_color="#4a5568", line_width=0.8)
        fig.add_vline(x=-lfc_thresh, line_dash="dash", line_color="#4a5568", line_width=0.8)
        fig.update_traces(marker_size=5, marker_opacity=0.75)
        fig.update_layout(
            paper_bgcolor="#0f1117", plot_bgcolor="#0f1117",
            margin=dict(l=0, r=0, t=20, b=0), height=450,
            legend_title="Direction",
            xaxis=dict(gridcolor="#1e2535", title="log2 fold change"),
            yaxis=dict(gridcolor="#1e2535", title="-log10(adj p-value)"),
        )
        st.plotly_chart(fig, use_container_width=True)

        # DEG table
        st.markdown("<div class='section-header'>Significant DEGs</div>",
                    unsafe_allow_html=True)
        if not sig_df.empty:
            display = sig_df[["symbol", "log2FoldChange", "padj"]].copy()
            display.columns = ["Gene", "log2FC", "Adj p-value"]
            display["log2FC"] = display["log2FC"].round(3)
            display["Adj p-value"] = display["Adj p-value"].apply(lambda x: f"{x:.2e}")
            st.dataframe(display.reset_index(drop=True), use_container_width=True, height=300)
        else:
            st.info("No significant DEGs at current thresholds.")

# ── PATHWAYS PAGE ─────────────────────────────────────────────
elif page == "Pathways":
    st.markdown(f"## {selected_ct} — pathway enrichment")

    db_choice = st.selectbox("Database", [
        "KEGG_2021_Human", "Reactome_2022",
        "GO_Biological_Process_2023", "MSigDB_Hallmark_2020"
    ])
    path_df2 = load_pathway(selected_ct, db_choice)

    if path_df2.empty:
        st.warning("No pathway results found.")
    else:
        top_n = st.slider("Show top N pathways", 5, 30, 15)
        top_p = path_df2.head(top_n).copy()
        top_p["-log10padj"] = -np.log10(top_p["adj_pval"].clip(lower=1e-50))
        top_p["n_genes"]    = top_p["genes"].apply(
            lambda x: len(x) if isinstance(x, list) else len(str(x).split(","))
        )
        top_p["term_short"] = top_p["term"].str[:60]

        fig = px.scatter(
            top_p.sort_values("-log10padj"),
            x="-log10padj", y="term_short",
            size="n_genes", color="combined_score",
            color_continuous_scale=["#1e2535", "#3b82f6", "#e05c2e"],
            template="plotly_dark",
            hover_data={"term": True, "adj_pval": ":.3f",
                        "combined_score": ":.1f", "term_short": False},
        )
        fig.update_layout(
            paper_bgcolor="#0f1117", plot_bgcolor="#0f1117",
            margin=dict(l=0, r=0, t=20, b=0),
            height=max(350, top_n * 28),
            coloraxis_showscale=True,
            coloraxis_colorbar=dict(title="Combined score"),
            xaxis=dict(gridcolor="#1e2535", title="-log10(adj p-value)"),
            yaxis=dict(gridcolor="#1e2535", title=""),
        )
        st.plotly_chart(fig, use_container_width=True)

        st.markdown("<div class='section-header'>Pathway table</div>",
                    unsafe_allow_html=True)
        show = top_p[["term", "adj_pval", "combined_score", "n_genes"]].copy()
        show.columns = ["Pathway", "Adj p-value", "Combined score", "Genes"]
        show["Adj p-value"] = show["Adj p-value"].apply(lambda x: f"{x:.3f}")
        show["Combined score"] = show["Combined score"].round(1)
        st.dataframe(show.reset_index(drop=True), use_container_width=True)

# ── DRUG TARGETS PAGE ─────────────────────────────────────────
elif page == "Drug targets":
    st.markdown(f"## {selected_ct} — drug target actionability")

    if drug_df.empty:
        st.warning("No drug target results found.")
    else:
        col1, col2, col3 = st.columns(3)
        col1.metric("Total associations", len(drug_df))
        col2.metric("FDA-approved drugs",
                    int((drug_df["phase"] >= 4).sum()))
        col3.metric("Unique genes targeted",
                    drug_df["gene"].nunique() if "gene" in drug_df.columns else "—")

        # Phase filter
        phase_filter = st.multiselect(
            "Filter by phase",
            options=[4, 3, 2, 1, 0],
            default=[4, 3],
            format_func=lambda x: {4:"FDA approved", 3:"Phase 3",
                                    2:"Phase 2", 1:"Phase 1", 0:"Preclinical"}[x]
        )
        filtered = drug_df[drug_df["phase"].isin(phase_filter)]

        # Actionability bar chart
        if not filtered.empty:
            st.markdown("<div class='section-header'>Top drug-gene pairs</div>",
                        unsafe_allow_html=True)
            top_drugs = filtered.head(20).copy()
            top_drugs["label"] = top_drugs["gene"] + " → " + top_drugs["drug"].str[:25]
            phase_colors = {4: "#68d391", 3: "#f6ad55", 2: "#63b3ed", 1: "#a0aec0", 0: "#4a5568"}
            top_drugs["color"] = top_drugs["phase"].map(phase_colors)

            fig = go.Figure()
            for phase, grp in top_drugs.groupby("phase"):
                fig.add_trace(go.Bar(
                    x=grp["actionability_score"],
                    y=grp["label"],
                    orientation="h",
                    name={4:"FDA approved",3:"Phase 3",2:"Phase 2",1:"Phase 1",0:"Preclinical"}[phase],
                    marker_color=phase_colors[phase],
                ))
            fig.update_layout(
                paper_bgcolor="#0f1117", plot_bgcolor="#0f1117",
                barmode="stack", height=max(300, len(top_drugs) * 28),
                margin=dict(l=0, r=0, t=10, b=0),
                legend=dict(bgcolor="#0f1117", font=dict(color="#a0aec0")),
                xaxis=dict(gridcolor="#1e2535", title="Actionability score"),
                yaxis=dict(gridcolor="#1e2535", title=""),
                template="plotly_dark",
            )
            st.plotly_chart(fig, use_container_width=True)

        # Full table
        st.markdown("<div class='section-header'>Full drug table</div>",
                    unsafe_allow_html=True)
        show = filtered[["gene", "drug", "phase", "mechanism"]].copy()
        show["phase"] = show["phase"].map({4:"FDA approved",3:"Phase 3",
                                           2:"Phase 2",1:"Phase 1",0:"Preclinical"})
        show.columns = ["Gene", "Drug", "Phase", "Mechanism"]
        st.dataframe(show.reset_index(drop=True), use_container_width=True, height=400)

# ── HYPOTHESES PAGE ───────────────────────────────────────────
elif page == "Hypotheses":
    st.markdown(f"## {selected_ct} — LLM biological hypotheses")
    st.markdown(
        "<div style='font-size:13px;color:#4a5568;margin-bottom:1.5rem'>"
        "Generated by Mistral 7B via OmicsOracle hypothesis agent — "
        "grounded in DEG + pathway + drug target evidence."
        "</div>",
        unsafe_allow_html=True
    )

    if not hyp or "hypotheses" not in hyp:
        st.warning("No hypotheses found. Run the pipeline first.")
    else:
        st.markdown(
            f"<div style='background:#161b27;border:1px solid #1e2535;"
            f"border-radius:10px;padding:1rem 1.25rem;margin-bottom:1.5rem;"
            f"font-size:14px;color:#a0aec0;line-height:1.7'>"
            f"<strong style='color:#e2e8f0'>Summary</strong><br>{hyp.get('summary','')}"
            f"</div>",
            unsafe_allow_html=True
        )

        for h in hyp["hypotheses"]:
            conf = h.get("confidence", "Low")
            conf_class = {"High": "conf-high", "Medium": "conf-medium", "Low": "conf-low"}.get(conf, "conf-low")

            genes_html = "".join(
                f"<span class='tag'>{g}</span>" for g in h.get("supporting_genes", [])
            )
            paths_html = "".join(
                f"<span class='tag tag-path'>{p[:45]}</span>"
                for p in h.get("supporting_pathways", [])
            )
            drug_html = f"<span class='tag tag-drug'>{h.get('drug_opportunity','')[:60]}</span>"

            st.markdown(f"""
<div class='hyp-card'>
  <div style='display:flex;align-items:center;gap:10px;margin-bottom:8px'>
    <span style='color:#4a5568;font-size:12px'>#{h['rank']}</span>
    <span class='hyp-title'>{h['title']}</span>
    <span class='{conf_class}' style='margin-left:auto'>{conf}</span>
  </div>
  <div class='hyp-body'>{h['mechanism']}</div>
  <div style='margin-bottom:6px'>{genes_html}{paths_html}</div>
  <div>{drug_html}</div>
  <div style='margin-top:10px;padding-top:10px;border-top:1px solid #1e2535;
              font-size:12px;color:#4a5568'>
    <strong style='color:#718096'>Validation:</strong> {h.get('validation_experiment','')}
  </div>
</div>
""", unsafe_allow_html=True)

        # Reasoning trace
        if "reasoning_trace" in hyp:
            with st.expander("Reasoning trace"):
                st.json(hyp["reasoning_trace"])