import os
import sqlite3
import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Marine Bioactives & Breast Cancer Explorer", page_icon="🧬", layout="wide")
DB_PATH = "marine_bioactives.db"

CHEMICAL_DB = {
    "Araguspongine C": {"formula": "C28H50N2O4", "mw": "478.71 g/mol", "class": "Macrocyclic bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "10459942", "smiles": "CC1CCCC2N1CCCC3(C(O2)CCCC4N3CCCC(O4)C)O", "mechanism": "Reported PI3K/AKT/mTOR inhibitor inducing autophagic cell death across HR+, HER2+, and TNBC models."},
    "Meso-araguspongine C": {"formula": "C28H50N2O4", "mw": "478.71 g/mol", "class": "Macrocyclic bis-1-oxaquinolizidine alkaloid (meso stereoisomer)", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "10459942", "smiles": "Stereoisomer of Araguspongine C (meso configuration)", "mechanism": "Potent antiproliferative activity in ER+ MCF-7 cells (IC50 = 0.44 µM)."},
    "Araguspongine A": {"formula": "C28H50N2O2", "mw": "446.71 g/mol", "class": "Bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "101552483", "smiles": "C1CCN2CCCC3C(O2)CCCC4N3CCCC(O4)CC1", "mechanism": "Cytotoxic activity against HR+ breast carcinoma (MCF-7)."},
    "Araguspongine E": {"formula": "C28H48N2O3", "mw": "460.69 g/mol", "class": "Bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "102124567", "smiles": "Oxo-substituted bis-1-oxaquinolizidine scaffold", "mechanism": "Moderate growth inhibition in MCF-7 cells (IC50 = 24.85 µM)."},
    "Araguspongine L": {"formula": "C28H50N2O3", "mw": "462.71 g/mol", "class": "Hydroxylated bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "102124568", "smiles": "Monohydroxy bis-1-oxaquinolizidine scaffold", "mechanism": "Moderate growth inhibition in MCF-7 cells (IC50 = 24.85 µM)."},
    "Araguspongine N": {"formula": "C28H50N2O3", "mw": "462.71 g/mol", "class": "Bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "102124569", "smiles": "Bis-1-oxaquinolizidine derivative", "mechanism": "Cytotoxic activity in MCF-7 breast cancer cells (IC50 = 7.36 µM)."},
    "Araguspongine O": {"formula": "C28H50N2O3", "mw": "462.71 g/mol", "class": "Bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "102124570", "smiles": "Bis-1-oxaquinolizidine derivative", "mechanism": "Cytotoxic activity in MCF-7 breast cancer cells (IC50 = 5.32 µM)."},
    "Araguspongine P": {"formula": "C28H50N2O4", "mw": "478.71 g/mol", "class": "Dihydroxy bis-1-oxaquinolizidine alkaloid", "source": "Marine sponge (Xestospongia sp.)", "pubchem_cid": "102124571", "smiles": "Dihydroxy bis-1-oxaquinolizidine scaffold", "mechanism": "Cytotoxic activity in MCF-7 breast cancer cells (IC50 = 5.68 µM)."},
    "Manzamine A": {"formula": "C36H44N4O", "mw": "548.76 g/mol", "class": "Beta-carboline polycyclic alkaloid", "source": "Marine sponge (Haliclona / Acanthostrongylophora sp.)", "pubchem_cid": "6436222", "smiles": "C1CC=CCCC2C3C=CC4(C3C5(C2N1)CC6=C(C5)C7=CC=CC=C7N6)CCN8C4CCC=CC8O", "mechanism": "Inhibits GSK-3β, CDK5, and autophagy flux; active in both HR+ (MCF-7) and TNBC (MDA-MB-231)."},
    "Renieramycin M": {"formula": "C31H34N2O9", "mw": "578.61 g/mol", "class": "Bistetrahydroisoquinolinequinone alkaloid", "source": "Blue marine sponge (Xestospongia sp.)", "pubchem_cid": "11618026", "smiles": "CC1=C(C(=O)C2=C(C1=O)C3C4CC5=C(C(=C(C(=C5OC)C)O)C)C(N4C(C(N3C)C#N)COC(=O)C(=CC)C)C2)OC", "mechanism": "Ultra-potent nanomolar cytotoxicity (3.8 nM in MDA-MB-231); induces p53-dependent apoptosis."},
    "Fascaplysin": {"formula": "C18H11ClN2O", "mw": "306.75 g/mol", "class": "Pentacyclic pyrido[1,2-a]benzimidazole alkaloid", "source": "Marine sponge (Fascaplysinopsis sp.)", "pubchem_cid": "73432", "smiles": "C1=CC=C2C(=C1)C3=C(N2)C4=[N+](C=CC5=CC=CC=C54)C3=O.[Cl-]", "mechanism": "Selective CDK4/cyclin D1 inhibitor causing G1 cell cycle arrest."},
    "4-Chlorofascaplysin": {"formula": "C18H10Cl2N2O", "mw": "341.19 g/mol", "class": "Halogenated fascaplysin derivative", "source": "Marine sponge / synthetic analog", "pubchem_cid": "101587112", "smiles": "Chlorinated pentacyclic pyrido-benzimidazole scaffold", "mechanism": "Enhanced potency in TNBC (MDA-MB-231, IC50 = 0.3 µM) via PI3K/AKT/mTOR suppression."},
    "Ptilomycalin A": {"formula": "C45H80N6O5", "mw": "785.15 g/mol", "class": "Pentacyclic guanidine alkaloid", "source": "Marine sponge (Ptilocaulis / Monanchora sp.)", "pubchem_cid": "6438509", "smiles": "Polycyclic guanidinium core linked to ω-hydroxyhexadecanoic acid ester", "mechanism": "Inhibits Na+/K+-ATPase and triggers apoptosis in TNBC (MDA-MB-231, IC50 = 4.3 µM)."},
    "Monanchomycalin C": {"formula": "C46H82N6O5", "mw": "799.18 g/mol", "class": "Pentacyclic guanidine alkaloid", "source": "Marine sponge (Monanchora pulchra)", "pubchem_cid": "71588154", "smiles": "Pentacyclic guanidinium spermidine-ester derivative", "mechanism": "Induces caspase-dependent apoptosis and JNK activation in MDA-MB-231 (IC50 = 8.2 µM)."},
    "Crambescidin 800": {"formula": "C45H80N6O6", "mw": "801.15 g/mol", "class": "Pentacyclic guanidine alkaloid", "source": "Marine sponge (Crambe crambe / Monanchora sp.)", "pubchem_cid": "6438661", "smiles": "Hydroxylated pentacyclic guanidinium spermidine alkaloid", "mechanism": "Potent sub-micromolar cytotoxicity against TNBC stem-like SUM159PT cells (IC50 = 0.59 µM)."},
    "7-Hydroxyneolamellarin A": {"formula": "C27H21NO6", "mw": "455.46 g/mol", "class": "Pyrrole-derived marine alkaloid", "source": "Marine sponge (Dendrilla nigra)", "pubchem_cid": "24756847", "smiles": "COC1=CC=C(C=C1)C2=C(C3=C(N2CC4=CC=C(C=C4)O)C=C(C(=C3)O)OC)C(=O)OC", "mechanism": "HIF-1α inhibitor active against HR+ T-47D cells (IC50 = 1.9 µM)."},
    "T11": {"formula": "Marine-derived lead", "mw": "N/A", "class": "Marine bioactive fraction/compound", "source": "Marine organism", "pubchem_cid": "N/A", "smiles": "Not specified in source dataset", "mechanism": "Sub-micromolar activity in TNBC SUM159PT cells (IC50 = 0.07 µM)."}
}

@st.cache_data
def load_from_db():
  db_valid = os.path.exists(DB_PATH) and os.path.getsize(DB_PATH) > 0

  if not db_valid and os.path.exists("build_database.py"):
    import build_database

    build_database.main(DB_PATH)
    db_valid = os.path.exists(DB_PATH) and os.path.getsize(DB_PATH) > 0

  if db_valid:
    try:
      conn = sqlite3.connect(DB_PATH)
      ic50_df = pd.read_sql("SELECT * FROM ic50_rankings;", conn)
      summary_df = pd.read_sql("SELECT * FROM cell_line_summary;", conn)
      multi_df = pd.read_sql(
          "SELECT * FROM multi_cell_line_compounds;", conn
      )
      bands_df = pd.read_sql("SELECT * FROM potency_bands;", conn)
      ppi_db_df = pd.read_sql("SELECT * FROM ppi_interactions;", conn)
      conn.close()
    except Exception:
      db_valid = False

  if not db_valid:
    if os.path.exists("ic50_ranking_by_cell_line.csv"):
      ic50_df = pd.read_csv("ic50_ranking_by_cell_line.csv")
    elif os.path.exists("ic50_ranking_by_cell_line-1.csv"):
      ic50_df = pd.read_csv("ic50_ranking_by_cell_line-1.csv")
    else:
      ic50_df = pd.DataFrame(
          [
              ("BT-474", 1, "Araguspongine C", "15.2 ± 2.1", 15.2),
              ("MCF-7", 1, "Meso-araguspongine C", "0.44 ± 0.05 µM", 0.44),
              ("MCF-7", 2, "Manzamine A", "2.86 ± 0.19 µM", 2.86),
              ("MDA-MB-231", 1, "Renieramycin M", "3.8 nM", 0.0038),
              ("T-47D", 2, "Fascaplysin", "5 ± 0.2 µM (T-47D)", 5.0),
          ],
          columns=[
              "cell_line",
              "rank_in_cell_line",
              "compound",
              "reported_text",
              "value_uM",
          ],
      )

    summary_df = (
        pd.read_csv("cell_line_summary.csv")
        if os.path.exists("cell_line_summary.csv")
        else pd.DataFrame()
    )
    multi_df = (
        pd.read_csv("multi_cell_line_compounds.csv")
        if os.path.exists("multi_cell_line_compounds.csv")
        else pd.DataFrame()
    )
    bands_df = (
        pd.read_csv("potency_bands.csv")
        if os.path.exists("potency_bands.csv")
        else pd.DataFrame()
    )
    ppi_db_df = pd.DataFrame(
        columns=["protein_a", "protein_b", "string_score", "context"]
    )

  subtype_map = {
      "MCF-7": "HR+",
      "T-47D": "HR+",
      "BT-474": "HER2+",
      "SKBR3": "HER2+",
      "MDA-MB-231": "TNBC",
      "SUM159PT": "TNBC",
      "MDA-MB-435": "Unclassified",
  }
  ic50_df["subtype"] = ic50_df["cell_line"].map(subtype_map).fillna("Other")
  ic50_df["pIC50"] = (6 - np.log10(ic50_df["value_uM"])).round(2)
  return ic50_df, summary_df, multi_df, bands_df, ppi_db_df


@st.cache_data
def fetch_string_ppi(genes, min_score=0.700, db_edges=None):
    url = "https://string-db.org/api/json/network"
    params = {"identifiers": "%0d".join(genes), "species": 9606, "required_score": int(min_score * 1000), "caller_identity": "marine_bc_app"}
    try:
        resp = requests.get(url, params=params, timeout=8)
        if resp.status_code == 200:
            edges = [(i["preferredName_A"], i["preferredName_B"], round(float(i["score"]), 3)) for i in resp.json() if i.get("preferredName_A") in genes and i.get("preferredName_B") in genes and float(i.get("score", 0)) >= min_score]
            if edges:
                return edges, "Live STRING v12 API"
    except Exception:
        pass
    if db_edges is not None:
        filtered = [(r["protein_a"], r["protein_b"], float(r["string_score"])) for _, r in db_edges.iterrows() if r["protein_a"] in genes and r["protein_b"] in genes and float(r["string_score"]) >= min_score]
        return filtered, "SQLite Database (ppi_interactions table)"
    return [], "No edges"

ic50_df, summary_df, multi_df, bands_df, ppi_db_df = load_from_db()

st.sidebar.title("Navigation")
page = st.sidebar.radio("Select Module:", [
    "1. Compound Search & Chemical Info",
    "2. Anticancer Activity Explorer",
    "3. Compare Selected Compounds",
    "4. Study Visualizations & Gallery",
    "5. Protein–Protein Interaction (PPI) Network",
    "6. SQLite Database Inspector"
])
st.sidebar.success(f"Connected: {DB_PATH}")

if page == "1. Compound Search & Chemical Info":
    st.title("Marine Compound Search & Chemical Profiles")
    all_compounds = sorted(ic50_df["compound"].unique())
    col1, col2 = st.columns([1, 2])
    with col1:
        q = st.text_input("Filter by name:", "")
        filtered_list = [c for c in all_compounds if q.lower() in c.lower()]
        selected = st.selectbox("Select Marine Compound:", filtered_list if filtered_list else all_compounds)
    chem = CHEMICAL_DB.get(selected, {})
    comp_act = ic50_df[ic50_df["compound"] == selected]
    with col2:
        st.subheader(selected)
        c1, c2, c3 = st.columns(3)
        c1.metric("Formula", chem.get("formula", "N/A"))
        c2.metric("Molecular Weight", chem.get("mw", "N/A"))
        c3.metric("Cell Lines Tested", len(comp_act))
        st.markdown(f"**Chemical Class:** {chem.get('class', 'Marine Natural Product')}")
        st.markdown(f"**Biological Source:** {chem.get('source', 'Marine organism')}")
        st.markdown(f"**SMILES:** `{chem.get('smiles', 'N/A')}`")
        st.markdown(f"**Reported Mechanism:** {chem.get('mechanism', 'N/A')}")
        if chem.get("pubchem_cid") and chem.get("pubchem_cid") != "N/A":
            st.markdown(f"[Open in PubChem (CID: {chem['pubchem_cid']})](https://pubchem.ncbi.nlm.nih.gov/compound/{chem['pubchem_cid']})")
    st.dataframe(comp_act[["cell_line", "subtype", "reported_text", "value_uM", "pIC50", "rank_in_cell_line"]], use_container_width=True)

elif page == "2. Anticancer Activity Explorer":
    st.title("Anticancer Activity Explorer")
    f1, f2, f3 = st.columns(3)
    subtypes = f1.multiselect("Subtype:", sorted(ic50_df["subtype"].unique()), default=sorted(ic50_df["subtype"].unique()))
    cell_lines = f2.multiselect("Cell Line:", sorted(ic50_df["cell_line"].unique()), default=sorted(ic50_df["cell_line"].unique()))
    max_ic50 = f3.slider("Max IC50 (µM):", 0.005, 50.0, 50.0, 0.5)
    filtered = ic50_df[(ic50_df["subtype"].isin(subtypes)) & (ic50_df["cell_line"].isin(cell_lines)) & (ic50_df["value_uM"] <= max_ic50)].sort_values(by="value_uM")
    st.dataframe(filtered[["compound", "cell_line", "subtype", "reported_text", "value_uM", "pIC50", "rank_in_cell_line"]], use_container_width=True)
    col_a, col_b = st.columns(2)
    with col_a:
        st.subheader("Cell Line Summary")
        st.dataframe(summary_df, use_container_width=True)
    with col_b:
        st.subheader("Potency Distribution Bands")
        st.dataframe(bands_df, use_container_width=True)

elif page == "3. Compare Selected Compounds":
    st.title("Head-to-Head Compound Comparison")
    all_compounds = sorted(ic50_df["compound"].unique())
    defaults = [c for c in ["Araguspongine C", "Manzamine A", "Renieramycin M", "Fascaplysin"] if c in all_compounds]
    selected_comps = st.multiselect("Select 2 or more compounds:", all_compounds, default=defaults)
    if len(selected_comps) >= 2:
        comp_df = ic50_df[ic50_df["compound"].isin(selected_comps)]
        pivot = comp_df.pivot_table(index="compound", columns="cell_line", values="value_uM", aggfunc="min")
        st.subheader("IC50 Matrix (µM — lower is more potent)")
        st.dataframe(pivot.style.format("{:.4f}", na_rep="—"), use_container_width=True)
        fig, ax = plt.subplots(figsize=(9, 4.5))
        for comp in selected_comps:
            sub = comp_df[comp_df["compound"] == comp]
            ax.bar(sub["cell_line"] + " (" + sub["compound"] + ")", sub["pIC50"], label=comp)
        ax.set_ylabel("pIC50 (6 - log10 IC50 µM)")
        plt.xticks(rotation=35, ha="right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.tight_layout()
        st.pyplot(fig)
    else:
        st.warning("Please select at least 2 compounds.")

elif page == "4. Study Visualizations & Gallery":
    st.title("Quantitative Summary Graphs")
    charts = [
        ("fig1_ic50_by_cell_line.png", "Figure 1: Reported IC50 / pIC50 Ranking Within Major Cell Lines"),
        ("fig2_araguspongine_c_cell_lines.png", "Figure 2: Pan-Subtype Activity of Araguspongine C"),
        ("fig3_cell_line_distribution.png", "Figure 3: Compound Screening Distribution by Cell Line and Subtype"),
        ("fig5_multi_cell_line_compounds.png", "Figure 4: Multi-Cell-Line Marine Compounds")
    ]
    found = False
    for fname, caption in charts:
        if os.path.exists(fname):
            found = True
            st.subheader(caption)
            st.image(fname, use_container_width=True)
    if not found:
        fig, ax = plt.subplots(figsize=(8, 4))
        counts = ic50_df["cell_line"].value_counts()
        ax.barh(counts.index, counts.values, color="#2980B9")
        ax.set_xlabel("Reported Assays")
        st.pyplot(fig)

elif page == "5. Protein–Protein Interaction (PPI) Network":
    st.title("Stage 7: Breast Cancer Protein–Protein Interaction (PPI) Network")
    st.warning("**Important Scientific Distinction:** STRING maps **functional and physical protein–protein associations**. It **does not** prove that marine compounds directly bind to these proteins.")
    default_genes = ["TP53", "BRCA1", "ESR1", "ERBB2", "EGFR", "PIK3CA", "AKT1", "MTOR", "CDK4", "BCL2", "MDM2"]
    col1, col2 = st.columns([2, 1])
    selected_genes = col1.multiselect("Select Breast Cancer Proteins:", sorted(list(set(default_genes + ["MYC", "CCND1", "RB1", "PTEN", "PARP1", "VEGFA", "CASP3"]))), default=default_genes)
    min_score = col2.slider("Minimum STRING Confidence:", 0.40, 0.99, 0.70, 0.05)
    if len(selected_genes) >= 3:
        edges, src = fetch_string_ppi(selected_genes, min_score=min_score, db_edges=ppi_db_df)
        st.caption(f"Source: **{src}** | Nodes: **{len(selected_genes)}** | Edges: **{len(edges)}**")
        G = nx.Graph()
        G.add_nodes_from(selected_genes)
        for u, v, w in edges:
            G.add_edge(u, v, weight=w)
        deg_c = nx.degree_centrality(G)
        bet_c = nx.betweenness_centrality(G)
        eig_c = nx.eigenvector_centrality(G, max_iter=500) if len(edges) > 0 else {n: 0 for n in G.nodes()}
        metrics_df = pd.DataFrame({"Protein": list(G.nodes()), "Degree": [G.degree(n) for n in G.nodes()], "Degree Centrality": [round(deg_c[n], 3) for n in G.nodes()], "Betweenness": [round(bet_c[n], 3) for n in G.nodes()], "Eigenvector": [round(eig_c[n], 3) for n in G.nodes()]}).sort_values(by=["Degree", "Eigenvector"], ascending=False)
        left, right = st.columns([1.4, 1])
        with left:
            fig, ax = plt.subplots(figsize=(8.5, 6.5))
            pos = nx.spring_layout(G, k=1.2, seed=42)
            nx.draw_networkx_edges(G, pos, width=[1.5 + (G[u][v]["weight"] - min_score) * 8 for u, v in G.edges()], alpha=0.45, edge_color="#7F8C8D", ax=ax)
            nx.draw_networkx_nodes(G, pos, node_color="#2980B9", node_size=[600 + G.degree(n) * 280 for n in G.nodes()], edgecolors="white", linewidths=2, ax=ax)
            nx.draw_networkx_labels(G, pos, font_size=9.5, font_weight="bold", font_color="white", ax=ax)
            nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}, font_size=7, ax=ax)
            ax.axis("off")
            st.pyplot(fig)
        with right:
            st.subheader("NetworkX Centrality Metrics")
            st.dataframe(metrics_df, use_container_width=True)
            edges_df = pd.DataFrame(edges, columns=["Protein A", "Protein B", "STRING Score"]).sort_values(by="STRING Score", ascending=False)
            st.download_button("Download Edge List (CSV)", data=edges_df.to_csv(index=False).encode("utf-8"), file_name="string_ppi_edges.csv", mime="text/csv")

elif page == "6. SQLite Database Inspector":
  st.title("SQLite Database Inspector (marine_bioactives.db)")
  st.markdown(
      "Inspect any table directly or run live SQL `SELECT` queries for"
      " examiners."
  )


  def ensure_db_populated(db_path):
    conn = sqlite3.connect(db_path)
    existing = pd.read_sql(
        "SELECT name FROM sqlite_master WHERE type='table';", conn
    )["name"].tolist()
    if existing:
      conn.close()
      return

    cur = conn.cursor()
    cur.executescript("""
            DROP TABLE IF EXISTS ic50_rankings;
            DROP TABLE IF EXISTS cell_line_summary;
            DROP TABLE IF EXISTS multi_cell_line_compounds;
            DROP TABLE IF EXISTS potency_bands;
            DROP TABLE IF EXISTS ppi_interactions;
        """)

    ic50_fallback = pd.DataFrame(
        [
            ("BT-474", 1, "Araguspongine C", "15.2 ± 2.1", 15.2),
            ("MCF-7", 1, "Meso-araguspongine C", "0.44 ± 0.05 µM", 0.44),
            ("MCF-7", 2, "Manzamine A", "2.86 ± 0.19 µM", 2.86),
            ("MCF-7", 3, "Araguspongine O", "5.32 ± 0.67 µM", 5.32),
            ("MCF-7", 4, "Araguspongine P", "5.68 ± 0.89 µM", 5.68),
            ("MCF-7", 5, "Araguspongine N", "7.36 ± 1.16 µM", 7.36),
            ("MCF-7", 6, "Araguspongine A", "7.82 ± 0.53 µM", 7.82),
            ("MCF-7", 7, "Araguspongine C", "8.5 ± 1.6", 8.5),
            ("MCF-7", 8, "Araguspongine E", "24.85 ± 1.2 µM", 24.85),
            ("MCF-7", 8, "Araguspongine L", "24.85 ± 0.91 µM", 24.85),
            ("MDA-MB-231", 1, "Renieramycin M", "3.8 nM", 0.0038),
            ("MDA-MB-231", 2, "4-Chlorofascaplysin", "0.3 µM", 0.3),
            ("MDA-MB-231", 3, "Ptilomycalin A", "4.3 µM", 4.3),
            ("MDA-MB-231", 4, "Manzamine A", "7.87 ± 0.30 µM", 7.87),
            ("MDA-MB-231", 5, "Monanchomycalin C", "8.2 µM", 8.2),
            ("MDA-MB-231", 6, "Araguspongine C", "10.1 ± 2.3", 10.1),
            ("MDA-MB-435", 1, "Renieramycin M", "6.3 ± 0.1 nM", 0.0063),
            ("SKBR3", 1, "Araguspongine C", "18.3 ± 2.5", 18.3),
            ("SUM159PT", 1, "T11", "0.07 ± 0.01 µM", 0.07),
            ("SUM159PT", 2, "Crambescidin 800", "0.59 ± 0.08 µM", 0.59),
            ("T-47D", 1, "7-Hydroxyneolamellarin A", "1.9 µM", 1.9),
            ("T-47D", 2, "Fascaplysin", "5 ± 0.2 µM (T-47D)", 5.0),
            ("T-47D", 3, "Araguspongine C", "46.1 ± 4.8 µM", 46.1),
        ],
        columns=[
            "cell_line",
            "rank_in_cell_line",
            "compound",
            "reported_text",
            "value_uM",
        ],
    )
    ic50_fallback.to_sql("ic50_rankings", conn, if_exists="replace", index=False)

    summary_fallback = pd.DataFrame(
        [
            ("MCF-7", "HR+", 9, 0.44, 7.36, 24.85),
            ("MDA-MB-231", "TNBC", 6, 0.0038, 6.085, 10.1),
            ("T-47D", "HR+", 3, 1.9, 5.0, 46.1),
            ("SUM159PT", "TNBC", 2, 0.07, 0.33, 0.59),
            ("BT-474", "HER2+", 1, 15.2, 15.2, 15.2),
            ("SKBR3", "HER2+", 1, 18.3, 18.3, 18.3),
            (
                "MDA-MB-435",
                "not classified in source sheet",
                1,
                0.0063,
                0.0063,
                0.0063,
            ),
        ],
        columns=[
            "cell_line",
            "subtype",
            "n_compounds",
            "min_uM",
            "median_uM",
            "max_uM",
        ],
    )
    summary_fallback.to_sql(
        "cell_line_summary", conn, if_exists="replace", index=False
    )

    multi_fallback = pd.DataFrame(
        [
            (
                "Araguspongine C",
                5,
                "BT-474, MCF-7, MDA-MB-231, SKBR3, T-47D",
                "HER2+, HR+, TNBC",
            ),
            (
                "Monanchoradin A",
                3,
                "MDA-MB-231, MDA-MB-435, MDA-MB-468",
                "TNBC, not classified",
            ),
            ("Fascaplysin", 2, "MCF-7, T-47D", "HR+"),
            ("Manzamine A", 2, "MCF-7, MDA-MB-231", "HR+, TNBC"),
            (
                "Renieramycin M",
                2,
                "MDA-MB-231, MDA-MB-435",
                "TNBC, not classified",
            ),
        ],
        columns=["compound", "n_cell_lines", "cell_lines", "subtypes"],
    )
    multi_fallback.to_sql(
        "multi_cell_line_compounds", conn, if_exists="replace", index=False
    )

    bands_fallback = pd.DataFrame(
        [
            ("BT-474", 0, 0, 0, 1),
            ("MCF-7", 0, 1, 6, 2),
            ("MDA-MB-231", 1, 1, 3, 1),
            ("MDA-MB-435", 1, 0, 0, 0),
            ("SKBR3", 0, 0, 0, 1),
            ("SUM159PT", 1, 1, 0, 0),
            ("T-47D", 0, 0, 2, 1),
        ],
        columns=["cell_line", "< 0.1 µM", "0.1–1 µM", "1–10 µM", "> 10 µM"],
    )
    bands_fallback.to_sql(
        "potency_bands", conn, if_exists="replace", index=False
    )

    ppi_fallback = pd.DataFrame(
        [
            ("TP53", "BRCA1", 0.999, "DNA damage & tumor suppression"),
            ("TP53", "MDM2", 0.999, "p53 autoregulation / ubiquitination"),
            ("TP53", "CDK4", 0.965, "Cell cycle checkpoint cross-talk"),
            ("TP53", "BCL2", 0.992, "Apoptosis regulation"),
            ("TP53", "AKT1", 0.986, "Survival vs. stress signaling"),
            ("TP53", "ESR1", 0.978, "Hormone & tumor suppressor signaling"),
            ("TP53", "MTOR", 0.912, "Stress & metabolic signaling"),
            ("BRCA1", "ESR1", 0.985, "DNA repair & ER transcriptional control"),
            ("BRCA1", "AKT1", 0.915, "Repair & survival pathway link"),
            ("BRCA1", "MDM2", 0.941, "Genome stability network"),
            ("ERBB2", "EGFR", 0.999, "RTK heterodimerization / signaling"),
            ("ERBB2", "PIK3CA", 0.994, "HER2 -> PI3K activation"),
            ("ERBB2", "AKT1", 0.991, "HER2 -> AKT survival cascade"),
            ("ERBB2", "ESR1", 0.976, "HER2 / ER cross-talk"),
            ("ERBB2", "MTOR", 0.934, "Growth & translational control"),
            ("EGFR", "PIK3CA", 0.992, "EGFR -> PI3K signaling"),
            ("EGFR", "AKT1", 0.995, "EGFR -> AKT survival axis"),
            ("EGFR", "MTOR", 0.945, "RTK -> mTOR axis"),
            ("PIK3CA", "AKT1", 0.999, "Core PI3K/AKT axis"),
            ("PIK3CA", "MTOR", 0.996, "PI3K/mTOR pathway"),
            ("PIK3CA", "ESR1", 0.918, "PI3K / ER pathway interaction"),
            ("AKT1", "MTOR", 0.999, "Central AKT/mTOR survival hub"),
            ("AKT1", "BCL2", 0.988, "Anti-apoptotic signaling"),
            ("AKT1", "ESR1", 0.984, "AKT phosphorylation / ER link"),
            ("AKT1", "MDM2", 0.993, "AKT-mediated MDM2 activation"),
            ("AKT1", "CDK4", 0.910, "Cell cycle & survival link"),
            ("MTOR", "BCL2", 0.905, "Autophagy & apoptosis balance"),
            ("ESR1", "BCL2", 0.974, "Estrogen-regulated survival"),
            ("ESR1", "CDK4", 0.948, "ER -> Cyclin D1/CDK4 cell cycle"),
            ("ESR1", "MDM2", 0.952, "ER & MDM2 interaction"),
            ("CDK4", "BCL2", 0.875, "Cell cycle & survival cross-talk"),
            ("CDK4", "MDM2", 0.920, "Cell cycle & p53 pathway link"),
        ],
        columns=["protein_a", "protein_b", "string_score", "context"],
    )
    ppi_fallback.to_sql(
        "ppi_interactions", conn, if_exists="replace", index=False
    )
    conn.commit()
    conn.close()


  ensure_db_populated(DB_PATH)
  conn = sqlite3.connect(DB_PATH)
  tables = pd.read_sql(
      "SELECT name FROM sqlite_master WHERE type='table';", conn
  )["name"].tolist()

  if tables:
    chosen_table = st.selectbox("Select Database Table:", tables)
    if chosen_table:
      st.dataframe(
          pd.read_sql(f'SELECT * FROM "{chosen_table}";', conn),
          use_container_width=True,
      )

  st.subheader("Run Custom SQL Query")
  default_sql = (
      "SELECT compound, cell_line, reported_text, value_uM FROM ic50_rankings"
      " ORDER BY value_uM ASC LIMIT 10;"
  )
  user_sql = st.text_area("SQL Query (SELECT only):", default_sql, height=90)
  if st.button("Execute SQL"):
    if user_sql.strip().lower().startswith("select"):
      try:
        st.dataframe(pd.read_sql(user_sql, conn), use_container_width=True)
      except Exception as e:
        st.error(f"Query error: {e}")
    else:
      st.error("Only SELECT queries are permitted in the web inspector.")
  conn.close()
    
