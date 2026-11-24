import os
import sqlite3
import pandas as pd

def main(db_path: str = "marine_bioactives.db") -> None:
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.executescript("""
        DROP TABLE IF EXISTS ic50_rankings;
        DROP TABLE IF EXISTS cell_line_summary;
        DROP TABLE IF EXISTS multi_cell_line_compounds;
        DROP TABLE IF EXISTS potency_bands;
        DROP TABLE IF EXISTS ppi_interactions;
    """)

    # 1. IC50 Rankings
    ic50_file = "ic50_ranking_by_cell_line.csv" if os.path.exists("ic50_ranking_by_cell_line.csv") else (
        "ic50_ranking_by_cell_line-1.csv" if os.path.exists("ic50_ranking_by_cell_line-1.csv") else None
    )
    if ic50_file:
        ic50_df = pd.read_csv(ic50_file)
    else:
        ic50_df = pd.DataFrame([
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
        ], columns=["cell_line", "rank_in_cell_line", "compound", "reported_text", "value_uM"])
    ic50_df.to_sql("ic50_rankings", conn, if_exists="replace", index=False)

    # 2. Cell Line Summary
    summary_file = "cell_line_summary.csv" if os.path.exists("cell_line_summary.csv") else (
        "cell_line_summary-1.csv" if os.path.exists("cell_line_summary-1.csv") else None
    )
    if summary_file:
        summary_df = pd.read_csv(summary_file)
    else:
        summary_df = pd.DataFrame([
            ("MCF-7", "HR+", 9, 0.44, 7.36, 24.85),
            ("MDA-MB-231", "TNBC", 6, 0.0038, 6.085, 10.1),
            ("T-47D", "HR+", 3, 1.9, 5.0, 46.1),
            ("SUM159PT", "TNBC", 2, 0.07, 0.33, 0.59),
            ("BT-474", "HER2+", 1, 15.2, 15.2, 15.2),
            ("SKBR3", "HER2+", 1, 18.3, 18.3, 18.3),
            ("MDA-MB-435", "not classified in source sheet", 1, 0.0063, 0.0063, 0.0063),
        ], columns=["cell_line", "subtype", "n_compounds", "min_uM", "median_uM", "max_uM"])
    summary_df.to_sql("cell_line_summary", conn, if_exists="replace", index=False)

    # 3. Multi-Cell-Line Compounds
    if os.path.exists("multi_cell_line_compounds.csv"):
        multi_df = pd.read_csv("multi_cell_line_compounds.csv")
    else:
        multi_df = pd.DataFrame([
            ("Araguspongine C", 5, "BT-474, MCF-7, MDA-MB-231, SKBR3, T-47D", "HER2+, HR+, TNBC"),
            ("Monanchoradin A", 3, "MDA-MB-231, MDA-MB-435, MDA-MB-468", "TNBC, not classified"),
            ("Fascaplysin", 2, "MCF-7, T-47D", "HR+"),
            ("Manzamine A", 2, "MCF-7, MDA-MB-231", "HR+, TNBC"),
            ("Renieramycin M", 2, "MDA-MB-231, MDA-MB-435", "TNBC, not classified"),
        ], columns=["compound", "n_cell_lines", "cell_lines", "subtypes"])
    multi_df.to_sql("multi_cell_line_compounds", conn, if_exists="replace", index=False)

    # 4. Potency Bands
    if os.path.exists("potency_bands.csv"):
        bands_df = pd.read_csv("potency_bands.csv")
    else:
        bands_df = pd.DataFrame([
            ("BT-474", 0, 0, 0, 1),
            ("MCF-7", 0, 1, 6, 2),
            ("MDA-MB-231", 1, 1, 3, 1),
            ("MDA-MB-435", 1, 0, 0, 0),
            ("SKBR3", 0, 0, 0, 1),
            ("SUM159PT", 1, 1, 0, 0),
            ("T-47D", 0, 0, 2, 1),
        ], columns=["cell_line", "< 0.1 µM", "0.1–1 µM", "1–10 µM", "> 10 µM"])
    bands_df.to_sql("potency_bands", conn, if_exists="replace", index=False)

    # 5. Stage 7 PPI Interactions
    ppi_edges = [
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
        ("CDK4", "MDM2", 0.920, "Cell cycle & p53 pathway link")
    ]
    ppi_df = pd.DataFrame(ppi_edges, columns=["protein_a", "protein_b", "string_score", "context"])
    ppi_df.to_sql("ppi_interactions", conn, if_exists="replace", index=False)

    conn.commit()
    conn.close()
    print(f"Database created successfully: {db_path}")

if __name__ == "__main__":
    main()
