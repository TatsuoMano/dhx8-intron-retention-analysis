#!/usr/bin/env python3
"""Fig. 4B event counts (written for the code release, 2026-10-20).

Recounts every number in the Fig. 4B table and in the Results text around it from the tool outputs:

  rMATS (JCEC, FDR < 0.05, no |IncLevelDifference| cutoff): SE, A3SS, A5SS, plus MXE and RI
      for DHX8 KD (rmats/dhx8_v2), GSE65973 and GSE126543.
  limma eBayes IR (retained = logFC > 0 and adj_p < 0.05): introns tested, retained introns,
      retained genes, released introns (logFC < 0 and adj_p < 0.05).
  MAJIQ deltapsi (deltapsi_annotated.tsv): significant LSV junctions are those with
      probability_changing >= 0.95 and |dpsi_mean| >= 0.1 (the filter of
      majiq/analysis_deltapsi_integration.py, which writes all_datasets_significant_events.tsv);
      de novo junctions (is_denovo) among them; CE candidates = control PSI < 0.05 and
      disease/KD PSI > 0.2 (raw PSI means); strong CE = CE candidate with disease/KD PSI > 0.5.

The Fig. 4B panel itself is drawn by 07_figure_panels/local/_fig4e_summary_v2.py, which has the
values written in; this script is the code behind those values. It only reads files.

Run on the server (read-only, streamed):
    ssh ncnp-linux 'python3 -' < 02_splicing/count_events_fig4b.py
"""
import pandas as pd

U = "/mnt/HDD2/unified_results"
RMATS = {"DHX8 KD": f"{U}/rmats/dhx8_v2", "GSE65973": f"{U}/rmats/gse65973", "GSE126543": f"{U}/rmats/gse126543"}
LIMMA = {"DHX8 KD": "/mnt/HDD2/TDP43_comparison/validation_20260226/limma_ir_results_dhx8_s2.tsv",
         "GSE65973": f"{U}/limma_ir_results_gse65973.tsv",
         "GSE126543": f"{U}/limma_ir_results_gse126543.tsv"}
# (directory, control-group column prefix, disease/KD column prefix)
MAJIQ = {"DHX8 KD": ("dhx8", "Control", "DHX8_KD"),
         "GSE65973": ("gse65973", "Control", "dNLS_TDP"),
         "GSE126543": ("gse126543", "TDP43_positive", "TDP43_negative")}  # TDP-43(+) = normal

print("== rMATS JCEC, FDR < 0.05 (column found by name; MXE has two extra coordinate columns) ==")
print(f"{'dataset':<10} {'SE':>6} {'A3SS':>6} {'A5SS':>6} {'MXE':>6} {'RI':>6}")
for ds, d in RMATS.items():
    n = {}
    for ev in ["SE", "A3SS", "A5SS", "MXE", "RI"]:
        df = pd.read_csv(f"{d}/{ev}.MATS.JCEC.txt", sep="\t")
        n[ev] = int((df["FDR"] < 0.05).sum())
    print(f"{ds:<10} {n['SE']:>6} {n['A3SS']:>6} {n['A5SS']:>6} {n['MXE']:>6} {n['RI']:>6}")

print("\n== limma eBayes IR (retained: logFC > 0 & adj_p < 0.05) ==")
for ds, f in LIMMA.items():
    df = pd.read_csv(f, sep="\t")
    sig = df["adj_p"] < 0.05
    ret = sig & (df["logFC"] > 0)
    rel = sig & (df["logFC"] < 0)
    n_genes = df.loc[ret, "gene_name"].nunique()
    agree = "" if "retained" not in df.columns else f"  (matches 'retained' column: {bool((df['retained'].astype(str).str.upper() == 'TRUE').sum() == ret.sum())})"
    print(f"{ds:<10} tested={len(df):,}  retained introns={int(ret.sum()):,}  retained genes={n_genes:,}  released introns={int(rel.sum()):,}{agree}")

print("\n== MAJIQ deltapsi_annotated.tsv (P_changing >= 0.95 & |dpsi_mean| >= 0.1; de novo; CE) ==")
for ds, (d, ctl, dis) in MAJIQ.items():
    df = pd.read_csv(f"{U}/majiq/{d}/deltapsi_annotated.tsv", sep="\t", comment="#")
    sig = df[(df["probability_changing"] >= 0.95) & (df["dpsi_mean"].abs() >= 0.1)]
    dn = sig[sig["is_denovo"] == True]
    c, k = dn[f"{ctl}_raw_psi_mean"], dn[f"{dis}_raw_psi_mean"]
    ce = (c < 0.05) & (k > 0.2)
    strong = ce & (k > 0.5)
    sign = 1 if MAJIQ[ds][0] != "gse126543" else -1   # dpsi_mean = group2 - group1; GSE126543 group1 = TDP43_negative
    genes = ", ".join(f"{g} (raw PSI difference {dk - dc:+.2f}, dpsi_mean {sign * dm:+.2f})" for g, dc, dk, dm in
                      zip(dn.loc[strong, "gene_name"], c[strong], k[strong], dn.loc[strong, "dpsi_mean"]))
    print(f"{ds:<10} significant junctions={len(sig):,}  de novo={len(dn)}  CE candidates={int(ce.sum())}  strong CE={int(strong.sum())}")
    if genes:
        print(f"           strong CE: {genes}")
    # strict reading of the Methods wording (P(|dPSI| > 0.2) > 0.95, no |dPSI| floor)
    sig2 = df[df["probability_changing"] > 0.95]
    dn2 = sig2[sig2["is_denovo"] == True]
    ce2 = (dn2[f"{ctl}_raw_psi_mean"] < 0.05) & (dn2[f"{dis}_raw_psi_mean"] > 0.2)
    print(f"           [Methods wording, P > 0.95 only] de novo={len(dn2)}  CE={int(ce2.sum())}  strong={int((ce2 & (dn2[f'{dis}_raw_psi_mean'] > 0.5)).sum())}")

# cross-check against the table the CE figure reads
alls = pd.read_csv(f"{U}/majiq/all_datasets_significant_events.tsv", sep="\t")
print("\n== cross-check: majiq/all_datasets_significant_events.tsv, is_denovo rows per dataset ==")
print(alls[alls["is_denovo"] == True].groupby("dataset").size().to_string())
