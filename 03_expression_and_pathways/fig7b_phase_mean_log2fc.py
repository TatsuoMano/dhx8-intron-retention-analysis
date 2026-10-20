#!/usr/bin/env python3
"""Per-phase mean log2FC of the curated core spliceosomal factors (Fig 7B labels).

Fig 7B prints, above each spliceosome-phase NES bar, the mean DESeq2 log2 fold
change (DHX8 KD vs control) of a curated subset of core factors.  The subset is
the 83-factor PHASE_MAP of
Analysis/2026-03-11_phase3_specificity/analysis_splicing_factor_phase.py
(27 E/A, 16 B, 6 Bact/NTC, 2 catalytic, 2 post-catalytic, 30 regulatory),
reproduced verbatim below.  It is a different, smaller set than the GMT gene
sets used for the phase GSEA (spliceosome_phase_genesets.gmt), whose all-member
means are also reported for comparison.

The script reads only the files given on the command line and writes to stdout.
Standard library only, so it can be streamed to the server read-only:

    ssh ncnp-linux 'python3 - /mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv \
        --gmt /mnt/HDD2/unified_results/gsea/data/spliceosome_phase_genesets.gmt' \
        < fig7b_phase_mean_log2fc.py

Each DESeq2 file may be TSV or CSV with columns log2FoldChange, padj and a gene
symbol column (gene / symbol / first column).  Gene symbols are matched
case-insensitively.  Output is in the Fig 7B display order
(E/A, B, Bact/NTC, Regulatory, Post-catalytic, Catalytic).

Result (2026-10-20): the printed labels +0.408/+0.369/+0.342/+0.260/-0.011/
+0.011 reproduce to 3 decimals from the current (v2) deseq2_dhx8_kd_vs_ctrl.tsv
(identical to sg_analysis/gsea_dhx8/data/deseq2_results.csv).  The v1 table
(_v1_corrupt_snapshot_20260706/downstream_results/deseq2_dhx8_kd_vs_ctrl.tsv)
gives +0.400/+0.361/+0.354/+0.261/-0.011/-0.013.  The all-member GMT means
(+0.30 to +0.38 in every phase) are a different quantity.
"""
import argparse
import csv
import math
import os
import statistics

CURATED = {
    "E/A": ["Snrnp70", "Snrpa", "Snrpc", "Prpf39", "Sf3a1", "Sf3a2", "Sf3a3",
            "Sf3b1", "Sf3b2", "Sf3b3", "Sf3b4", "Sf3b5", "Sf3b6", "Snrpa1",
            "Snrpb2", "Snrpd1", "Snrpd2", "Snrpd3", "Snrpe", "Snrpf", "Snrpg",
            "Snrpb", "Ddx46", "Ddx39b", "Snrnp25", "Snrnp35", "Snrnp48"],
    "B": ["Prpf3", "Prpf4", "Prpf6", "Prpf31", "Prpf8", "Snrnp200", "Snrnp40",
          "Eftud2", "Snrnp27", "Lsm2", "Lsm3", "Lsm4", "Lsm5", "Lsm6", "Lsm7",
          "Lsm8"],
    "Bact/NTC": ["Prpf19", "Cdc5l", "Bcas2", "Plrg1", "Cwc15", "Cwc27"],
    "Regulatory": ["Srsf1", "Srsf2", "Srsf3", "Srsf4", "Srsf5", "Srsf6",
                   "Srsf7", "Srsf9", "Srsf10", "Srsf11", "Hnrnpa1",
                   "Hnrnpa2b1", "Hnrnpc", "Hnrnpf", "Hnrnph1", "Hnrnpk",
                   "Hnrnpl", "Hnrnpm", "Hnrnpu", "Ptbp1", "Ptbp2", "Ddx5",
                   "Ddx17", "Dhx36", "Rbm25", "Rbm27", "Rbm15", "Rbm17",
                   "Sfpq", "Nono"],
    "Post-catalytic": ["Dhx8", "Dhx15"],
    "Catalytic": ["Dhx16", "Dhx38"],
}
DISPLAY = ["E/A", "B", "Bact/NTC", "Regulatory", "Post-catalytic", "Catalytic"]
GMT_NAME = {
    "E/A": "SPLICEOSOME_PHASE1_EA_COMPLEX",
    "B": "SPLICEOSOME_PHASE2_B_COMPLEX",
    "Bact/NTC": "SPLICEOSOME_PHASE3_BACT_NTC",
    "Regulatory": "SPLICEOSOME_PHASE6_REGULATORY",
    "Post-catalytic": "SPLICEOSOME_PHASE5_POST_CATALYTIC",
    "Catalytic": "SPLICEOSOME_PHASE4_CATALYTIC",
}
# values hard-coded in 07_figure_panels/local/_fig7a_phase_response_v4.py
PRINTED = [0.408, 0.369, 0.342, 0.260, -0.011, 0.011]


def load_deseq2(path):
    delim = "," if path.endswith(".csv") else "\t"
    lfc, padj = {}, {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter=delim)
        head = next(rd)
        low = [h.strip().strip('"').lower() for h in head]
        i_l = low.index("log2foldchange")
        i_p = low.index("padj") if "padj" in low else None
        i_g = low.index("gene") if "gene" in low else (
            low.index("symbol") if "symbol" in low else 0)
        for row in rd:
            if len(row) < len(head):  # R write.csv with unnamed rowname column
                row = [""] * (len(head) - len(row)) + row
            g = row[i_g].strip().strip('"')
            if not g or g == "NA":
                continue
            try:
                v = float(row[i_l])
            except ValueError:
                continue
            key = g.lower()
            if key in lfc:  # keep the first occurrence
                continue
            lfc[key] = v
            try:
                padj[key] = float(row[i_p]) if i_p is not None else math.nan
            except ValueError:
                padj[key] = math.nan
    return lfc, padj


def load_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                sets[p[0]] = [g for g in p[2:] if g]
    return sets


def summarise(lfc, padj, genes):
    vals = [(g, lfc[g.lower()], padj.get(g.lower(), math.nan))
            for g in genes if g.lower() in lfc]
    missing = [g for g in genes if g.lower() not in lfc]
    xs = [v for _, v, _ in vals]
    mean = statistics.fmean(xs) if xs else math.nan
    med = statistics.median(xs) if xs else math.nan
    return vals, missing, mean, med


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("deseq2", nargs="+")
    ap.add_argument("--gmt", default=None)
    ap.add_argument("--genes", action="store_true",
                    help="print per-gene log2FC for the small phases")
    a = ap.parse_args()
    gmt = load_gmt(a.gmt) if a.gmt and os.path.exists(a.gmt) else None
    print("printed in Fig 7B: " + "  ".join(
        f"{p}={v:+.3f}" for p, v in zip(DISPLAY, PRINTED)))
    for path in a.deseq2:
        lfc, padj = load_deseq2(path)
        st = os.stat(path)
        print(f"\n=== {path}  (genes={len(lfc)}, bytes={st.st_size}) ===")
        print(f"{'phase':15s} {'n':>5s} {'curated mean':>13s} {'median':>8s}"
              f" {'printed':>8s} {'diff':>7s}" + (
                  f" {'GMT n':>7s} {'GMT mean':>9s}" if gmt else ""))
        for p, pv in zip(DISPLAY, PRINTED):
            vals, miss, m, md = summarise(lfc, padj, CURATED[p])
            line = (f"{p:15s} {len(vals):2d}/{len(CURATED[p]):<2d} {m:+13.4f}"
                    f" {md:+8.4f} {pv:+8.3f} {m - pv:+7.3f}")
            if gmt:
                gs = gmt.get(GMT_NAME[p], [])
                gv, _, gm, _ = summarise(lfc, padj, gs)
                line += f" {len(gv):3d}/{len(gs):<3d} {gm:+9.4f}"
            print(line + (f"  missing={miss}" if miss else ""))
            if a.genes and len(CURATED[p]) <= 6:
                for g, v, q in vals:
                    print(f"      {g:10s} log2FC={v:+.4f} padj={q:.3g}")


if __name__ == "__main__":
    main()
