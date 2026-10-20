#!/usr/bin/env python3
"""Supplementary Table S2: spliceosome phase classification with DHX8 KD log2FC and padj.

RECONSTRUCTED FILE (code release, 2026-10-20). The table was first written on 2026-06-03 by
/tmp/build_s3.py on the server (the table was then numbered S3). That file no longer exists; its
body was recovered from the session log ~/.claude-mem/logs/claude-mem-2026-06-03.log (22:45 JST).
The logic below is that body unchanged (phase labels, the direction rule, the sort order and the
number formats); only the hard-coded /tmp output paths were replaced by command-line arguments
and the xlsx step was made optional.

One row per gene-phase membership of the six phase gene sets in spliceosome_phase_genesets.gmt
(the SPLICEOSOME_ALL_CORE union is skipped), joined with the DESeq2 result of DHX8 KD vs control.
direction = UP / DOWN when padj < 0.05 (sign of log2FC), ns otherwise, NA when the gene has no
DESeq2 log2FC.

Fig7B_subset (column added 2026-10-20, revision 3): "yes" when the gene-phase row belongs to the
curated 83-factor subset whose per-phase mean log2FC is printed above the NES bars of Fig 7B, "no"
otherwise. The subset is PHASE_MAP of
Analysis/2026-03-11_phase3_specificity/analysis_splicing_factor_phase.py (27 E/A, 16 B, 6 Bact/NTC,
2 catalytic, 2 post-catalytic, 30 regulatory), the same list as CURATED in
fig7b_phase_mean_log2fc.py. It is copied below because the script is streamed to the server and
cannot import local files. Every subset gene is a member of the GMT gene set of the same phase; the
script stops if that is not so.

Run on the server (read-only, streamed; the table goes to stdout):
    ssh ncnp-linux 'python3 - --stdout' < 03_expression_and_pathways/build_table_s2.py > Supplementary_Table_S2.tsv
Options: --gmt PATH  --deseq PATH  --out PATH (tsv)  --xlsx PATH
"""
import argparse
import csv
import io
import sys
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--gmt", default="/mnt/HDD2/unified_results/gsea/spliceosome_phase_genesets.gmt")
ap.add_argument("--deseq", default="/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv")
ap.add_argument("--out", default=None)
ap.add_argument("--xlsx", default=None)
ap.add_argument("--stdout", action="store_true")
a = ap.parse_args()

# Fig 7B curated subset (PHASE_MAP of analysis_splicing_factor_phase.py), by GMT phase number
FIG7B_SUBSET = {
    1: ["Snrnp70", "Snrpa", "Snrpc", "Prpf39", "Sf3a1", "Sf3a2", "Sf3a3",
        "Sf3b1", "Sf3b2", "Sf3b3", "Sf3b4", "Sf3b5", "Sf3b6", "Snrpa1",
        "Snrpb2", "Snrpd1", "Snrpd2", "Snrpd3", "Snrpe", "Snrpf", "Snrpg",
        "Snrpb", "Ddx46", "Ddx39b", "Snrnp25", "Snrnp35", "Snrnp48"],
    2: ["Prpf3", "Prpf4", "Prpf6", "Prpf31", "Prpf8", "Snrnp200", "Snrnp40",
        "Eftud2", "Snrnp27", "Lsm2", "Lsm3", "Lsm4", "Lsm5", "Lsm6", "Lsm7",
        "Lsm8"],
    3: ["Prpf19", "Cdc5l", "Bcas2", "Plrg1", "Cwc15", "Cwc27"],
    4: ["Dhx16", "Dhx38"],
    5: ["Dhx8", "Dhx15"],
    6: ["Srsf1", "Srsf2", "Srsf3", "Srsf4", "Srsf5", "Srsf6", "Srsf7",
        "Srsf9", "Srsf10", "Srsf11", "Hnrnpa1", "Hnrnpa2b1", "Hnrnpc",
        "Hnrnpf", "Hnrnph1", "Hnrnpk", "Hnrnpl", "Hnrnpm", "Hnrnpu", "Ptbp1",
        "Ptbp2", "Ddx5", "Ddx17", "Dhx36", "Rbm25", "Rbm27", "Rbm15", "Rbm17",
        "Sfpq", "Nono"],
}
SUBSET_KEYS = {(g, pn) for pn, gs in FIG7B_SUBSET.items() for g in gs}
assert len(SUBSET_KEYS) == 83

# DESeq2 log2FC/padj by gene
lfc = {}
with open(a.deseq) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        g = row["gene"]
        try:
            lfc[g] = (float(row["log2FoldChange"]), row["padj"])
        except Exception:
            lfc[g] = (None, row["padj"])

phase_clean = {
    "SPLICEOSOME_PHASE1_EA_COMPLEX": (1, "E/A complex"),
    "SPLICEOSOME_PHASE2_B_COMPLEX": (2, "B complex"),
    "SPLICEOSOME_PHASE3_BACT_NTC": (3, "Bact/NTC complex"),
    "SPLICEOSOME_PHASE4_CATALYTIC": (4, "Catalytic (C/C*)"),
    "SPLICEOSOME_PHASE5_POST_CATALYTIC": (5, "Post-catalytic (P/ILS)"),
    "SPLICEOSOME_PHASE6_REGULATORY": (6, "Regulatory"),
}
rows = []
with open(a.gmt) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        name = p[0]
        if name not in phase_clean:
            continue  # skip ALL_CORE union
        pnum, plabel = phase_clean[name]
        desc = p[1]
        for g in p[2:]:
            g = g.strip()
            if not g:
                continue
            l2, padj = lfc.get(g, (None, ""))
            try:
                padjf = float(padj)
            except Exception:
                padjf = None
            if l2 is None:
                direction = "NA"
            elif padjf is not None and padjf < 0.05:
                direction = "UP" if l2 > 0 else "DOWN"
            else:
                direction = "ns"
            rows.append((g, pnum, plabel, desc, l2, padj, direction))

rows.sort(key=lambda x: (x[1], x[0].lower()))
missing = SUBSET_KEYS - {(r[0], r[1]) for r in rows}
if missing:
    sys.exit(f"Fig 7B subset genes not in the GMT phase of the same number: {sorted(missing)}")
buf = io.StringIO()
w = csv.writer(buf, delimiter="\t", lineterminator="\r\n")  # csv default terminator, as in the original
w.writerow(["gene_symbol", "phase_number", "phase", "phase_description",
            "DHX8_KD_log2FC", "DHX8_KD_padj", "direction", "Fig7B_subset"])
for g, pn, pl, desc, l2, padj, d in rows:
    w.writerow([g, pn, pl, desc, "" if l2 is None else f"{l2:.4f}", padj, d,
                "yes" if (g, pn) in SUBSET_KEYS else "no"])
text = buf.getvalue()

if a.out:
    with open(a.out, "w", newline="") as fh:
        fh.write(text)
if a.stdout:
    sys.stdout.write(text)

uniq = len(set(r[0] for r in rows))
log = sys.stderr
print(f"rows(memberships)={len(rows)}  unique_genes={uniq}", file=log)
for pn in sorted(set(r[1] for r in rows)):
    c = [r for r in rows if r[1] == pn]
    print(f"  phase{pn} {c[0][2]}: {len(c)} genes", file=log)
print("direction:", dict(Counter(r[6] for r in rows)), file=log)
sub = [r for r in rows if (r[0], r[1]) in SUBSET_KEYS]
print(f"Fig7B_subset=yes: {len(sub)} rows; direction:", dict(Counter(r[6] for r in sub)), file=log)
if a.xlsx:
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Table S2"
    for line in text.splitlines():
        ws.append(line.split("\t"))
    wb.save(a.xlsx)
    print("xlsx written:", a.xlsx, file=log)
