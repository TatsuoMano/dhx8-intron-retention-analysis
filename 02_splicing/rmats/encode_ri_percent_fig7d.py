#!/usr/bin/env python3
"""Fig. 7D ENCODE rMATS RI proportion (written for the code release, 2026-10-20).

For each ENCODE knockdown, counts rMATS JCEC events with FDR < 0.05 and |IncLevelDifference| > 0.1
and reports RI% = RI / (SE + A3SS + A5SS + RI) x 100. MXE is excluded from the denominator (it is
printed for reference, with the percentage that including it would give).

Columns are taken by name. The summary block at the end of run_encode_splicing.sh (meant to write
encode_rbp_kd/splicing/rmats_summary.tsv) filters with the fixed awk columns $22 (as FDR) and $23
(as IncLevelDifference) in every event file. In the SE, A3SS, A5SS and RI files FDR is column 20,
column 22 is the comma-separated IncLevel2 string and IncLevelDifference is column 23; in the MXE
file FDR is column 22 and IncLevelDifference column 25 ($23 is IncLevel1). That block therefore
does not count significant events, and it also adds MXE to the denominator. It is not the source
of the Fig. 7D values (rmats_summary.tsv is not on the server; RBFOX2 and RBM39 were run by
run_remaining.sh, which has no summary block). run_encode_splicing.sh is kept unchanged as it was
run; this script is the computation behind the values written in
07_figure_panels/local/_fig7b_encode_v2.py.

Result of the 2026-10-20 run: the SE and RI counts and the RI% of all six factors reproduce the
panel. The panel and the text first printed 35.0 (PRPF8) and 24.4 (EFTUD2); the computed values are
1306/3726 = 35.05% and 546/2233 = 24.45%, and the panel (_fig7b_encode_v2.py) and the Results were
corrected to 35.1 and 24.5 on 2026-10-20. The last column (A3SS and A5SS events with
|IncLevelDifference| exactly 0.1 kept) gives the earlier 35.0 and 24.4; it is printed only as a
candidate explanation of the earlier values, not as a recovered record of how they were computed.

Run on the server (read-only, streamed):
    ssh ncnp-linux 'python3 -' < 02_splicing/rmats/encode_ri_percent_fig7d.py
"""
import pandas as pd

BASE = "/mnt/HDD2/unified_results/encode_rbp_kd/splicing"
FACTORS = ["PRPF8", "EFTUD2", "SNRNP200", "PRPF6", "RBFOX2", "RBM39"]
# values written in _fig7b_encode_v2.py (SE, RI, RI%)
PANEL = {"PRPF8": (1756, 1306, 35.1), "EFTUD2": (1173, 546, 24.5), "SNRNP200": (1900, 729, 20.9),
         "PRPF6": (737, 308, 22.4), "RBFOX2": (648, 250, 22.1), "RBM39": (5016, 1601, 19.3)}


def n_sig(path, ge=False):
    df = pd.read_csv(path, sep="\t")
    a = df["IncLevelDifference"].abs()
    return int(((df["FDR"] < 0.05) & ((a >= 0.1) if ge else (a > 0.1))).sum())


print(f"{'factor':<9} {'SE':>5} {'A3SS':>5} {'A5SS':>5} {'RI':>5} {'MXE':>5} {'RI%':>6} {'RI%+MXE':>8}  panel(SE,RI,RI%)   match  RI%(A3SS/A5SS >= 0.1)")
for f in FACTORS:
    d = f"{BASE}/rmats_{f}"
    n = {ev: n_sig(f"{d}/{ev}.MATS.JCEC.txt") for ev in ["SE", "A3SS", "A5SS", "RI", "MXE"]}
    den = n["SE"] + n["A3SS"] + n["A5SS"] + n["RI"]
    pct = 100 * n["RI"] / den
    pct_mxe = 100 * n["RI"] / (den + n["MXE"])
    # diagnostic only: A3SS and A5SS counted with |IncLevelDifference| >= 0.1 (boundary events kept)
    ge = n["SE"] + n_sig(f"{d}/A3SS.MATS.JCEC.txt", True) + n_sig(f"{d}/A5SS.MATS.JCEC.txt", True) + n["RI"]
    pct_ge = 100 * n["RI"] / ge
    p = PANEL[f]
    ok = (n["SE"], n["RI"], round(pct, 1)) == p
    print(f"{f:<9} {n['SE']:>5} {n['A3SS']:>5} {n['A5SS']:>5} {n['RI']:>5} {n['MXE']:>5} {pct:>6.2f} {pct_mxe:>8.2f}  {str(p):<18} {'OK' if ok else 'DIFF':<5}  {pct_ge:.2f}")
