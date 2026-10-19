#!/usr/bin/env python3
"""RECONSTRUCTED FILE (code release, 2026-10-19).

Brown-Forsythe test of equal variances across the four genotype groups of each 2 x 2 Drosophila
eye experiment (Methods, "Drosophila experiments": "Brown-Forsythe P of approximately 1 x 10^-5"
for the second RNAi line). The original value was computed interactively on 2026-05-29 and no
script was kept. This script recomputes it with scipy.stats.levene(center="median"), the
Brown-Forsythe form, from the same Prism file read by _fig8a_eye_v2.py.
Recomputed on 2026-10-19: line #1 P = 0.51, line #2 P = 9.9e-06.

Usage: python3 brown_forsythe_eye.py [path/to/Fig8a.2 260514_analysis.pzfx]
"""
import sys
import xml.etree.ElementTree as ET
import numpy as np
from scipy import stats

PRISM = sys.argv[1] if len(sys.argv) > 1 else (
    "/Users/tatsuomano/Dropbox/Research/DHX8-Sakai/Manuscript/Figures/"
    "figures_for_final/Fig08/Fig8a.2 260514_analysis.pzfx")


def _tag(t): return t.split("}")[-1]
def _txt(e): return (e.text or "").strip()


root = ET.parse(PRISM).getroot()
pig = []
for tbl in root:
    if _tag(tbl.tag) == "Table" and next((_txt(k) for k in tbl if _tag(k.tag) == "Title"), "") == "pigmantation":
        for k in tbl:
            if _tag(k.tag) == "YColumn":
                for sub in k:
                    if _tag(sub.tag) == "Subcolumn":
                        pig.append(np.array([float(_txt(d)) for d in sub
                                             if _tag(d.tag) == "d" and _txt(d) != ""]))
# column order (as in _fig8a_eye_v2.py):
# 0 GFP#1-T  1 GFP#2-T  2 dDHX8#1-T  3 dDHX8#2-T  4 GFP#1+T  5 GFP#2+T  6 dDHX8#1+T  7 dDHX8#2+T
for name, idx in [("RNAi line #1", [0, 2, 4, 6]), ("RNAi line #2", [1, 3, 5, 7])]:
    groups = [pig[i] for i in idx]
    res = stats.levene(*groups, center="median")
    print(f"{name}: Brown-Forsythe W = {res.statistic:.3f}, P = {res.pvalue:.3g}")
