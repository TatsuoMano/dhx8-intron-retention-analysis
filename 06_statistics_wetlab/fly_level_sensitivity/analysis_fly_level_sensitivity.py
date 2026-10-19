#!/usr/bin/env python3
"""
Fly-level sensitivity analysis of the Drosophila compound-eye 2 x 2 experiment
(TDP-43 x dDHX8-RNAi), run separately for RNAi line #1 (R2) and line #2 (R3).

The reported tests treat 8 eyes per genotype as independent. The 8 eyes come
from 4 flies (both eyes of each fly), and which two eyes share a fly is not
recorded. Here eyes are averaged into fly means (n = 4 per genotype) under
(a) the stored-order pairing (eyes 1+2, 3+4, 5+6, 7+8),
(b) every possible pairing: 105 per genotype, 105^4 = 121,550,625 joint
    pairings per line. These are enumerated exhaustively, because each test
    depends on a pairing only through the per-genotype variance of the fly means
    (cell means do not change with pairing). A random sample of 20,000 joint
    pairings (seed 20261019) is also reported, as requested,
(c) a coordinate-wise search for the worst case (largest P), one genotype at a
    time, from several starts, checked against the exhaustive maximum.

Tests (all as in Analysis/2026-05-29_drosophila_eye_twoway/analysis_drosophila_eye_twoway.py):
  wj_interaction      single-df Welch-James interaction contrast [1,-1,-1,1]
  anova_interaction   conventional two-way ANOVA interaction (pooled MSE)
  kdT_vs_kdnoT        Welch t, dDHX8-RNAi + TDP-43 vs dDHX8-RNAi alone
                      (the value reported as the simple-effect contrast, 6.5e-6)
  kdT_vs_ctrlT        Welch t, dDHX8-RNAi + TDP-43 vs TDP-43 alone
  simple_both_max     max of the two above (both must hold for "versus each")
  tdp_alone           Welch t, GFP-RNAi -TDP vs +TDP
  kd_alone            Welch t, GFP-RNAi vs dDHX8-RNAi, no TDP-43

Data: Fig8a.2 260514_analysis.pzfx ('pigmantation' table), cross-checked
against 260508-2.xlsx (pigmentation = %Area(row 2k) / %Area(row 2k-1)).
"""
import os
import sys
import json
import itertools
import xml.etree.ElementTree as ET

import numpy as np
from scipy import stats
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = "/Users/tatsuomano/Dropbox/Research/DHX8-Sakai"
PRISM = os.path.join(PROJ, "Manuscript/Figures/figures_for_final/Fig08/Fig8a.2 260514_analysis.pzfx")
PRISM_COPY = os.path.join(PROJ, "Manuscript/Collaborators/Fujikake/260514_analysis.pzfx")
XLSX = os.path.join(PROJ, "Manuscript/Collaborators/Fujikake/260508-2.xlsx")
SEED = 20261019
N_RANDOM = 20000
N_STARTS = 30

CELLS = ["ctrl_noT", "ctrl_T", "kd_noT", "kd_T"]   # order of the WJ contrasts
TESTS = ["wj_interaction", "anova_interaction", "kdT_vs_kdnoT", "kdT_vs_ctrlT",
         "simple_both_max", "tdp_alone", "kd_alone"]
# which cells each test depends on
DEPENDS = {"wj_interaction": CELLS, "anova_interaction": CELLS,
           "kdT_vs_kdnoT": ["kd_noT", "kd_T"], "kdT_vs_ctrlT": ["ctrl_T", "kd_T"],
           "simple_both_max": ["ctrl_T", "kd_noT", "kd_T"],
           "tdp_alone": ["ctrl_noT", "ctrl_T"], "kd_alone": ["ctrl_noT", "kd_noT"]}


# ----------------------------------------------------------------- data
def _tag(t): return t.split("}")[-1]
def _txt(e): return (e.text or "").strip()


def read_prism(path):
    root = ET.parse(path).getroot()
    cols = []
    for tbl in root:
        if _tag(tbl.tag) == "Table" and next((_txt(k) for k in tbl if _tag(k.tag) == "Title"), "") == "pigmantation":
            for k in tbl:
                if _tag(k.tag) == "YColumn":
                    for sub in k:
                        if _tag(sub.tag) == "Subcolumn":
                            cols.append(np.array([float(_txt(d)) for d in sub
                                                  if _tag(d.tag) == "d" and _txt(d) != ""]))
    return cols


def read_xlsx(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    out = {}
    for sheet in ["GFPR2", "GFPR3", "dDHX8R2", "dDHX8R3",
                  "TDP-43_GFPR2", "TDP-43_GFPR3", "TDP-43_dDHX8R2", "TDP-43_dDHX8R3"]:
        rows = [r for r in wb[sheet].iter_rows(min_row=2, max_row=17, values_only=True)]
        area = np.array([r[5] for r in rows], dtype=float)          # %Area
        out[sheet] = area[1::2] / area[0::2]                       # even row / odd row
    return out


pig = read_prism(PRISM)
pig_copy = read_prism(PRISM_COPY)
xl = read_xlsx(XLSX)
# Prism column order (verified against the xlsx below):
# 0 GFP#1-T 1 GFP#2-T 2 dDHX8#1-T 3 dDHX8#2-T 4 GFP#1+T 5 GFP#2+T 6 dDHX8#1+T 7 dDHX8#2+T
PRISM_TO_SHEET = ["GFPR2", "GFPR3", "dDHX8R2", "dDHX8R3",
                  "TDP-43_GFPR2", "TDP-43_GFPR3", "TDP-43_dDHX8R2", "TDP-43_dDHX8R3"]
data_check = {
    "prism_files_identical_values": bool(all(np.allclose(a, b, atol=0) for a, b in zip(pig, pig_copy))),
    "prism_vs_xlsx_max_abs_diff": float(max(np.max(np.abs(pig[i] - xl[s])) for i, s in enumerate(PRISM_TO_SHEET))),
    "prism_vs_xlsx_same_eye_order": bool(all(np.allclose(pig[i], xl[s], atol=1e-6) for i, s in enumerate(PRISM_TO_SHEET))),
}

LINES = {
    "line1_R2": dict(ctrl_noT=pig[0], kd_noT=pig[2], ctrl_T=pig[4], kd_T=pig[6]),
    "line2_R3": dict(ctrl_noT=pig[1], kd_noT=pig[3], ctrl_T=pig[5], kd_T=pig[7]),
}
SHEETS = {"line1_R2": dict(ctrl_noT="GFPR2", kd_noT="dDHX8R2", ctrl_T="TDP-43_GFPR2", kd_T="TDP-43_dDHX8R2"),
          "line2_R3": dict(ctrl_noT="GFPR3", kd_noT="dDHX8R3", ctrl_T="TDP-43_GFPR3", kd_T="TDP-43_dDHX8R3")}


# ----------------------------------------------------------------- vectorised tests
# Every test is a function of the cell means m[c], cell variances v[c] (ddof=1)
# and the per-cell n. Arrays broadcast, so the same code serves one data set or
# millions of pairings.
def welch_t_p(m1, v1, m2, v2, n):
    s1, s2 = v1 / n, v2 / n
    t = (m1 - m2) / np.sqrt(s1 + s2)
    df = (s1 + s2) ** 2 / (s1 ** 2 / (n - 1) + s2 ** 2 / (n - 1))
    return 2 * stats.t.sf(np.abs(t), df)


def wj_interaction_p(m, v, n):
    co = {"ctrl_noT": 1, "ctrl_T": -1, "kd_noT": -1, "kd_T": 1}
    L = sum(co[c] * m[c] for c in CELLS)
    vi = [co[c] ** 2 * v[c] / n for c in CELLS]
    V = sum(vi)
    df = V ** 2 / sum(x ** 2 / (n - 1) for x in vi)
    return 2 * stats.t.sf(np.abs(L / np.sqrt(V)), df)


def anova_interaction_p(m, v, n):
    L = m["ctrl_noT"] - m["ctrl_T"] - m["kd_noT"] + m["kd_T"]
    mse = sum(v[c] for c in CELLS) / 4.0            # balanced design
    F = L ** 2 / (mse * 4.0 / n)
    return stats.f.sf(F, 1, 4 * (n - 1))


def all_tests(m, v, n):
    out = {
        "wj_interaction": wj_interaction_p(m, v, n),
        "anova_interaction": anova_interaction_p(m, v, n),
        "kdT_vs_kdnoT": welch_t_p(m["kd_T"], v["kd_T"], m["kd_noT"], v["kd_noT"], n),
        "kdT_vs_ctrlT": welch_t_p(m["kd_T"], v["kd_T"], m["ctrl_T"], v["ctrl_T"], n),
        "tdp_alone": welch_t_p(m["ctrl_noT"], v["ctrl_noT"], m["ctrl_T"], v["ctrl_T"], n),
        "kd_alone": welch_t_p(m["ctrl_noT"], v["ctrl_noT"], m["kd_noT"], v["kd_noT"], n),
    }
    out["simple_both_max"] = np.maximum(out["kdT_vs_kdnoT"], out["kdT_vs_ctrlT"])
    return out


def tests_on_data(cells):
    n = len(cells["kd_T"])
    m = {c: cells[c].mean() for c in CELLS}
    v = {c: cells[c].var(ddof=1) for c in CELLS}
    return {k: float(x) for k, x in all_tests(m, v, n).items()}


# independent cross-check of the eye-level conventional ANOVA via OLS
def ols_interaction_p(cells):
    y, X = [], []
    for c in CELLS:
        T = 1.0 if c.endswith("_T") else 0.0
        K = 1.0 if c.startswith("kd") else 0.0
        for val in cells[c]:
            y.append(val); X.append([1.0, T, K, T * K])
    y, X = np.array(y), np.array(X)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    dfe = len(y) - X.shape[1]
    s2 = res @ res / dfe
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[3, 3])
    return float(2 * stats.t.sf(abs(beta[3] / se), dfe))


# ----------------------------------------------------------------- pairings
def perfect_matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        rest = items[1:i] + items[i + 1:]
        for m in perfect_matchings(rest):
            yield [(a, b)] + m


PAIRINGS = list(perfect_matchings(list(range(8))))
assert len(PAIRINGS) == 105
STORED = [(0, 1), (2, 3), (4, 5), (6, 7)]
STORED_IDX = PAIRINGS.index(STORED)


def fly_means(x, pairing):
    return np.array([(x[a] + x[b]) / 2.0 for a, b in pairing])


def pairing_label(p):
    return " ".join(f"{a + 1}+{b + 1}" for a, b in p)


def summarise(p):
    p = np.asarray(p, dtype=np.float64)
    q = np.quantile(p, [0.05, 0.5, 0.95])
    return {"min": float(p.min()), "p5": float(q[0]), "median": float(q[1]),
            "p95": float(q[2]), "max": float(p.max()),
            "frac_lt_0.05": float(np.mean(p < 0.05)), "n": int(p.size)}


def test_p(t, m, v, n):
    if t == "wj_interaction":
        return wj_interaction_p(m, v, n)
    if t == "anova_interaction":
        return anova_interaction_p(m, v, n)
    if t == "kdT_vs_kdnoT":
        return welch_t_p(m["kd_T"], v["kd_T"], m["kd_noT"], v["kd_noT"], n)
    if t == "kdT_vs_ctrlT":
        return welch_t_p(m["kd_T"], v["kd_T"], m["ctrl_T"], v["ctrl_T"], n)
    if t == "simple_both_max":
        return np.maximum(test_p("kdT_vs_kdnoT", m, v, n), test_p("kdT_vs_ctrlT", m, v, n))
    if t == "tdp_alone":
        return welch_t_p(m["ctrl_noT"], v["ctrl_noT"], m["ctrl_T"], v["ctrl_T"], n)
    if t == "kd_alone":
        return welch_t_p(m["ctrl_noT"], v["ctrl_noT"], m["kd_noT"], v["kd_noT"], n)
    raise KeyError(t)


def exhaustive_p(t, dep, VAR, MEAN):
    """P for every combination of the 105 pairings of the cells in `dep`
    (row-major over dep). Chunked over the first cell; stored as float64
    for up to 3 cells and float32 for 4 cells (121.6 M values)."""
    k = len(dep)
    block = 105 ** (k - 1)
    out = np.empty(105 ** k, dtype=np.float32 if k == 4 else np.float64)
    rest = dep[1:]
    grids = np.meshgrid(*[np.arange(105)] * len(rest), indexing="ij") if rest else []
    gi = {c: g.ravel() for c, g in zip(rest, grids)}
    for i0 in range(105):
        v = {c: (VAR[c][gi[c]] if c in gi else (VAR[c][i0] if c == dep[0] else VAR[c][0])) for c in CELLS}
        out[i0 * block:(i0 + 1) * block] = test_p(t, MEAN, v, 4)
    def decode(j):
        idx = np.unravel_index(j, (105,) * k)
        return {c: int(i) for c, i in zip(dep, idx)}
    return out, decode


def analyse_line(name, cells):
    res = {}
    n_eye = len(cells["kd_T"])
    # per-genotype fly means / variances for all 105 pairings (cell means are pairing-invariant)
    FM = {c: np.array([fly_means(cells[c], p) for p in PAIRINGS]) for c in CELLS}   # 105 x 4
    MEAN = {c: float(cells[c].mean()) for c in CELLS}
    VAR = {c: FM[c].var(axis=1, ddof=1) for c in CELLS}                           # 105
    assert all(np.allclose(FM[c].mean(axis=1), MEAN[c]) for c in CELLS)

    # eye level
    res["eye_level"] = tests_on_data(cells)
    res["eye_level"]["anova_interaction_ols_check"] = ols_interaction_p(cells)
    res["eye_level"]["n_per_genotype"] = n_eye
    res["cell_means"] = MEAN
    res["eye_sd"] = {c: float(cells[c].std(ddof=1)) for c in CELLS}

    # (a) stored order
    fm_stored = {c: FM[c][STORED_IDX] for c in CELLS}
    res["stored_order"] = tests_on_data(fm_stored)
    res["stored_order"]["fly_means"] = {c: [round(float(x), 6) for x in fm_stored[c]] for c in CELLS}
    res["stored_order"]["fly_sd"] = {c: float(fm_stored[c].std(ddof=1)) for c in CELLS}

    # (b1) random sample of joint pairings
    rng = np.random.default_rng(SEED)
    idx = {c: rng.integers(0, 105, N_RANDOM) for c in CELLS}
    m = {c: np.full(N_RANDOM, MEAN[c]) for c in CELLS}
    v = {c: VAR[c][idx[c]] for c in CELLS}
    pr = all_tests(m, v, 4)
    res["random_sample"] = {"seed": SEED, "n": N_RANDOM,
                            "tests": {k: summarise(pr[k]) for k in TESTS}}

    # (b2) exhaustive: each test enumerated over the cells it depends on
    #      (uniform over 105^4 joint pairings == uniform over 105^k for k cells)
    exh = {}
    for t in TESTS:
        dep = DEPENDS[t]
        p, decode = exhaustive_p(t, dep, VAR, MEAN)
        s = summarise(p)
        s["n_joint_pairings"] = 105 ** 4
        s["enumerated_over_cells"] = dep
        s["argmax_pairing"] = {c: pairing_label(PAIRINGS[k]) for c, k in decode(int(np.argmax(p))).items()}
        s["argmin_pairing"] = {c: pairing_label(PAIRINGS[k]) for c, k in decode(int(np.argmin(p))).items()}
        exh[t] = s
        del p
    res["exhaustive"] = exh

    # (c) coordinate-wise search for the worst case, several starts
    srch = {}
    rng2 = np.random.default_rng(SEED + 1)
    for t in ["wj_interaction", "anova_interaction", "simple_both_max", "kdT_vs_kdnoT"]:
        dep = DEPENDS[t]
        starts = [{c: STORED_IDX for c in CELLS}]
        starts += [{c: int(rng2.integers(0, 105)) for c in CELLS} for _ in range(N_STARTS - 1)]
        best_vals = []
        best = (-1.0, None)
        for st in starts:
            cur = dict(st)
            def pval(state):
                mm = {c: MEAN[c] for c in CELLS}
                vv = {c: VAR[c][state[c]] for c in CELLS}
                return float(test_p(t, mm, vv, 4))
            curp = pval(cur)
            changed, sweeps = True, 0
            while changed and sweeps < 50:
                changed, sweeps = False, sweeps + 1
                for c in dep:
                    vv = {cc: (VAR[cc] if cc == c else VAR[cc][cur[cc]]) for cc in CELLS}
                    ps = test_p(t, MEAN, vv, 4)
                    k = int(np.argmax(ps))
                    if ps[k] > curp * (1 + 1e-12):
                        cur[c], curp, changed = k, float(ps[k]), True
            best_vals.append(curp)
            if curp > best[0]:
                best = (curp, {c: pairing_label(PAIRINGS[cur[c]]) for c in dep})
        srch[t] = {"n_starts": len(starts), "max_found": best[0], "pairing": best[1],
                   "n_starts_reaching_exhaustive_max": int(np.sum(np.isclose(best_vals, exh[t]["max"], rtol=1e-9))),
                   "matches_exhaustive_max": bool(np.isclose(best[0], exh[t]["max"], rtol=1e-9))}
    res["coordinate_search"] = srch

    # fly-SD range per genotype across the 105 pairings
    res["fly_sd_range"] = {c: [float(np.sqrt(VAR[c].min())), float(np.sqrt(VAR[c].max()))] for c in CELLS}
    res["_VAR"] = VAR
    res["_MEAN"] = MEAN
    return res


# ----------------------------------------------------------------- run
REPORTED = {"line1_R2": {"wj_interaction": 7.3e-5, "anova_interaction": 4.2e-5,
                         "kdT_vs_kdnoT": 6.5e-6, "tdp_alone": 0.063, "kd_alone": 0.18},
            "line2_R3": {"wj_interaction": 1.3e-3}}

results = {"data": {"prism": PRISM, "xlsx_check": XLSX, **data_check,
                    "eye_order_note": "eye k = xlsx rows 2k-1 / 2k; Prism stores eyes in the same order",
                    "per_eye_values": {ln: {c: [float(x) for x in LINES[ln][c]] for c in CELLS} for ln in LINES},
                    "xlsx_sheet": SHEETS},
           "settings": {"seed": SEED, "n_random": N_RANDOM, "n_search_starts": N_STARTS,
                        "pairings_per_genotype": 105, "joint_pairings_per_line": 105 ** 4},
           "lines": {}}
internal = {}
for ln, cells in LINES.items():
    r = analyse_line(ln, cells)
    internal[ln] = {"VAR": r.pop("_VAR"), "MEAN": r.pop("_MEAN")}
    r["reproduced_vs_reported"] = {k: {"reported": v, "computed": r["eye_level"][k],
                                       "match_2sf": float(f"{r['eye_level'][k]:.2g}") == float(f"{v:.2g}")}
                                   for k, v in REPORTED[ln].items()}
    results["lines"][ln] = r

with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(results, fh, indent=1)

# ----------------------------------------------------------------- console summary
def f(x): return f"{x:.2g}"
for ln, r in results["lines"].items():
    print(f"\n=== {ln}")
    for k, d in r["reproduced_vs_reported"].items():
        print(f"  reproduce {k:18s} reported {d['reported']:<8g} computed {d['computed']:.3e} match={d['match_2sf']}")
    print(f"  eye ANOVA OLS check {r['eye_level']['anova_interaction_ols_check']:.3e}")
    print("  test               eye       stored    exh.min   exh.p5    exh.med   exh.p95   exh.max   frac<.05 | rand med  rand max  rand frac | search max")
    for t in TESTS:
        e, s, x, rs = r["eye_level"][t], r["stored_order"][t], r["exhaustive"][t], r["random_sample"]["tests"][t]
        cs = r["coordinate_search"].get(t, {}).get("max_found", float("nan"))
        print(f"  {t:18s} {f(e):9s} {f(s):9s} {f(x['min']):9s} {f(x['p5']):9s} {f(x['median']):9s} {f(x['p95']):9s} {f(x['max']):9s} {x['frac_lt_0.05']:.4f}   | {f(rs['median']):9s} {f(rs['max']):9s} {rs['frac_lt_0.05']:.4f}    | {f(cs)}")
    print("  worst-case pairing (WJ):", r["exhaustive"]["wj_interaction"]["argmax_pairing"])

# ----------------------------------------------------------------- figure (optional)
try:
    sys.path.insert(0, os.path.expanduser("~/bin"))
    import figprov
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figprov.apply_arial() if hasattr(figprov, "apply_arial") else None
    plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.4), constrained_layout=True)
    for ax, (ln, title) in zip(axes, [("line1_R2", "dDHX8-RNAi line #1"), ("line2_R3", "dDHX8-RNAi line #2")]):
        V, M = internal[ln]["VAR"], internal[ln]["MEAN"]
        pp, _ = exhaustive_p("wj_interaction", CELLS, V, M)
        lp = -np.log10(pp)
        cnt, edges = np.histogram(lp, bins=120)
        ax.stairs(cnt / lp.size, edges, fill=True, color="#4C72B0")
        del lp, pp
        r = results["lines"][ln]
        ax.axvline(-np.log10(0.05), color="black", lw=0.8, ls="--")
        ax.axvline(-np.log10(r["stored_order"]["wj_interaction"]), color="#DD8452", lw=1.0)
        ax.axvline(-np.log10(r["eye_level"]["wj_interaction"]), color="#55A868", lw=1.0)
        ax.set_xlabel("−log10 P, Welch–James interaction (fly means, n = 4)")
        ax.set_ylabel("Fraction of joint pairings")
        ax.set_title(title, fontsize=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    figprov.save_figure(fig, os.path.join(HERE, "figures", "fly_level_wj_interaction_distribution"),
                        inputs=[PRISM], script=__file__,
                        note="All 105^4 joint eye-to-fly pairings per line; dashed P=0.05, orange stored order, green eye level")
    plt.close(fig)
except Exception as exc:  # figure is optional
    print("figure skipped:", exc)
