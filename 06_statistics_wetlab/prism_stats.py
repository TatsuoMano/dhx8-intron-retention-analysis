#!/usr/bin/env python3
"""Statistics for the Prism-derived graphs, recomputed in Python with the tests named in the manuscript legends,
and compared with the P-values printed on the current Illustrator figures.

  Fig 1C  two-way ANOVA (siRNA x genotype) + Šídák multiple comparisons (4 simple effects)   [legend]
  Fig 1E  "pairwise P-values" (test not named): one-way ANOVA + Šídák on the pairs bracketed in the figure,
          and Tukey HSD, both reported so the test used by Prism can be identified
  Fig 8B  one-way ANOVA + Dunnett vs Cont.                                                  [legend]
  S2B     test not named in the legend: paired and Welch t-tests T=1 vs T=20 min per group are reported
  S7A/B   unpaired t-test                                                                    [legend]
Writes panels/prism_stats.tsv and returns a dict used by prism_panels.py for the annotations.
"""
import os, sys, itertools, math
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prism_panels import read_pzfx, oneway, FFF  # noqa: E402

# P-values printed on the current figures (pdftotext of MainFigures_TM2.pdf / SuppFigures_TM2.pdf, 2026-09-16)
FIG = {"1C": {"WT: siCont vs siDhx8": "<0.0001", "KO: siCont vs siDhx8": "<0.0001", "siCont: WT vs KO": "0.0017", "siDhx8: WT vs KO": "<0.0005 (printed)"},
       "8B_TDP-43": {"dDHX8 #1": "0.8540", "dDHX8 #2": "0.3097"}, "8B_TDP-35": {"dDHX8 #1": "0.0478", "dDHX8 #2": "0.0327"},
       "S2B": {"siCont": "0.0003", "siDhx8": "0.9154"}, "S7A": {"siCont vs siDHX8": "<0.0001"}, "S7B": {"siCont vs siDhx8": "<0.0001"},
       "1E": "0.8050 0.9846 0.5973 <0.0001 0.0012 >0.9999 0.0631 0.0158 0.3232 0.9991 0.0009 0.9789 0.9936 >0.9999 (14 brackets)"}


def sidak(p, k):
    return 1 - (1 - p) ** k


def fmt(p):
    return "<0.0001" if p < 0.00005 else (">0.9999" if p >= 0.99995 else "%.4f" % p)


def fig1c():
    t = read_pzfx(os.path.join(FFF, "Fig01", "Fig1c.pzfx"))["Data 3"]
    sg = np.array([[s[r] for s in t["cols"][0][1]] for r in range(4)])
    rows = []
    for r, (si, ge) in enumerate([("siCont", "WT"), ("siDhx8", "WT"), ("siCont", "KO"), ("siDhx8", "KO")]):
        for v in sg[r]:
            if not np.isnan(v): rows.append((si, ge, v))
    df = pd.DataFrame(rows, columns=["si", "geno", "y"])
    m = smf.ols("y ~ C(si) * C(geno)", df).fit(); an = sm.stats.anova_lm(m, typ=2)
    mse, dfr = m.mse_resid, m.df_resid
    g = df.groupby(["si", "geno"])["y"].agg(["mean", "count"])
    def cmp(a, b):
        d = g.loc[a, "mean"] - g.loc[b, "mean"]; se = math.sqrt(mse * (1 / g.loc[a, "count"] + 1 / g.loc[b, "count"]))
        return 2 * stats.t.sf(abs(d / se), dfr)
    pairs = {"WT: siCont vs siDhx8": (("siCont", "WT"), ("siDhx8", "WT")), "KO: siCont vs siDhx8": (("siCont", "KO"), ("siDhx8", "KO")),
             "siCont: WT vs KO": (("siCont", "WT"), ("siCont", "KO")), "siDhx8: WT vs KO": (("siDhx8", "WT"), ("siDhx8", "KO"))}
    out = {k: sidak(cmp(*v), 4) for k, v in pairs.items()}
    return out, {"anova": an.to_dict()}


def fig1e():
    g = oneway(read_pzfx(os.path.join(FFF, "Fig01", "Fig1f. P-EIF2a定量_240112 .pzfx"))["P-EIF2α"])
    names = ["AS-/WT/-", "AS-/WT/+", "AS-/KO/-", "AS-/KO/+", "AS+/WT/-", "AS+/WT/+", "AS+/KO/-", "AS+/KO/+"]
    df = pd.DataFrame([(n, v) for (_, vals), n in zip(g, names) for v in vals], columns=["grp", "y"])
    m = smf.ols("y ~ C(grp)", df).fit(); mse, dfr = m.mse_resid, m.df_resid
    means = df.groupby("grp")["y"].agg(["mean", "count"])
    def raw(a, b):
        d = means.loc[a, "mean"] - means.loc[b, "mean"]; se = math.sqrt(mse * (1 / means.loc[a, "count"] + 1 / means.loc[b, "count"]))
        return 2 * stats.t.sf(abs(d / se), dfr)
    pairs = list(itertools.combinations(names, 2))
    tuk = pairwise_tukeyhsd(df["y"], df["grp"])
    tk = {(row[0], row[1]): row[3] for row in tuk._results_table.data[1:]}
    res = {}
    for a, b in pairs:
        res["%s vs %s" % (a, b)] = {"sidak28": sidak(raw(a, b), 28), "sidak14": sidak(raw(a, b), 14), "tukey": tk.get((a, b), tk.get((b, a)))}
    return res


def fig8b():
    T = read_pzfx(os.path.join(FFF, "Fig08", "Fig8b4. Fujikake-_定量データ_250716.pzfx")); out = {}
    for key in ("TDP-43", "TDP-35"):
        g = oneway(T[key]); ctrl = g[0][1]; others = [v for _, v in g[1:]]
        d = stats.dunnett(*others, control=ctrl, random_state=20260917)   # fixed seed: the Dunnett p is Monte-Carlo integrated
        out[key] = {g[i + 1][0]: float(p) for i, p in enumerate(d.pvalue)}
    return out


def s2b():
    t = read_pzfx(os.path.join(FFF, "SupFig02", "SupFig02a. 2SGs定量_230609IC.pzfx"))["SGs number"]; out = {}
    for name, subs in t["cols"]:
        a = np.array(subs); first, last = a[:, 0], a[:, -1]; ok = ~np.isnan(first) & ~np.isnan(last)
        out[name] = {"paired_t_T1_vs_T20": stats.ttest_rel(first[ok], last[ok]).pvalue, "welch_t_T1_vs_T20": stats.ttest_ind(first[ok], last[ok], equal_var=False).pvalue,
                     "student_t_T1_vs_T20": stats.ttest_ind(first[ok], last[ok]).pvalue}
    return out


def s7():
    g = oneway(read_pzfx(os.path.join(FFF, "SupFig12", "SupFig12a._DHX8KD_tdp43_actb_gapdh まとめ_250306.pzfx"))["TDP-43/ActB qPCR (DHX8 #1のみ)"])
    a = stats.ttest_ind(g[0][1], g[1][1]).pvalue
    g2 = oneway(read_pzfx(os.path.join(FFF, "SupFig12", "SupFig12b&c_231222 定量.pzfx"))["T=0 TDP-43"])
    b = stats.ttest_ind(g2[0][1], g2[1][1]).pvalue
    return {"S7A": a, "S7B": b}


def compute():
    R = {}
    R["1C"], R["1C_anova"] = fig1c(); R["1E"] = fig1e(); R["8B"] = fig8b(); R["S2B"] = s2b(); R.update(s7())
    return R


if __name__ == "__main__":
    R = compute(); lines = ["panel\tcomparison\ttest\tP (python)\tP printed on current figure"]
    for k, p in R["1C"].items(): lines.append("Fig 1C\t%s\ttwo-way ANOVA + Šídák (4)\t%s\t%s" % (k, fmt(p), FIG["1C"][k]))
    for k, d in R["1E"].items(): lines.append("Fig 1E\t%s\tŠídák(14) / Šídák(28) / Tukey\t%s / %s / %s\t(figure: %s)" % (k, fmt(d["sidak14"]), fmt(d["sidak28"]), fmt(d["tukey"]), FIG["1E"]))
    for key, d in R["8B"].items():
        for k, p in d.items(): lines.append("Fig 8B %s\t%s vs Cont.\tone-way ANOVA + Dunnett\t%s\t%s" % (key, k, fmt(p), FIG["8B_" + key][k]))
    for k, d in R["S2B"].items():
        for tname, p in d.items(): lines.append("Fig S2B\t%s\t%s\t%s\t%s" % (k, tname, fmt(p), FIG["S2B"]["siCont" if "Cont" in k else "siDhx8"]))
    lines.append("Fig S7A\tsiCont vs siDHX8\tunpaired t\t%s\t%s" % (fmt(R["S7A"]), FIG["S7A"]["siCont vs siDHX8"]))
    lines.append("Fig S7B\tsiCont vs siDhx8\tunpaired t\t%s\t%s" % (fmt(R["S7B"]), FIG["S7B"]["siCont vs siDhx8"]))
    os.makedirs("panels", exist_ok=True); open("panels/prism_stats.tsv", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
