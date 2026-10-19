#!/usr/bin/env python3
"""
Drosophila compound-eye pigmented-area-fraction reanalysis as a 2x2 factorial
design (TDP-43 [-/+] x dDHX8 RNAi [control/KD]), with two independent RNAi lines.

Source data : Fig8a.2 260514_analysis.pzfx (Fujikake), 'pigmantation' table, n=8/group.
Design       : main  = RNAi line #1 (R2) ; supp = RNAi line #2 (R3).
Stats        : two-way ANOVA via single-df Welch-James (unequal-variance) contrasts;
               simple effects by Welch t-tests.

Outputs (this dir / figures/):
  fig8a_eye_pigmented_area_R2_main.{pdf,png}   single panel, main figure (line #1)
  supp_eye_pigmented_area_R3.{pdf,png}         single panel, supplement (line #2)
  eye_twoway_review.{pdf,png}                  2x2 review grid (bars + interaction plots)
  eye_twoway_summary_stats.tsv                 descriptive + test table
"""
import os
import xml.etree.ElementTree as ET
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False

PRISM = ("/Users/tatsuomano/Dropbox/Research/DHX8-Sakai/Manuscript/Figures/"
         "figures_for_final/Fig08/Fig8a.2 260514_analysis.pzfx")
HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = "/Users/tatsuomano/Dropbox/Research/DHX8-Sakai/Analysis/2026-06-26_figure_fontfix"

# ----------------------------------------------------------------- parse Prism
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
# column order verified against 260508-2.xlsx sheet values:
# 0 GFP#1-T  1 GFP#2-T  2 dDHX8#1-T  3 dDHX8#2-T  4 GFP#1+T  5 GFP#2+T  6 dDHX8#1+T  7 dDHX8#2+T
SETS = {
    "main_R2": dict(title="RNAi line #1 (R2) — main",
                    ctrl_noT=pig[0], kd_noT=pig[2], ctrl_T=pig[4], kd_T=pig[6]),
    "supp_R3": dict(title="RNAi line #2 (R3) — supplement",
                    ctrl_noT=pig[1], kd_noT=pig[3], ctrl_T=pig[5], kd_T=pig[7]),
}

# ----------------------------------------------------------------- statistics
def welch_t(a, b):
    se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    t = (a.mean() - b.mean()) / se
    df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / (
        (a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1)
        + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
    return 2 * stats.t.sf(abs(t), df)

def welch_james(c):
    """unequal-variance two-way ANOVA effects (each 1 df) for a 2x2."""
    cells = [c["ctrl_noT"], c["ctrl_T"], c["kd_noT"], c["kd_T"]]
    contr = {"TDP main": [-.5, .5, -.5, .5],
             "DHX8 main": [-.5, -.5, .5, .5],
             "interaction": [1, -1, -1, 1]}
    out = {}
    for name, co in contr.items():
        L = sum(ci * x.mean() for ci, x in zip(co, cells))
        v = [(ci ** 2) * (x.var(ddof=1) / len(x)) for ci, x in zip(co, cells)]
        V = sum(v)
        df = V ** 2 / sum(vi ** 2 / (len(x) - 1) for vi, x in zip(v, cells))
        out[name] = (L ** 2 / V, 2 * stats.t.sf(abs(L / np.sqrt(V)), df), df)
    return out

def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

# ----------------------------------------------------------------- plotting
CTRL, KD = "#B0B0B0", "#D55E00"
RNG = np.random.RandomState(7)

def bar_panel(ax, c, show_pint=True):
    order = ["ctrl_noT", "kd_noT", "ctrl_T", "kd_T"]
    xs = [0, 0.85, 2.05, 2.90]
    cols = [CTRL, KD, CTRL, KD]
    for x, k, col in zip(xs, order, cols):
        d = c[k]
        ax.bar(x, d.mean(), width=0.7, color=col, edgecolor="black", lw=0.8, zorder=2)
        ax.errorbar(x, d.mean(), yerr=d.std(ddof=1), fmt="none", ecolor="black",
                    capsize=3, lw=1.0, zorder=3)
        ax.scatter(x + RNG.uniform(-0.15, 0.15, len(d)), d, s=13, color="black",
                   alpha=0.5, zorder=4, edgecolor="none")
    p_dhx8_noT = welch_t(c["ctrl_noT"], c["kd_noT"])
    p_dhx8_T = welch_t(c["ctrl_T"], c["kd_T"])
    p_tdp_kd = welch_t(c["kd_noT"], c["kd_T"])

    def br(x1, x2, y, p, h=0.016):
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.1, c="black")
        ax.text((x1 + x2) / 2, y + h + 0.002, stars(p), ha="center", va="bottom", fontsize=8)

    br(0, 0.85, 1.05, p_dhx8_noT)
    br(2.05, 2.90, 1.05, p_dhx8_T)
    br(0.85, 2.90, 1.18, p_tdp_kd)
    if show_pint:
        pint = welch_james(c)["interaction"][1]
        ax.text(1.45, 1.33, f"interaction  p = {pint:.1e}",
                ha="center", va="center", fontsize=7, style="italic")
    ax.set_xticks(xs)
    ax.set_xticklabels(["Ctrl\nRNAi", "dDHX8\nRNAi", "Ctrl\nRNAi", "dDHX8\nRNAi"], fontsize=7)
    for xc, lab in [(0.42, "− TDP-43"), (2.47, "+ TDP-43")]:
        ax.text(xc, -0.265, lab, ha="center", va="top",
                transform=ax.get_xaxis_transform(), fontsize=8, fontweight="bold")
    ax.plot([-0.2, 1.05], [-0.205, -0.205], color="black", lw=0.9,
            transform=ax.get_xaxis_transform(), clip_on=False)
    ax.plot([1.85, 3.10], [-0.205, -0.205], color="black", lw=0.9,
            transform=ax.get_xaxis_transform(), clip_on=False)
    ax.set_ylim(0, 1.42)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_ylabel("Pigmented area fraction", fontsize=8); ax.tick_params(labelsize=6.5)
    ax.set_xlim(-0.55, 3.45)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

def interaction_panel(ax, c):
    x = [0, 1]
    noT = [c["ctrl_noT"].mean(), c["kd_noT"].mean()]
    noT_e = [c["ctrl_noT"].std(ddof=1), c["kd_noT"].std(ddof=1)]
    T = [c["ctrl_T"].mean(), c["kd_T"].mean()]
    T_e = [c["ctrl_T"].std(ddof=1), c["kd_T"].std(ddof=1)]
    ax.errorbar(x, noT, yerr=noT_e, marker="o", ms=7, color="#7a7a7a", capsize=3,
                lw=2.0, label="− TDP-43")
    ax.errorbar(x, T, yerr=T_e, marker="s", ms=7, color=KD, capsize=3,
                lw=2.0, label="+ TDP-43")
    ax.set_xticks(x)
    ax.set_xticklabels(["Control\nRNAi", "dDHX8\nRNAi"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0.75, 1.02)
    ax.set_ylabel("Pigmented area fraction", fontsize=8); ax.tick_params(labelsize=6.5)
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    pint = welch_james(c)["interaction"][1]
    ax.text(0.98, 0.97, f"interaction\np = {pint:.1e}", transform=ax.transAxes,
            fontsize=7, style="italic", ha="right", va="top")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

# single panels (publication-ready): bar chart + standalone interaction plot
PANEL_LABEL = {"main_R2": "dDHX8-RNAi line #1", "supp_R3": "dDHX8-RNAi line #2"}
for key, fname in [("main_R2", "fig8a_eye_pigmented_area_R2_main"),
                   ("supp_R3", "supp_eye_pigmented_area_R3")]:
    fig, ax = plt.subplots(figsize=(2.4, 2.6))
    bar_panel(ax, SETS[key])
    ax.set_title(PANEL_LABEL[key], fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{fname}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(2.4, 2.4))
    interaction_panel(ax, SETS[key])
    ax.set_title(PANEL_LABEL[key], fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{fname}_interaction.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

# 2x2 review grid
fig, axes = plt.subplots(2, 2, figsize=(8.4, 8.2))
for j, key in enumerate(["main_R2", "supp_R3"]):
    bar_panel(axes[0, j], SETS[key])
    axes[0, j].set_title(SETS[key]["title"], fontsize=8)
    interaction_panel(axes[1, j], SETS[key])
    axes[1, j].set_title("interaction plot (mean ± SD)", fontsize=8)
fig.suptitle("Drosophila eye: 2×2 (TDP-43 × dDHX8 RNAi), n=8/group", fontsize=12, y=1.0)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIGDIR, f"eye_twoway_review.{ext}"), dpi=300, bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------- stats table
lines = ["set\tgroup\tmean\tSD\tn"]
for key, c in SETS.items():
    for k in ["ctrl_noT", "kd_noT", "ctrl_T", "kd_T"]:
        d = c[k]
        lines.append(f"{key}\t{k}\t{d.mean():.4f}\t{d.std(ddof=1):.4f}\t{len(d)}")
lines.append("")
lines.append("set\teffect/comparison\tstat\tp")
for key, c in SETS.items():
    wj = welch_james(c)
    for e in ["TDP main", "DHX8 main", "interaction"]:
        F, p, df = wj[e]
        lines.append(f"{key}\tANOVA {e} (Welch-James)\tF(1,{df:.1f})={F:.2f}\t{p:.3e}")
    lines.append(f"{key}\tsimple: TDP-43 alone (ctrl -T vs +T)\t-\t{welch_t(c['ctrl_noT'], c['ctrl_T']):.3e}")
    lines.append(f"{key}\tsimple: DHX8 alone (-TDP ctrl vs KD)\t-\t{welch_t(c['ctrl_noT'], c['kd_noT']):.3e}")
    lines.append(f"{key}\tsimple: KD+TDP vs TDP alone\t-\t{welch_t(c['ctrl_T'], c['kd_T']):.3e}")
    lines.append(f"{key}\tsimple: KD+TDP vs DHX8 alone\t-\t{welch_t(c['kd_noT'], c['kd_T']):.3e}")
with open(os.path.join(HERE, "eye_twoway_summary_stats.tsv"), "w") as fh:
    fh.write("\n".join(lines) + "\n")

print("\n".join(lines))
print("\nFigures written to", FIGDIR)
