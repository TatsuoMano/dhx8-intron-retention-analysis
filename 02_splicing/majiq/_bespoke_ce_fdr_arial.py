#!/usr/bin/env python3
"""Faithful Arial vector-PDF re-creation of fig4d_ce_combined_fdr.
Panel A: cryptic-exon candidate scatter (is_denovo events; ctrl PSI<0.05 & disease PSI>0.2
= CE; PSI>0.5 = strong), 2 rows (DHX8 KD, GSE126543). Panel B: mirror-method empirical FDR
eFDR(t)=#{dPSI<-t}/#{dPSI>t} over all de-novo events. Uses is_denovo -> counts match the
manuscript exactly (DHX8 22/12/1, GSE126543 160/24/9)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
plt.rcParams["font.family"] = "Arial"; plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = False

TSV = "/mnt/HDD2/unified_results/majiq/all_datasets_significant_events.tsv"
OUT = "/mnt/HDD2/unified_results/postfix_figures_pdf/fig4d_ce_combined_fdr"
df = pd.read_csv(TSV, sep="\t")

DS = [("dhx8", "DHX8 KD", "Control_raw_psi_mean", "DHX8_KD_raw_psi_mean", "#e07b54"),
      ("gse126543", "GSE126543\n(TDP-43+/-)", "TDP43_positive_raw_psi_mean",
       "TDP43_negative_raw_psi_mean", "#8e6fb0")]
C_STRONG, C_WEAK, C_NS = "#c0392b", "#e8a598", "#c8c8c8"

fig = plt.figure(figsize=(13, 6.8))
gs = GridSpec(1, 3, width_ratios=[2.05, 0.06, 1.0], wspace=0.25)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 2])

XLO, XHI = -0.8, 1.05
rng = np.random.default_rng(0)
yrows = {"dhx8": 1.0, "gse126543": 0.0}
# label ladders: top row labels go ABOVE its band, bottom row labels go BELOW
LADDER = {"dhx8": dict(y=1.70, va="bottom", anchor=+0.34),
          "gse126543": dict(y=-0.70, va="top", anchor=-0.34)}
axA.axvspan(0.2, XHI, color="#f6e4e4", zorder=0)
axA.axvspan(XLO, 0.2, color="#e9edf3", zorder=0)
axA.text(0.92, 1.40, "CE direction →", color="#c0392b", fontsize=10, ha="right", fontweight="bold")

def place_labels(ax, rows, yc, lad):
    """Spread strong-CE gene labels on an evenly-spaced ladder with thin leader
    lines so they neither overlap each other nor the dPSI points."""
    rows = sorted(rows, key=lambda t: t[0])          # by dPSI
    n = len(rows)
    if n == 0:
        return
    xs = [t[0] for t in rows]
    lo = max(XLO + 0.04, min(xs) - 0.10)
    hi = min(XHI - 0.04, max(xs) + 0.10)
    need = 0.085 * (n - 1)                            # min spacing for vertical text
    if hi - lo < need:                               # widen if the cluster is too tight
        mid = (lo + hi) / 2
        lo, hi = max(XLO + 0.04, mid - need / 2), min(XHI - 0.04, mid + need / 2)
    lab_x = np.linspace(lo, hi, n) if n > 1 else [np.mean(xs)]
    for (xp, gene, psi), lx in zip(rows, lab_x):
        ax.annotate(f"{gene} ({psi:.2f})", xy=(xp, yc + lad["anchor"]),
                    xytext=(lx, lad["y"]), fontsize=6.5, color="#7a1f12",
                    rotation=90, ha="center", va=lad["va"], zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.55",
                                    shrinkA=1, shrinkB=2))

for ds, label, cc, dc, _ in DS:
    s = df[(df.dataset == ds) & (df.is_denovo == True)].copy()
    s["dpsi"] = s[dc] - s[cc]
    yc = yrows[ds]
    y = yc + rng.uniform(-0.32, 0.32, len(s))
    is_ce = (s[cc] < 0.05) & (s[dc] > 0.2)
    is_strong = is_ce & (s[dc] > 0.5)
    is_weak = is_ce & ~is_strong
    ns = ~is_ce
    axA.scatter(s["dpsi"][ns], y[ns.values], s=12, c=C_NS, alpha=0.5, linewidths=0, zorder=2)
    axA.scatter(s["dpsi"][is_weak], y[is_weak.values], s=42, c=C_WEAK, marker="D",
                edgecolors="#a04030", linewidths=0.4, zorder=3)
    axA.scatter(s["dpsi"][is_strong], y[is_strong.values], s=70, c=C_STRONG, marker="D",
                edgecolors="k", linewidths=0.6, zorder=4)
    axA.text(-0.75, yc + 0.44, f"n={len(s)}  CE={int(is_ce.sum())}  (strong={int(is_strong.sum())})",
             fontsize=8, color="0.3")
    # strong-CE gene labels on a leader-line ladder (no overlap with points)
    place_labels(axA, [(r["dpsi"], r["gene_name"], r[dc]) for _, r in s[is_strong].iterrows()],
                 yc, LADDER[ds])
axA.axvline(0, color="0.4", lw=0.8, ls="--", zorder=1)
axA.set_yticks([1.0, 0.0]); axA.set_yticklabels(["DHX8 KD", "GSE126543\n(TDP-43+/-)"], fontsize=9)
axA.set_ylim(-1.45, 2.15); axA.set_xlim(XLO, XHI)
axA.set_xlabel("dPSI (disease/KD − control/normal)", fontsize=10)
axA.set_title("A. Cryptic Exon Candidates (ctrl PSI<0.05, disease PSI>0.2)", fontsize=10, fontweight="bold")
from matplotlib.lines import Line2D
axA.legend(handles=[Line2D([], [], marker="D", ls="", mfc=C_STRONG, mec="k", label="Strong CE (PSI>0.5)"),
                    Line2D([], [], marker="D", ls="", mfc=C_WEAK, mec="#a04030", label="Weak CE (PSI 0.2-0.5)"),
                    Line2D([], [], marker="o", ls="", mfc=C_NS, mec="none", label="Not CE candidate")],
           fontsize=7, loc="lower left", framealpha=0.9, ncol=1)

# Panel B: mirror-method eFDR
thr = np.linspace(0.2, 0.85, 80)
for ds, label, cc, dc, col in DS:
    s = df[(df.dataset == ds) & (df.is_denovo == True)].copy()
    d = (s[dc] - s[cc]).values
    efdr = [ (np.sum(d < -t) / np.sum(d > t)) if np.sum(d > t) > 0 else np.nan for t in thr ]
    axB.plot(thr, np.clip(efdr, 0, 1.05), color=col, lw=1.6, label=label.replace("\n", " "))
axB.axhline(0.1, color="0.6", ls="--", lw=0.7); axB.text(0.83, 0.11, "FDR=0.1", fontsize=6.5, color="0.5", ha="right")
axB.axhline(0.05, color="0.7", ls="--", lw=0.7); axB.text(0.83, 0.058, "FDR=0.05", fontsize=6.5, color="0.5", ha="right")
axB.set_xlabel("dPSI threshold", fontsize=10); axB.set_ylabel("Empirical FDR (mirror method)", fontsize=9)
axB.set_title("B. Empirical FDR (all denovo events)\n#{dPSI < −t} / #{dPSI > t}", fontsize=10, fontweight="bold")
axB.set_ylim(-0.03, 1.05); axB.legend(fontsize=8, loc="upper left")
fig.tight_layout()
fig.savefig(OUT + ".pdf"); fig.savefig(OUT + ".png", dpi=200)
print("WROTE", OUT + ".pdf")
