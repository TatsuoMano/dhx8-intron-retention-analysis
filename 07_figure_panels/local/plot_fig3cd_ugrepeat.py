#!/usr/bin/env python3
"""Fig3C (k-mer volcano) + Fig3D (top-hexamer table) rebuilt for the corrected
UG-repeat / TDP-43 result. Foreground = strict DEG-UP (log2FC>1 & padj<0.05, =
manuscript Fig3A definition), background = MEF-expressed 3'UTRs (Transite k-mer
TSMA). poly-U (TIA1) is NOT enriched and is marked as such. Arial, pdf.fonttype=42."""
import os, sys, csv, math
sys.path.insert(0, os.path.expanduser("~/bin"))
import figprov
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"; plt.rcParams["font.sans-serif"] = ["Arial","Helvetica","DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = False

WD = os.path.dirname(os.path.abspath(__file__))
KMER_TSV = os.path.join(WD, "fig3c_kmer_strict_mef.tsv")

rows = list(csv.DictReader(open(KMER_TSV), delimiter="\t"))
kmer = np.array([r["kmer"] for r in rows])
enr  = np.array([float(r["enrichment"]) for r in rows])
adjp = np.array([float(r["adj_p_value"]) for r in rows])
log2e = np.log2(np.clip(enr, 1e-6, None))
neglp = -np.log10(np.clip(adjp, 1e-300, None))

def is_rep(k, a, b):  # (ab)3 in either phase
    return k in (a + b + a + b + a + b, b + a + b + a + b + a)
cls = {}
for k in kmer:
    if is_rep(k, "U", "G"): cls[k] = "UG"
    elif is_rep(k, "U", "C"): cls[k] = "UC"
    elif is_rep(k, "G", "C"): cls[k] = "GC"
    else: cls[k] = "other"
cl = np.array([cls[k] for k in kmer])

# ---------- Fig3C: volcano ----------
fig, ax = plt.subplots(figsize=(2.55, 2.35))  # ~65x60 mm placed size
sig = adjp < 0.05
ax.scatter(log2e[cl == "other"], neglp[cl == "other"], s=5, c="0.78", edgecolors="none", rasterized=True, zorder=1)
palette = {"UG": "#c0392b", "UC": "#e08214", "GC": "#7f7f7f"}
for c, lab in [("GC", "GC-rich (low-complexity)"), ("UC", "UC-repeat"), ("UG", "UG-repeat")]:
    m = cl == c
    ax.scatter(log2e[m], neglp[m], s=16, c=palette[c], edgecolors="none", zorder=3, label=lab)
# label top UG-repeats
for k in ("GUGUGU", "UGUGUG"):
    i = np.where(kmer == k)[0][0]
    ax.annotate(k, (log2e[i], neglp[i]), textcoords="offset points", xytext=(4, -1),
                fontsize=5.6, color="#c0392b", ha="left", va="top", fontweight="bold")
ax.axhline(-math.log10(0.05), ls="--", lw=0.6, color="0.5")
ax.text(ax.get_xlim()[0] if False else -1.02, -math.log10(0.05) + 0.4, "adj. P = 0.05", fontsize=5.0, color="0.5")
ax.axvline(0, ls="-", lw=0.5, color="0.8", zorder=0)
ax.set_xlabel("log2 enrichment (up-regulated 3'UTR / expressed)", fontsize=6.5)
ax.set_ylabel("−log10 adjusted P", fontsize=6.5)
ax.tick_params(labelsize=6)
ax.set_title("Hexamer motif enrichment", fontsize=7, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
leg = ax.legend(fontsize=5.4, frameon=False, loc="upper left", handletextpad=0.2, borderpad=0.1, labelspacing=0.2)
fig.subplots_adjust(left=0.16, right=0.97, top=0.90, bottom=0.17)
figprov.save_figure(fig, "fig3c_kmer_volcano_ugrepeat.pdf", inputs=[KMER_TSV], script=__file__,
    note="Fig3C k-mer TSMA volcano, strict DEG-UP (log2FC>1) fg vs MEF-expressed 3'UTR bg. "
         "UG-repeat (UGUGUG/GUGUGU) top-significant; poly-U(TIA1) not enriched (enr 0.91, ns). "
         "Robust to background (all-3'UTR) and to 3'UTR length-matching (KS p=1, UG-repeat adj_p<1e-16).")
plt.close(fig)
print("WROTE fig3c_kmer_volcano_ugrepeat")

# ---------- Fig3D: top-hexamer table ----------
kd = {r["kmer"]: r for r in rows}
def row(k):
    r = kd[k]; return (k, float(r["enrichment"]), float(r["adj_p_value"]))
table = [
    # hexamer, class, associated RBP, enr, adjP
    ("GUGUGU", "UG-repeat", "TDP-43 (UG-repeat); CELF"),
    ("UGUGUG", "UG-repeat", "TDP-43 (UG-repeat); CELF"),
    ("UCUCUC", "UC-repeat", "PTBP1 (pyrimidine)"),
    ("CUCUCU", "UC-repeat", "PTBP1 (pyrimidine)"),
    ("GCGCGC", "GC-rich",   "GC-rich, low-complexity"),
]
def fmtp(p):
    if p >= 0.05: return f"{p:.2f}"
    e = int(math.floor(math.log10(p))); m = p / 10**e
    return f"{m:.1f}e{e:d}"  # e.g. 3.3e-15 (Arial lacks Unicode super/subscripts)

fig, ax = plt.subplots(figsize=(3.35, 1.75))
ax.axis("off")
cols = ["Hexamer", "Class", "Associated RBP", "Enrich.", "Adj. P"]
cw = [0.14, 0.16, 0.42, 0.11, 0.17]
x0 = [sum(cw[:i]) for i in range(len(cw))]
yh = 0.135; ytop = 0.90
# header
for j, c in enumerate(cols):
    ax.text(x0[j] + 0.004, ytop, c, fontsize=6.4, fontweight="bold", va="center", ha="left")
ax.plot([0, 1], [ytop - yh * 0.55] * 2, lw=0.8, color="0.3")
rowcolors = {"UG-repeat": "#c0392b", "UC-repeat": "#e08214", "GC-rich": "#7f7f7f", "poly-U": "#2166ac"}
for i, (k, klass, rbp) in enumerate(table):
    y = ytop - (i + 1) * yh
    _, e, p = row(k)
    col = rowcolors[klass]
    cells = [k, klass, rbp, f"{e:.2f}", fmtp(p)]
    for j, txt in enumerate(cells):
        w = "bold" if (j == 0 and klass in ("UG-repeat",)) else "normal"
        ax.text(x0[j] + 0.004, y, txt, fontsize=6.0, va="center", ha="left",
                color=col if j <= 1 else "0.1", fontweight=w)
ax.plot([0, 1], [ytop - (len(table) + 0.4) * yh] * 2, lw=0.6, color="0.6")
ax.set_title("Top enriched hexamer motifs (Transite TSMA)", fontsize=7, fontweight="bold", pad=2)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.04)
figprov.save_figure(fig, "fig3d_motif_table_ugrepeat.pdf", inputs=[KMER_TSV], script=__file__,
    note="Fig3D top-hexamer table (matches legend 'hexamer motifs with associated RBPs'). "
         "UG-repeat = TDP-43 canonical UG-repeat site (in vivo CLIP; CELF in Transite in-vitro DB). "
         "poly-U/TIA1 row shown as not-enriched. Enrichment/adj_p from strict-fg MEF-bg k-mer TSMA.")
plt.close(fig)
print("WROTE fig3d_motif_table_ugrepeat")
