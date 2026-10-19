#!/usr/bin/env python3
"""Fig6E Biancon as TWO separate single venns, geometry EXACTLY matched to the placed
v2 Khong panel (fig6a_venn_khong_v2.pdf). The v2 Khong script is server-side/unavailable,
so its output geometry was reverse-engineered from the PDF:
  page 136.8 x 115.2 pt (1.90 x 1.60 in); Tf font sizes 8.0 (subset) / 6.5 (set) / 8.5 (title).
Reproduced with figsize (1.90,1.60) + those fonts + PLAIN savefig (no bbox='tight'), so the
page size and font sizes match the v2 Khong 1:1 -> identical text size when placed. v2 values."""
import os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
plt.rcParams["font.family"] = "Arial"; plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
OUTD = os.path.dirname(os.path.abspath(__file__))
IRC = "#ec9a87"; BLUE = "#7bb3d9"
FIGSIZE = (1.90, 1.60)            # = 136.8 x 115.2 pt, matches v2 Khong page
FS_SUBSET, FS_SET, FS_TITLE = 8.0, 6.5, 8.5   # matches v2 Khong Tf sizes

def single(name, ir_only, ov, sg_only, sg_label, or_text):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    v = venn2(subsets=(ir_only, sg_only, ov), set_labels=("IR genes\n(1,149)", sg_label),
              set_colors=(IRC, BLUE), alpha=0.65, ax=ax)
    for tid in ("10", "01", "11"):
        t = v.get_label_by_id(tid)
        if t: t.set_fontsize(FS_SUBSET)
    for sl in v.set_labels:
        if sl: sl.set_fontsize(FS_SET)
    ax.set_title(or_text, fontsize=FS_TITLE, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTD, name + ".pdf"))
    fig.savefig(os.path.join(OUTD, name + ".png"), dpi=300)
    plt.close(fig)

single("fig6e_biancon_high_khongmatch",    999, 150, 765,
       "SG RNA meta\n(Biancon 2025, high)",     "OR = 2.14  ***")
single("fig6e_biancon_medhigh_khongmatch",  669, 480, 3247,
       "SG RNA meta\n(Biancon 2025, med-high)", "OR = 1.86  ***")
print("done")
