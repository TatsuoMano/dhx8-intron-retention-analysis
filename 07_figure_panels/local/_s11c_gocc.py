#!/usr/bin/env python3
"""SuppS11 panel C: GO:CC enrichment of IR genes (g:Profiler) — faithful bespoke
re-creation from displayed values, authored 1:1 (~92x46mm), Arial, fonts 7-8pt."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.expanduser("~/bin"))
import figprov  # provenance-emitting save
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
# bottom -> top. v2 (siCont_2-corrected) g:Profiler g:SCS -log10 adjusted P.
# P-body dropped below significance in v2 (not returned by g:Profiler); other granule terms retained.
terms=["cytoplasmic stress granule","cytoplasmic ribonucleoprotein granule","ribonucleoprotein granule"]
vals=[2.74,3.87,10.36]
fig,ax=plt.subplots(figsize=(3.62,1.95))  # keep original ~2:1 aspect (Illustrator frame fits to height); labels fit via smaller font + left margin
ax.barh(range(len(terms)),vals,color="#d9685f",edgecolor="none",height=0.66)
for i,v in enumerate(vals):
    ax.text(v+0.12,i,f"{v:.1f}",va="center",ha="left",fontsize=7)
ax.axvline(-np.log10(0.05),ls="--",lw=0.7,color="0.5")  # FDR<0.05
ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms,fontsize=6.5)
ax.set_xlabel("-log10 (adjusted P)",fontsize=8); ax.tick_params(axis="x",labelsize=7)
ax.set_xlim(0,11.5); ax.set_title("GO:CC enrichment of IR genes",fontsize=8,fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
fig.subplots_adjust(left=0.52, right=0.96, top=0.82, bottom=0.24)  # reserve left space for long labels; no tight-crop so aspect stays fixed
figprov.save_figure(fig, "sfig_s11c_gocc.pdf", inputs=[], script=__file__,
    note="SuppFig S11C: GO:CC enrichment of IR genes (g:Profiler), v2 values; fixed 3.62x1.95 aspect + left margin so labels are not clipped by the Illustrator frame",
    savefig_kwargs={"bbox_inches": None})
print("WROTE sfig_s11c_gocc")
