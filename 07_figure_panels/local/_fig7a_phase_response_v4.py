#!/usr/bin/env python3
"""Fig7A DATA-EXACT (post-strandfix). NES from fgsea_spliceosome_phase_dhx8_s2,
per-phase mean log2FC from GMT x post-strandfix deseq2, key factors from deseq2."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.expanduser("~/bin"))
import figprov  # provenance-emitting save
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
RED,GREEN="#d9534f","#5cb85c"
fig,(axL,axR)=plt.subplots(1,2,figsize=(6.4,3.2),gridspec_kw={"width_ratios":[1.15,1]})

phases=["E/A\nComplex","B\nComplex","Bact/\nNTC","Regulatory","Post-\ncatalytic","Catalytic"]
nes=[2.494,2.254,2.072,2.224,1.937,1.554]
lfc=[0.408,0.369,0.342,0.260,-0.011,0.011]
cols=["#d9534f","#e8857a","#f0a35e","#f5c24f","#5b9bd5","#7fb77f"]
x=np.arange(len(phases))
axL.bar(x,nes,width=0.74,color=cols,edgecolor="none")
for xi,h,l in zip(x,nes,lfc):
    axL.text(xi,h+0.04,f"LFC\n{l:+.2f}",ha="center",va="bottom",fontsize=5.4,fontweight="bold",linespacing=0.9,color=("#222" if l>=0 else "#b03030"))
axL.axhline(1.30,ls="--",lw=0.7,color="0.55")
axL.text(5.45,1.33,"FDR<0.05",ha="right",va="bottom",fontsize=6,color="0.45")
axL.set_xticks(x); axL.set_xticklabels(phases,fontsize=6.5)
axL.set_ylabel("GSEA NES",fontsize=8); axL.set_ylim(0,2.75); axL.tick_params(axis="y",labelsize=7)
axL.set_title("Spliceosome phase-specific\ncompensatory response",fontsize=8,fontweight="bold")
axL.spines[["top","right"]].set_visible(False)

# Fig7C: ALL canonical spliceosomal ATP-dependent RNA helicases + Prpf3 factor,
# grouped & coloured by family. Bar direction (left/right of 0) = down/up.
import matplotlib.patches as mpatches
C_DEAH="#2b7bba"  # DEAH-box (DHX / AQR) = DHX8's own family
C_DEAD="#e8873b"  # DEAD-box (DDX)
C_SKI ="#5aa469"  # Ski2-like (SNRNP200/Brr2)
C_FAC ="#9a9a9a"  # non-helicase tri-snRNP factor
# label, log2FC, colour  (bottom -> top, grouped by family with gaps)
fac =["Ddx46\n(Prp5)","Ddx23\n(Prp28)", "Snrnp200\n(Brr2)", "Dhx16\n(Prp2)","Dhx8\n(Prp22)","Aqr\n(IBP160)","Dhx15\n(Prp43)","Dhx38\n(Prp16)", "Prpf3"]
val =[0.447,0.501,                       0.395,             -0.641,-0.553,0.389,0.531,0.663,                             1.646]
bcol=[C_DEAD,C_DEAD,                     C_SKI,             C_DEAH,C_DEAH,C_DEAH,C_DEAH,C_DEAH,                          C_FAC]
y   =np.array([0.0,1.0,                  2.4,               3.8,4.8,5.8,6.8,7.8,                                        9.2])
axR.barh(y,val,height=0.7,color=bcol,edgecolor="none")
for yi,v in zip(y,val):
    axR.text(v+(0.04 if v>0 else -0.04),yi,f"{v:+.2f}",va="center",ha="left" if v>0 else "right",fontsize=5.6,fontweight="bold")
axR.axvline(0,color="0.3",lw=0.7)
axR.set_yticks(y); axR.set_yticklabels(fac,fontsize=5.2)
axR.set_xlabel("log2 Fold Change",fontsize=7.5); axR.tick_params(axis="x",labelsize=6.5)
axR.set_xlim(-1.0,2.0); axR.set_ylim(-0.7,9.9)
axR.set_title("Spliceosomal RNA helicases\n(and Prpf3) in DHX8 KD",fontsize=8,fontweight="bold")
axR.spines[["top","right"]].set_visible(False)
_leg=[mpatches.Patch(color=C_DEAH,label="DEAH-box (DHX/AQR)"),
      mpatches.Patch(color=C_DEAD,label="DEAD-box (DDX)"),
      mpatches.Patch(color=C_SKI,label="Ski2-like (Brr2)"),
      mpatches.Patch(color=C_FAC,label="tri-snRNP factor")]
axR.legend(handles=_leg,fontsize=4.6,loc="lower right",frameon=False,
           handlelength=0.9,handletextpad=0.35,labelspacing=0.2,borderaxespad=0.2)

fig.tight_layout(w_pad=1.5)
figprov.save_figure(fig, "fig7a_phase_response_v4.pdf", inputs=[], script=__file__,
    note="Fig7A/B/C: per-phase mean log2FC (curated) + helicase family LFC; v2 post-strandfix deseq2 values")
print("WROTE fig7a_phase_response_v4 (data-exact: NES=GMT fgsea, LFC=curated mean_log2FC)")
