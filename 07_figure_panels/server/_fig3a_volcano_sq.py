#!/usr/bin/env python3
"""Fig3A volcano, SQUARE (~84x84mm), authored 1:1, fonts >=7pt, dev-title removed
(counts kept in legend). red UP / blue DOWN / gray NS; dashed thresholds."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
TSV="/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv"
OUT="/mnt/HDD2/unified_results/postfix_figures_pdf/fig3a_volcano_sq"
LFC_TH,FDR_TH,YCAP=1.0,0.05,300.0
C_UP,C_DN,C_NS="#d9534f","#4f81bd","#bdbdbd"
df=pd.read_csv(TSV,sep="\t").dropna(subset=["padj","log2FoldChange"])
y=(-np.log10(df["padj"].clip(lower=10**(-YCAP)))).clip(upper=YCAP)
x=df["log2FoldChange"].values
up=(df["padj"]<FDR_TH)&(df["log2FoldChange"]>LFC_TH)
dn=(df["padj"]<FDR_TH)&(df["log2FoldChange"]<-LFC_TH)
ns=~(up|dn); n_up,n_dn=int(up.sum()),int(dn.sum())
fig,ax=plt.subplots(figsize=(3.35,3.35))   # ~85x85mm square, 1:1
ax.scatter(x[ns],y[ns],s=5,c=C_NS,alpha=0.5,linewidths=0,rasterized=True)
ax.scatter(x[dn],y[dn],s=6,c=C_DN,alpha=0.7,linewidths=0,rasterized=True)
ax.scatter(x[up],y[up],s=6,c=C_UP,alpha=0.7,linewidths=0,rasterized=True)
ax.axvline(LFC_TH,ls="--",lw=0.7,color="0.4"); ax.axvline(-LFC_TH,ls="--",lw=0.7,color="0.4")
ax.axhline(-np.log10(FDR_TH),ls="--",lw=0.7,color="0.4")
ax.set_xlabel("log2 fold change",fontsize=8)
ax.set_ylabel("-log10 adjusted P",fontsize=8)
ax.tick_params(labelsize=7)
h=[plt.Line2D([],[],marker="o",ls="",color=C_UP,label=f"UP: {n_up}"),
   plt.Line2D([],[],marker="o",ls="",color=C_DN,label=f"DOWN: {n_dn}")]
ax.legend(handles=h,fontsize=7.5,loc="upper right",framealpha=0.9,handletextpad=0.3,borderpad=0.3)
ax.grid(True,lw=0.4,color="0.92"); ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(OUT+".pdf",dpi=600); fig.savefig(OUT+".png",dpi=600)
print("WROTE",OUT,"| UP",n_up,"DOWN",n_dn)
