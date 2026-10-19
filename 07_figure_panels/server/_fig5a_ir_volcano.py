#!/usr/bin/env python3
"""Fig5A IR volcano, authored 1:1 at placed size (~85x64mm, 4:3 ratio kept),
Arial vector PDF, fonts >=7pt, dev-title removed (counts in legend)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
TSV="/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8/data/limma_ir_results_dhx8.tsv"
OUT="/mnt/HDD2/unified_results/postfix_figures_pdf/fig5a_ir_volcano_v1"
df=pd.read_csv(TSV,sep="\t")
lfccol=next(c for c in ["logFC","log2FoldChange","logfc"] if c in df.columns)
pcol=next(c for c in ["adj.P.Val","padj","FDR","adj_p","adj.P"] if c in df.columns)
df=df.dropna(subset=[lfccol,pcol])
y=-np.log10(df[pcol].clip(lower=1e-300))
up=(df[pcol]<0.05)&(df[lfccol]>0); dn=(df[pcol]<0.05)&(df[lfccol]<0); ns=~(up|dn)
nret,nrel=int(up.sum()),int(dn.sum())
C_UP,C_DN,C_NS="#d9534f","#4f81bd","#bdbdbd"
fig,ax=plt.subplots(figsize=(3.35,2.5))   # ~85x64mm, 4:3 kept
ax.scatter(df[lfccol][ns],y[ns],s=4,c=C_NS,alpha=0.5,linewidths=0,rasterized=True)
ax.scatter(df[lfccol][dn],y[dn],s=6,c=C_DN,alpha=0.8,linewidths=0,rasterized=True,label=f"released: {nrel}")
ax.scatter(df[lfccol][up],y[up],s=6,c=C_UP,alpha=0.8,linewidths=0,rasterized=True,label=f"retained: {nret}")
ax.axhline(-np.log10(0.05),ls="--",lw=0.7,color="0.4")
ax.axvline(0,ls="-",lw=0.5,color="0.5")
ax.set_xlabel("logFC (intron retention)",fontsize=8)
ax.set_ylabel("-log10 adjusted P",fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=7.5,loc="lower left",framealpha=0.9,handletextpad=0.3,borderpad=0.4,markerscale=1.4)
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT+".pdf",dpi=600); fig.savefig(OUT+".png",dpi=600)
print(f"WROTE {OUT} | {len(df)} introns | retained {nret} released {nrel} | cols {lfccol}/{pcol}")
