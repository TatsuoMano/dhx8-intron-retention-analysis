#!/usr/bin/env python3
"""Fig5B v5: d-value+sig moved INTO each panel title (removes in-plot d-text that
collided with the legend); legend kept top-right with y-headroom clearing the peak."""
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
B="/mnt/HDD2/unified_results"; OUT=f"{B}/postfix_figures_pdf/fig5_rna_quantity_v7"
deseq=pd.read_csv(f"{B}/deseq2_dhx8_kd_vs_ctrl.tsv",sep="\t").dropna(subset=["log2FoldChange","gene"])
lfc=deseq.groupby("gene")["log2FoldChange"].mean()
def rmats_genes(fn):
    df=pd.read_csv(f"{B}/rmats/dhx8_v2/{fn}",sep="\t"); sig=df[df["FDR"]<0.05]
    return set(sig["geneSymbol"].astype(str).str.strip("\"'"))
ir=set(pd.read_csv("/mnt/HDD2/TDP43_comparison/validation_20260226/limma_retained_dhx8_s2.tsv",sep="\t")["gene_name"].dropna().unique())
se=rmats_genes("SE.MATS.JCEC.txt"); ass=rmats_genes("A3SS.MATS.JCEC.txt")|rmats_genes("A5SS.MATS.JCEC.txt")
def cohend(a,b): return (np.mean(a)-np.mean(b))/np.sqrt((np.var(a,ddof=1)+np.var(b,ddof=1))/2)
panels=[("IR",ir,"#d9534f"),("SE",se,"#4f81bd"),("A3SS/A5SS",ass,"#8e6fb0")]
fig,axes=plt.subplots(1,4,figsize=(6.7,2.46))
summary=[]
for i,(name,gset,col) in enumerate(panels):
    aff=lfc[lfc.index.isin(gset)].values; non=lfc[~lfc.index.isin(gset)].values
    d=cohend(aff,non); p=min(1.0,stats.ttest_ind(aff,non,equal_var=False)[1]*3)
    sig="***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
    summary.append((name,len(aff),np.mean(aff),np.std(aff,ddof=1),d,sig))
    ax=axes[i]; xs=np.linspace(-4,4,300)
    kn=gaussian_kde(non)(xs); ka=gaussian_kde(aff)(xs)
    ax.plot(xs,kn,color="0.55",lw=1.1,label=f"non-{name} ({len(non)})")
    ax.plot(xs,ka,color=col,lw=1.5,label=f"{name} ({len(aff)})")
    ax.axvline(np.median(non),color="0.5",ls="--",lw=0.7); ax.axvline(np.median(aff),color=col,ls="--",lw=0.8)
    ax.set_title(f"{name}   d={d:+.3f} {sig}",fontsize=7.5,fontweight="bold",
                 color="black" if sig!="ns" else "0.45")
    ax.set_xlabel("log2 FC",fontsize=7.5); ax.set_ylabel("Density" if i==0 else "",fontsize=7.5)
    ax.legend(fontsize=6,loc="upper right",handlelength=1.0,handletextpad=0.4,borderpad=0.3,labelspacing=0.25)
    ax.set_xlim(-4,4); ax.set_ylim(0, max(kn.max(),ka.max())*1.35); ax.tick_params(labelsize=6.5)
ax=axes[3]
for i,(name,n,m,sd,d,sig) in enumerate(summary):
    col=panels[i][2]; aff=lfc[lfc.index.isin(panels[i][1])].values
    k=min(len(aff),400)
    ax.scatter(np.full(k,i)+np.random.default_rng(0).uniform(-0.12,0.12,k),
               np.random.default_rng(i).choice(aff,k,replace=False),s=2,c=col,alpha=0.22)
    ax.errorbar(i,m,yerr=sd,fmt="o",color=col,capsize=3,ms=5,zorder=5)
    ax.text(i+0.17,m,f"{m:+.3f}",fontsize=6,color=col,va="center",ha="left",fontweight="bold")
    ax.text(i,3.75,sig,ha="center",fontsize=8,color=col if sig!="ns" else "0.5",fontweight="bold")
ax.axhline(0,color="0.5",ls="--",lw=0.7)
ax.set_xticks(range(3)); ax.set_xticklabels([p[0] for p in panels],fontsize=7)
ax.set_ylabel("log2 Fold Change",fontsize=7.5); ax.set_title("Mean ± SD",fontsize=8.5,fontweight="bold")
ax.set_ylim(-3.5,4.2); ax.tick_params(labelsize=6.5)
fig.tight_layout(w_pad=0.6)
fig.savefig(OUT+".pdf"); fig.savefig(OUT+".png",dpi=200)
print("WROTE",OUT)
print("=== SUMMARY (name, n, mean_lfc, sd, cohend, sig) ===")
for row in summary:
    print(f"{row[0]:>10}  n={row[1]:>5}  mean={row[2]:+.4f}  sd={row[3]:.4f}  d={row[4]:+.4f}  {row[5]}")
