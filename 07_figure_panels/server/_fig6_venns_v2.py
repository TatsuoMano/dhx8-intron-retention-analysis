#!/usr/bin/env python3
"""Fig6 venns v2 (siCont_2-corrected, data-driven).
Reads all Venn subset counts + odds ratios directly from
sg_analysis/gsea_dhx8/sg_rna_ir_overlap.tsv (regenerated 2026-07-07 after the
corrupt siCont_2 sample was restored), so no counts are hardcoded.
Layout/style/output-names preserved from the pre-correction version.
Biancon med-high panel: nudge colliding overlap/ir_only labels apart."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
OUTD="/mnt/HDD2/unified_results/postfix_figures_pdf"
TSV="/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8/sg_rna_ir_overlap.tsv"
IRC="#ec9a87"; BLUE="#7bb3d9"; PURPLE=GREEN=TAN=BLUE
F_SUB,F_SET,F_TITLE=8,6.5,8.5

df=pd.read_csv(TSV,sep="\t").set_index("gene_set")
def row(name):
    r=df.loc[name]
    ir=int(r["ir_in_expressed"]); sg=int(r["sg_in_expressed"]); ov=int(r["overlap"])
    return ir-ov, ov, sg-ov, float(r["odds_ratio"]), int(ir), float(r["adj_p"])
def stars(p):
    return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"

def draw(ax,ir_only,ov,sg_only,ir_label,sg_label,sg_color,or_text):
    v=venn2(subsets=(ir_only,sg_only,ov),set_labels=(ir_label,sg_label),
            set_colors=(IRC,sg_color),alpha=0.65,ax=ax)
    for tid in ("10","01","11"):
        t=v.get_label_by_id(tid)
        if t: t.set_fontsize(F_SUB)
    for sl in v.set_labels:
        if sl: sl.set_fontsize(F_SET)
    ax.set_title(or_text,fontsize=F_TITLE,fontweight="bold")
    return v
def nudge(v,tid,dx,dy):
    t=v.get_label_by_id(tid)
    if t: x,y=t.get_position(); t.set_position((x+dx,y+dy))
def single(name,setname,sg_label,sg_color):
    ir_only,ov,sg_only,orv,irtot,adjp=row(setname)
    fig,ax=plt.subplots(figsize=(1.9,1.6))
    draw(ax,ir_only,ov,sg_only,f"IR genes\n({irtot:,})",sg_label,sg_color,
         f"OR = {orv:.2f}  {stars(adjp)}")
    fig.tight_layout(); fig.savefig(f"{OUTD}/{name}_v2.pdf"); fig.savefig(f"{OUTD}/{name}_v2.png",dpi=300)
    print("WROTE",name+"_v2","OR=%.3f"%orv,"subsets=",(ir_only,ov,sg_only),"IR=",irtot)

single("fig6a_venn_khong",      "SG RNA Enriched (Khong 2017)",       "SG RNA Enriched\n(Khong 2017)",   BLUE)
single("fig6c_venn_moon",       "SG-proximal (Moon 2023 sorbitol)",   "SG-proximal Sorb\n(Ren 2023)",    TAN)
single("fig6d_venn_namkoong",   "RG Enriched C1-2 (Namkoong 2018 AS)","RG Enriched AS\n(Namkoong 2018)", BLUE)
single("sfig_venn_namkoong_hs", "RG Enriched C1-2 (Namkoong 2018 HS)","RG Enriched HS\n(Namkoong 2018)", BLUE)

# Biancon two-panel (high + med-high)
hi=row("SG RNA High Consensus (Biancon 2025)")
mh=row("SG RNA Med-High Consensus (Biancon 2025)")
fig,axes=plt.subplots(1,2,figsize=(3.8,1.65))
fig.suptitle("Biancon et al. 2025 (RNA) — IR × SG Meta-analysis",fontsize=8,fontweight="bold")
draw(axes[0],hi[0],hi[1],hi[2],f"IR genes\n({hi[4]:,})",
     f"High consensus\n(>=8 datasets, {hi[2]+hi[1]:,})",PURPLE,f"OR = {hi[3]:.2f}  {stars(hi[5])}")
vR=draw(axes[1],mh[0],mh[1],mh[2],f"IR genes\n({mh[4]:,})",
     f"Med-high consensus\n(>=4 datasets, {mh[2]+mh[1]:,})",PURPLE,f"OR = {mh[3]:.2f}  {stars(mh[5])}")
nudge(vR,"10",-0.13,0.0)   # ir_only -> left
nudge(vR,"11",+0.06,0.0)   # overlap -> right (into lens)
fig.tight_layout(rect=[0,0,1,0.9])
fig.savefig(f"{OUTD}/fig6b_venn_biancon_v2.pdf"); fig.savefig(f"{OUTD}/fig6b_venn_biancon_v2.png",dpi=300)
print("WROTE fig6b_venn_biancon_v2 high_OR=%.3f mh_OR=%.3f"%(hi[3],mh[3]))
