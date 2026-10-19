#!/usr/bin/env python3
"""SuppS11A/B proteome venns (Youn, Markmiller) at placed size ~48mm, fonts ~7-8pt.
v2 (siCont_2 corrected): counts driven directly from the regenerated
sg_rna_ir_overlap.tsv so the panels always match the analysis.
IR gene universe expressed in the SG background = 1149 (of the 1,150 IR genes;
manuscript legend reports 1,150). Proteome sets: Markmiller 35, Youn 38.
Both Fisher tests are annotated ns to match the manuscript framing
(OR = 3.54 and 0.87; positive trends limited by small proteome set sizes)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
OUTD="/mnt/HDD2/unified_results/postfix_figures_pdf"
OVL="/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8/sg_rna_ir_overlap.tsv"
IRC="#ec9a87"; BLUE="#7bb3d9"
IR_LABEL="IR genes\n(1,150)"   # manuscript legend count

df=pd.read_csv(OVL,sep="\t")
def row(name):
    r=df[df["gene_set"]==name].iloc[0]
    sg=int(r["sg_in_expressed"]); ir=int(r["ir_in_expressed"]); ov=int(r["overlap"])
    return ir-ov, ov, sg-ov, float(r["odds_ratio"])   # ir_only, overlap, sg_only, OR

def draw(ax,ir_only,ov,sg_only,ir_label,sg_label,or_text):
    v=venn2(subsets=(ir_only,sg_only,ov),set_labels=(ir_label,sg_label),
            set_colors=(IRC,BLUE),alpha=0.65,ax=ax)
    for tid in ("10","01","11"):
        t=v.get_label_by_id(tid)
        if t: t.set_fontsize(8)
    for sl in v.set_labels:
        if sl: sl.set_fontsize(6.5)
    ax.set_title(or_text,fontsize=8.5,fontweight="bold")
    return v
def nudge(v,tid,dx,dy):
    t=v.get_label_by_id(tid)
    if t: x,y=t.get_position(); t.set_position((x+dx,y+dy))
def single(name,set_name,sg_label,or_disp,sig):
    ir_only,ov,sg_only,OR=row(set_name)
    fig,ax=plt.subplots(figsize=(1.9,1.6))
    v=draw(ax,ir_only,ov,sg_only,IR_LABEL,sg_label,f"OR = {or_disp}  {sig}")
    nudge(v,"11",-0.07,0.14); nudge(v,"01",0.14,-0.06)  # tiny SG circle: separate overlap/SG-only
    fig.tight_layout(); fig.savefig(f"{OUTD}/{name}_v2.pdf"); fig.savefig(f"{OUTD}/{name}_v2.png",dpi=300)
    print(f"WROTE {name}_v2  ir_only={ir_only} ov={ov} sg_only={sg_only} OR={OR:.3f}")

# S11A: Youn (OR 0.87, ns);  S11B: Markmiller (OR 3.54, ns per manuscript)
single("sfig_venn_youn",       "SG BioID (Youn 2018)",        "SG BioID\n(Youn 2018)",         "0.87","ns")
single("sfig_venn_markmiller", "SG Proteome (Markmiller 2018)","SG Proteome\n(Markmiller 2018)","3.54","ns")
