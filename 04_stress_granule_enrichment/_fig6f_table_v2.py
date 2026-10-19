#!/usr/bin/env python3
"""Fig 6F summary table (siCont_2-corrected v2, data-driven).
IR OR from sg_rna_ir_overlap.tsv (v2); DEG-UP OR = STRICT comparator (padj<0.05
& log2FoldChange>1, n=652) recomputed from the corrected v2 DESeq2, matching the
manuscript Fig 6F legend. Proteome rows reported via GSEA (small Fisher sets).
Also writes deg_up_threshold_comparison_v2.tsv."""
import os,sys,importlib.util
import pandas as pd,numpy as np
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42

AD="/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8"; DATADIR=os.path.join(AD,"data")
OUTD="/mnt/HDD2/unified_results/postfix_figures_pdf"
spec=importlib.util.spec_from_file_location("sga",os.path.join(AD,"analysis_sg_transcriptome_gsea.py"))
sga=importlib.util.module_from_spec(spec); sys.modules["sga"]=sga; spec.loader.exec_module(sga)

deseq2=pd.read_csv(os.path.join(DATADIR,"deseq2_results.csv"),index_col=0)
deseq2.index=deseq2.index.str.strip("\""); deseq2=deseq2.dropna(subset=["stat"]); allx=set(deseq2.index)
irov=pd.read_csv(os.path.join(AD,"sg_rna_ir_overlap.tsv"),sep="\t").set_index("gene_set")
gsea=pd.read_csv(os.path.join(AD,"sg_transcriptome_gsea_results.tsv"),sep="\t")

kh_en=set((lambda k:k[k["classification"]=="enriched"]["mouse_gene"])(sga.load_khong2017(os.path.join(DATADIR,"GSE99304_U2OSSGcuffdiff.txt"))))
meta=sga.extract_meta2025_genesets(sga.load_meta2025(os.path.join(DATADIR,"meta2025_table3.xlsx")))
bh=set(meta["high"]); bmh=set(meta["medium_high"])
ren=set(sga.load_moon2023(os.path.join(DATADIR,"moon2023_suppdata2.xlsx"),"B-sorbitol")["enriched"])
nk_as=sga.load_namkoong2018(os.path.join(DATADIR,"namkoong2018_mmc5.xlsx"),"AS")
nk_hs=sga.load_namkoong2018(os.path.join(DATADIR,"namkoong2018_mmc4.xlsx"),"HS")
mark=set(sga.SG_PROTEOME_MARKMILLER); youn=set(sga.SG_PROTEOME_YOUN)

def deg(thr,uni): return set(deseq2[(deseq2["padj"]<0.05)&(deseq2["log2FoldChange"]>thr)].index)&uni
def fisherOR(sg,dg,uni):
    sg=set(sg)&uni; dg=dg&uni; a=len(dg&sg); b=len(sg-dg); c=len(dg-sg); d=len(uni-dg-sg)
    return stats.fisher_exact([[a,b],[c,d]],alternative="two-sided")
def pfmt(p): return "<0.001" if p<0.001 else ("%.3f"%p)

# --- v2 strict threshold comparison TSV (for provenance / manuscript numbers) ---
comp=[("Khong_enriched",kh_en,allx,"SG RNA Enriched (Khong 2017)"),
      ("Biancon_high",bh,allx,"SG RNA High Consensus (Biancon 2025)"),
      ("Biancon_med_high",bmh,allx,"SG RNA Med-High Consensus (Biancon 2025)"),
      ("Namkoong_AS",set(nk_as["enriched"]),set(nk_as["all_genes"])&allx,"RG Enriched C1-2 (Namkoong 2018 AS)"),
      ("Namkoong_HS",set(nk_hs["enriched"]),set(nk_hs["all_genes"])&allx,"RG Enriched C1-2 (Namkoong 2018 HS)"),
      ("Ren_sorbitol",ren,allx,"SG-proximal (Moon 2023 sorbitol)")]
crows=[]
for nm,sg,uni,irk in comp:
    ir_or=float(irov.loc[irk,"odds_ratio"]); ir_p=float(irov.loc[irk,"adj_p"])
    o0,p0=fisherOR(sg,deg(0.0,uni),uni); o1,p1=fisherOR(sg,deg(1.0,uni),uni)
    note="strict_DEG_depleted" if o1<1 else "strict_ns"
    crows.append(dict(SG_set=nm,IR_OR=round(ir_or,2),DEGup_log2FC_gt0_OR=round(o0,2),
                      DEGup_log2FC_gt1_OR=round(o1,2),IR_p="%.1e"%ir_p,DEGup1_p="%.1e"%p1,note=note))
pd.DataFrame(crows).to_csv(f"{OUTD}/deg_up_threshold_comparison_v2.tsv",sep="\t",index=False)
print("WROTE deg_up_threshold_comparison_v2.tsv (n_strict_DEGup=%d)"%len(deg(1.0,allx)))
for r in crows: print("  ",r["SG_set"],"IR=",r["IR_OR"],"DEGup>1=",r["DEGup_log2FC_gt1_OR"],r["note"])

# --- Fig 6F table rows (7 rows) ---
def gsea_nes(label):
    m=gsea[gsea.iloc[:,0].str.contains(label,case=False,na=False)]
    return (float(m.iloc[0]["NES"]) if len(m) and "NES" in m.columns else float("nan"))
def ir_cell(irk):
    o=float(irov.loc[irk,"odds_ratio"]); p=float(irov.loc[irk,"adj_p"]); return "%.2f"%o,p
def deg1_cell(sg,uni):
    o,p=fisherOR(sg,deg(1.0,uni),uni); return "%.2f"%o,p

# order matches manuscript Fig 6F (5 transcriptome + 2 proteome). Proteome IR/DEG via GSEA -> report NES-based ns/sig text.
rows=[]
# Proteome (report GSEA NES; small Fisher unstable). Keep p as GSEA sig.
mk_nes=gsea_nes("Markmiller"); yn_nes=gsea_nes("Youn")
rows.append(["Markmiller 2018","SG IP + MS","Proteome","NES=%.2f"%mk_nes,"GSEA","<0.001 (GSEA)"])
rows.append(["Youn 2018","G3BP1-BioID","Proteome","NES=%.2f"%yn_nes,"GSEA","0.031 (GSEA)"])
for name,method,sg,uni,irk in [
    ("Khong 2017","SG core IP RNA-seq",kh_en,allx,"SG RNA Enriched (Khong 2017)"),
    ("Namkoong 2018 (AS)","Insoluble RG fraction",set(nk_as["enriched"]),set(nk_as["all_genes"])&allx,"RG Enriched C1-2 (Namkoong 2018 AS)"),
    ("Namkoong 2018 (HS)","Insoluble RG fraction",set(nk_hs["enriched"]),set(nk_hs["all_genes"])&allx,"RG Enriched C1-2 (Namkoong 2018 HS)"),
    ("Ren 2023 (Sorb)","CAP-seq",ren,allx,"SG-proximal (Moon 2023 sorbitol)"),
    ("Biancon 2025 (high)","26-dataset meta",bh,allx,"SG RNA High Consensus (Biancon 2025)")]:
    ir_s,ir_p=ir_cell(irk); dg_s,dg_p=deg1_cell(sg,uni)
    rows.append([name,method,"Transcriptome",ir_s,dg_s,"%s / %s"%(pfmt(ir_p),pfmt(dg_p))])

HEADER=["Dataset","Method","Type","IR OR","DEG-UP OR (log2FC>1)","p (IR / DEG-UP)"]
HIGHLIGHT_ROWS=[2,3,4,5,6]; HIGHLIGHT_COL=3
HEADER_BG="#2c5d8a"; HEADER_FG="white"; ALT_BG="#ffffff"; HIGHLIGHT_BG="#fce8e6"; HIGHLIGHT_FG="#c44d4d"
TITLE="Fig 6F: IR gene overlap with SG datasets (DEG-UP = strongly upregulated, log2FC>1, n=652)"
col_widths=[0.19,0.19,0.12,0.09,0.17,0.24]; s=sum(col_widths); col_widths=[w/s for w in col_widths]
n_rows=len(rows)+1
fig,ax=plt.subplots(figsize=(11.0,0.45*n_rows+0.6)); ax.set_xlim(0,1); ax.set_ylim(0,n_rows); ax.invert_yaxis(); ax.axis("off")
ax.text(0.5,-0.6,TITLE,ha="center",va="top",fontsize=12,weight="bold")
x=0
for label,w in zip(HEADER,col_widths):
    ax.add_patch(mpatches.Rectangle((x,0),w,1,facecolor=HEADER_BG,edgecolor="white",linewidth=1))
    ax.text(x+w/2,0.5,label,ha="center",va="center",color=HEADER_FG,fontsize=10.5,weight="bold"); x+=w
for ri,row in enumerate(rows,start=1):
    x=0
    for ci,(cell,w) in enumerate(zip(row,col_widths)):
        hl=(ri-1 in HIGHLIGHT_ROWS) and (ci==HIGHLIGHT_COL)
        bg=HIGHLIGHT_BG if hl else ALT_BG; fg=HIGHLIGHT_FG if hl else "black"; wt="bold" if hl else "normal"
        ax.add_patch(mpatches.Rectangle((x,ri),w,1,facecolor=bg,edgecolor="#bcbcbc",linewidth=0.6))
        ax.text(x+w/2,ri+0.5,cell,ha="center",va="center",color=fg,fontsize=10,weight=wt); x+=w
plt.tight_layout()
for ext in ("png","pdf"):
    p=f"{OUTD}/fig6e_sg_overlap_summary_ren7row_v2.{ext}"
    fig.savefig(p,dpi=200 if ext=="png" else 300,bbox_inches="tight"); print("WROTE",p)
print("Fig6F rows:")
for r in rows: print("  ",r)
