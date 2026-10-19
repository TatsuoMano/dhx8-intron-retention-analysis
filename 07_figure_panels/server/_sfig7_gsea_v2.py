#!/usr/bin/env python3
"""SuppFig S7 (GO:BP + Hallmark) — Fig3B recipe: 10 up + 10 down, sig_only(padj<0.05),
full names (no shorten), plain black labels, authored 1:1, fonts ~7pt."""
import re, csv, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
B="/mnt/HDD2/unified_results/gsea"; OUTD="/mnt/HDD2/unified_results/postfix_figures_pdf"
UP,DN="#e07b5a","#5b9bd5"
ACR={"Mrna":"mRNA","Rna":"RNA","Snrnp":"snRNP","Trna":"tRNA","Ncam":"NCAM","Ecm":"ECM",
     "Tcr":"TCR","Dna":"DNA","Atp":"ATP","Nf":"NF","Kb":"κB","Kappab":"κB","Dvl":"DVL",
     "Orc1":"ORC1","Auf1":"AUF1","Hnrnp":"hnRNP","P53":"P53","G2m":"G2M","Myc":"MYC",
     "Tnfa":"TNFα","Il2":"IL2","Stat5":"STAT5","Kras":"KRAS","Uv":"UV","Dna4":"DNA"}
def clean(p):
    s=re.sub(r"^(GOBP|HALLMARK|GO|REACTOME|KEGG|WP|BIOCARTA|PID|NABA|ST|SA|SIG)_","",p).replace("_"," ").title()
    s=" ".join(ACR.get(w,w) for w in s.split())
    return s.replace("NF κB","NF-κB")
def load(tsv):
    rows=[]
    for d in csv.DictReader(open(f"{B}/{tsv}"),delimiter="\t"):
        try: nes=float(d.get("NES") or "nan"); pa=float(d.get("padj") or "nan")
        except: continue
        if nes!=nes: continue
        rows.append((d.get("pathway") or list(d.values())[0],nes,pa))
    return rows
def panel(tsv,out,W_mm):
    rows=load(tsv); sig=[r for r in rows if r[2]==r[2] and r[2]<0.05]
    up=sorted([r for r in sig if r[1]>0],key=lambda x:-x[1])[:10]
    dn=sorted([r for r in sig if r[1]<0],key=lambda x:x[1])[:10]
    sel=sorted(up+dn,key=lambda x:x[1])
    fig,ax=plt.subplots(figsize=(W_mm/25.4,84/25.4))
    nes=[r[1] for r in sel]; col=[UP if v>0 else DN for v in nes]
    ax.barh(range(len(sel)),nes,color=col,edgecolor="none",height=0.72)
    ax.set_yticks(range(len(sel))); ax.set_yticklabels([clean(r[0]) for r in sel],fontsize=7)
    ax.axvline(0,color="0.3",lw=0.7); ax.set_xlabel("NES",fontsize=8); ax.tick_params(axis="x",labelsize=7)
    ax.spines[["top","right"]].set_visible(False); ax.set_ylim(-0.6,len(sel)-0.4)
    fig.tight_layout(); fig.savefig(out+".pdf"); fig.savefig(out+".png",dpi=300)
    print("WROTE",out,"| sig",len(sig),"| maxname",max(len(clean(r[0])) for r in sel))
panel("fgsea_gobp_dhx8.tsv",     f"{OUTD}/sfig_gsea_gobp_v2",     120)
panel("fgsea_hallmark_dhx8.tsv", f"{OUTD}/sfig_gsea_hallmark_v2", 82)
