#!/usr/bin/env python3
"""Fig3B: C2:CP GSEA, 10 up + 10 down, sig_only(padj<0.05), full names (no shorten,
single line, left margin auto-fit), acronym-corrected, authored 1:1 (~108x59mm)."""
import re, csv, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
plt.rcParams["axes.unicode_minus"]=False
UP,DN="#e07b5a","#5b9bd5"
SPLICE=re.compile(r"splic|spliceosom|snrnp|mrna.process|pre.?mrna|capped.intron|transesterif",re.I)
TSV="/mnt/HDD2/unified_results/gsea/fgsea_c2cp_dhx8.tsv"
OUT="/mnt/HDD2/unified_results/postfix_figures_pdf/fig3b_gsea_c2cp_v5_plain"
ACR={"Mrna":"mRNA","Rna":"RNA","Snrnp":"snRNP","Trna":"tRNA","Ncam":"NCAM","Ncam1":"NCAM1",
     "Ecm":"ECM","Tcr":"TCR","Dvl":"DVL","Orc1":"ORC1","Auf1":"AUF1","Hnrnp":"hnRNP",
     "Nf":"NF","Kb":"κB","Kappab":"κB","Fceri":"FcεRI","Dna":"DNA","Atp":"ATP","D0":"D0"}
def clean(p):
    s=re.sub(r"^(GOBP|HALLMARK|GO|REACTOME|KEGG|WP|BIOCARTA|PID|NABA|ST|SA|SIG)_","",p).replace("_"," ").title()
    s=" ".join(ACR.get(w,w) for w in s.split())
    return s.replace("NF κB","NF-κB")
rows=[]
for d in csv.DictReader(open(TSV),delimiter="\t"):
    try: nes=float(d.get("NES") or "nan"); pa=float(d.get("padj") or "nan")
    except: continue
    if nes!=nes: continue
    rows.append((d.get("pathway") or list(d.values())[0],nes,pa))
sig=[r for r in rows if r[2]==r[2] and r[2]<0.05]
up=sorted([r for r in sig if r[1]>0],key=lambda x:-x[1])[:10]
dn=sorted([r for r in sig if r[1]<0],key=lambda x:x[1])[:10]
sel=sorted(up+dn,key=lambda x:x[1])
W_mm,H_mm=108,84
fig,ax=plt.subplots(figsize=(W_mm/25.4,H_mm/25.4))
nes=[r[1] for r in sel]; col=[UP if v>0 else DN for v in nes]
ax.barh(range(len(sel)),nes,color=col,edgecolor="none",height=0.66)
ax.set_yticks(range(len(sel))); ax.set_yticklabels([clean(r[0]) for r in sel],fontsize=7)
for tick,r in zip(ax.get_yticklabels(),sel):
    if SPLICE.search(r[0]): pass
ax.axvline(0,color="0.3",lw=0.7)
ax.set_xlabel("NES",fontsize=8); ax.tick_params(axis="x",labelsize=7)
ax.spines[["top","right"]].set_visible(False)
ax.set_ylim(-0.6,len(sel)-0.4)
fig.tight_layout()
fig.savefig(OUT+".pdf"); fig.savefig(OUT+".png",dpi=300)
print("WROTE",OUT,"| max name len:",max(len(clean(r[0])) for r in sel))
