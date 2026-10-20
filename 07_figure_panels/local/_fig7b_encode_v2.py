#!/usr/bin/env python3
"""Fig7B ENCODE table, authored 1:1 at placed size (~152x63mm, ratio 2.4 kept),
fonts bumped so they render ~8-9pt at the placed size."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"]="Arial"; plt.rcParams["font.sans-serif"]=["Arial"]
plt.rcParams["pdf.fonttype"]=42; plt.rcParams["ps.fonttype"]=42
OUT="fig7b_encode_comparison_v2"
# Values: rMATS JCEC events with FDR < 0.05 and |IncLevelDifference| > 0.1 in the ENCODE knockdowns
# (ncnp-linux:/mnt/HDD2/unified_results/encode_rbp_kd/splicing/rmats_<factor>/); RI% = RI / (SE + A3SS + A5SS + RI),
# MXE excluded from the denominator. Computed by 02_splicing/rmats/encode_ri_percent_fig7d.py in the code release
# (run 2026-10-20, streamed read-only to ncnp-linux). All SE and RI counts reproduced; PRPF8 and EFTUD2 RI% corrected
# 2026-10-20 from 35.0% and 24.4% (earlier hard-coded values) to 35.1% (1306/3726 = 35.05%) and 24.5% (546/2233 = 24.45%).
cols=["Factor","Complex","Cell","SE","RI","RI%"]
rows=[["PRPF8","C*/U5 core","HepG2","1,756","1,306","35.1%"],
      ["EFTUD2","C*/U5","HepG2","1,173","546","24.5%"],
      ["SNRNP200","C*/U5","HepG2","1,900","729","20.9%"],
      ["PRPF6","tri-snRNP","HepG2","737","308","22.4%"],
      ["RBFOX2","RBP (control)","K562","648","250","22.1%"],
      ["RBM39","RBP (control)","K562","5,016","1,601","19.3%"]]
fig,ax=plt.subplots(figsize=(6.0,2.5)); ax.axis("off")
ax.set_title("ENCODE C* complex factor KD: rMATS RI proportion",fontsize=10,fontweight="bold",pad=8)
tbl=ax.table(cellText=rows,colLabels=cols,loc="center",cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1,1.5)
for (r,c),cell in tbl.get_celld().items():
    cell.set_edgecolor("black"); cell.set_linewidth(1.0)
    if r==0:
        cell.set_facecolor("#2b6cb0"); cell.set_text_props(color="white",fontweight="bold")
    elif c==5 and r in (1,2,3):
        cell.set_facecolor("#fbe0e6"); cell.set_text_props(fontweight="bold")
fig.tight_layout()
fig.savefig(OUT+".pdf"); fig.savefig(OUT+".png",dpi=300)
print("WROTE",OUT)
