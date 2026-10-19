#!/usr/bin/env python3
"""Faithful Arial vector-PDF re-creation of Fig 4B (= fig4e_integrated_summary_3ds):
Integrated splicing analysis, 3 tools x 3 datasets (DHX8 KD + GSE65973 + GSE126543).
All values fixed/known. Changes vs prior _3ds PNG/PDF:
  - descriptive 2-line dataset column headers (TK review #41/#42)
  - column order swapped to Event, Tool (was Tool, Event)
Style: blue header, pink DHX8 KD column (Tool text plain black)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"; plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["ps.fonttype"] = 42
OUT = "fig4e_integrated_summary_3ds_v2"

# Event, Tool, then the three datasets with descriptive 2-line headers
cols = ["Event", "Tool",
        "DHX8 KD\n(this study, MEF)",
        "GSE65973\n(dNLS-TDP-43,\nmouse cortex)",
        "GSE126543\n(FANS-sorted TDP-43±\nnuclei, human ALS)"]
rows = [   # Event, Tool, DHX8 KD (v2, corrected-siCont_2), GSE65973, GSE126543
    ["SE",           "rMATS",        "4,204", "204", "667"],
    ["A3SS",         "rMATS",        "599",   "42",  "170"],
    ["A5SS",         "rMATS",        "458",   "36",  "86"],
    ["IR (introns)", "limma eBayes", "1,258", "0",   "0"],
    ["IR (genes)",   "limma eBayes", "1,150", "0",   "0"],
    ["CE (strong)",  "MAJIQ",        "1",     "0",   "9"],
]
widths = [0.14, 0.16, 0.20, 0.25, 0.25]   # Event/Tool narrow; GSE columns wide for descriptors

fig, ax = plt.subplots(figsize=(6.2, 2.85)); ax.axis("off")
ax.set_title("Integrated splicing analysis: 3 tools x 3 datasets",
             fontsize=11, fontweight="bold", pad=10)
tbl = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.55)
for (r, c), cell in tbl.get_celld().items():
    cell.set_width(widths[c])
    cell.set_edgecolor("black"); cell.set_linewidth(0.8)
    if r == 0:                                       # header (taller for 3-line labels)
        cell.set_height(cell.get_height() * 2.4)
        cell.set_facecolor("#2b6cb0"); cell.set_text_props(color="white", fontweight="bold")
    else:
        if c == 2:                                   # DHX8 KD column highlight (pink + bold)
            cell.set_facecolor("#fbe0e6"); cell.set_text_props(fontweight="bold")
fig.tight_layout()
fig.savefig(OUT + ".pdf"); fig.savefig(OUT + ".png", dpi=200)
print("WROTE", OUT + ".pdf")
