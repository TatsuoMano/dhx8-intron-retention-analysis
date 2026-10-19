"""Regenerate Fig 6F (a.k.a. fig6e_sg_overlap_summary) 7-row summary table.

v2 (siCont_2-corrected, 2026-07-07): a corrupt siCont_2 RNA-seq sample was
restored and all analyses re-run. This table now reflects the corrected v2
numbers:
  - IR OR from sg_analysis/gsea_dhx8/sg_rna_ir_overlap.tsv (Khong 1.86,
    Namkoong AS 1.07, HS 1.21, Ren/Sorb 2.17, Biancon high 2.14).
  - DEG-UP OR is now the STRICT comparator (padj<0.05 & log2FC>1, n=652),
    matching the manuscript Fig 6F legend. The strongly upregulated genes are
    DEPLETED among SG-enriched IR datasets (Khong 0.57, Biancon high 0.46,
    Biancon med-high 0.35, Ren 0.42), contrasting the IR-gene enrichment.
  - Proteome rows (Markmiller/Youn) reported via GSEA NES (small Fisher sets).

The reproducible, data-driven generator lives on ncnp-linux as
/mnt/HDD2/unified_results (fig6f_table_v2.py). It reads the v2 TSVs directly and
writes fig6e_sg_overlap_summary_ren7row_v2.{pdf,png} plus
deg_up_threshold_comparison_v2.tsv. The arrays below are the frozen v2 values
for local provenance.
"""
import os
import matplotlib
matplotlib.rcParams["font.family"] = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT = os.path.dirname(os.path.abspath(__file__))

# Header + 7 data rows
# p-value column now shows "IR p / DEG-UP p" uniformly across all rows.
# IR p-values are from the manuscript (Fig 6 / Suppl Fig S8 legends);
# DEG-UP p-values are from Analysis/2026-03-30_sg_transcriptome_reanalysis/sg_deg_overlap.tsv.
# --- v2 (siCont_2-corrected, 2026-07-07) ---
# IR OR from sg_analysis/gsea_dhx8/sg_rna_ir_overlap.tsv; DEG-UP OR is now the
# STRICT comparator (padj<0.05 & log2FC>1, n=652) recomputed from the corrected
# v2 DESeq2, matching the manuscript Fig 6F legend. Proteome rows via GSEA (NES).
# Reproducible generator: /mnt/HDD2/unified_results (fig6f_table_v2.py) writing
# fig6e_sg_overlap_summary_ren7row_v2.{pdf,png} + deg_up_threshold_comparison_v2.tsv.
HEADER = ["Dataset", "Method", "Type", "IR enrichment", "DEG-UP OR (log2FC>1)", "p (IR / DEG-UP)"]
# Transcriptome (main-figure focus) first, then Proteome (supplementary).
# IR enrichment metric is stated per cell: Fisher OR for transcriptome, GSEA NES for proteome.
ROWS = [
    ["Khong 2017",          "SG core IP RNA-seq",    "Transcriptome", "OR = 1.86",  "0.57", "<0.001 / <0.001"],
    ["Namkoong 2018 (AS)",  "Insoluble RG fraction", "Transcriptome", "OR = 1.07",  "1.27", "0.379 / 0.176"],
    ["Namkoong 2018 (HS)",  "Insoluble RG fraction", "Transcriptome", "OR = 1.21",  "1.00", "0.014 / 1.000"],
    ["Ren 2023 (Sorb)",     "CAP-seq",               "Transcriptome", "OR = 2.17",  "0.42", "<0.001 / <0.001"],
    ["Biancon 2025 (high)", "26-dataset meta",       "Transcriptome", "OR = 2.14",  "0.46", "<0.001 / <0.001"],
    ["Markmiller 2018",     "SG IP + MS",            "Proteome",      "NES = 1.83", "GSEA", "<0.001 (GSEA)"],
    ["Youn 2018",           "G3BP1-BioID",           "Proteome",      "NES = 1.39", "GSEA", "0.031 (GSEA)"],
]

# No color highlight (plain table, per request)
HIGHLIGHT_ROWS = []
HIGHLIGHT_COL  = 3                  # "IR enrichment"

HEADER_BG  = "#2c5d8a"     # blue
HEADER_FG  = "white"
ALT_BG     = "#ffffff"
HIGHLIGHT_BG = "#fce8e6"   # light salmon
HIGHLIGHT_FG = "#c44d4d"   # dark red
TITLE = "Fig 6F: IR gene overlap with SG datasets (DEG-UP = strongly upregulated, log2FC>1, n=652)"


def main():
    n_rows = len(ROWS) + 1   # +1 for header
    n_cols = len(HEADER)

    # Per-column widths (proportional). "DEG-UP OR (log2FC>1)" header is wide,
    # so give the IR OR and DEG-UP OR columns more room; p-value column holds
    # "<0.001 / <0.001".
    col_widths = [0.18, 0.18, 0.11, 0.12, 0.18, 0.23]
    # normalize to sum to 1
    s = sum(col_widths)
    col_widths = [w / s for w in col_widths]

    fig_w, fig_h = 12.0, 0.45 * n_rows + 0.6
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, 1); ax.set_ylim(0, n_rows)
    ax.invert_yaxis()
    ax.axis("off")

    # Title
    ax.text(0.5, -0.6, TITLE, ha="center", va="top",
            fontsize=12, weight="bold")

    # Header row
    x = 0
    for ci, (label, w) in enumerate(zip(HEADER, col_widths)):
        ax.add_patch(mpatches.Rectangle((x, 0), w, 1, facecolor=HEADER_BG,
                                        edgecolor="white", linewidth=1))
        ax.text(x + w / 2, 0.5, label, ha="center", va="center",
                color=HEADER_FG, fontsize=10, weight="bold")
        x += w

    # Data rows
    for ri, row in enumerate(ROWS, start=1):
        x = 0
        for ci, (cell, w) in enumerate(zip(row, col_widths)):
            is_highlight = (ri - 1 in HIGHLIGHT_ROWS) and (ci == HIGHLIGHT_COL)
            bg = HIGHLIGHT_BG if is_highlight else ALT_BG
            fg = HIGHLIGHT_FG if is_highlight else "black"
            weight = "bold" if is_highlight else "normal"
            ax.add_patch(mpatches.Rectangle((x, ri), w, 1, facecolor=bg,
                                            edgecolor="#bcbcbc", linewidth=0.6))
            ax.text(x + w / 2, ri + 0.5, cell, ha="center", va="center",
                    color=fg, fontsize=10.5, weight=weight)
            x += w

    plt.tight_layout()
    for ext in ("png", "pdf"):
        path = os.path.join(OUT, f"fig6e_sg_overlap_summary_ren7row.{ext}")
        fig.savefig(path, dpi=200 if ext == "png" else 300, bbox_inches="tight")
        print(f"  wrote {path}")
    plt.close(fig)

    print("\nInject into the manuscript with:")
    print("  cd ~/Dropbox/Research/DHX8-Sakai/Manuscript")
    print("  mscript-fig-update \\")
    print("    --src ../Analysis/2026-05-14_relabel_moon_to_ren/fig6e_sg_overlap_summary_ren7row.png \\")
    print("    --target Figures/figures_for_final/Fig06/fig6e_sg_overlap_summary.png \\")
    print("    --note 'ren_label_7rows_with_namkoong_as'")


if __name__ == "__main__":
    main()
