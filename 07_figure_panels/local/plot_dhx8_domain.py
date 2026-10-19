"""
Suppl Fig S13 — Human DHX8/Prp22 domain architecture.

Positions:
  - UniProt Q14562 (length 1220 aa) verified ranges for N-term/disordered,
    S1 motif, RecA1 (Helicase ATP-binding), RecA2 (Helicase C-terminal),
    DEAH-box motif, Walker A mutagenesis site.
  - WH / HB (Ratchet) / OB-fold / Sec63 boundaries derived from the
    human DHX8 helicase-core crystal structure of Felisberto-Rodrigues
    et al. 2019 (PDB 6HYT/6HYS/6HYU, Biochem J 476:2521) and the
    homologous yeast Prp43 architecture (Tauchert et al. 2017, eLife
    6:e21510).
  - R620A is in RecA1, between motif I (GKT/Walker A, 588-595) and
    motif II (DEAH, 685-688).

Output: sfig_dhx8_domain.png and sfig_dhx8_domain.pdf in the same folder.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
LENGTH = 1220  # DHX8 / UniProt Q14562

# (start, end, label, color)
# Boundaries: S1, RecA1, RecA2 and the DEAH motif are from UniProt Q14562.
# Disordered regions (153-266, 337-372) are also UniProt-annotated.
# WH/HB/OB/Sec63 boundaries are approximate values consistent with the
# yeast Prp22 / DHX8 crystal structures (Tauchert 2017; Hamann 2019, PDB 6HYT).
DOMAINS = [
    (1,    152, "N-term\n(disorder)",  "#cfd8dc"),
    (153,  264, "Disordered",          "#eceff1"),
    (265,  336, "S1",                  "#90caf9"),
    (337,  574, "Linker",              "#eceff1"),
    (575,  738, "RecA1\n(DExxH)",      "#4fc3f7"),
    (739,  755, "",                    "#eceff1"),
    (756,  936, "RecA2\n(HA2)",        "#81c784"),
    (937, 1010, "WH",                  "#ffe082"),
    (1011,1090, "HB\n(Ratchet)",       "#ffb74d"),
    (1091,1170, "OB-fold",             "#e57373"),
    (1171,1220, "CTD\n(Sec63)",        "#ba68c8"),
]

# DEAH-box helicase motifs (approximate positions for human DHX8;
# motif I/GKT and II/DEAH are from UniProt, others are placed at
# canonical DEAH-box motif spacings).
MOTIFS = [
    ( 591, "I (GKT)"),     # Walker A, UniProt motif 588-595
    ( 615, "Ia"),
    ( 645, "Ib"),
    ( 687, "II (DEAH)"),   # UniProt motif 685-688
    ( 715, "III (SAT)"),
    ( 778, "IV"),
    ( 845, "V"),
    ( 898, "VI (QRxGR)"),
]

# Mutation
R620A_POS = 620


def main():
    fig, ax = plt.subplots(figsize=(13, 3.6), dpi=200)
    y0, h = 0.0, 1.0

    # draw domain blocks
    for s, e, lab, color in DOMAINS:
        ax.add_patch(mpatches.Rectangle((s, y0), e - s + 1, h,
                                        facecolor=color, edgecolor="#37474f",
                                        linewidth=0.6))
        if lab:
            ax.text((s + e) / 2, y0 + h / 2, lab,
                    ha="center", va="center", fontsize=8.5,
                    color="#1a237e" if color != "#eceff1" else "#90a4ae")

    # helicase motif triangles + labels (alternating staggered rows below the bar)
    for i, (pos, lab) in enumerate(MOTIFS):
        ax.plot(pos, y0 - 0.05, marker="^", markersize=8,
                color="#212121", clip_on=False)
        # alternate two rows of labels so adjacent motifs do not collide
        y_label = y0 - 0.45 if i % 2 == 0 else y0 - 0.85
        # small connector tick
        ax.plot([pos, pos], [y0 - 0.10, y_label + 0.08], color="#9e9e9e",
                lw=0.4, clip_on=False)
        ax.text(pos, y_label, lab, ha="center", va="top",
                fontsize=7.5, color="#212121")

    # R620A red star (above the bar)
    ax.plot(R620A_POS, y0 + h + 0.05, marker="*", markersize=18,
            color="#d32f2f", clip_on=False)
    ax.annotate("R620A\n(helicase-deficient)",
                xy=(R620A_POS, y0 + h + 0.1),
                xytext=(R620A_POS, y0 + h + 0.85),
                ha="center", va="bottom", fontsize=10, color="#d32f2f",
                weight="bold",
                arrowprops=dict(arrowstyle="-", color="#d32f2f",
                                lw=0.8, shrinkA=0, shrinkB=2))

    # axis cosmetics
    ax.set_xlim(-15, LENGTH + 15)
    ax.set_ylim(-1.6, 2.6)
    ax.set_yticks([])

    # custom x-axis below the motif labels
    ax.set_xticks([])  # remove default ticks (they collide with motifs)
    for x in [1, 200, 400, 600, 800, 1000, 1220]:
        ax.plot([x, x], [-1.18, -1.10], color="#37474f", lw=0.6,
                clip_on=False)
        ax.text(x, -1.30, str(x), ha="center", va="top", fontsize=8,
                color="#37474f")
    ax.plot([0, LENGTH], [-1.10, -1.10], color="#37474f", lw=0.6,
            clip_on=False)
    ax.text(LENGTH / 2, -1.55, "Amino acid position (UniProt Q14562, 1220 aa)",
            ha="center", va="top", fontsize=9.5, color="#37474f")

    ax.set_title("Human DHX8 / Prp22 domain architecture",
                 fontsize=12, weight="bold", pad=14)

    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = os.path.join(OUT_DIR, f"sfig_dhx8_domain.{ext}")
        plt.savefig(out, dpi=300, bbox_inches="tight")
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
