#!/usr/bin/env python3
"""S12 (motif & sequence-feature negative result) re-plotted in the ORIGINAL v1
composite style (1x3: RBP-motif Cohen's d forest / GC-content boxplot / length-
matching histogram), updated with v2 (corrected-siCont_2) data. Reproduces the
look of Figures/figures_for_final/SupFig10/sfig_motif_features.png; k-mer volcano
and G-run/complexity sub-panels are NOT included (reported in text only).

Data (dumped by analysis_motif_enrichment.py on ncnp-linux):
  motif_enrichment_results.tsv       - RBP motif Cohen's d (panel A)
  s12_per_intron_features_v2.tsv     - per-intron GC content by group (panel B)
  s12_length_arrays_v2.tsv           - per-intron length by group (panel C)
  s12_ks_v2.txt                      - length-matching KS (D, p)
Arial, pdf.fonttype=42, literal Unicode. Sized 1:1 for placement (~180x58 mm).
"""
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.expanduser("~/bin"))
import figprov  # provenance-emitting save (script/data hashes, git, env -> .prov.json)
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["axes.unicode_minus"] = False

C_RET = "#dd7a6b"   # retained (coral)
C_CTRL = "#b7bcc0"  # control (gray)
C_BAR = "#8f9296"   # forest bars (neutral gray)

# ---- load ----
res = pd.read_csv("motif_enrichment_results.tsv", sep="\t")
feat = pd.read_csv("s12_per_intron_features_v2.tsv", sep="\t")
length = pd.read_csv("s12_length_arrays_v2.tsv", sep="\t")
ks = dict(x.split("=") for x in open("s12_ks_v2.txt").read().strip().split("\t"))
ks_D = float(ks["D"]); ks_p = float(ks["p"])

LAB = {"TDP43": "TDP-43", "FUS": "FUS", "MBNL1": "MBNL1", "SF1": "SF1",
       "SRSF1": "SRSF1", "SRSF2": "SRSF2", "SRSF3": "SRSF3", "SRSF6": "SRSF6",
       "U2AF65": "U2AF65", "hnRNPA1": "hnRNP A1", "hnRNPA2B1": "hnRNP A2B1",
       "hnRNPC": "hnRNP C", "hnRNPF": "hnRNP F", "hnRNPH": "hnRNP H"}

m = res[res["feature"].str.startswith("motif_")].copy()
m["rbp"] = m["feature"].str.replace("motif_", "", regex=False).str.replace("_perkb", "", regex=False)
m = m.sort_values("cohens_d")
n_sig = int(m["significant_bonf"].astype(str).str.lower().isin(["true", "1"]).sum())

gc_ret = feat.loc[feat["category"] == "retained", "gc_content"].dropna().values
gc_ctrl = feat.loc[feat["category"] != "retained", "gc_content"].dropna().values
gc_row = res[res["feature"] == "gc_content"].iloc[0]
gc_p, gc_d = float(gc_row["mwu_p"]), float(gc_row["cohens_d"])

log_ret = np.log10(length.loc[length["category"] == "retained", "intron_length"].values)
log_ctrl = np.log10(length.loc[length["category"] == "control", "intron_length"].values)
n_ret, n_ctrl = len(log_ret), len(log_ctrl)

# ---- figure ----
fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.3))
FS_T, FS_L, FS_TK = 8, 7, 6.5

# (A) RBP motif Cohen's d forest
ax = axes[0]
y = np.arange(len(m))
ax.barh(y, m["cohens_d"].values, color=C_BAR, edgecolor="black", linewidth=0.4, height=0.7)
ax.set_yticks(y); ax.set_yticklabels([LAB.get(r, r) for r in m["rbp"]], fontsize=FS_TK)
ax.axvline(0, color="black", lw=0.7)
# neutral reference at Cohen's d = ±0.2 (recognised small-effect threshold), gray/dotted
# so it reads as a magnitude reference, not a significance cutoff.
ax.axvline(-0.2, color="0.6", ls=":", lw=0.7)
ax.axvline(0.2, color="0.6", ls=":", lw=0.7)
ax.set_xlim(-0.22, 0.22)
ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
ax.set_xlabel("Cohen's d", fontsize=FS_L)
ax.set_title("RBP motif enrichment\n(retained vs matched control)", fontsize=FS_T, fontweight="bold")
ax.tick_params(axis="x", labelsize=FS_TK)
ax.text(0.80, 0.50, f"All Bonferroni\nn.s. ({len(m)-n_sig}/{len(m)})",
        transform=ax.transAxes, ha="center", va="center", fontsize=FS_TK, color="#2e7d32")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (B) GC content boxplot
ax = axes[1]
bp = ax.boxplot([gc_ret, gc_ctrl], widths=0.6, patch_artist=True, showfliers=False,
                medianprops=dict(color="#e08214", lw=1.1),
                whiskerprops=dict(color="black", lw=0.7),
                capprops=dict(color="black", lw=0.7),
                boxprops=dict(edgecolor="black", lw=0.6))
bp["boxes"][0].set_facecolor(C_RET); bp["boxes"][1].set_facecolor(C_CTRL)
ax.set_xticks([1, 2])
ax.set_xticklabels([f"Retained\n(n={n_ret:,})", f"Control\n(n={n_ctrl:,})"], fontsize=FS_TK)
ax.set_ylabel("GC content", fontsize=FS_L)
ax.set_title(f"GC content\nMWU P = {gc_p:.2g}, d = {gc_d:+.3f}", fontsize=FS_T, fontweight="bold")
ax.tick_params(axis="y", labelsize=FS_TK)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (C) Length-matching histogram
ax = axes[2]
bins = np.linspace(min(log_ret.min(), log_ctrl.min()), max(log_ret.max(), log_ctrl.max()), 45)
ax.hist(log_ctrl, bins=bins, density=True, color=C_CTRL, alpha=0.75,
        label=f"Control (n={n_ctrl:,})", edgecolor="none")
ax.hist(log_ret, bins=bins, density=True, color=C_RET, alpha=0.6,
        label=f"Retained (n={n_ret:,})", edgecolor="none")
p_txt = "P > 0.999" if ks_p > 0.999 else f"P = {ks_p:.3g}"
ax.set_xlabel("log10(intron length)", fontsize=FS_L)
ax.set_ylabel("Density", fontsize=FS_L)
ax.set_title(f"Length matching\nKS D = {ks_D:.3f}, {p_txt}", fontsize=FS_T, fontweight="bold")
ax.tick_params(labelsize=FS_TK)
ax.legend(fontsize=FS_TK - 0.5, frameon=False, loc="upper left", handlelength=1.0)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

fig.tight_layout(w_pad=1.4)
figprov.save_figure(
    fig, "sfig_motif_features_v2.pdf",
    inputs=["motif_enrichment_results.tsv", "s12_per_intron_features_v2.tsv",
            "s12_length_arrays_v2.tsv", "s12_ks_v2.txt"],
    script=__file__,
    note="S12 (SupFig10) v1-style composite: RBP motif Cohen's d forest / GC boxplot / "
         "length-matching histogram; v2 corrected-siCont_2 data",
)
print(f"WROTE sfig_motif_features_v2  | RBP {len(m)-n_sig}/{len(m)} ns, "
      f"GC MWU P={gc_p:.3g} d={gc_d:+.3f}, KS D={ks_D:.4f} p={ks_p:.4f}, "
      f"n_ret={n_ret} n_ctrl={n_ctrl}")
