#!/usr/bin/env python3
"""
Permutation test (R1 reviewer concern): SG enrichment in IR genes —
controls for BOTH expression (baseMean) AND transcript length jointly.

Rationale:
  - Khong et al. 2017 showed SG-enriched transcripts tend to be long.
  - Longer genes with more introns are more likely to have IR events.
  - Hence transcript length is a candidate confounder for the IR↔SG association.

Strategy:
  - Compute per-gene length = sum of exon lengths of the longest transcript
    from mm10 ncbiRefSeq GTF.
  - For each of 10,000 permutations, sample a null set of len(IR) genes
    matched JOINTLY to the IR gene (log10 baseMean) × (log10 length) distribution
    using a 2-D binning scheme (n_bins x n_bins), with neighbor-cell fallback
    when a cell is undersupplied.
  - Compare OR(null ∩ SG) distribution to observed OR(IR ∩ SG).
  - Report observed OR, mean null OR, 95th percentile, empirical p-value.
  - Repeat for: (1) Khong 2017 SG-enriched, (2) Biancon 2025 High Consensus.

Outputs (alongside original expression-only baseline for comparison):
  - sg_permutation_test_lenmatched.tsv
  - sg_permutation_test_lenmatched.png
  - sg_permutation_diagnostics_lenmatched.png  (matching quality QC)
"""

import os
import re
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, ks_2samp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASEDIR      = "/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8"
IR_FILE      = f"{BASEDIR}/data/limma_retained_dhx8.tsv"
DESEQ_FILE   = "/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv"
KHONG_FILE   = f"{BASEDIR}/data/GSE99304_U2OSSGcuffdiff.txt"
BIANCON_FILE = f"{BASEDIR}/data/meta2025_table3.xlsx"
GTF_FILE     = "/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf"

OUTDIR       = "/mnt/HDD2/unified_results/postfix_figures"
GENE_LEN_TSV = f"{OUTDIR}/mm10_gene_length.tsv"
OUT_TSV      = f"{OUTDIR}/sg_permutation_test_lenmatched.tsv"
OUT_PNG      = f"{OUTDIR}/sg_permutation_test_lenmatched.png"
QC_PNG       = f"{OUTDIR}/sg_permutation_diagnostics_lenmatched.png"

SEED       = 42
N_PERMS    = 10_000
N_BINS_BM  = 10   # 10 x 10 = 100 joint cells
N_BINS_LEN = 10


# ── Gene symbol utilities ──────────────────────────────────────────────────────
def human_to_mouse(symbol):
    if not isinstance(symbol, str) or len(symbol.strip()) == 0:
        return None
    s = symbol.strip()
    return s[0].upper() + s[1:].lower()


# ── 1. Build per-gene length table from GTF (longest transcript) ──────────────
def build_gene_length_table(gtf_path, out_tsv):
    """For each gene, length = sum of exon lengths of its longest transcript."""
    if os.path.exists(out_tsv):
        print(f"  Loading cached gene length table from {out_tsv}")
        return pd.read_csv(out_tsv, sep="\t")

    print(f"  Parsing GTF: {gtf_path}")
    gn_re = re.compile(r'gene_name "([^"]+)"')
    tx_re = re.compile(r'transcript_id "([^"]+)"')

    # tx_len[(gene, tx)] = sum of exon widths
    tx_len  = {}
    tx_gene = {}

    t0 = time.time()
    n_lines = 0
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            start = int(parts[3])
            end   = int(parts[4])
            attrs = parts[8]
            mg = gn_re.search(attrs)
            mt = tx_re.search(attrs)
            if not (mg and mt):
                continue
            gene = mg.group(1)
            tx   = mt.group(1)
            key  = (gene, tx)
            tx_len[key]  = tx_len.get(key, 0) + (end - start + 1)
            tx_gene[tx]  = gene
            n_lines += 1

    print(f"  Parsed {n_lines:,} exon lines in {time.time()-t0:.1f}s, {len(tx_len):,} transcripts")

    # Per gene: longest transcript length
    gene_max = {}
    for (gene, tx), L in tx_len.items():
        if L > gene_max.get(gene, -1):
            gene_max[gene] = L

    df = pd.DataFrame({
        "gene":   list(gene_max.keys()),
        "length": list(gene_max.values()),
    })
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"  Wrote {len(df):,} genes → {out_tsv}")
    return df


# ── 2. OR helper ──────────────────────────────────────────────────────────────
def or_from_sets(test_set, sg_set, universe_set):
    ts  = test_set & universe_set
    sg  = sg_set & universe_set
    not_ts = universe_set - ts
    a = len(ts & sg)
    b = len(ts - sg)
    c = len(not_ts & sg)
    d = len(not_ts - sg)
    if b == 0 or c == 0:
        return np.nan
    return (a * d) / (b * c)


# ── 3. Joint (baseMean × length) matched sampling ─────────────────────────────
def joint_matched_sample(ir_idx, bg_idx, log_bm, log_len, n_bins_bm, n_bins_len, rng):
    """
    ir_idx, bg_idx: numpy index arrays into shared log_bm, log_len arrays.
    Returns numpy array of background indices forming a length+expression-matched null set.
    """
    # Bin edges from IR distribution (percentile-based, in each axis)
    bm_edges  = np.percentile(log_bm[ir_idx],  np.linspace(0, 100, n_bins_bm  + 1))
    len_edges = np.percentile(log_len[ir_idx], np.linspace(0, 100, n_bins_len + 1))
    bm_edges  = np.unique(bm_edges)
    len_edges = np.unique(len_edges)
    nb = len(bm_edges) - 1
    nl = len(len_edges) - 1

    def assign(arr_bm, arr_len):
        b = np.clip(np.digitize(arr_bm,  bm_edges[1:]),  0, nb - 1)
        l = np.clip(np.digitize(arr_len, len_edges[1:]), 0, nl - 1)
        return b, l

    ir_b, ir_l = assign(log_bm[ir_idx],  log_len[ir_idx])
    bg_b, bg_l = assign(log_bm[bg_idx],  log_len[bg_idx])

    # For each (b, l) cell, count how many IR genes fall in it
    sampled = []
    for b in range(nb):
        for l in range(nl):
            n_needed = int(((ir_b == b) & (ir_l == l)).sum())
            if n_needed == 0:
                continue
            cand_local = np.where((bg_b == b) & (bg_l == l))[0]
            # neighbor fallback: expand 1 cell, then 2, then global
            r = 1
            while len(cand_local) < n_needed and r <= max(nb, nl):
                cand_local = np.where(
                    (np.abs(bg_b.astype(int) - b) <= r) &
                    (np.abs(bg_l.astype(int) - l) <= r)
                )[0]
                r += 1
            if len(cand_local) < n_needed:
                cand_local = np.arange(len(bg_idx))
            n_draw = min(n_needed, len(cand_local))
            chosen = rng.choice(cand_local, size=n_draw, replace=False)
            sampled.extend(bg_idx[chosen].tolist())
    return np.array(sampled, dtype=int)


# ─────────────────────────────────────────────────────────────────────────────
# Load IR genes
print("Loading IR genes …")
ir_df    = pd.read_csv(IR_FILE, sep="\t")
ir_genes = set(ir_df["gene_name"].dropna().unique())
print(f"  Total IR genes: {len(ir_genes)}")

# Load DESeq2 universe
print("Loading DESeq2 results …")
deseq_df = pd.read_csv(DESEQ_FILE, sep="\t").dropna(subset=["baseMean", "gene"])
deseq_df = deseq_df[deseq_df["baseMean"] > 0].copy()
universe   = set(deseq_df["gene"].unique())
ir_in_univ = ir_genes & universe
print(f"  Universe: {len(universe)}; IR in universe: {len(ir_in_univ)}")

# Load gene length table
print("Building/loading gene length table from GTF …")
glen_df = build_gene_length_table(GTF_FILE, GENE_LEN_TSV)

# Merge length onto DESeq2 universe
deseq_df = deseq_df.merge(glen_df, on="gene", how="inner")
deseq_df = deseq_df[deseq_df["length"] > 0].copy()
universe = set(deseq_df["gene"].unique())
ir_in_univ = ir_genes & universe
print(f"  After length annotation: universe={len(universe)}, IR in universe={len(ir_in_univ)}")

# Build aligned arrays
gene_arr  = deseq_df["gene"].values
log_bm    = np.log10(deseq_df["baseMean"].values + 1.0)
log_len   = np.log10(deseq_df["length"].values   + 1.0)
gene2idx  = {g: i for i, g in enumerate(gene_arr)}
ir_idx    = np.array([gene2idx[g] for g in ir_in_univ if g in gene2idx], dtype=int)
bg_idx    = np.arange(len(gene_arr))

# ── Load SG sets ──────────────────────────────────────────────────────────────
print("Loading Khong 2017 SG-enriched genes …")
khong_raw = pd.read_csv(KHONG_FILE, sep="\t")
khong_ok  = khong_raw[khong_raw["status"] == "OK"].copy()
khong_ok["mouse_gene"] = khong_ok["gene"].apply(human_to_mouse)
khong_ok = khong_ok.dropna(subset=["mouse_gene"]).drop_duplicates(subset="mouse_gene", keep="first")
khong_enriched_all = khong_ok[khong_ok["classification"] == "enriched"]["mouse_gene"].tolist()
khong_sg_set = set(khong_enriched_all) & universe
print(f"  Khong SG-enriched in universe: {len(khong_sg_set)}")

print("Loading Biancon 2025 meta-analysis …")
biancon_df = pd.read_excel(BIANCON_FILE, engine="openpyxl")
gene_info = biancon_df.groupby("gene_symbol").agg(max_consensus=("consensus", "max")).reset_index()
gene_info["mouse_gene"] = gene_info["gene_symbol"].apply(human_to_mouse)
gene_info = gene_info.dropna(subset=["mouse_gene"])
high_sg_all = gene_info[gene_info["max_consensus"] >= 8]["mouse_gene"].tolist()
biancon_high_set = set(high_sg_all) & universe
print(f"  Biancon high consensus (>=8) in universe: {len(biancon_high_set)}")

# ── Run permutations ──────────────────────────────────────────────────────────
sg_targets = [
    ("SG RNA Enriched (Khong 2017)",        khong_sg_set),
    ("SG RNA High Consensus (Biancon 2025)", biancon_high_set),
]

rng = np.random.default_rng(SEED)
results_summary = []
null_or_dict    = {}
qc_first_null   = None

for sg_name, sg_set in sg_targets:
    print(f"\n{'='*70}")
    print(f"SG set: {sg_name}  ({len(sg_set)} genes)")

    obs_or = or_from_sets(ir_in_univ, sg_set, universe)
    a = len(ir_in_univ & sg_set)
    b = len(ir_in_univ - sg_set)
    c = len((universe - ir_in_univ) & sg_set)
    d = len(universe - ir_in_univ - sg_set)
    _, fisher_p = fisher_exact([[a, b], [c, d]], alternative="greater")
    print(f"  Contingency a={a} b={b} c={c} d={d}")
    print(f"  Observed OR: {obs_or:.4f}   Fisher exact p (one-sided): {fisher_p:.2e}")

    print(f"  Running {N_PERMS:,} length+expression-matched permutations …", flush=True)
    null_ors = []
    for i in range(N_PERMS):
        null_local = joint_matched_sample(
            ir_idx, bg_idx, log_bm, log_len,
            N_BINS_BM, N_BINS_LEN, rng
        )
        null_genes = set(gene_arr[null_local].tolist())
        nor = or_from_sets(null_genes, sg_set, universe)
        if not np.isnan(nor):
            null_ors.append(nor)
        if (i + 1) % 1000 == 0:
            print(f"    {i+1}/{N_PERMS} done …", flush=True)
        if i == 0 and qc_first_null is None:
            qc_first_null = null_local.copy()

    null_ors  = np.array(null_ors)
    mean_null = float(np.mean(null_ors))
    p95_null  = float(np.percentile(null_ors, 95))
    n_exceed  = int(np.sum(null_ors >= obs_or))
    emp_p     = float(n_exceed / len(null_ors))

    print(f"\n  *** RESULTS (length + expression matched) ***")
    print(f"    Observed OR:             {obs_or:.4f}")
    print(f"    Mean null OR:            {mean_null:.4f}")
    print(f"    95th percentile null OR: {p95_null:.4f}")
    print(f"    Empirical p-value:       {emp_p:.5f}  ({n_exceed}/{len(null_ors)})")
    print(f"    Fisher exact p (ref):    {fisher_p:.2e}")

    results_summary.append({
        "sg_set":            sg_name,
        "matching":          "length+expression (10x10 bins)",
        "observed_or":       obs_or,
        "mean_null_or":      mean_null,
        "p95_null_or":       p95_null,
        "n_permutations":    len(null_ors),
        "n_exceed":          n_exceed,
        "empirical_pvalue":  emp_p,
        "fisher_p":          fisher_p,
        "overlap":           a,
        "ir_in_univ":        len(ir_in_univ),
        "sg_in_univ":        len(sg_set),
    })
    null_or_dict[sg_name] = (null_ors, obs_or)

# ── Save TSV ──────────────────────────────────────────────────────────────────
res_df = pd.DataFrame(results_summary)
os.makedirs(OUTDIR, exist_ok=True)
res_df.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.5f")
print(f"\nResults TSV → {OUT_TSV}")
print(res_df[["sg_set", "observed_or", "mean_null_or", "p95_null_or", "empirical_pvalue", "fisher_p"]].to_string(index=False))

# ── Plot main: null OR distributions ──────────────────────────────────────────
n_plots = len(null_or_dict)
fig, axes = plt.subplots(1, n_plots, figsize=(6.5 * n_plots, 5.5))
if n_plots == 1:
    axes = [axes]

colors = ["#1565C0", "#B71C1C"]

for ax, (sg_name, (null_ors, obs_or)), col in zip(axes, null_or_dict.items(), colors):
    emp_p = float(np.mean(null_ors >= obs_or))
    p95   = float(np.percentile(null_ors, 95))
    mean  = float(np.mean(null_ors))

    ax.hist(null_ors, bins=60, color=col, alpha=0.65, edgecolor="white", linewidth=0.2)
    ax.axvline(obs_or, color="black",   linewidth=2.0, linestyle="-",
               label=f"Observed OR = {obs_or:.2f}")
    ax.axvline(mean,   color="#555555", linewidth=1.2, linestyle=":",
               label=f"Mean null = {mean:.2f}")
    ax.axvline(p95,    color="#FF8F00", linewidth=1.5, linestyle="--",
               label=f"95th pctile = {p95:.2f}")

    ax.set_xlabel("Odds Ratio (null permutations)", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    short_name = sg_name.replace("SG RNA ", "")
    p_text = f"p = {emp_p:.4f}" if emp_p > 0 else f"p < {1/len(null_ors):.5f}"
    ax.set_title(f"{short_name}\nEmpirical {p_text}", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.85, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.suptitle(
    "SG Enrichment Permutation Test — joint baseMean × transcript length matching\n"
    f"(10×10 percentile bins, n={N_PERMS:,} permutations)",
    fontsize=12, fontweight="bold", y=1.01
)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"\nMain figure → {OUT_PNG}")

# ── QC: matching quality ─────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ir_lbm  = log_bm[ir_idx]
ir_llen = log_len[ir_idx]
nl_lbm  = log_bm[qc_first_null]
nl_llen = log_len[qc_first_null]

ks_bm,  p_bm  = ks_2samp(ir_lbm,  nl_lbm)
ks_len, p_len = ks_2samp(ir_llen, nl_llen)

ax = axes[0]
ax.hist(ir_lbm,  bins=40, alpha=0.55, color="#B71C1C", label=f"IR genes (n={len(ir_lbm)})", density=True)
ax.hist(nl_lbm,  bins=40, alpha=0.55, color="#1565C0", label=f"Matched null #1 (n={len(nl_lbm)})", density=True)
ax.set_xlabel("log10(baseMean + 1)")
ax.set_ylabel("Density")
ax.set_title(f"Expression matching\nKS D={ks_bm:.3f}, p={p_bm:.3f}")
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

ax = axes[1]
ax.hist(ir_llen, bins=40, alpha=0.55, color="#B71C1C", label="IR genes",          density=True)
ax.hist(nl_llen, bins=40, alpha=0.55, color="#1565C0", label="Matched null #1", density=True)
ax.set_xlabel("log10(transcript length + 1)")
ax.set_ylabel("Density")
ax.set_title(f"Length matching\nKS D={ks_len:.3f}, p={p_len:.3f}")
ax.legend(fontsize=9)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(QC_PNG, dpi=150, bbox_inches="tight")
print(f"QC figure  → {QC_PNG}")
print(f"Matching diagnostics: KS(baseMean) D={ks_bm:.4f} p={p_bm:.4f}; KS(length) D={ks_len:.4f} p={p_len:.4f}")
print("Done.")
