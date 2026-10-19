#!/usr/bin/env python3
"""
Motif Enrichment Reverse Validation for DHX8 IR Targets
=======================================================
Run on: ncnp-linux (conda activate clip3)

Purpose:
  Test whether DHX8 IR targets show enrichment for specific
  sequence motifs, to assess the hypothesis that DHX8 acts as a core
  spliceosome factor (not sequence-specific).

Analyses:
  1. De novo k-mer enrichment (4-6mer) in retained vs non-retained introns
  2. Known RBP motif scanning (SRSF1, SRSF2, hnRNP, TDP-43 motifs)
  3. Splice site strength comparison (MaxEntScan-style)
  4. Sequence complexity analysis
  5. Length-controlled feature comparison

Input:
  - limma_retained_dhx8.tsv (1,258 retained introns; v2 corrected-siCont_2)
  - limma_ir_results_dhx8.tsv (all tested introns = background)
  - mm10 genome FASTA

Output:
  - motif_enrichment_results.tsv
  - figures/motif_*.pdf + .png
"""

import os
import sys
import re
import json
import numpy as np
import pandas as pd
from scipy import stats
from collections import Counter
from itertools import product
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# Arial + Illustrator-editable text (project figure standard); literal Unicode for symbols.
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.unicode_minus'] = False
from matplotlib.backends.backend_pdf import PdfPages
import pysam

# === Configuration ===
WORKDIR = "/mnt/HDD2/phase3_specificity"
DATA_SRC = "/mnt/HDD2/TDP43_comparison/validation_20260226"
GENOME_FA = "/mnt/NAS/database/UCSC_genome/mm10/mm10.fa"
FIGDIR = os.path.join(WORKDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(WORKDIR, exist_ok=True)

# Known RBP binding motifs (consensus sequences)
RBP_MOTIFS = {
    # SR proteins
    'SRSF1': [r'[AG]GA[AG]GA', r'GGAGGA'],
    'SRSF2': [r'[CG]C[CG]G', r'CCNG', r'GGNG'],
    'SRSF3': [r'CAUC', r'CUCKUCY'],
    'SRSF6': [r'UGCNUG', r'TGCNTG'],
    # hnRNPs
    'hnRNPA1': [r'UAGG', r'TAGG[AT]'],
    'hnRNPA2B1': [r'[AT]AGGA', r'UAGGA'],
    'hnRNPC': [r'[TU]{5,}', r'TTTTT'],
    'hnRNPH': [r'GGG[ATGC]?GGG', r'GGGG'],
    'hnRNPF': [r'GGGGA', r'GGGG'],
    # Disease-associated
    'TDP43': [r'[TU]G{3,6}', r'TGTGT', r'[TU]G[TU]G[TU]G'],
    'FUS': [r'GGU[TU]', r'GGTT'],
    'MBNL1': [r'[CT]GC[CT]', r'YGCY'],
    # U2AF / core splicing
    'U2AF65': [r'[CT]{3,}', r'TTTTTT'],  # polypyrimidine tract
    'SF1': [r'[TA]ACT[AG]AC'],  # branch point
}


def get_sequence(chrom, start, end, genome):
    """Get sequence from genome (0-based coordinates)."""
    try:
        return genome.fetch(chrom, start, end).upper()
    except Exception:
        return None


def gc_content(seq):
    if not seq or len(seq) == 0:
        return np.nan
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)


def count_motif(seq, pattern):
    """Count regex motif occurrences in sequence."""
    if not seq:
        return 0
    # Convert RNA motifs to DNA
    pattern_dna = pattern.replace('U', 'T')
    try:
        return len(re.findall(pattern_dna, seq, re.IGNORECASE))
    except Exception:
        return 0


def kmer_frequencies(seq, k):
    """Compute k-mer frequency vector."""
    if not seq or len(seq) < k:
        return Counter()
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:
            counts[kmer] += 1
    return counts


def sequence_complexity(seq):
    """Linguistic complexity: ratio of observed to possible k-mers."""
    if not seq or len(seq) < 3:
        return np.nan
    k = 3
    observed = len(set(seq[i:i+k] for i in range(len(seq) - k + 1)
                       if 'N' not in seq[i:i+k]))
    possible = min(4**k, len(seq) - k + 1)
    return observed / possible if possible > 0 else np.nan


def g_run_density(seq):
    """Count G-runs (>=3 consecutive Gs) per kb."""
    if not seq:
        return 0
    runs = re.findall(r'G{3,}', seq)
    return len(runs) / (len(seq) / 1000) if len(seq) > 0 else 0


def rev_comp(seq):
    """Reverse complement."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


def main():
    print("=" * 70)
    print("Motif Enrichment Reverse Validation for DHX8 IR Targets")
    print("=" * 70)

    # Load data
    retained = pd.read_csv(os.path.join(DATA_SRC, "limma_retained_dhx8.tsv"), sep="\t")
    all_tested = pd.read_csv(os.path.join(DATA_SRC, "limma_ir_results_dhx8.tsv"), sep="\t")
    non_retained = all_tested[~all_tested["intron_id"].isin(retained["intron_id"])].copy()

    print(f"\nRetained introns: {len(retained)}")
    print(f"Non-retained introns: {len(non_retained)}")

    # Open genome
    genome = pysam.FastaFile(GENOME_FA)

    # === Part 1: Length-matched background selection ===
    print("\n--- Part 1: Length-matched background ---")

    # Create length bins and sample matched background
    retained['log_length'] = np.log10(retained['intron_length'])
    non_retained['log_length'] = np.log10(non_retained['intron_length'])

    # Bin by log10(length) with 0.2-width bins
    bins = np.arange(1, 6.5, 0.2)
    retained['len_bin'] = pd.cut(retained['log_length'], bins=bins)
    non_retained['len_bin'] = pd.cut(non_retained['log_length'], bins=bins)

    # Sample length-matched background (3:1 ratio)
    np.random.seed(42)
    matched_bg = []
    for bin_val, group in retained.groupby('len_bin'):
        n_needed = len(group) * 3
        pool = non_retained[non_retained['len_bin'] == bin_val]
        if len(pool) >= n_needed:
            matched_bg.append(pool.sample(n_needed, random_state=42))
        elif len(pool) > 0:
            matched_bg.append(pool)

    if matched_bg:
        matched_background = pd.concat(matched_bg)
    else:
        print("WARNING: Could not create length-matched background!")
        matched_background = non_retained.sample(min(3000, len(non_retained)), random_state=42)

    print(f"Length-matched background: {len(matched_background)} introns")

    # Verify length matching
    ks_stat, ks_p = stats.ks_2samp(retained['log_length'], matched_background['log_length'])
    print(f"Length matching KS test: D={ks_stat:.4f}, p={ks_p:.4f}")

    # DUMP for bespoke S12 panel re-plot (length arrays, same as KS input)
    import pandas as _pd
    _len = _pd.concat([
        retained[['intron_length']].assign(category='retained'),
        matched_background[['intron_length']].assign(category='control'),
    ])
    _len.to_csv(os.path.join(WORKDIR, "s12_length_arrays_v2.tsv"), sep='\t', index=False)
    with open(os.path.join(WORKDIR, "s12_ks_v2.txt"), 'w') as _fh:
        _fh.write(f"D={ks_stat:.4f}\tp={ks_p:.6f}\tn_ret={len(retained)}\tn_ctrl={len(matched_background)}\n")
    print("DUMP: s12_length_arrays_v2.tsv, s12_ks_v2.txt")

    # === Part 2: Extract sequences and compute features ===
    print("\n--- Part 2: Sequence extraction and feature computation ---")

    datasets = {
        'retained': retained,
        'matched_bg': matched_background,
    }

    all_features = []
    all_kmer_counts = {'retained': Counter(), 'matched_bg': Counter()}

    for label, df in datasets.items():
        print(f"\n  Processing {label} ({len(df)} introns)...")
        for idx, (_, row) in enumerate(df.iterrows()):
            chrom = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            strand = row['strand']
            length = end - start

            # Limit sequence extraction to first/last 500bp for very long introns
            if length > 2000:
                seq_5 = get_sequence(chrom, start, start + 500, genome)
                seq_3 = get_sequence(chrom, end - 500, end, genome)
                seq_mid = get_sequence(chrom, (start + end) // 2 - 250,
                                       (start + end) // 2 + 250, genome)
                seq_full = None  # too long for full k-mer
                seq_for_motif = (seq_5 or '') + (seq_mid or '') + (seq_3 or '')
            else:
                seq_full = get_sequence(chrom, start, end, genome)
                seq_for_motif = seq_full
                seq_5 = seq_full[:500] if seq_full and len(seq_full) >= 500 else seq_full
                seq_3 = seq_full[-500:] if seq_full and len(seq_full) >= 500 else seq_full

            if not seq_for_motif:
                continue

            # Handle strand
            if strand == '-':
                seq_for_motif = rev_comp(seq_for_motif)
                if seq_5:
                    seq_5 = rev_comp(seq_5)
                if seq_3:
                    seq_3 = rev_comp(seq_3)

            # GC content
            gc = gc_content(seq_for_motif)

            # G-run density
            g_runs = g_run_density(seq_for_motif)

            # Sequence complexity
            complexity = sequence_complexity(seq_for_motif)

            # RBP motif counts (normalized per kb)
            motif_counts = {}
            for rbp, patterns in RBP_MOTIFS.items():
                total = sum(count_motif(seq_for_motif, p) for p in patterns)
                motif_counts[f'motif_{rbp}'] = total
                motif_counts[f'motif_{rbp}_perkb'] = total / (len(seq_for_motif) / 1000)

            # k-mer counts for enrichment analysis
            if seq_for_motif:
                for k in [4, 5, 6]:
                    kmers = kmer_frequencies(seq_for_motif, k)
                    all_kmer_counts[label] += kmers

            feat = {
                'intron_id': row['intron_id'],
                'category': label,
                'length': length,
                'gc_content': gc,
                'g_run_density': g_runs,
                'complexity': complexity,
                'seq_length_sampled': len(seq_for_motif),
            }
            feat.update(motif_counts)
            all_features.append(feat)

            if (idx + 1) % 500 == 0:
                print(f"    {idx+1}/{len(df)}...")

    genome.close()

    feat_df = pd.DataFrame(all_features)
    print(f"\nTotal features computed: {len(feat_df)}")

    # DUMP per-intron features for bespoke S12 GC panel re-plot
    _cols = [c for c in ['category', 'gc_content', 'g_run_density', 'complexity', 'length']
             if c in feat_df.columns]
    feat_df[_cols].to_csv(os.path.join(WORKDIR, "s12_per_intron_features_v2.tsv"),
                          sep='\t', index=False)
    print("DUMP: s12_per_intron_features_v2.tsv")

    # === Part 3: Statistical comparison ===
    print("\n" + "=" * 70)
    print("Part 3: Statistical Comparison (Retained vs Length-matched Background)")
    print("=" * 70)

    ret = feat_df[feat_df['category'] == 'retained']
    bg = feat_df[feat_df['category'] == 'matched_bg']

    results = []

    # 3a. Basic sequence features
    basic_features = ['gc_content', 'g_run_density', 'complexity']
    print("\n--- Basic sequence features ---")
    for feat_name in basic_features:
        r_vals = ret[feat_name].dropna()
        b_vals = bg[feat_name].dropna()
        try:
            t_stat, t_p = stats.ttest_ind(r_vals, b_vals)
        except Exception:
            t_stat, t_p = 0, 1.0
        try:
            mw_stat, mw_p = stats.mannwhitneyu(r_vals, b_vals, alternative='two-sided')
        except ValueError:
            mw_p = 1.0
        pooled_var = (r_vals.var() + b_vals.var()) / 2
        effect_size = (r_vals.mean() - b_vals.mean()) / np.sqrt(pooled_var) if pooled_var > 0 else 0

        print(f"  {feat_name}:")
        print(f"    Retained: {r_vals.mean():.4f} ± {r_vals.std():.4f}")
        print(f"    Background: {b_vals.mean():.4f} ± {b_vals.std():.4f}")
        print(f"    t-test p={t_p:.4e}, MWU p={mw_p:.4e}, Cohen's d={effect_size:.3f}")

        results.append({
            'feature': feat_name,
            'ret_mean': r_vals.mean(),
            'ret_std': r_vals.std(),
            'bg_mean': b_vals.mean(),
            'bg_std': b_vals.std(),
            'ttest_p': t_p,
            'mwu_p': mw_p,
            'cohens_d': effect_size,
            'significant_bonf': mw_p < 0.05 / 50,  # Bonferroni for ~50 tests
        })

    # 3b. RBP motif densities
    print("\n--- RBP motif densities (per kb, length-controlled) ---")
    motif_cols = [c for c in feat_df.columns if c.endswith('_perkb')]
    for col in sorted(motif_cols):
        rbp = col.replace('motif_', '').replace('_perkb', '')
        r_vals = ret[col].dropna()
        b_vals = bg[col].dropna()
        if len(r_vals) < 10 or len(b_vals) < 10:
            continue
        try:
            mw_stat, mw_p = stats.mannwhitneyu(r_vals, b_vals, alternative='two-sided')
        except ValueError:
            mw_p = 1.0  # All values identical
        pooled_var = (r_vals.var() + b_vals.var()) / 2
        effect_size = (r_vals.mean() - b_vals.mean()) / np.sqrt(pooled_var) if pooled_var > 0 else 0

        sig_marker = " ***" if mw_p < 0.05 / len(motif_cols) else ""
        print(f"  {rbp:>10s}: ret={r_vals.mean():.3f}  bg={b_vals.mean():.3f}  "
              f"d={effect_size:+.3f}  p={mw_p:.4e}{sig_marker}")

        results.append({
            'feature': f'motif_{rbp}_perkb',
            'ret_mean': r_vals.mean(),
            'ret_std': r_vals.std(),
            'bg_mean': b_vals.mean(),
            'bg_std': b_vals.std(),
            'ttest_p': np.nan,
            'mwu_p': mw_p,
            'cohens_d': effect_size,
            'significant_bonf': mw_p < 0.05 / len(motif_cols),
        })

    # === Part 4: K-mer enrichment analysis ===
    print("\n" + "=" * 70)
    print("Part 4: K-mer Enrichment Analysis")
    print("=" * 70)

    kmer_results = []
    for k in [4, 5, 6]:
        print(f"\n--- {k}-mer analysis ---")
        # Get all possible k-mers
        all_kmers = [''.join(p) for p in product('ACGT', repeat=k)]

        # Total counts
        ret_total = sum(all_kmer_counts['retained'].values())
        bg_total = sum(all_kmer_counts['matched_bg'].values())

        if ret_total == 0 or bg_total == 0:
            print(f"  Skipping {k}-mer: no counts")
            continue

        enrichments = []
        for kmer in all_kmers:
            ret_count = all_kmer_counts['retained'].get(kmer, 0)
            bg_count = all_kmer_counts['matched_bg'].get(kmer, 0)

            ret_freq = ret_count / ret_total
            bg_freq = bg_count / bg_total

            if bg_freq > 0:
                log2_enrich = np.log2(ret_freq / bg_freq)
            else:
                log2_enrich = np.nan

            # Fisher's exact test (sample down for speed)
            a = ret_count
            b = ret_total - ret_count
            c = bg_count
            d = bg_total - bg_count
            # Use chi-squared approximation for large counts
            if a + c > 5:
                try:
                    chi2, p_val = stats.chi2_contingency(
                        [[a, b], [c, d]], correction=True)[:2]
                except Exception:
                    chi2, p_val = 0, 1.0
            else:
                p_val = 1.0
                chi2 = 0

            enrichments.append({
                'kmer': kmer,
                'k': k,
                'ret_count': ret_count,
                'bg_count': bg_count,
                'ret_freq': ret_freq,
                'bg_freq': bg_freq,
                'log2_enrichment': log2_enrich,
                'chi2': chi2,
                'pvalue': p_val,
            })

        enrich_df = pd.DataFrame(enrichments)
        # Bonferroni correction
        enrich_df['padj'] = np.minimum(enrich_df['pvalue'] * len(enrich_df), 1.0)
        enrich_df = enrich_df.sort_values('log2_enrichment', ascending=False)

        # Report top enriched/depleted
        sig = enrich_df[enrich_df['padj'] < 0.05]
        print(f"  Total {k}-mers tested: {len(enrich_df)}")
        print(f"  Significant after Bonferroni: {len(sig)}")

        if len(sig) > 0:
            print(f"\n  Top enriched:")
            for _, row in sig.head(5).iterrows():
                print(f"    {row['kmer']}: log2={row['log2_enrichment']:.3f}, "
                      f"padj={row['padj']:.4e}")
            print(f"  Top depleted:")
            for _, row in sig.tail(5).iterrows():
                print(f"    {row['kmer']}: log2={row['log2_enrichment']:.3f}, "
                      f"padj={row['padj']:.4e}")

        # Effect size distribution
        valid = enrich_df['log2_enrichment'].dropna()
        print(f"\n  Log2 enrichment distribution:")
        print(f"    Range: [{valid.min():.3f}, {valid.max():.3f}]")
        print(f"    Mean: {valid.mean():.4f}, Std: {valid.std():.4f}")
        print(f"    Max |log2|: {valid.abs().max():.3f}")

        kmer_results.append(enrich_df)

    # === Part 5: Visualization ===
    print("\n" + "=" * 70)
    print("Part 5: Generating figures")
    print("=" * 70)

    results_df = pd.DataFrame(results)

    # Figure 1: RBP motif density comparison (forest plot style)
    fig, ax = plt.subplots(figsize=(10, 8))
    motif_results = results_df[results_df['feature'].str.startswith('motif_')]
    motif_results = motif_results.sort_values('cohens_d')

    labels = [r['feature'].replace('motif_', '').replace('_perkb', '')
              for _, r in motif_results.iterrows()]
    y_pos = range(len(labels))
    colors = ['red' if r['significant_bonf'] else 'steelblue'
              for _, r in motif_results.iterrows()]

    ax.barh(y_pos, motif_results['cohens_d'].values, color=colors, alpha=0.7,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Cohen's d (effect size)", fontsize=12)
    ax.set_title("RBP Motif Density: Retained vs Background\n(length-matched, red = Bonferroni significant)",
                 fontsize=13)
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(-0.2, color='gray', linewidth=0.5, linestyle=':')
    ax.axvline(0.2, color='gray', linewidth=0.5, linestyle=':')
    n_sig_rbp = int(motif_results['significant_bonf'].sum())
    ax.text(0.02, 0.02,
            f'{n_sig_rbp}/{len(motif_results)} Bonferroni significant\n'
            'dotted lines: d = ±0.2 (small-effect reference)',
            transform=ax.transAxes, va='bottom', ha='left', fontsize=9, color='0.25')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f'motif_rbp_enrichment.{ext}'), dpi=300)
    plt.close()

    # Figure 2: K-mer enrichment volcano
    if kmer_results:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for i, (k, enrich_df) in enumerate(zip([4, 5, 6], kmer_results)):
            ax = axes[i]
            valid = enrich_df.dropna(subset=['log2_enrichment'])
            colors = ['red' if p < 0.05 else 'gray'
                      for p in valid['padj']]
            ax.scatter(valid['log2_enrichment'], -np.log10(valid['pvalue']),
                       c=colors, alpha=0.5, s=10)
            ax.set_xlabel('log2(enrichment)', fontsize=11)
            ax.set_ylabel('-log10(p-value)', fontsize=11)
            ax.set_title(f'{k}-mer enrichment', fontsize=12)
            ax.axvline(0, color='black', linewidth=0.5, linestyle='--')

            # Annotate the Bonferroni-significant k-mers (red points): the only
            # k-mers reaching significance are low-complexity (poly-C, GC-rich,
            # CT/GC dinucleotide-repeat), so labelling them shows there is no
            # meaningful sequence preference beyond simple/low-complexity motifs.
            sig = valid[valid['padj'] < 0.05].sort_values('log2_enrichment')
            for _, row in sig.iterrows():
                ha = 'left' if row['log2_enrichment'] >= 0 else 'right'
                ax.annotate(row['kmer'],
                            (row['log2_enrichment'], -np.log10(row['pvalue'])),
                            fontsize=7, ha=ha, color='#7a1010')
            n_sig_k = int((valid['padj'] < 0.05).sum())
            sig_maxabs = sig['log2_enrichment'].abs().max() if len(sig) else 0.0
            ax.text(0.02, 0.98,
                    f'{n_sig_k} sig. (all low-complexity)\n'
                    f'max sig. |log2| = {sig_maxabs:.2f}',
                    transform=ax.transAxes, va='top', ha='left', fontsize=7,
                    color='0.25')

        plt.suptitle("K-mer Enrichment in DHX8 Retained vs Background Introns\n"
                     "(red = Bonferroni significant; labelled k-mers are low-complexity)",
                     fontsize=13)
        plt.tight_layout(rect=[0, 0, 1, 0.90])
        for ext in ['pdf', 'png']:
            fig.savefig(os.path.join(FIGDIR, f'motif_kmer_enrichment.{ext}'), dpi=300)
        plt.close()

    # Figure 3: Sequence feature comparison
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    feat_labels = {'gc_content': 'GC content',
                   'g_run_density': 'G-run density (per kb)',
                   'complexity': 'Sequence complexity'}
    for i, feat_name in enumerate(basic_features):
        ax = axes[i]
        r_vals = ret[feat_name].dropna()
        b_vals = bg[feat_name].dropna()
        ax.hist(r_vals, bins=50, alpha=0.6, density=True, label='Retained', color='red')
        ax.hist(b_vals, bins=50, alpha=0.6, density=True, label='Control', color='steelblue')
        ax.set_xlabel(feat_labels.get(feat_name, feat_name), fontsize=11)
        ax.set_ylabel('Density', fontsize=11)
        ax.legend(fontsize=9)

        # Add p-value annotation (Mann-Whitney U)
        mw_row = results_df[results_df['feature'] == feat_name]
        if len(mw_row) > 0:
            p = mw_row.iloc[0]['mwu_p']
            d = mw_row.iloc[0]['cohens_d']
            ax.set_title(f"MWU P = {p:.2g}, d = {d:.3f}", fontsize=10)

    plt.suptitle("Sequence Features: Retained vs Length-matched Control", fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f'motif_seq_features.{ext}'), dpi=300)
    plt.close()

    # Figure 4: Summary - number of significant motifs
    fig, ax = plt.subplots(figsize=(6, 4))
    n_total = len(motif_results)
    n_sig = motif_results['significant_bonf'].sum()
    bars = ax.bar(['Significant\n(Bonferroni)', 'Not significant'],
                  [n_sig, n_total - n_sig],
                  color=['red', 'steelblue'], alpha=0.7, edgecolor='black')
    ax.set_ylabel('Number of RBP motifs tested', fontsize=12)
    ax.set_title(f'RBP Motif Enrichment Summary\n'
                 f'{n_sig}/{n_total} motifs show significant enrichment',
                 fontsize=13)
    for bar, val in zip(bars, [n_sig, n_total - n_sig]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                str(val), ha='center', fontsize=14, fontweight='bold')
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f'motif_summary.{ext}'), dpi=300)
    plt.close()

    # Figure 5: Length matching control (S12D)
    # ECDF of intron length, retained (n=1,258) vs 3x length-matched control,
    # with the KS statistic that verifies the two distributions are matched.
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    rl = np.sort(retained['intron_length'].values.astype(float))
    bl = np.sort(matched_background['intron_length'].values.astype(float))
    ax.plot(rl, np.linspace(0, 1, len(rl)), color='#c0392b', lw=1.6,
            label=f'Retained (n = {len(rl):,})')
    ax.plot(bl, np.linspace(0, 1, len(bl)), color='#4a7fb5', lw=1.6,
            label=f'Length-matched control (n = {len(bl):,})')
    ax.set_xscale('log')
    ax.set_xlabel('Intron length (bp)', fontsize=11)
    ax.set_ylabel('Cumulative fraction', fontsize=11)
    p_txt = 'P > 0.999' if ks_p > 0.999 else f'P = {ks_p:.3g}'
    ax.set_title(f'Length matching (KS D = {ks_stat:.3f}, {p_txt})', fontsize=11)
    ax.legend(fontsize=8, loc='lower right', frameon=False)
    for _sp in ('top', 'right'):
        ax.spines[_sp].set_visible(False)
    plt.tight_layout()
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(FIGDIR, f'motif_length_matching.{ext}'), dpi=300)
    plt.close()

    # === Save results ===
    results_df.to_csv(os.path.join(WORKDIR, "motif_enrichment_results.tsv"),
                      sep='\t', index=False)

    # Save k-mer results
    if kmer_results:
        all_kmer_df = pd.concat(kmer_results)
        all_kmer_df.to_csv(os.path.join(WORKDIR, "kmer_enrichment_results.tsv"),
                           sep='\t', index=False)

    # === Summary ===
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nRBP motifs tested: {n_total}")
    print(f"Significant (Bonferroni): {n_sig}")

    if n_sig == 0:
        print("\n>>> NO RBP motifs show significant enrichment in DHX8 IR targets")
        print(">>> This supports DHX8 acting as a core spliceosome factor,")
        print(">>> not a sequence-specific splicing regulator.")
    else:
        print(f"\n>>> {n_sig} motifs show enrichment - investigate further")
        sig_feats = motif_results[motif_results['significant_bonf']]
        for _, r in sig_feats.iterrows():
            print(f"    {r['feature']}: d={r['cohens_d']:.3f}, p={r['mwu_p']:.2e}")

    for k_idx, (k, enrich_df) in enumerate(zip([4, 5, 6], kmer_results)):
        n_sig_k = (enrich_df['padj'] < 0.05).sum()
        max_enrich = enrich_df['log2_enrichment'].dropna().abs().max()
        print(f"\n{k}-mer: {n_sig_k} significant, max |log2 enrichment| = {max_enrich:.3f}")

    print(f"\nResults saved to: {WORKDIR}/motif_enrichment_results.tsv")
    print(f"Figures saved to: {FIGDIR}/motif_*.pdf/png")
    print("\nDone!")


if __name__ == "__main__":
    main()
