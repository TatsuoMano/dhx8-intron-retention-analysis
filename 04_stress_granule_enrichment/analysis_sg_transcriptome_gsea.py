#!/usr/bin/env python3
"""
SG Transcriptome GSEA — Do IR mRNAs go to stress granules or stay in the nucleus?

Analyses:
  1. Multi-source pre-ranked GSEA (DEG rank list vs SG RNA gene sets)
  2. IR gene × SG transcriptome overlap (Fisher's exact)
  3. SG enrichment ratio distribution (IR vs all, Khong 2017)
  4. Transcript length confound check
  5. RNA vs Protein level comparison
  6. Per-dataset overlap (26 Biancon meta-analysis datasets)
  7. DEG × SG overlap (Fisher's exact, UP/DOWN separately)
  8. Venn diagrams (IR∩SG, DEG-UP∩SG for each dataset)

Data sources:
  - Khong et al. 2017 Mol Cell (SG core IP + RNA-seq, GEO: GSE99304)
  - Namkoong et al. 2018 Mol Cell (Insoluble RG fraction, 3 stress conditions)
  - Moon et al. 2023 Nat Commun (CAP-seq photocatalytic proximity labeling)
  - Biancon et al. 2025 RNA (meta-analysis of 26 SG transcriptome datasets)
  - Existing SG proteome sets (Markmiller 2018, Youn 2018)

Input:
  - data/deseq2_results.csv
  - data/limma_retained_dhx8.tsv (significant IR introns)
  - data/limma_ir_results_dhx8.tsv (all tested introns)
  - data/GSE99304_U2OSSGcuffdiff.txt (Khong 2017)
  - data/namkoong2018_mmc3.xlsx (Namkoong 2018 ER stress)
  - data/namkoong2018_mmc4.xlsx (Namkoong 2018 heat shock)
  - data/namkoong2018_mmc5.xlsx (Namkoong 2018 arsenite)
  - data/moon2023_suppdata2.xlsx (Moon 2023 CAP-seq)
  - data/meta2025_table3.xlsx (Biancon 2025 meta-analysis)
  - data/meta2025_table1.xlsx (Biancon 2025 dataset metadata)

Output:
  - sg_transcriptome_gsea_results.tsv
  - sg_rna_ir_overlap.tsv
  - sg_rna_enrichment_distribution.tsv
  - sg_per_dataset_overlap.tsv
  - figures/sg_transcriptome_gsea_enrichment.{pdf,png}
  - figures/sg_rna_ir_overlap.{pdf,png}
  - figures/sg_rna_distribution.{pdf,png}
  - figures/sg_rna_length_control.{pdf,png}
  - figures/sg_per_dataset_forest.{pdf,png}
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# === Configuration ===
BASEDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(BASEDIR, "data")
FIGDIR = os.path.join(BASEDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)


# =============================================
# Human → Mouse gene symbol conversion
# =============================================
def human_to_mouse_symbol(human_symbol):
    """Convert human gene symbol to mouse (Title case: GAPDH → Gapdh)."""
    if not isinstance(human_symbol, str):
        return None
    s = human_symbol.strip()
    if len(s) == 0:
        return None
    # Standard conversion: first letter uppercase, rest lowercase
    # Handles most cases (GAPDH→Gapdh, TP53→Tp53, etc.)
    return s[0].upper() + s[1:].lower()


def batch_human_to_mouse(human_genes):
    """Convert a list of human gene symbols to mouse symbols."""
    mouse_genes = []
    for g in human_genes:
        m = human_to_mouse_symbol(g)
        if m:
            mouse_genes.append(m)
    return mouse_genes


# =============================================
# Data loading functions
# =============================================
def load_khong2017(filepath):
    """Load Khong et al. 2017 SG transcriptome data from GEO cuffdiff output."""
    df = pd.read_csv(filepath, sep='\t')
    # Filter to OK status (testable genes)
    df_ok = df[df['status'] == 'OK'].copy()
    # Ensure numeric columns
    for col in ['SG(FPKM)', 'Total(FPKM)', 'log2(fold_change)', 'test_stat', 'p_value', 'q_value']:
        df_ok[col] = pd.to_numeric(df_ok[col], errors='coerce')
    # Convert human → mouse gene symbols
    df_ok['mouse_gene'] = df_ok['gene'].apply(human_to_mouse_symbol)
    df_ok = df_ok.dropna(subset=['mouse_gene'])
    # Remove duplicates (keep first occurrence)
    df_ok = df_ok.drop_duplicates(subset='mouse_gene', keep='first')
    return df_ok


def load_meta2025(filepath):
    """Load Biancon et al. 2025 meta-analysis data."""
    import openpyxl
    wb = openpyxl.load_workbook(filepath)
    ws = wb[wb.sheetnames[0]]
    header = [cell.value for cell in ws[1]]
    data = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        data.append(dict(zip(header, row)))
    df = pd.DataFrame(data)
    return df


def load_namkoong2018(filepath, stress_label):
    """Load Namkoong et al. 2018 stress-induced RNA granulation data.

    Args:
        filepath: path to mmc3/4/5 xlsx (Table S2/S3/S4)
        stress_label: 'ER', 'HS', or 'AS'

    Returns:
        dict with 'enriched' and 'depleted' mouse gene lists, and full dataframe
    """
    import openpyxl
    wb = openpyxl.load_workbook(filepath, read_only=True)
    ws = wb[wb.sheetnames[0]]
    header_row = list(ws.iter_rows(min_row=1, max_row=1, values_only=True))[0]
    data = []
    for row in ws.iter_rows(min_row=2, values_only=True):
        gene = row[0]
        if not gene:
            continue
        del_rg = row[3]
        cluster = row[4]
        if del_rg is not None and cluster is not None:
            data.append({
                'gene': str(gene).strip(),
                'del_RG': float(del_rg),
                'cluster': int(cluster),
            })
    wb.close()

    df = pd.DataFrame(data)
    # Namkoong data is already in mouse gene symbols
    # Enriched = Clusters 1-2 (highest del_RG, positive)
    # Depleted = Clusters 5-6 (lowest del_RG, negative/low)
    enriched = df[df['cluster'].isin([1, 2])]['gene'].tolist()
    depleted = df[df['cluster'].isin([5, 6])]['gene'].tolist()
    all_genes = df['gene'].unique().tolist()
    return {
        'enriched': enriched,
        'depleted': depleted,
        'all_genes': all_genes,
        'df': df,
        'stress': stress_label,
    }


def load_moon2023(filepath, sheet_name):
    """Load Moon et al. 2023 CAP-seq SG proximity data.

    Args:
        filepath: path to suppdata2 xlsx
        sheet_name: 'A-arsenite' or 'B-sorbitol'

    Returns:
        dict with 'enriched' (SG-proximal) and 'depleted' (SG-excluded) gene lists
    """
    import openpyxl
    wb = openpyxl.load_workbook(filepath, read_only=True)
    ws = wb[sheet_name]
    enriched_human = []
    depleted_human = []
    all_genes = []
    for row in ws.iter_rows(min_row=3, values_only=True):
        gene = row[0]
        if not gene:
            continue
        gene = str(gene).strip()
        all_genes.append(gene)
        sg_proximal = row[3]
        sg_excluded = row[4]
        log2fc = row[5]
        if sg_proximal == 1:
            enriched_human.append(gene)
        elif sg_excluded == 1:
            depleted_human.append(gene)
    wb.close()

    enriched_mouse = batch_human_to_mouse(enriched_human)
    depleted_mouse = batch_human_to_mouse(depleted_human)
    return {
        'enriched': enriched_mouse,
        'depleted': depleted_mouse,
        'enriched_human': enriched_human,
        'depleted_human': depleted_human,
        'n_total': len(all_genes),
    }


def extract_meta2025_genesets(df):
    """Extract gene sets from meta-analysis by consensus level."""
    # Get per-gene maximum consensus score
    gene_info = df.groupby('gene_symbol').agg(
        max_consensus=('consensus', 'max'),
        consensus_class=('consensus_class', 'first'),
        n_datasets=('dataset_ID', 'nunique')
    ).reset_index()

    # Convert to mouse symbols
    gene_info['mouse_gene'] = gene_info['gene_symbol'].apply(human_to_mouse_symbol)
    gene_info = gene_info.dropna(subset=['mouse_gene'])

    sets = {}
    # High consensus (detected in many datasets)
    high = gene_info[gene_info['max_consensus'] >= 8]
    sets['high'] = high['mouse_gene'].tolist()
    # Medium-high (consensus >= 4)
    med_high = gene_info[gene_info['max_consensus'] >= 4]
    sets['medium_high'] = med_high['mouse_gene'].tolist()
    # All SG-enriched (any consensus)
    sets['all'] = gene_info['mouse_gene'].tolist()

    return sets


# =============================================
# SG Proteome gene sets (from existing analysis)
# =============================================
SG_PROTEOME_MARKMILLER = [
    "G3bp1", "G3bp2", "Tia1", "Tial1", "Eif3a", "Eif3b",
    "Eif3c", "Eif3d", "Eif3e", "Eif3f", "Eif3g", "Eif3h",
    "Eif3i", "Eif3j1", "Eif3k", "Eif3l", "Eif3m", "Eif4g1",
    "Eif4g2", "Eif4a1", "Eif4b", "Eif4e", "Eif4h", "Atxn2",
    "Atxn2l", "Caprin1", "Ubap2l", "Fmr1", "Fxr1", "Fxr2",
    "Pabpc1", "Pabpc4", "Larp1", "Larp4", "Larp4b",
]

SG_PROTEOME_YOUN = [
    "Tardbp", "Fus", "Hnrnpa1", "Hnrnpa2b1", "Matr3",
    "Ewsr1", "Taf15", "Hnrnpd", "Hnrnpdl", "Khdrbs1", "Ybx1",
    "Ybx3", "Ddx3x", "Ddx1", "Elavl1", "Elavl2", "Elavl3",
    "Elavl4", "Zfp36", "Zfp36l1", "Zfp36l2", "Cpeb1", "Cpeb2",
    "Dis3l2", "Lsm14a", "Lsm14b", "Dcp1a", "Dcp1b", "Dcp2",
    "Xrn1", "Exosc2", "Exosc3", "Exosc4", "Stau1", "Stau2",
    "Ankhd1", "Prrc2a", "Prrc2b", "Prrc2c",
]


# =============================================
# GSEA engine (adapted from analysis_gsea_literature.py)
# =============================================
def compute_es(rank_values, hits, N, n_hit):
    """Compute enrichment score using weighted Kolmogorov-Smirnov."""
    hit_set = set(hits)
    hit_weights = np.abs(rank_values[hits])
    total_hit_weight = hit_weights.sum()
    if total_hit_weight == 0:
        total_hit_weight = 1

    miss_penalty = 1.0 / (N - n_hit)
    running = np.zeros(N)
    max_es = -np.inf
    min_es = np.inf
    max_idx = 0
    min_idx = 0
    cumsum = 0.0

    for i in range(N):
        if i in hit_set:
            cumsum += np.abs(rank_values[i]) / total_hit_weight
        else:
            cumsum -= miss_penalty
        running[i] = cumsum
        if cumsum > max_es:
            max_es = cumsum
            max_idx = i
        if cumsum < min_es:
            min_es = cumsum
            min_idx = i

    if abs(max_es) >= abs(min_es):
        return max_es, max_idx, running
    else:
        return min_es, min_idx, running


def bh_correction(pvals):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    order = np.argsort(pvals)
    adj = np.zeros(n)
    for i, rank in enumerate(range(1, n + 1)):
        adj[order[i]] = pvals[order[i]] * n / rank
    adj[order[-1]] = min(adj[order[-1]], 1.0)
    for i in range(n - 2, -1, -1):
        adj[order[i]] = min(adj[order[i]], adj[order[i + 1]], 1.0)
    return adj


def preranked_gsea(ranking, gene_sets, n_perm=10000, min_size=10):
    """Run pre-ranked GSEA with permutation testing."""
    results = []
    ranking = ranking.sort_values(ascending=False)
    gene_list = ranking.index.tolist()
    rank_values = ranking.values
    N = len(gene_list)
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}

    for gs_name, gs_info in gene_sets.items():
        gs_genes = gs_info["genes"]
        hits = [gene_to_idx[g] for g in gs_genes if g in gene_to_idx]
        n_hit = len(hits)

        if n_hit < min_size:
            print(f"    SKIP {gs_info['label_short']}: only {n_hit} genes matched (min={min_size})")
            continue

        hits = sorted(hits)
        es, es_idx, running = compute_es(rank_values, hits, N, n_hit)

        # Permutation test
        perm_es = np.zeros(n_perm)
        rng = np.random.default_rng(42)
        for p in range(n_perm):
            perm_hits = sorted(rng.choice(N, n_hit, replace=False))
            perm_es[p], _, _ = compute_es(rank_values, perm_hits, N, n_hit)

        if es >= 0:
            pval = np.mean(perm_es >= es)
        else:
            pval = np.mean(perm_es <= es)

        if es >= 0:
            pos_mean = np.mean(perm_es[perm_es >= 0]) if np.any(perm_es >= 0) else 1
            nes = es / pos_mean if pos_mean > 0 else 0
        else:
            neg_mean = np.abs(np.mean(perm_es[perm_es < 0])) if np.any(perm_es < 0) else 1
            nes = es / neg_mean if neg_mean > 0 else 0

        if es >= 0:
            le_genes = [gene_list[h] for h in hits if h <= es_idx]
        else:
            le_genes = [gene_list[h] for h in hits if h >= es_idx]

        results.append({
            "gene_set": gs_name,
            "label_short": gs_info["label_short"],
            "citation": gs_info["citation"],
            "level": gs_info.get("level", "RNA"),
            "n_genes": n_hit,
            "n_total": len(gs_genes),
            "ES": es,
            "NES": nes,
            "pvalue": max(pval, 1 / n_perm),
            "leading_edge": ",".join(le_genes[:20]),
            "n_leading_edge": len(le_genes),
            "running_es": running,
        })

    # BH correction
    if results:
        pvals = np.array([r["pvalue"] for r in results])
        adj_p = bh_correction(pvals)
        for i, r in enumerate(results):
            r["adj_p"] = adj_p[i]

    return sorted(results, key=lambda x: -abs(x["NES"]))


# =============================================
# Analysis functions
# =============================================
def analysis1_gsea(ranking, gene_sets):
    """Analysis 1: Multi-source pre-ranked GSEA."""
    print("\n" + "=" * 60)
    print("Analysis 1: Pre-ranked GSEA (SG RNA + Protein gene sets)")
    print("=" * 60)
    results = preranked_gsea(ranking, gene_sets, n_perm=10000)

    print(f"\n  {'Gene Set':<45s} {'Level':>5s} {'n':>5s} {'NES':>7s} "
          f"{'pval':>10s} {'FDR':>10s}")
    print("  " + "-" * 95)
    for r in results:
        sig = ("***" if r["adj_p"] < 0.001 else
               "**" if r["adj_p"] < 0.01 else
               "*" if r["adj_p"] < 0.05 else " ")
        print(f"  {r['label_short']:<45s} {r['level']:>5s} {r['n_genes']:>5d} "
              f"{r['NES']:>7.2f} {r['pvalue']:>10.4f} {r['adj_p']:>10.4f} {sig}")

    return results


def analysis2_overlap(ir_genes, sg_gene_sets, all_expressed):
    """Analysis 2: IR genes × SG transcriptome overlap (Fisher's exact)."""
    print("\n" + "=" * 60)
    print("Analysis 2: IR gene × SG RNA overlap (Fisher's exact)")
    print("=" * 60)

    results = []
    for name, info in sg_gene_sets.items():
        # Use dataset-specific universe if available (e.g. Namkoong intersection)
        universe = info.get("universe", all_expressed)
        sg_genes = set(info["genes"]) & universe
        ir_set = ir_genes & universe

        a = len(ir_set & sg_genes)
        b = len(sg_genes - ir_set)
        c = len(ir_set - sg_genes)
        d = len(universe - ir_set - sg_genes)

        odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        overlap_genes = sorted(ir_set & sg_genes)

        bg_label = f" (universe={len(universe)})" if "universe" in info else ""
        results.append({
            "gene_set": info["label_short"],
            "level": info.get("level", "RNA"),
            "citation": info["citation"],
            "sg_in_expressed": len(sg_genes),
            "ir_in_expressed": len(ir_set),
            "background": len(universe),
            "overlap": a,
            "odds_ratio": odds,
            "pvalue": pval,
            "overlap_genes": ",".join(overlap_genes[:30]),
        })

        print(f"\n  {info['label_short']}:{bg_label}")
        print(f"    Background: {len(universe)}, SG genes: {len(sg_genes)}")
        print(f"    Overlap with IR ({len(ir_set)}): {a}")
        print(f"    OR={odds:.3f}, p={pval:.2e}")
        if a > 0 and a <= 20:
            print(f"    Genes: {', '.join(overlap_genes)}")

    # BH correction across all tests
    pvals = np.array([r["pvalue"] for r in results])
    adj_p = bh_correction(pvals)
    for i, r in enumerate(results):
        r["adj_p"] = adj_p[i]

    return results


def analysis3_distribution(khong_data, ir_genes, all_expressed):
    """Analysis 3: SG enrichment ratio distribution (IR vs all)."""
    print("\n" + "=" * 60)
    print("Analysis 3: SG enrichment ratio distribution (Khong 2017)")
    print("=" * 60)

    # Get enrichment values for mouse genes in our expression data
    khong_expr = khong_data[khong_data['mouse_gene'].isin(all_expressed)].copy()
    khong_expr['is_IR'] = khong_expr['mouse_gene'].isin(ir_genes)

    # Filter out infinite values for continuous analysis
    khong_finite = khong_expr[np.isfinite(khong_expr['log2(fold_change)'])].copy()
    ir_enrich = khong_finite[khong_finite['is_IR']]['log2(fold_change)']
    non_ir_enrich = khong_finite[~khong_finite['is_IR']]['log2(fold_change)']

    print(f"\n  Genes matched (Khong → expressed): {len(khong_expr)}")
    print(f"  Genes with finite enrichment: {len(khong_finite)}")
    print(f"  IR genes matched (finite): {len(ir_enrich)}")
    print(f"  Non-IR genes (finite): {len(non_ir_enrich)}")

    print(f"\n  SG enrichment (log2 SG/Total):")
    print(f"    IR genes:     mean={ir_enrich.mean():.3f}, median={ir_enrich.median():.3f}")
    print(f"    Non-IR genes: mean={non_ir_enrich.mean():.3f}, median={non_ir_enrich.median():.3f}")

    # KS test
    ks_stat, ks_p = stats.ks_2samp(ir_enrich, non_ir_enrich)
    print(f"    KS test: D={ks_stat:.3f}, p={ks_p:.2e}")

    # Mann-Whitney
    mw_stat, mw_p = stats.mannwhitneyu(ir_enrich, non_ir_enrich, alternative='two-sided')
    print(f"    Mann-Whitney: U={mw_stat:.0f}, p={mw_p:.2e}")

    # Classification distribution (uses all genes, not just finite)
    ir_all = khong_expr[khong_expr['is_IR']]
    ir_class = ir_all['classification'].value_counts()
    all_class = khong_expr['classification'].value_counts()
    print(f"\n  SG classification of IR genes:")
    for cls in ['enriched', 'depleted', 'neither']:
        ir_n = ir_class.get(cls, 0)
        all_n = all_class.get(cls, 0)
        ir_pct = 100 * ir_n / max(1, len(ir_all))
        all_pct = 100 * all_n / max(1, len(khong_expr))
        print(f"    {cls:>10s}: IR={ir_n:>4d} ({ir_pct:5.1f}%)  All={all_n:>5d} ({all_pct:5.1f}%)")

    return khong_expr


def analysis4_length_confound(khong_data, ir_genes, ir_data, deseq2, all_expressed):
    """Analysis 4: Transcript length confound check."""
    print("\n" + "=" * 60)
    print("Analysis 4: Transcript length confound check")
    print("=" * 60)

    # Get intron lengths as proxy for gene length
    gene_lengths = ir_data.groupby('gene_name')['intron_length'].sum().to_dict()

    khong_expr = khong_data[khong_data['mouse_gene'].isin(all_expressed)].copy()
    khong_expr['is_IR'] = khong_expr['mouse_gene'].isin(ir_genes)

    # Add gene length info where available
    khong_expr['total_intron_length'] = khong_expr['mouse_gene'].map(gene_lengths)
    has_length = khong_expr.dropna(subset=['total_intron_length'])

    if len(has_length) < 50:
        print("  Insufficient genes with length data for confound analysis")
        return None

    # Correlation: gene length vs SG enrichment
    rho_all, rho_p = stats.spearmanr(
        has_length['total_intron_length'],
        has_length['log2(fold_change)']
    )
    print(f"\n  Total intron length vs SG enrichment (all genes with both):")
    print(f"    n={len(has_length)}, Spearman rho={rho_all:.3f}, p={rho_p:.2e}")

    # Length-binned Fisher's exact test
    # Bin genes by total intron length into quartiles
    has_length['length_quartile'] = pd.qcut(
        has_length['total_intron_length'], 4,
        labels=['Q1_short', 'Q2', 'Q3', 'Q4_long']
    )

    print(f"\n  Length-binned overlap (IR × SG-enriched):")
    print(f"  {'Quartile':<12s} {'n_genes':>8s} {'n_IR':>6s} {'n_SG_enr':>9s} "
          f"{'overlap':>8s} {'OR':>8s} {'p':>12s}")
    print("  " + "-" * 70)

    binned_results = []
    for q in ['Q1_short', 'Q2', 'Q3', 'Q4_long']:
        subset = has_length[has_length['length_quartile'] == q]
        n_total = len(subset)
        n_ir = subset['is_IR'].sum()
        n_sg = (subset['classification'] == 'enriched').sum()
        n_both = ((subset['is_IR']) & (subset['classification'] == 'enriched')).sum()

        a = n_both
        b = n_sg - n_both
        c = n_ir - n_both
        d = n_total - n_ir - n_sg + n_both

        if min(a + b, c + d, a + c, b + d) > 0:
            odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
        else:
            odds, pval = float('nan'), float('nan')

        print(f"  {q:<12s} {n_total:>8d} {n_ir:>6d} {n_sg:>9d} "
              f"{n_both:>8d} {odds:>8.2f} {pval:>12.2e}")
        binned_results.append({
            'quartile': q, 'n_genes': n_total, 'n_ir': n_ir,
            'n_sg_enriched': n_sg, 'overlap': n_both,
            'odds_ratio': odds, 'pvalue': pval
        })

    return binned_results


def extract_per_dataset_genesets(meta_df):
    """Extract per-dataset gene sets from Biancon 2025 meta-analysis."""
    datasets = {}
    for ds_id in meta_df['dataset_ID'].unique():
        sub = meta_df[meta_df['dataset_ID'] == ds_id]
        genes_human = sub['gene_symbol'].dropna().unique().tolist()
        genes_mouse = batch_human_to_mouse(genes_human)
        stress = sub['stress_type'].iloc[0] if 'stress_type' in sub.columns else 'unknown'
        system = sub['system'].iloc[0] if 'system' in sub.columns else 'unknown'
        technique = sub['technique'].iloc[0] if 'technique' in sub.columns else 'unknown'
        datasets[ds_id] = {
            'genes_mouse': genes_mouse,
            'n_human': len(genes_human),
            'stress': stress,
            'system': system,
            'technique': technique,
        }
    return datasets


def load_meta2025_table1(filepath):
    """Load Biancon 2025 Table 1 (dataset metadata)."""
    import openpyxl
    wb = openpyxl.load_workbook(filepath)
    ws = wb[wb.sheetnames[0]]
    header = [cell.value for cell in ws[1]]
    data = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        d = dict(zip(header, row))
        if d.get('Dataset ID'):
            data[d['Dataset ID']] = {
                'reference': d.get('Reference', ''),
                'technique_full': d.get('Technique', ''),
                'system': d.get('System', ''),
                'stress': d.get('Stress Type', ''),
            }
    return data


def analysis6_per_dataset_overlap(ir_genes, meta_df, all_expressed, meta_table1=None):
    """Analysis 6: Per-dataset IR × SG overlap across all 26 meta-analysis datasets."""
    print("\n" + "=" * 60)
    print("Analysis 6: Per-dataset IR × SG RNA overlap (26 datasets)")
    print("=" * 60)

    per_ds = extract_per_dataset_genesets(meta_df)
    ir_set = ir_genes & all_expressed

    results = []
    for ds_id, ds_info in sorted(per_ds.items()):
        sg_genes = set(ds_info['genes_mouse']) & all_expressed
        if len(sg_genes) < 10:
            continue

        a = len(ir_set & sg_genes)
        b = len(sg_genes - ir_set)
        c = len(ir_set - sg_genes)
        d = len(all_expressed - ir_set - sg_genes)

        odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')

        ref = ''
        if meta_table1 and ds_id in meta_table1:
            ref = meta_table1[ds_id].get('reference', '')

        results.append({
            'dataset_id': ds_id,
            'reference': ref,
            'stress': ds_info['stress'],
            'system': ds_info['system'],
            'technique': ds_info['technique'],
            'n_sg_expressed': len(sg_genes),
            'n_ir': len(ir_set),
            'overlap': a,
            'odds_ratio': odds,
            'pvalue': pval,
            'log2_or': np.log2(odds) if odds > 0 else 0,
        })

    # BH correction
    pvals = np.array([r['pvalue'] for r in results])
    adj_p = bh_correction(pvals)
    for i, r in enumerate(results):
        r['adj_p'] = adj_p[i]

    # Sort by OR
    results.sort(key=lambda x: -x['odds_ratio'])

    print(f"\n  {'Dataset':<50s} {'Ref':<25s} {'n_SG':>6s} {'Ovlp':>6s} "
          f"{'OR':>7s} {'FDR':>10s}")
    print("  " + "-" * 110)
    n_sig = 0
    for r in results:
        sig = "***" if r['adj_p'] < 0.001 else "**" if r['adj_p'] < 0.01 else "*" if r['adj_p'] < 0.05 else " "
        if r['adj_p'] < 0.05:
            n_sig += 1
        print(f"  {r['dataset_id']:<50s} {r['reference']:<25s} {r['n_sg_expressed']:>6d} "
              f"{r['overlap']:>6d} {r['odds_ratio']:>7.2f} {r['adj_p']:>10.4f} {sig}")

    print(f"\n  Summary: {n_sig}/{len(results)} datasets show significant overlap (FDR<0.05)")
    sig_pos = sum(1 for r in results if r['adj_p'] < 0.05 and r['odds_ratio'] > 1)
    print(f"  Of significant: {sig_pos} positive (OR>1), {n_sig - sig_pos} negative")

    return results


def analysis5_rna_vs_protein(gsea_results):
    """Analysis 5: RNA vs Protein level comparison."""
    print("\n" + "=" * 60)
    print("Analysis 5: RNA vs Protein level GSEA comparison")
    print("=" * 60)

    rna_results = [r for r in gsea_results if r['level'] == 'RNA']
    protein_results = [r for r in gsea_results if r['level'] == 'Protein']

    print(f"\n  {'Level':<8s} {'Gene Set':<45s} {'NES':>7s} {'FDR':>10s}")
    print("  " + "-" * 75)
    for r in rna_results:
        sig = "***" if r["adj_p"] < 0.001 else "**" if r["adj_p"] < 0.01 else "*" if r["adj_p"] < 0.05 else "ns"
        print(f"  {'RNA':<8s} {r['label_short']:<45s} {r['NES']:>7.2f} {r['adj_p']:>10.4f} {sig}")
    print()
    for r in protein_results:
        sig = "***" if r["adj_p"] < 0.001 else "**" if r["adj_p"] < 0.01 else "*" if r["adj_p"] < 0.05 else "ns"
        print(f"  {'Protein':<8s} {r['label_short']:<45s} {r['NES']:>7.2f} {r['adj_p']:>10.4f} {sig}")


def analysis7_deg_overlap(deseq2, sg_gene_sets, all_expressed):
    """Analysis 7: DEG (UP/DOWN) × SG gene set overlap (Fisher's exact)."""
    print("\n" + "=" * 60)
    print("Analysis 7: DEG × SG gene set overlap (Fisher's exact)")
    print("=" * 60)

    deg_up_all = set(deseq2[(deseq2['padj'] < 0.05) & (deseq2['log2FoldChange'] > 0)].index) & all_expressed
    deg_down_all = set(deseq2[(deseq2['padj'] < 0.05) & (deseq2['log2FoldChange'] < 0)].index) & all_expressed
    print(f"  DEG UP: {len(deg_up_all)}, DEG DOWN: {len(deg_down_all)}")

    results = []
    for name, info in sg_gene_sets.items():
        # Use dataset-specific universe if available
        universe = info.get("universe", all_expressed)
        sg_genes = set(info["genes"]) & universe
        for deg_label, deg_set_all in [("DEG_UP", deg_up_all), ("DEG_DOWN", deg_down_all)]:
            deg_set = deg_set_all & universe
            a = len(deg_set & sg_genes)
            b = len(sg_genes - deg_set)
            c = len(deg_set - sg_genes)
            d = len(universe - deg_set - sg_genes)
            odds, pval = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            results.append({
                "gene_set": info["label_short"],
                "deg_direction": deg_label,
                "level": info.get("level", "RNA"),
                "citation": info["citation"],
                "sg_in_expressed": len(sg_genes),
                "deg_in_expressed": len(deg_set),
                "background": len(universe),
                "overlap": a,
                "odds_ratio": odds,
                "pvalue": pval,
            })

    pvals = np.array([r["pvalue"] for r in results])
    adj_p = bh_correction(pvals)
    for i, r in enumerate(results):
        r["adj_p"] = adj_p[i]

    print(f"\n  {'Gene Set':<45s} {'Dir':>8s} {'SG':>6s} {'DEG':>6s} "
          f"{'Ovlp':>6s} {'OR':>7s} {'FDR':>10s}")
    print("  " + "-" * 95)
    for r in results:
        sig = "***" if r['adj_p'] < 0.001 else "**" if r['adj_p'] < 0.01 else "*" if r['adj_p'] < 0.05 else " "
        print(f"  {r['gene_set']:<45s} {r['deg_direction']:>8s} {r['sg_in_expressed']:>6d} "
              f"{r['deg_in_expressed']:>6d} {r['overlap']:>6d} {r['odds_ratio']:>7.2f} "
              f"{r['adj_p']:>10.4f} {sig}")

    return results, deg_up_all, deg_down_all


def plot_venn_diagrams(ir_genes, deg_up, deg_down, sg_gene_sets, all_expressed, output_base):
    """Generate Venn diagrams for IR∩SG and DEG-UP∩SG for each enriched dataset."""
    from matplotlib_venn import venn2

    # Select only enriched/proximal sets (skip depleted/excluded)
    enriched_sets = {}
    for name, info in sg_gene_sets.items():
        label = info["label_short"]
        if any(kw in label.lower() for kw in ["depleted", "excluded"]):
            continue
        enriched_sets[label] = set(info["genes"]) & all_expressed

    ir_set = ir_genes & all_expressed
    deg_up_set = deg_up & all_expressed

    n_sets = len(enriched_sets)
    n_cols = min(4, n_sets)
    n_rows = (n_sets + n_cols - 1) // n_cols * 2  # 2 rows per set row (IR + DEG)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    for idx, (label, sg_set) in enumerate(enriched_sets.items()):
        row_base = (idx // n_cols) * 2
        col = idx % n_cols

        # IR ∩ SG Venn
        ax = axes[row_base, col]
        v_ir = venn2([ir_set, sg_set], set_labels=('IR genes', 'SG RNA'),
                     ax=ax, set_colors=('#e65100', '#1565c0'), alpha=0.5)
        ir_only = len(ir_set - sg_set)
        sg_only = len(sg_set - ir_set)
        both = len(ir_set & sg_set)
        odds_ir, p_ir = stats.fisher_exact([
            [both, sg_only], [ir_only, len(all_expressed) - len(ir_set) - sg_only]
        ], alternative="two-sided")
        sig_ir = "***" if p_ir < 0.001 else "**" if p_ir < 0.01 else "*" if p_ir < 0.05 else "ns"
        short_label = label.replace("(", "\n(")
        ax.set_title(f"IR x {short_label}\nOR={odds_ir:.2f} {sig_ir}", fontsize=8, fontweight='bold')

        # DEG-UP ∩ SG Venn
        ax = axes[row_base + 1, col]
        v_deg = venn2([deg_up_set, sg_set], set_labels=('DEG UP', 'SG RNA'),
                      ax=ax, set_colors=('#c62828', '#1565c0'), alpha=0.5)
        deg_only = len(deg_up_set - sg_set)
        sg_only2 = len(sg_set - deg_up_set)
        both2 = len(deg_up_set & sg_set)
        odds_deg, p_deg = stats.fisher_exact([
            [both2, sg_only2], [deg_only, len(all_expressed) - len(deg_up_set) - sg_only2]
        ], alternative="two-sided")
        sig_deg = "***" if p_deg < 0.001 else "**" if p_deg < 0.01 else "*" if p_deg < 0.05 else "ns"
        ax.set_title(f"DEG-UP x {short_label}\nOR={odds_deg:.2f} {sig_deg}", fontsize=8, fontweight='bold')

    # Hide unused axes
    for i in range(n_sets, n_rows * n_cols // 2):
        row_base = (i // n_cols) * 2
        col = i % n_cols
        if row_base < n_rows and col < n_cols:
            axes[row_base, col].set_visible(False)
            axes[row_base + 1, col].set_visible(False)

    plt.suptitle("IR genes and DEG-UP genes overlap with SG-enriched RNA sets\n"
                 "(within MEF-expressed gene background)",
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f"{output_base}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_base}.png", dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.pdf / .png")


# =============================================
# Plotting functions
# =============================================
def plot_gsea_barplot(results, output_base):
    """Bar plot of NES for all gene sets, colored by RNA vs Protein level."""
    fig, ax = plt.subplots(figsize=(12, 7))

    names = []
    nes_vals = []
    fdr_vals = []
    levels = []
    for r in reversed(results):
        names.append(r["label_short"])
        nes_vals.append(r["NES"])
        fdr_vals.append(r["adj_p"])
        levels.append(r["level"])

    y_pos = np.arange(len(names))
    colors = []
    for n, fdr, level in zip(nes_vals, fdr_vals, levels):
        if fdr >= 0.05:
            colors.append("#9e9e9e")
        elif level == "RNA":
            if n > 0:
                colors.append("#e65100")  # orange-red for RNA up
            else:
                colors.append("#01579b")  # dark blue for RNA down
        else:  # Protein
            if n > 0:
                colors.append("#c62828")  # red for protein up
            else:
                colors.append("#1565c0")  # blue for protein down

    # Use two-column layout: bars on left, stats table on right
    bars = ax.barh(y_pos, nes_vals, color=colors, edgecolor="none", height=0.55, alpha=0.9)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=7.5)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
    ax.set_title("Pre-ranked GSEA: SG Transcriptome vs Proteome\n"
                 "(Ranking: DESeq2 Wald statistic, 10,000 permutations)",
                 fontsize=12, fontweight='bold')
    ax.axvline(0, color="black", linewidth=0.5)

    # Place annotations inside bars or at bar edge using compact format
    for i, (n, fdr, r) in enumerate(zip(nes_vals, fdr_vals, reversed(results))):
        fdr_str = "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else "ns"
        label = f"{n:+.2f} {fdr_str} (n={r['n_genes']})"
        x_offset = 0.06
        if n >= 0:
            ax.text(n + x_offset, i, label, va="center", ha="left", fontsize=6.5,
                    fontfamily='monospace')
        else:
            ax.text(n - x_offset, i, label, va="center", ha="right", fontsize=6.5,
                    fontfamily='monospace')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e65100', label='RNA enriched/proximal (sig)'),
        Patch(facecolor='#01579b', label='RNA depleted/excluded (sig)'),
        Patch(facecolor='#c62828', label='Protein level (sig)'),
        Patch(facecolor='#9e9e9e', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7)

    ax.text(0.98, 0.02, "*** FDR<0.001  ** FDR<0.01  * FDR<0.05  ns: not significant",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="grey")

    ax.set_xlim(ax.get_xlim()[0] - 0.3, ax.get_xlim()[1] + 1.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{output_base}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_base}.png", dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.pdf / .png")


def plot_overlap(overlap_results, output_base):
    """Bar plot of overlap statistics."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Odds ratios
    ax = axes[0]
    names = [r['gene_set'] for r in overlap_results]
    ors = [r['odds_ratio'] for r in overlap_results]
    pvals = [r['adj_p'] for r in overlap_results]
    levels = [r['level'] for r in overlap_results]
    overlaps = [r['overlap'] for r in overlap_results]

    y_pos = np.arange(len(names))
    colors = ['#e65100' if l == 'RNA' else '#c62828' for l in levels]
    alphas = [0.9 if p < 0.05 else 0.4 for p in pvals]

    for i, (y, o, c, a) in enumerate(zip(y_pos, ors, colors, alphas)):
        ax.barh(y, np.log2(o) if o > 0 else 0, color=c, alpha=a, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("log2(Odds Ratio)", fontsize=10)
    ax.set_title("IR Gene × SG Gene Set Overlap\n(Fisher's exact test)", fontsize=11, fontweight='bold')
    ax.axvline(0, color="black", linewidth=0.5)

    for i, (o, p, n_ov) in enumerate(zip(ors, pvals, overlaps)):
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        x = np.log2(o) if o > 0 else 0
        ax.text(max(x, 0) + 0.1, i, f"OR={o:.2f} {sig} (n={n_ov})",
                va="center", ha="left", fontsize=7.5)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel B: Overlap counts
    ax = axes[1]
    rna_results = [r for r in overlap_results if r['level'] == 'RNA']
    prot_results = [r for r in overlap_results if r['level'] == 'Protein']

    categories = ['RNA level\n(SG transcriptome)', 'Protein level\n(SG proteome)']
    rna_overlaps = [r['overlap'] for r in rna_results]
    prot_overlaps = [r['overlap'] for r in prot_results]

    x = np.arange(len(rna_results))
    width = 0.35
    if rna_results:
        ax.bar(x, rna_overlaps, width, color='#e65100', alpha=0.8,
               label='RNA sets')
    if prot_results:
        ax.bar(x[:len(prot_results)] + width, prot_overlaps, width,
               color='#c62828', alpha=0.8, label='Protein sets')

    all_results = rna_results + prot_results
    all_names = [r['gene_set'] for r in all_results]
    ax.set_xticks(np.arange(max(len(rna_results), len(prot_results))) + width / 2)
    tick_labels = []
    for i in range(max(len(rna_results), len(prot_results))):
        parts = []
        if i < len(rna_results):
            parts.append(rna_results[i]['gene_set'][:25])
        if i < len(prot_results):
            parts.append(prot_results[i]['gene_set'][:25])
        tick_labels.append('\n'.join(parts))
    ax.set_xticklabels(tick_labels, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel("Number of overlapping genes", fontsize=10)
    ax.set_title("Overlap counts: IR genes × SG gene sets", fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(f"{output_base}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_base}.png", dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.pdf / .png")


def plot_distribution(khong_expr, output_base):
    """SG enrichment ratio distribution (IR vs all)."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Filter to finite values for continuous plots
    finite = khong_expr[np.isfinite(khong_expr['log2(fold_change)'])].copy()
    ir_data = finite[finite['is_IR']]
    non_ir_data = finite[~finite['is_IR']]

    # Panel A: Histogram
    ax = axes[0]
    bins = np.linspace(-6, 6, 60)
    ax.hist(non_ir_data['log2(fold_change)'].clip(-6, 6), bins=bins, alpha=0.5,
            label=f"Non-IR (n={len(non_ir_data)})", color="#757575", density=True)
    ax.hist(ir_data['log2(fold_change)'].clip(-6, 6), bins=bins, alpha=0.7,
            label=f"IR (n={len(ir_data)})", color="#e65100", density=True)
    ax.set_xlabel("log2(SG/Total) enrichment")
    ax.set_ylabel("Density")
    ax.set_title("SG enrichment distribution\n(Khong et al. 2017)")
    ax.legend(fontsize=8)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")

    # Panel B: CDF
    ax = axes[1]
    ir_sorted = np.sort(ir_data['log2(fold_change)'].values)
    non_ir_sorted = np.sort(non_ir_data['log2(fold_change)'].values)
    ax.plot(non_ir_sorted, np.arange(len(non_ir_sorted)) / len(non_ir_sorted),
            color="#757575", label="Non-IR", linewidth=1.5)
    ax.plot(ir_sorted, np.arange(len(ir_sorted)) / len(ir_sorted),
            color="#e65100", label="IR", linewidth=1.5)

    ks_stat, ks_p = stats.ks_2samp(ir_data['log2(fold_change)'],
                                    non_ir_data['log2(fold_change)'])
    ax.set_xlabel("log2(SG/Total) enrichment")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title(f"CDF comparison\nKS D={ks_stat:.3f}, p={ks_p:.2e}")
    ax.legend(fontsize=8)
    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")

    # Panel C: Classification proportions (use all data, not just finite)
    ax = axes[2]
    ir_all = khong_expr[khong_expr['is_IR']]
    non_ir_all = khong_expr[~khong_expr['is_IR']]
    classes = ['enriched', 'neither', 'depleted']
    ir_counts = [ir_all['classification'].eq(c).sum() for c in classes]
    non_ir_counts = [non_ir_all['classification'].eq(c).sum() for c in classes]
    ir_pcts = [100 * c / max(1, len(ir_all)) for c in ir_counts]
    non_ir_pcts = [100 * c / max(1, len(non_ir_all)) for c in non_ir_counts]

    x = np.arange(len(classes))
    width = 0.35
    ax.bar(x - width / 2, ir_pcts, width, label='IR genes', color='#e65100', alpha=0.8)
    ax.bar(x + width / 2, non_ir_pcts, width, label='Non-IR', color='#757575', alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(['SG-enriched', 'Neither', 'SG-depleted'])
    ax.set_ylabel("Percentage of genes (%)")
    ax.set_title("SG classification\n(Khong 2017 categories)")
    ax.legend(fontsize=8)

    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(f"{output_base}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_base}.png", dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.pdf / .png")


def plot_length_control(khong_data, ir_genes, ir_all_data, all_expressed, output_base):
    """Length confound check plots."""
    khong_expr = khong_data[khong_data['mouse_gene'].isin(all_expressed)].copy()
    khong_expr['is_IR'] = khong_expr['mouse_gene'].isin(ir_genes)

    gene_lengths = ir_all_data.groupby('gene_name')['intron_length'].sum()
    khong_expr['total_intron_length'] = khong_expr['mouse_gene'].map(gene_lengths)
    has_length = khong_expr.dropna(subset=['total_intron_length'])

    if len(has_length) < 50:
        print("  Insufficient data for length control plot")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Length vs SG enrichment scatter
    ax = axes[0]
    non_ir = has_length[~has_length['is_IR']]
    ir = has_length[has_length['is_IR']]
    ax.scatter(np.log10(non_ir['total_intron_length'] + 1),
               non_ir['log2(fold_change)'],
               s=3, alpha=0.1, c='#757575', label='Non-IR', rasterized=True)
    ax.scatter(np.log10(ir['total_intron_length'] + 1),
               ir['log2(fold_change)'],
               s=8, alpha=0.5, c='#e65100', label='IR', rasterized=True)
    rho, rho_p = stats.spearmanr(has_length['total_intron_length'],
                                  has_length['log2(fold_change)'])
    ax.set_xlabel("log10(total intron length)")
    ax.set_ylabel("log2(SG/Total) enrichment")
    ax.set_title(f"Length vs SG enrichment\nSpearman ρ={rho:.3f}, p={rho_p:.2e}")
    ax.legend(fontsize=8, markerscale=3)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")

    # Panel B: Length-binned OR
    has_length_copy = has_length.copy()
    has_length_copy['length_quartile'] = pd.qcut(
        has_length_copy['total_intron_length'], 4,
        labels=['Q1\n(short)', 'Q2', 'Q3', 'Q4\n(long)']
    )
    ax = axes[1]
    quartiles = ['Q1\n(short)', 'Q2', 'Q3', 'Q4\n(long)']
    ors = []
    for q in quartiles:
        subset = has_length_copy[has_length_copy['length_quartile'] == q]
        n_ir = subset['is_IR'].sum()
        n_sg = (subset['classification'] == 'enriched').sum()
        n_both = ((subset['is_IR']) & (subset['classification'] == 'enriched')).sum()
        n_total = len(subset)
        a, b = n_both, n_sg - n_both
        c, d = n_ir - n_both, n_total - n_ir - n_sg + n_both
        if min(a + b, c + d, a + c, b + d) > 0:
            odds, _ = stats.fisher_exact([[a, b], [c, d]])
        else:
            odds = float('nan')
        ors.append(odds)

    x_pos = np.arange(len(quartiles))
    colors = ['#e65100' if o > 1 else '#01579b' for o in ors]
    ax.bar(x_pos, [np.log2(o) if o > 0 and not np.isnan(o) else 0 for o in ors],
           color=colors, alpha=0.8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(quartiles)
    ax.set_ylabel("log2(Odds Ratio)")
    ax.set_xlabel("Gene length quartile")
    ax.set_title("Length-controlled overlap\n(IR × SG-enriched)")
    ax.axhline(0, color="black", linewidth=0.5)
    for i, o in enumerate(ors):
        if not np.isnan(o):
            ax.text(i, np.log2(o) + 0.05 if o > 1 else np.log2(o) - 0.1,
                    f"OR={o:.2f}", ha='center', fontsize=8)

    # Panel C: IR gene length distribution vs non-IR
    ax = axes[2]
    ir_lengths = has_length[has_length['is_IR']]['total_intron_length']
    non_ir_lengths = has_length[~has_length['is_IR']]['total_intron_length']
    ax.hist(np.log10(non_ir_lengths + 1), bins=50, alpha=0.5,
            label=f"Non-IR (n={len(non_ir_lengths)})", color="#757575", density=True)
    ax.hist(np.log10(ir_lengths + 1), bins=50, alpha=0.7,
            label=f"IR (n={len(ir_lengths)})", color="#e65100", density=True)
    mw_stat, mw_p = stats.mannwhitneyu(ir_lengths, non_ir_lengths, alternative='two-sided')
    ax.set_xlabel("log10(total intron length)")
    ax.set_ylabel("Density")
    ax.set_title(f"Gene length distribution\nMW p={mw_p:.2e}")
    ax.legend(fontsize=8)

    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(f"{output_base}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_base}.png", dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.pdf / .png")


def plot_per_dataset_forest(per_ds_results, output_base):
    """Forest plot of OR across all 26 meta-analysis datasets."""
    # Sort by OR
    results = sorted(per_ds_results, key=lambda x: x['odds_ratio'])
    n = len(results)

    fig, ax = plt.subplots(figsize=(12, max(6, n * 0.35)))

    y_pos = np.arange(n)
    log2_ors = [r['log2_or'] for r in results]
    colors = []
    for r in results:
        if r['adj_p'] >= 0.05:
            colors.append('#9e9e9e')
        elif r['odds_ratio'] > 1:
            colors.append('#e65100')
        else:
            colors.append('#01579b')

    bars = ax.barh(y_pos, log2_ors, color=colors, edgecolor='none', height=0.6, alpha=0.85)

    # Labels
    labels = []
    for r in results:
        ref_short = r['reference'].split(',')[0] if r['reference'] else r['dataset_id'][:20]
        stress_short = r['stress'][:10]
        labels.append(f"{ref_short} ({r['system'][:8]}, {stress_short})")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)

    # Annotations
    for i, r in enumerate(results):
        sig = "***" if r['adj_p'] < 0.001 else "**" if r['adj_p'] < 0.01 else "*" if r['adj_p'] < 0.05 else "ns"
        label = f"OR={r['odds_ratio']:.2f} {sig} (n={r['overlap']})"
        x = r['log2_or']
        offset = 0.05
        if x >= 0:
            ax.text(x + offset, i, label, va='center', ha='left', fontsize=6.5)
        else:
            ax.text(x - offset, i, label, va='center', ha='right', fontsize=6.5)

    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('log2(Odds Ratio)', fontsize=11)
    ax.set_title('IR Gene × SG RNA Overlap: Per-Dataset Reproducibility\n'
                 '(Biancon 2025 meta-analysis, 26 individual datasets)',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Summary text
    n_sig = sum(1 for r in results if r['adj_p'] < 0.05 and r['odds_ratio'] > 1)
    n_total = len(results)
    ax.text(0.98, 0.02,
            f"{n_sig}/{n_total} datasets: significant positive overlap (FDR<0.05)\n"
            "*** FDR<0.001  ** FDR<0.01  * FDR<0.05  ns: not significant",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=7, color='grey')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e65100', label='Significant (OR>1)'),
        Patch(facecolor='#01579b', label='Significant (OR<1)'),
        Patch(facecolor='#9e9e9e', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=8)

    plt.tight_layout()
    fig.savefig(f"{output_base}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_base}.png", dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_base}.pdf / .png")


# =============================================
# Main
# =============================================
def main():
    print("=" * 70)
    print("SG Transcriptome GSEA: Do IR mRNAs go to SGs or stay nuclear?")
    print("=" * 70)

    # --- Load input data ---
    print("\n--- Loading input data ---")

    deseq2 = pd.read_csv(os.path.join(DATADIR, "deseq2_results.csv"), index_col=0)
    deseq2.index = deseq2.index.str.strip('"')
    deseq2 = deseq2.dropna(subset=["stat"])
    ranking = deseq2["stat"]
    all_expressed = set(deseq2.index)
    print(f"  DESeq2 genes: {len(deseq2)}")

    retained = pd.read_csv(os.path.join(DATADIR, "limma_retained_dhx8.tsv"), sep="\t")
    ir_genes = set(retained["gene_name"].unique())
    print(f"  Retained introns: {len(retained)} ({len(ir_genes)} genes)")

    all_ir = pd.read_csv(os.path.join(DATADIR, "limma_ir_results_dhx8.tsv"), sep="\t")
    print(f"  All tested introns: {len(all_ir)}")

    # --- Load SG transcriptome data ---
    print("\n--- Loading SG transcriptome data ---")

    # Khong 2017
    khong_path = os.path.join(DATADIR, "GSE99304_U2OSSGcuffdiff.txt")
    khong = load_khong2017(khong_path)
    khong_enriched = set(khong[khong['classification'] == 'enriched']['mouse_gene'])
    khong_depleted = set(khong[khong['classification'] == 'depleted']['mouse_gene'])
    print(f"  Khong 2017: {len(khong)} genes (OK status, human→mouse)")
    print(f"    Enriched: {len(khong_enriched)}, Depleted: {len(khong_depleted)}")
    print(f"    Enriched in expressed: {len(khong_enriched & all_expressed)}")
    print(f"    Depleted in expressed: {len(khong_depleted & all_expressed)}")

    # Meta-analysis 2025
    meta_path = os.path.join(DATADIR, "meta2025_table3.xlsx")
    meta_df = load_meta2025(meta_path)
    meta_sets = extract_meta2025_genesets(meta_df)
    print(f"  Meta-analysis 2025: {len(meta_sets['all'])} total genes")
    print(f"    High consensus (>=8): {len(meta_sets['high'])}")
    print(f"    Medium-high (>=4): {len(meta_sets['medium_high'])}")
    print(f"    High in expressed: {len(set(meta_sets['high']) & all_expressed)}")
    print(f"    Medium-high in expressed: {len(set(meta_sets['medium_high']) & all_expressed)}")

    # --- Load Namkoong 2018 data (3 stress conditions) ---
    print("\n--- Loading Namkoong 2018 data ---")
    namkoong = {}
    namkoong_files = {
        'ER': os.path.join(DATADIR, "namkoong2018_mmc3.xlsx"),
        'HS': os.path.join(DATADIR, "namkoong2018_mmc4.xlsx"),
        'AS': os.path.join(DATADIR, "namkoong2018_mmc5.xlsx"),
    }
    for stress, fpath in namkoong_files.items():
        if os.path.exists(fpath):
            nk = load_namkoong2018(fpath, stress)
            namkoong[stress] = nk
            print(f"  Namkoong 2018 {stress}: {len(nk['enriched'])} enriched (C1-2), "
                  f"{len(nk['depleted'])} depleted (C5-6), "
                  f"enriched in expressed: {len(set(nk['enriched']) & all_expressed)}")
        else:
            print(f"  Namkoong 2018 {stress}: file not found ({fpath})")

    # --- Load Moon 2023 CAP-seq data ---
    print("\n--- Loading Moon 2023 CAP-seq data ---")
    moon = {}
    moon_path = os.path.join(DATADIR, "moon2023_suppdata2.xlsx")
    if os.path.exists(moon_path):
        for sheet, label in [('A-arsenite', 'arsenite'), ('B-sorbitol', 'sorbitol')]:
            mn = load_moon2023(moon_path, sheet)
            moon[label] = mn
            print(f"  Moon 2023 {label}: {len(mn['enriched'])} SG-proximal (mouse), "
                  f"{len(mn['depleted'])} SG-excluded (mouse), "
                  f"proximal in expressed: {len(set(mn['enriched']) & all_expressed)}")
    else:
        print(f"  Moon 2023: file not found ({moon_path})")

    # --- Define all gene sets ---
    gene_sets = {
        # RNA level - Khong 2017
        "SG RNA Enriched\n(Khong 2017)": {
            "genes": list(khong_enriched),
            "citation": "Khong et al. Mol Cell 2017; PMID:29129640",
            "label_short": "SG RNA Enriched (Khong 2017)",
            "level": "RNA",
        },
        "SG RNA Depleted\n(Khong 2017)": {
            "genes": list(khong_depleted),
            "citation": "Khong et al. Mol Cell 2017; PMID:29129640",
            "label_short": "SG RNA Depleted (Khong 2017)",
            "level": "RNA",
        },
        # RNA level - Meta-analysis 2025
        "SG RNA High Consensus\n(Biancon 2025)": {
            "genes": meta_sets['high'],
            "citation": "Biancon et al. RNA 2025; PMID:40086831",
            "label_short": "SG RNA High Consensus (Biancon 2025)",
            "level": "RNA",
        },
        "SG RNA Med-High Consensus\n(Biancon 2025)": {
            "genes": meta_sets['medium_high'],
            "citation": "Biancon et al. RNA 2025; PMID:40086831",
            "label_short": "SG RNA Med-High Consensus (Biancon 2025)",
            "level": "RNA",
        },
        # Protein level (existing)
        "SG Proteome\n(Markmiller 2018)": {
            "genes": SG_PROTEOME_MARKMILLER,
            "citation": "Markmiller et al. Cell 2018; PMID:29551268",
            "label_short": "SG Proteome (Markmiller 2018)",
            "level": "Protein",
        },
        "SG BioID\n(Youn 2018)": {
            "genes": SG_PROTEOME_YOUN,
            "citation": "Youn et al. Mol Cell 2018; PMID:30078726",
            "label_short": "SG BioID (Youn 2018)",
            "level": "Protein",
        },
    }

    # Add Namkoong 2018 gene sets (3 stress conditions)
    # Use intersection background: only genes measured by Namkoong AND expressed in MEF
    for stress, nk in namkoong.items():
        namkoong_universe = set(nk['all_genes']) & all_expressed
        gene_sets[f"RG Enriched C1-2\n(Namkoong 2018 {stress})"] = {
            "genes": nk['enriched'],
            "citation": "Namkoong et al. Mol Cell 2018; PMID:29576526",
            "label_short": f"RG Enriched C1-2 (Namkoong 2018 {stress})",
            "level": "RNA",
            "universe": namkoong_universe,
        }
        gene_sets[f"RG Depleted C5-6\n(Namkoong 2018 {stress})"] = {
            "genes": nk['depleted'],
            "citation": "Namkoong et al. Mol Cell 2018; PMID:29576526",
            "label_short": f"RG Depleted C5-6 (Namkoong 2018 {stress})",
            "level": "RNA",
            "universe": namkoong_universe,
        }

    # Add Moon 2023 gene sets (2 stress conditions)
    for stress, mn in moon.items():
        gene_sets[f"SG-proximal\n(Moon 2023 {stress})"] = {
            "genes": mn['enriched'],
            "citation": "Moon et al. Nat Commun 2023; PMID:37968266",
            "label_short": f"SG-proximal (Moon 2023 {stress})",
            "level": "RNA",
        }
        if len(mn['depleted']) >= 10:
            gene_sets[f"SG-excluded\n(Moon 2023 {stress})"] = {
                "genes": mn['depleted'],
                "citation": "Moon et al. Nat Commun 2023; PMID:37968266",
                "label_short": f"SG-excluded (Moon 2023 {stress})",
                "level": "RNA",
            }

    # === Analysis 1: Pre-ranked GSEA ===
    gsea_results = analysis1_gsea(ranking, gene_sets)

    # Save GSEA results
    gsea_df = pd.DataFrame([{
        "gene_set": r["label_short"],
        "level": r["level"],
        "citation": r["citation"],
        "n_matched": r["n_genes"],
        "n_total": r["n_total"],
        "ES": r["ES"],
        "NES": r["NES"],
        "pvalue": r["pvalue"],
        "FDR": r["adj_p"],
        "n_leading_edge": r["n_leading_edge"],
        "leading_edge": r["leading_edge"],
    } for r in gsea_results])
    gsea_path = os.path.join(BASEDIR, "sg_transcriptome_gsea_results.tsv")
    gsea_df.to_csv(gsea_path, sep="\t", index=False)
    print(f"\n  Results saved: {gsea_path}")

    # === Analysis 2: IR × SG overlap ===
    overlap_results = analysis2_overlap(ir_genes, gene_sets, all_expressed)

    overlap_df = pd.DataFrame([{
        "gene_set": r["gene_set"],
        "level": r["level"],
        "citation": r["citation"],
        "sg_in_expressed": r["sg_in_expressed"],
        "ir_in_expressed": r["ir_in_expressed"],
        "overlap": r["overlap"],
        "odds_ratio": r["odds_ratio"],
        "pvalue": r["pvalue"],
        "adj_p": r["adj_p"],
        "overlap_genes": r["overlap_genes"],
    } for r in overlap_results])
    overlap_path = os.path.join(BASEDIR, "sg_rna_ir_overlap.tsv")
    overlap_df.to_csv(overlap_path, sep="\t", index=False)
    print(f"\n  Results saved: {overlap_path}")

    # === Analysis 3: SG enrichment distribution ===
    khong_expr = analysis3_distribution(khong, ir_genes, all_expressed)

    # Save enrichment distribution data
    if khong_expr is not None:
        dist_df = khong_expr[['mouse_gene', 'log2(fold_change)', 'classification',
                              'is_IR', 'SG(FPKM)', 'Total(FPKM)']].copy()
        dist_path = os.path.join(BASEDIR, "sg_rna_enrichment_distribution.tsv")
        dist_df.to_csv(dist_path, sep="\t", index=False)
        print(f"\n  Results saved: {dist_path}")

    # === Analysis 4: Length confound ===
    binned_results = analysis4_length_confound(khong, ir_genes, all_ir, deseq2, all_expressed)

    # === Analysis 5: RNA vs Protein comparison ===
    analysis5_rna_vs_protein(gsea_results)

    # === Analysis 6: Per-dataset overlap (26 datasets) ===
    meta_table1_path = os.path.join(DATADIR, "meta2025_table1.xlsx")
    meta_table1 = None
    if os.path.exists(meta_table1_path):
        meta_table1 = load_meta2025_table1(meta_table1_path)
    per_ds_results = analysis6_per_dataset_overlap(ir_genes, meta_df, all_expressed, meta_table1)

    # Save per-dataset results
    per_ds_df = pd.DataFrame([{
        'dataset_id': r['dataset_id'],
        'reference': r['reference'],
        'stress': r['stress'],
        'system': r['system'],
        'technique': r['technique'],
        'n_sg_expressed': r['n_sg_expressed'],
        'n_ir': r['n_ir'],
        'overlap': r['overlap'],
        'odds_ratio': r['odds_ratio'],
        'log2_or': r['log2_or'],
        'pvalue': r['pvalue'],
        'adj_p': r['adj_p'],
    } for r in per_ds_results])
    per_ds_path = os.path.join(BASEDIR, "sg_per_dataset_overlap.tsv")
    per_ds_df.to_csv(per_ds_path, sep="\t", index=False)
    print(f"\n  Results saved: {per_ds_path}")

    # === Analysis 7: DEG × SG overlap ===
    deg_overlap_results, deg_up, deg_down = analysis7_deg_overlap(deseq2, gene_sets, all_expressed)

    # Save DEG overlap results
    deg_overlap_df = pd.DataFrame(deg_overlap_results)
    deg_overlap_path = os.path.join(BASEDIR, "sg_deg_overlap.tsv")
    deg_overlap_df.to_csv(deg_overlap_path, sep="\t", index=False)
    print(f"\n  Results saved: {deg_overlap_path}")

    # === Generate figures ===
    print("\n" + "=" * 60)
    print("Generating figures")
    print("=" * 60)

    plot_gsea_barplot(gsea_results,
                      os.path.join(FIGDIR, "sg_transcriptome_gsea_enrichment"))
    plot_overlap(overlap_results,
                 os.path.join(FIGDIR, "sg_rna_ir_overlap"))
    if khong_expr is not None:
        plot_distribution(khong_expr,
                          os.path.join(FIGDIR, "sg_rna_distribution"))
    plot_length_control(khong, ir_genes, all_ir, all_expressed,
                        os.path.join(FIGDIR, "sg_rna_length_control"))
    plot_per_dataset_forest(per_ds_results,
                            os.path.join(FIGDIR, "sg_per_dataset_forest"))

    # Venn diagrams (IR∩SG and DEG-UP∩SG)
    plot_venn_diagrams(ir_genes, deg_up, deg_down, gene_sets, all_expressed,
                       os.path.join(FIGDIR, "sg_venn_diagrams"))

    # === Summary ===
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("\n  Key question: Do IR mRNAs go to SGs or stay in the nucleus?")
    print()
    print("  GSEA results (RNA level):")
    for r in gsea_results:
        if r['level'] == 'RNA':
            sig = "***" if r["adj_p"] < 0.001 else "**" if r["adj_p"] < 0.01 else "*" if r["adj_p"] < 0.05 else "ns"
            direction = "UP in KD" if r['NES'] > 0 else "DOWN in KD"
            print(f"    {r['label_short']}: NES={r['NES']:.2f}, FDR={r['adj_p']:.4f} {sig} ({direction})")

    print("\n  Overlap (IR × SG RNA):")
    for r in overlap_results:
        if r['level'] == 'RNA':
            sig = "***" if r["adj_p"] < 0.001 else "**" if r["adj_p"] < 0.01 else "*" if r["adj_p"] < 0.05 else "ns"
            print(f"    {r['gene_set']}: OR={r['odds_ratio']:.2f}, FDR={r['adj_p']:.4f} {sig} (n={r['overlap']})")

    print("\n  Per-dataset reproducibility (26 datasets):")
    n_ds_sig = sum(1 for r in per_ds_results if r['adj_p'] < 0.05 and r['odds_ratio'] > 1)
    n_ds_total = len(per_ds_results)
    median_or = np.median([r['odds_ratio'] for r in per_ds_results])
    print(f"    Significant positive overlap: {n_ds_sig}/{n_ds_total} datasets (FDR<0.05)")
    print(f"    Median OR across all datasets: {median_or:.2f}")

    print("\n  Interpretation:")
    rna_gsea = [r for r in gsea_results if r['level'] == 'RNA']
    sig_rna = [r for r in rna_gsea if r['adj_p'] < 0.05]
    rna_overlap = [r for r in overlap_results if r['level'] == 'RNA']
    sig_overlap = [r for r in rna_overlap if r['adj_p'] < 0.05]

    if not sig_rna and not sig_overlap:
        print("    → IR mRNAs show NO significant association with SG transcriptome")
        print("    → Consistent with nuclear detention model")
    elif sig_overlap and all(r['odds_ratio'] > 1 for r in sig_overlap):
        print("    → IR genes ARE enriched in normal SG transcriptome")
        print("    → Nuclear detention of these mRNAs may deplete SGs of their usual RNA content")
    else:
        print("    → Mixed results - see detailed output above")

    print(f"\n{'='*70}")
    print("Analysis complete.")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
