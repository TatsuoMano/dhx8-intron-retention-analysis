#!/usr/bin/env python3
"""Fig. 5B statistics with the 95% confidence interval of Cohen's d (written for the code release,
2026-10-20).

Recomputes the Fig. 5B comparison exactly as 07_figure_panels/server/_fig5_rnaquant_v7_v2.py does
(same inputs, gene-level mean DESeq2 log2FC, affected vs non-affected genes, Welch t-test with
Bonferroni x3, d = mean difference / sqrt((var_a + var_b) / 2)) and adds two 95% CIs for d:

  analytic   d +/- 1.96 * SE, SE = sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
             (large-sample normal approximation, Hedges and Olkin 1985)
  bootstrap  percentile interval of d over 10,000 resamples of genes within each group (seed 0)

Run on the server (read-only, streamed):
    ssh ncnp-linux 'python3 -' < 02_splicing/ir_rna_level_effect_size.py
"""
import numpy as np
import pandas as pd
from scipy import stats

B = "/mnt/HDD2/unified_results"
deseq = pd.read_csv(f"{B}/deseq2_dhx8_kd_vs_ctrl.tsv", sep="\t").dropna(subset=["log2FoldChange", "gene"])
lfc = deseq.groupby("gene")["log2FoldChange"].mean()


def rmats_genes(fn):
    df = pd.read_csv(f"{B}/rmats/dhx8_v2/{fn}", sep="\t")
    sig = df[df["FDR"] < 0.05]
    return set(sig["geneSymbol"].astype(str).str.strip("\"'"))


ir = set(pd.read_csv("/mnt/HDD2/TDP43_comparison/validation_20260226/limma_retained_dhx8_s2.tsv",
                     sep="\t")["gene_name"].dropna().unique())
se = rmats_genes("SE.MATS.JCEC.txt")
ass = rmats_genes("A3SS.MATS.JCEC.txt") | rmats_genes("A5SS.MATS.JCEC.txt")


def cohend(a, b):
    return (np.mean(a) - np.mean(b)) / np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)


def boot_ci(a, b, n_boot=10000, seed=0, chunk=500):
    rng = np.random.default_rng(seed)
    out = []
    for start in range(0, n_boot, chunk):
        m = min(chunk, n_boot - start)
        ra = a[rng.integers(0, len(a), size=(m, len(a)))]
        rb = b[rng.integers(0, len(b), size=(m, len(b)))]
        d = (ra.mean(1) - rb.mean(1)) / np.sqrt((ra.var(1, ddof=1) + rb.var(1, ddof=1)) / 2)
        out.append(d)
    out = np.concatenate(out)
    return np.percentile(out, [2.5, 97.5])


print(f"{'group':<10} {'n_aff':>6} {'n_non':>6} {'mean':>8} {'d':>8} {'P_bonf':>10}  analytic 95% CI   bootstrap 95% CI")
for name, gset in [("IR", ir), ("SE", se), ("A3SS/A5SS", ass)]:
    aff = lfc[lfc.index.isin(gset)].values
    non = lfc[~lfc.index.isin(gset)].values
    d = cohend(aff, non)
    p = min(1.0, stats.ttest_ind(aff, non, equal_var=False)[1] * 3)
    n1, n2 = len(aff), len(non)
    se_d = np.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))
    lo, hi = d - 1.96 * se_d, d + 1.96 * se_d
    blo, bhi = boot_ci(aff, non)
    print(f"{name:<10} {n1:>6} {n2:>6} {np.mean(aff):>+8.3f} {d:>+8.3f} {p:>10.2e}  {lo:.3f}-{hi:.3f}      {blo:.3f}-{bhi:.3f}")
