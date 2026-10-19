#!/usr/bin/env python3
"""MAJIQ deltapsi integration analysis across 4 datasets."""

import pandas as pd
import numpy as np
from pathlib import Path
import json

BASE = Path("/mnt/HDD2/unified_results/majiq")

DATASETS = {
    "dhx8": {"dir": "dhx8", "grp1": "Control", "grp2": "DHX8_KD", "species": "mouse", "desc": "DHX8 KD (Sakai)"},
    "gse126543": {"dir": "gse126543", "grp1": "TDP43_negative", "grp2": "TDP43_positive", "species": "human", "desc": "TDP-43 pathology (Brown 2022, iPSC MN)"},
    "gse27218": {"dir": "gse27218", "grp1": "Control", "grp2": "TDP43_KD", "species": "human", "desc": "TDP-43 KD (Tollervey 2011, HEK293)"},
    "gse65973": {"dir": "gse65973", "grp1": "Control", "grp2": "dNLS_TDP", "species": "mouse", "desc": "TDP-43 dNLS (Fratta 2018, mouse brain)"},
}

# Significance thresholds
PROB_CHANGING_THRESH = 0.95  # P(|dPSI| > 0.2)
DPSI_THRESH = 0.1  # |dPSI| minimum

def load_deltapsi(dataset_key):
    info = DATASETS[dataset_key]
    f = BASE / info["dir"] / "deltapsi_annotated.tsv"
    df = pd.read_csv(f, sep="\t", comment="#")
    df["dataset"] = dataset_key
    return df

print("=" * 70)
print("MAJIQ DeltaPsi Integration Analysis")
print("=" * 70)

# Load all datasets
all_dfs = {}
for key in DATASETS:
    df = load_deltapsi(key)
    all_dfs[key] = df
    print(f"\n{DATASETS[key]['desc']} ({key}):")
    print(f"  Total events: {len(df):,}")
    
    # Significant events
    sig = df[(df["probability_changing"] >= PROB_CHANGING_THRESH) & (df["dpsi_mean"].abs() >= DPSI_THRESH)]
    sig_genes = sig["gene_name"].nunique()
    print(f"  Significant (P_changing≥{PROB_CHANGING_THRESH}, |dPSI|≥{DPSI_THRESH}): {len(sig):,} events in {sig_genes} genes")
    
    # By event type
    if len(sig) > 0:
        print(f"  Event types: {sig['event_type'].value_counts().to_dict()}")
        
        # De novo (unannotated) events
        denovo = sig[sig["event_denovo"] == True]
        print(f"  De novo (unannotated) events: {len(denovo)} ({100*len(denovo)/len(sig):.1f}%)")
        
        # Intron retention
        ir = sig[sig["is_intron"] == True]
        print(f"  Intron retention events: {len(ir)}")
        
        # Top 10 genes by |dpsi|
        top = sig.nlargest(10, "dpsi_mean", keep="first")[["gene_name", "event_type", "dpsi_mean", "probability_changing"]].copy()
        top["dpsi_mean"] = top["dpsi_mean"].round(3)
        print(f"\n  Top 10 inclusion events (dPSI > 0):")
        for _, r in top.iterrows():
            print(f"    {r['gene_name']:15s} {r['event_type']:3s} dPSI={r['dpsi_mean']:+.3f} P_change={r['probability_changing']:.4f}")

# ================================================================
# Cross-dataset comparison for DHX8 vs TDP-43 datasets
# ================================================================
print("\n" + "=" * 70)
print("Cross-Dataset Gene Overlap Analysis")
print("=" * 70)

sig_genes_per_dataset = {}
for key in DATASETS:
    df = all_dfs[key]
    sig = df[(df["probability_changing"] >= PROB_CHANGING_THRESH) & (df["dpsi_mean"].abs() >= DPSI_THRESH)]
    genes = set(sig["gene_name"].dropna().unique())
    sig_genes_per_dataset[key] = genes
    print(f"  {key}: {len(genes)} significant genes")

# Pairwise overlaps
print("\nPairwise gene overlaps:")
keys = list(DATASETS.keys())
for i in range(len(keys)):
    for j in range(i+1, len(keys)):
        k1, k2 = keys[i], keys[j]
        # For cross-species, use gene name matching (uppercase)
        g1 = {g.upper() for g in sig_genes_per_dataset[k1]}
        g2 = {g.upper() for g in sig_genes_per_dataset[k2]}
        overlap = g1 & g2
        print(f"  {k1} ∩ {k2}: {len(overlap)} genes")
        if len(overlap) > 0 and len(overlap) <= 30:
            print(f"    {sorted(overlap)}")

# ================================================================
# DHX8-specific vs TDP-43 comparison
# ================================================================
print("\n" + "=" * 70)
print("DHX8-specific splicing events (not in TDP-43 datasets)")
print("=" * 70)

dhx8_genes = {g.upper() for g in sig_genes_per_dataset["dhx8"]}
tdp43_genes = set()
for key in ["gse126543", "gse27218", "gse65973"]:
    tdp43_genes |= {g.upper() for g in sig_genes_per_dataset[key]}

dhx8_only = dhx8_genes - tdp43_genes
shared = dhx8_genes & tdp43_genes
tdp43_only = tdp43_genes - dhx8_genes

print(f"  DHX8-only genes: {len(dhx8_only)}")
print(f"  Shared with TDP-43: {len(shared)}")
print(f"  TDP-43-only genes: {len(tdp43_only)}")

if shared:
    print(f"\n  Shared genes: {sorted(shared)[:50]}")

# ================================================================
# Known TDP-43 targets check
# ================================================================
print("\n" + "=" * 70)
print("Known TDP-43 target gene check across datasets")
print("=" * 70)

KNOWN_TDP43_TARGETS = [
    "STMN2", "UNC13A", "KALRN", "AGRN", "SLC39A11", "KCNQ2",
    "PFKP", "UNC13B", "SORT1", "CAMK2B", "RAPGEF6", "EP400",
    "Stmn2", "Unc13a", "Kalrn", "Agrn", "Pfkp",  # mouse orthologs
]

for target in ["STMN2", "UNC13A", "KALRN", "AGRN", "SORT1", "PFKP", "CAMK2B"]:
    print(f"\n  {target}:")
    for key in DATASETS:
        df = all_dfs[key]
        # case-insensitive match
        target_events = df[df["gene_name"].str.upper() == target.upper()]
        if len(target_events) > 0:
            sig_events = target_events[
                (target_events["probability_changing"] >= PROB_CHANGING_THRESH) & 
                (target_events["dpsi_mean"].abs() >= DPSI_THRESH)
            ]
            max_dpsi = target_events["dpsi_mean"].abs().max()
            max_prob = target_events["probability_changing"].max()
            print(f"    {key:12s}: {len(target_events):3d} events, {len(sig_events)} sig, max|dPSI|={max_dpsi:.3f}, maxP={max_prob:.4f}")
        else:
            print(f"    {key:12s}: not found")

# ================================================================
# DHX8 top events detail
# ================================================================
print("\n" + "=" * 70)
print("DHX8 KD: Detailed significant events summary")
print("=" * 70)

dhx8_df = all_dfs["dhx8"]
dhx8_sig = dhx8_df[
    (dhx8_df["probability_changing"] >= PROB_CHANGING_THRESH) & 
    (dhx8_df["dpsi_mean"].abs() >= DPSI_THRESH)
].copy()

# Group by gene, summarize
gene_summary = dhx8_sig.groupby("gene_name").agg(
    n_events=("dpsi_mean", "size"),
    max_abs_dpsi=("dpsi_mean", lambda x: x.abs().max()),
    mean_dpsi=("dpsi_mean", "mean"),
    n_ir=("is_intron", "sum"),
    n_denovo=("event_denovo", "sum"),
    event_types=("event_type", lambda x: ",".join(sorted(x.unique()))),
).sort_values("max_abs_dpsi", ascending=False)

print(f"\nTop 30 DHX8-affected genes by |dPSI|:")
print(f"{'Gene':<15s} {'#Events':>7s} {'max|dPSI|':>10s} {'mean_dPSI':>10s} {'#IR':>4s} {'#denovo':>7s} {'types'}")
for gene, row in gene_summary.head(30).iterrows():
    print(f"{gene:<15s} {row['n_events']:>7d} {row['max_abs_dpsi']:>10.3f} {row['mean_dpsi']:>10.3f} {int(row['n_ir']):>4d} {int(row['n_denovo']):>7d} {row['event_types']}")

# Save summary tables
gene_summary.to_csv(BASE / "dhx8_deltapsi_gene_summary.tsv", sep="\t")
print(f"\nSaved: dhx8_deltapsi_gene_summary.tsv ({len(gene_summary)} genes)")

# All datasets significant genes
all_sig = []
for key in DATASETS:
    df = all_dfs[key]
    sig = df[(df["probability_changing"] >= PROB_CHANGING_THRESH) & (df["dpsi_mean"].abs() >= DPSI_THRESH)].copy()
    sig["dataset"] = key
    sig["species"] = DATASETS[key]["species"]
    all_sig.append(sig)
all_sig_df = pd.concat(all_sig, ignore_index=True)
all_sig_df.to_csv(BASE / "all_datasets_significant_events.tsv", sep="\t", index=False)
print(f"Saved: all_datasets_significant_events.tsv ({len(all_sig_df)} events)")

