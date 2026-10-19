#!/usr/bin/env python3
"""
Step 2: Per-intron IR ratio calculation + retained intron identification

Reads featureCounts output, computes IR ratio per intron using flanking exon density,
performs Welch's t-test + BH correction, and identifies retained introns.
Applies the same pipeline to all 3 datasets (DHX8, GSE27218, GSE65973).

Usage:
    python3 step2_identify_retained_introns.py

Runs on: ncnp-linux
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# === Configuration ===
WORKDIR = Path("/mnt/HDD2/TDP43_comparison/event_level_20260225")
INTRON_META = WORKDIR / "intron_metadata.tsv"

# Dataset configurations: (name, intron_counts, exon_counts, sample_info)
DATASETS = {
    "dhx8": {
        "intron_counts": WORKDIR / "dhx8_intron_counts.txt",
        "exon_counts": WORKDIR / "dhx8_exon_counts.txt",
        # Columns in featureCounts output correspond to BAM files in order
        # siCont_1, siCont_2, siCont_3, siDHX8_1, siDHX8_2, siDHX8_3
        "control_pattern": "siCont",
        "kd_pattern": "siDHX8",
    },
    "gse27218": {
        "intron_counts": WORKDIR / "gse27218_intron_counts.txt",
        "exon_counts": WORKDIR / "gse27218_exon_counts.txt",
        # SRR107061-SRR107071 = control (11), SRR107072-SRR107081 = tdp43_kd (10)
        "control_samples": [f"SRR10706{i}" for i in range(1, 10)] + ["SRR107070", "SRR107071"],
        "kd_samples": [f"SRR10707{i}" for i in range(2, 10)] + ["SRR107080", "SRR107081"],
    },
    "gse65973": {
        "intron_counts": WORKDIR / "gse65973_intron_counts.txt",
        "exon_counts": WORKDIR / "gse65973_exon_counts.txt",
        # SRR1805871-SRR1805874 = control (4), SRR1805867-SRR1805870 = bigenic (4)
        "control_samples": [f"SRR180587{i}" for i in range(1, 5)],
        "kd_samples": [f"SRR180586{i}" for i in range(7, 10)] + ["SRR1805870"],
    },
}

# Filter thresholds
MIN_FLANKING_EXON_COUNT = 10  # Minimum mean flanking exon count in at least one group
ADJ_P_THRESHOLD = 0.05


def load_featurecounts(filepath):
    """
    Load featureCounts output file.
    Returns DataFrame with GeneID as index and sample columns.
    featureCounts output has: GeneID, Chr, Start, End, Strand, Length, then count columns.
    First line is a comment starting with #, second line is header.
    """
    df = pd.read_csv(filepath, sep="\t", comment="#")
    df = df.set_index("Geneid")
    # Drop annotation columns, keep only count columns
    count_cols = [c for c in df.columns if c not in ("Chr", "Start", "End", "Strand", "Length")]
    lengths = df["Length"].copy()
    counts = df[count_cols].copy()
    # Clean column names: featureCounts uses full BAM paths
    counts.columns = [Path(c).stem.replace("_Aligned.sortedByCoord.out", "") for c in counts.columns]
    return counts, lengths


def classify_samples(counts, dataset_config):
    """Classify sample columns into control and KD groups."""
    cols = counts.columns.tolist()
    if "control_pattern" in dataset_config:
        ctrl_cols = [c for c in cols if dataset_config["control_pattern"] in c]
        kd_cols = [c for c in cols if dataset_config["kd_pattern"] in c]
    else:
        ctrl_samples = dataset_config["control_samples"]
        kd_samples = dataset_config["kd_samples"]
        ctrl_cols = [c for c in cols if any(s in c for s in ctrl_samples)]
        kd_cols = [c for c in cols if any(s in c for s in kd_samples)]
    return ctrl_cols, kd_cols


def load_intron_metadata(filepath):
    """Load intron metadata with flanking exon info."""
    return pd.read_csv(filepath, sep="\t")


def compute_ir_ratios(intron_counts, intron_lengths, exon_counts, exon_lengths,
                      metadata, ctrl_cols, kd_cols):
    """
    Compute IR ratio for each intron in each sample.

    IR ratio = (intron_reads / intron_length) / (mean_flanking_exon_reads / mean_flanking_exon_length)

    Returns DataFrame: intron_id x samples with IR ratios.
    """
    # Build mapping: intron_id -> list of flanking exon_ids
    intron_to_exons = {}
    for _, row in metadata.iterrows():
        iid = row["intron_id"]
        if row["flanking_exon_ids"] != "NA" and pd.notna(row["flanking_exon_ids"]):
            intron_to_exons[iid] = row["flanking_exon_ids"].split(",")
        else:
            intron_to_exons[iid] = []

    all_samples = ctrl_cols + kd_cols
    ir_ratios = pd.DataFrame(index=intron_counts.index, columns=all_samples, dtype=float)
    flanking_mean_counts = pd.DataFrame(index=intron_counts.index, columns=all_samples, dtype=float)

    for intron_id in intron_counts.index:
        if intron_id not in intron_to_exons:
            continue

        exon_ids = intron_to_exons.get(intron_id, [])
        intron_len = intron_lengths.loc[intron_id]

        if len(exon_ids) > 0:
            # Get flanking exon counts and lengths
            valid_exons = [eid for eid in exon_ids if eid in exon_counts.index]
            if len(valid_exons) == 0:
                continue

            for sample in all_samples:
                if sample not in intron_counts.columns:
                    continue
                intron_read = intron_counts.loc[intron_id, sample]
                intron_density = intron_read / intron_len if intron_len > 0 else 0

                exon_densities = []
                exon_raw_counts = []
                for eid in valid_exons:
                    if sample in exon_counts.columns:
                        ec = exon_counts.loc[eid, sample]
                        el = exon_lengths.loc[eid]
                        exon_raw_counts.append(ec)
                        if el > 0:
                            exon_densities.append(ec / el)

                if len(exon_densities) > 0:
                    mean_exon_density = np.mean(exon_densities)
                    flanking_mean_counts.loc[intron_id, sample] = np.mean(exon_raw_counts)
                    if mean_exon_density > 0:
                        ir_ratios.loc[intron_id, sample] = intron_density / mean_exon_density
                    else:
                        ir_ratios.loc[intron_id, sample] = np.nan
                else:
                    ir_ratios.loc[intron_id, sample] = np.nan
        else:
            # No flanking exons available
            for sample in all_samples:
                ir_ratios.loc[intron_id, sample] = np.nan

    return ir_ratios, flanking_mean_counts


def identify_retained_introns(ir_ratios, flanking_mean_counts, ctrl_cols, kd_cols):
    """
    Identify retained introns using Welch's t-test + BH correction.

    Filters:
    - At least one group has mean flanking exon count > MIN_FLANKING_EXON_COUNT
    - Welch's t-test on IR ratios (KD vs Control)
    - BH correction
    - Retained: delta_IR > 0 AND adj_p < 0.05
    """
    results = []

    for intron_id in ir_ratios.index:
        ctrl_ir = ir_ratios.loc[intron_id, ctrl_cols].astype(float).dropna()
        kd_ir = ir_ratios.loc[intron_id, kd_cols].astype(float).dropna()

        if len(ctrl_ir) < 2 or len(kd_ir) < 2:
            continue

        # Flanking exon count filter
        ctrl_flank = flanking_mean_counts.loc[intron_id, ctrl_cols].astype(float).dropna()
        kd_flank = flanking_mean_counts.loc[intron_id, kd_cols].astype(float).dropna()
        max_mean_flank = max(
            ctrl_flank.mean() if len(ctrl_flank) > 0 else 0,
            kd_flank.mean() if len(kd_flank) > 0 else 0
        )
        if max_mean_flank < MIN_FLANKING_EXON_COUNT:
            continue

        # Welch's t-test
        t_stat, p_value = stats.ttest_ind(kd_ir.values, ctrl_ir.values, equal_var=False)

        delta_ir = kd_ir.mean() - ctrl_ir.mean()

        results.append({
            "intron_id": intron_id,
            "ir_ratio_ctrl": ctrl_ir.mean(),
            "ir_ratio_kd": kd_ir.mean(),
            "delta_ir": delta_ir,
            "t_stat": t_stat,
            "pvalue": p_value,
            "n_ctrl": len(ctrl_ir),
            "n_kd": len(kd_ir),
            "mean_flanking_count": max_mean_flank,
        })

    df = pd.DataFrame(results)
    if len(df) == 0:
        return df

    # Drop rows with NaN p-values (from t-test with zero variance)
    n_before = len(df)
    df = df.dropna(subset=["pvalue"]).reset_index(drop=True)
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(f"    Dropped {n_dropped} introns with NaN p-values", file=sys.stderr)

    if len(df) == 0:
        return df

    # Manual BH correction
    pvals = df["pvalue"].values
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    adj_pvals = np.empty(n)
    adj_pvals[sorted_idx[-1]] = sorted_pvals[-1]
    for i in range(n - 2, -1, -1):
        adj_pvals[sorted_idx[i]] = min(
            adj_pvals[sorted_idx[i + 1]],
            sorted_pvals[i] * n / (i + 1)
        )
    adj_pvals = np.clip(adj_pvals, 0, 1)
    df["adj_p"] = adj_pvals

    # Parse intron_id into components
    df["chr"] = df["intron_id"].apply(lambda x: x.split(":")[0])
    df["coord"] = df["intron_id"].apply(lambda x: x.split(":")[1])
    df["start"] = df["coord"].apply(lambda x: int(x.split("-")[0]))
    df["end"] = df["coord"].apply(lambda x: int(x.split("-")[1]))
    df["strand"] = df["intron_id"].apply(lambda x: x.split(":")[2])
    df["gene_name"] = df["intron_id"].apply(lambda x: x.split(":")[3])
    df["intron_length"] = df["end"] - df["start"] + 1

    # Mark retained: use nominal p for primary criterion (robust for small n)
    # adj_p kept as reference for well-powered datasets
    df["retained_strict"] = (df["delta_ir"] > 0) & (df["adj_p"] < ADJ_P_THRESHOLD)
    df["retained"] = (df["delta_ir"] > 0) & (df["pvalue"] < ADJ_P_THRESHOLD)

    return df


def process_dataset(name, config, metadata):
    """Process a single dataset end-to-end."""
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Processing dataset: {name}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)

    # Load counts
    print(f"  Loading intron counts from {config['intron_counts']}...", file=sys.stderr)
    intron_counts, intron_lengths = load_featurecounts(config["intron_counts"])
    print(f"    Introns: {len(intron_counts)}, Samples: {len(intron_counts.columns)}", file=sys.stderr)
    print(f"    Sample names: {list(intron_counts.columns)}", file=sys.stderr)

    print(f"  Loading exon counts from {config['exon_counts']}...", file=sys.stderr)
    exon_counts, exon_lengths = load_featurecounts(config["exon_counts"])
    print(f"    Exons: {len(exon_counts)}, Samples: {len(exon_counts.columns)}", file=sys.stderr)

    # Classify samples
    ctrl_cols, kd_cols = classify_samples(intron_counts, config)
    print(f"  Control samples ({len(ctrl_cols)}): {ctrl_cols}", file=sys.stderr)
    print(f"  KD samples ({len(kd_cols)}): {kd_cols}", file=sys.stderr)

    if len(ctrl_cols) == 0 or len(kd_cols) == 0:
        print(f"  ERROR: Could not classify samples! Skipping.", file=sys.stderr)
        return None

    # Compute IR ratios
    print(f"  Computing IR ratios...", file=sys.stderr)
    ir_ratios, flanking_mean_counts = compute_ir_ratios(
        intron_counts, intron_lengths, exon_counts, exon_lengths,
        metadata, ctrl_cols, kd_cols
    )

    # Count valid ratios
    n_valid = ir_ratios.notna().any(axis=1).sum()
    print(f"    Introns with valid IR ratios: {n_valid}", file=sys.stderr)

    # Identify retained introns
    print(f"  Running statistical tests...", file=sys.stderr)
    results = identify_retained_introns(ir_ratios, flanking_mean_counts, ctrl_cols, kd_cols)

    if len(results) == 0:
        print(f"  WARNING: No introns passed filters!", file=sys.stderr)
        return results

    n_tested = len(results)
    n_retained = results["retained"].sum()
    n_retained_strict = results["retained_strict"].sum()
    n_sig_adj = (results["adj_p"] < ADJ_P_THRESHOLD).sum()
    n_sig_nom = (results["pvalue"] < ADJ_P_THRESHOLD).sum()

    print(f"  Results:", file=sys.stderr)
    print(f"    Introns tested: {n_tested}", file=sys.stderr)
    print(f"    Significant (nominal p < {ADJ_P_THRESHOLD}): {n_sig_nom}", file=sys.stderr)
    print(f"    Significant (adj_p < {ADJ_P_THRESHOLD}): {n_sig_adj}", file=sys.stderr)
    print(f"    Retained (nominal p, delta_IR > 0): {n_retained}", file=sys.stderr)
    print(f"    Retained strict (adj_p, delta_IR > 0): {n_retained_strict}", file=sys.stderr)

    # Save full results
    out_all = WORKDIR / f"ir_results_all_{name}.tsv"
    results.to_csv(out_all, sep="\t", index=False)
    print(f"    Full results: {out_all}", file=sys.stderr)

    # Save retained introns only
    retained = results[results["retained"]].copy()
    out_retained = WORKDIR / f"retained_introns_{name}.tsv"
    retained.to_csv(out_retained, sep="\t", index=False)
    print(f"    Retained introns: {out_retained}", file=sys.stderr)

    return results


def main():
    # Load intron metadata
    print("Loading intron metadata...", file=sys.stderr)
    metadata = load_intron_metadata(INTRON_META)
    print(f"  Introns in metadata: {len(metadata)}", file=sys.stderr)

    # Process all datasets
    all_results = {}
    for name, config in DATASETS.items():
        results = process_dataset(name, config, metadata)
        if results is not None:
            all_results[name] = results

    # Summary
    print(f"\n{'='*60}", file=sys.stderr)
    print("Summary across datasets:", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    for name, results in all_results.items():
        n_retained = results["retained"].sum()
        n_genes = results[results["retained"]]["gene_name"].nunique()
        print(f"  {name}: {n_retained} retained introns in {n_genes} genes", file=sys.stderr)

    # Sanity check for DHX8
    if "dhx8" in all_results:
        dhx8_genes = set(all_results["dhx8"][all_results["dhx8"]["retained"]]["gene_name"])
        print(f"\n  DHX8 sanity check: {len(dhx8_genes)} genes with retained introns", file=sys.stderr)
        print(f"  (Compare with 4,552 transcripts from previous analysis)", file=sys.stderr)
        # Check CDCA5
        if "Cdca5" in dhx8_genes:
            print(f"  Positive control: Cdca5 DETECTED among retained introns", file=sys.stderr)
        else:
            print(f"  WARNING: Cdca5 NOT detected among retained introns", file=sys.stderr)

    print(f"\n=== Step 2 complete ===", file=sys.stderr)


if __name__ == "__main__":
    main()
