#!/usr/bin/env python3
"""
Step 1: Non-redundant intron SAF creation + flanking exon SAF

Parses mm10.ncbiRefSeq.introns.gtf, deduplicates by coordinate (chr:start-end:strand),
outputs SAF format for featureCounts per-intron counting.
Also extracts flanking exons from the exon GTF for IR ratio normalization.

Usage:
    python3 step1_per_intron_count.py

Runs on: ncnp-linux
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

# === Configuration ===
INTRON_GTF = "/mnt/HDD2/Sakai_RNA-seq/output/mm10.ncbiRefSeq.introns.gtf"
EXON_GTF = "/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf"
OUTDIR = Path("/mnt/HDD2/TDP43_comparison/event_level_20260225")

INTRON_SAF = OUTDIR / "introns.saf"
FLANKING_EXON_SAF = OUTDIR / "flanking_exons.saf"
INTRON_META = OUTDIR / "intron_metadata.tsv"


def parse_gtf_attribute(attr_str, key):
    """Extract value for a given key from GTF attribute string."""
    m = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return m.group(1) if m else None


def load_introns(gtf_path):
    """Load all introns from GTF, return list of (chr, start, end, strand, gene_name, transcript_id)."""
    introns = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]
            gene_name = parse_gtf_attribute(attrs, "gene_name") or "unknown"
            transcript_id = parse_gtf_attribute(attrs, "transcript_id") or "unknown"
            introns.append((chrom, start, end, strand, gene_name, transcript_id))
    return introns


def deduplicate_introns(introns):
    """
    Deduplicate introns by coordinate (chr:start-end:strand).
    For each unique coordinate, pick one representative gene_name
    and collect all transcript_ids.
    """
    coord_to_info = {}  # (chr, start, end, strand) -> (gene_name, set of transcript_ids)
    for chrom, start, end, strand, gene_name, tid in introns:
        key = (chrom, start, end, strand)
        if key not in coord_to_info:
            coord_to_info[key] = (gene_name, {tid})
        else:
            existing_gene, tids = coord_to_info[key]
            tids.add(tid)
            # Keep gene_name consistent (prefer non-unknown)
            if existing_gene == "unknown" and gene_name != "unknown":
                coord_to_info[key] = (gene_name, tids)
    return coord_to_info


def load_exons(gtf_path):
    """
    Load all exons from GTF.
    Returns dict: transcript_id -> sorted list of (chr, start, end, strand, gene_name)
    """
    tx_exons = defaultdict(list)
    print(f"Loading exons from {gtf_path}...", file=sys.stderr)
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]
            tid = parse_gtf_attribute(attrs, "transcript_id")
            gene_name = parse_gtf_attribute(attrs, "gene_name") or "unknown"
            if tid:
                tx_exons[tid].append((chrom, start, end, strand, gene_name))
    # Sort exons by start position within each transcript
    for tid in tx_exons:
        tx_exons[tid].sort(key=lambda x: x[1])
    return tx_exons


def find_flanking_exons(intron_coords, tx_exons, introns_raw):
    """
    For each non-redundant intron, find the flanking (upstream + downstream) exons.
    Returns dict: intron_coord -> list of (chr, start, end, strand, exon_label)
    """
    # Build mapping: intron_coord -> set of transcript_ids
    coord_to_tids = defaultdict(set)
    for chrom, start, end, strand, gene_name, tid in introns_raw:
        coord_to_tids[(chrom, start, end, strand)].add(tid)

    # For each transcript, build intron -> flanking exon mapping
    # Intron between exon[i] and exon[i+1]: flanking = exon[i] and exon[i+1]
    intron_flanking = {}  # intron_coord -> set of (chr, start, end, strand)

    for tid, exons in tx_exons.items():
        if len(exons) < 2:
            continue
        for i in range(len(exons) - 1):
            e1 = exons[i]
            e2 = exons[i + 1]
            # Intron is between e1.end+1 and e2.start-1
            intron_start = e1[2] + 1  # e1.end + 1
            intron_end = e2[1] - 1    # e2.start - 1
            if intron_end <= intron_start:
                continue
            chrom = e1[0]
            strand = e1[3]
            ikey = (chrom, intron_start, intron_end, strand)

            if ikey in intron_coords:
                if ikey not in intron_flanking:
                    intron_flanking[ikey] = set()
                # Add upstream exon
                intron_flanking[ikey].add((chrom, e1[1], e1[2], strand))
                # Add downstream exon
                intron_flanking[ikey].add((chrom, e2[1], e2[2], strand))

    return intron_flanking


def write_intron_saf(intron_coords, outpath):
    """Write non-redundant introns in SAF format."""
    with open(outpath, "w") as f:
        f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for (chrom, start, end, strand), (gene_name, tids) in sorted(intron_coords.items()):
            intron_id = f"{chrom}:{start}-{end}:{strand}:{gene_name}"
            f.write(f"{intron_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
    print(f"Wrote {len(intron_coords)} non-redundant introns to {outpath}", file=sys.stderr)


def write_flanking_exon_saf(intron_coords, intron_flanking, outpath):
    """
    Write flanking exons in SAF format.
    Each flanking exon is labeled with its parent intron_id for later matching.
    We create combined "flanking region" entries per intron (upstream + downstream merged).
    To keep it simple for featureCounts, we write individual exon entries tagged by intron.
    """
    exon_entries = {}  # (chr, start, end, strand) -> set of intron_ids
    for icoord in intron_coords:
        gene_name = intron_coords[icoord][0]
        intron_id = f"{icoord[0]}:{icoord[1]}-{icoord[2]}:{icoord[3]}:{gene_name}"
        if icoord in intron_flanking:
            for echrom, estart, eend, estrand in intron_flanking[icoord]:
                ekey = (echrom, estart, eend, estrand)
                if ekey not in exon_entries:
                    exon_entries[ekey] = set()
                exon_entries[ekey].add(intron_id)

    # For featureCounts, we need unique GeneIDs.
    # Strategy: use exon coordinate as GeneID, then map back to introns via metadata.
    with open(outpath, "w") as f:
        f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for (chrom, start, end, strand) in sorted(exon_entries.keys()):
            exon_id = f"{chrom}:{start}-{end}:{strand}"
            f.write(f"{exon_id}\t{chrom}\t{start}\t{end}\t{strand}\n")

    print(f"Wrote {len(exon_entries)} unique flanking exon regions to {outpath}", file=sys.stderr)
    return exon_entries


def write_intron_metadata(intron_coords, intron_flanking, outpath):
    """
    Write metadata linking each intron to its flanking exons.
    Columns: intron_id, chr, start, end, strand, gene_name, intron_length,
             upstream_exon_id, downstream_exon_id
    """
    with open(outpath, "w") as f:
        f.write("intron_id\tchr\tstart\tend\tstrand\tgene_name\tintron_length\tflanking_exon_ids\n")
        for (chrom, start, end, strand), (gene_name, tids) in sorted(intron_coords.items()):
            intron_id = f"{chrom}:{start}-{end}:{strand}:{gene_name}"
            length = end - start + 1
            if (chrom, start, end, strand) in intron_flanking:
                exon_ids = []
                for echrom, estart, eend, estrand in sorted(intron_flanking[(chrom, start, end, strand)]):
                    exon_ids.append(f"{echrom}:{estart}-{eend}:{estrand}")
                flanking = ",".join(exon_ids)
            else:
                flanking = "NA"
            f.write(f"{intron_id}\t{chrom}\t{start}\t{end}\t{strand}\t{gene_name}\t{length}\t{flanking}\n")
    print(f"Wrote intron metadata to {outpath}", file=sys.stderr)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Load and deduplicate introns
    print("Loading intron GTF...", file=sys.stderr)
    introns_raw = load_introns(INTRON_GTF)
    print(f"  Total intron entries: {len(introns_raw)}", file=sys.stderr)

    intron_coords = deduplicate_introns(introns_raw)
    print(f"  Non-redundant introns: {len(intron_coords)}", file=sys.stderr)

    # Write intron SAF
    write_intron_saf(intron_coords, INTRON_SAF)

    # Load exons and find flanking exons
    tx_exons = load_exons(EXON_GTF)
    print(f"  Transcripts with exons: {len(tx_exons)}", file=sys.stderr)

    intron_flanking = find_flanking_exons(intron_coords, tx_exons, introns_raw)
    n_with_flanking = sum(1 for v in intron_flanking.values() if v)
    print(f"  Introns with flanking exons: {n_with_flanking}", file=sys.stderr)

    # Write flanking exon SAF
    write_flanking_exon_saf(intron_coords, intron_flanking, FLANKING_EXON_SAF)

    # Write metadata
    write_intron_metadata(intron_coords, intron_flanking, INTRON_META)

    print("\n=== Step 1 complete ===", file=sys.stderr)
    print(f"  Intron SAF: {INTRON_SAF}", file=sys.stderr)
    print(f"  Flanking exon SAF: {FLANKING_EXON_SAF}", file=sys.stderr)
    print(f"  Metadata: {INTRON_META}", file=sys.stderr)
    print("\nNext: run featureCounts with run_featurecounts.sh", file=sys.stderr)


if __name__ == "__main__":
    main()
