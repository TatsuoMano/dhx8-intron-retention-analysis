#!/usr/bin/env python3
"""
Step 1 (hg38): Generate intron SAF, flanking exon SAF, and intron metadata
from hg38.ncbiRefSeq.gtf exon annotations.

Introns are inferred from gaps between consecutive exons within each transcript.
Adapted from mm10 step1_per_intron_count.py with intron inference from exon GTF.

Output: /mnt/HDD2/hg38_ir_annotations/
  - hg38_introns.saf
  - hg38_flanking_exons.saf
  - hg38_intron_metadata.tsv
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

EXON_GTF = "/mnt/NAS/database/genes/hg38.ncbiRefSeq.gtf"
OUTDIR = Path("/mnt/HDD2/hg38_ir_annotations")

INTRON_SAF = OUTDIR / "hg38_introns.saf"
FLANKING_EXON_SAF = OUTDIR / "hg38_flanking_exons.saf"
INTRON_META = OUTDIR / "hg38_intron_metadata.tsv"


def parse_gtf_attribute(attr_str, key):
    m = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return m.group(1) if m else None


def load_exons(gtf_path):
    tx_exons = defaultdict(list)
    print(f"Loading exons from {gtf_path}...", file=sys.stderr)
    n = 0
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom = cols[0]
            if "_" in chrom or chrom.startswith("chrUn"):
                continue
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = cols[8]
            tid = parse_gtf_attribute(attrs, "transcript_id")
            gene_name = parse_gtf_attribute(attrs, "gene_name") or "unknown"
            if tid:
                tx_exons[tid].append((chrom, start, end, strand, gene_name))
                n += 1
    for tid in tx_exons:
        tx_exons[tid].sort(key=lambda x: x[1])
    print(f"  Loaded {n} exon entries from {len(tx_exons)} transcripts", file=sys.stderr)
    return tx_exons


def infer_introns_and_flanking(tx_exons):
    intron_coords = {}
    intron_flanking = defaultdict(set)

    for tid, exons in tx_exons.items():
        if len(exons) < 2:
            continue
        for i in range(len(exons) - 1):
            e1 = exons[i]
            e2 = exons[i + 1]
            if e1[0] != e2[0]:
                continue
            intron_start = e1[2] + 1
            intron_end = e2[1] - 1
            if intron_end <= intron_start:
                continue
            chrom = e1[0]
            strand = e1[3]
            gene_name = e1[4]
            ikey = (chrom, intron_start, intron_end, strand)

            if ikey not in intron_coords:
                intron_coords[ikey] = (gene_name, {tid})
            else:
                existing_gene, tids = intron_coords[ikey]
                tids.add(tid)
                if existing_gene == "unknown" and gene_name != "unknown":
                    intron_coords[ikey] = (gene_name, tids)

            intron_flanking[ikey].add((chrom, e1[1], e1[2], strand))
            intron_flanking[ikey].add((chrom, e2[1], e2[2], strand))

    print(f"  Inferred {len(intron_coords)} non-redundant introns", file=sys.stderr)
    n_with_flanking = sum(1 for v in intron_flanking.values() if v)
    print(f"  Introns with flanking exons: {n_with_flanking}", file=sys.stderr)
    return intron_coords, intron_flanking


def write_intron_saf(intron_coords, outpath):
    with open(outpath, "w") as f:
        f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for (chrom, start, end, strand), (gene_name, tids) in sorted(intron_coords.items()):
            intron_id = f"{chrom}:{start}-{end}:{strand}:{gene_name}"
            f.write(f"{intron_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
    print(f"Wrote {len(intron_coords)} introns to {outpath}", file=sys.stderr)


def write_flanking_exon_saf(intron_coords, intron_flanking, outpath):
    exon_entries = set()
    for icoord in intron_coords:
        if icoord in intron_flanking:
            for echrom, estart, eend, estrand in intron_flanking[icoord]:
                exon_entries.add((echrom, estart, eend, estrand))

    with open(outpath, "w") as f:
        f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for (chrom, start, end, strand) in sorted(exon_entries):
            exon_id = f"{chrom}:{start}-{end}:{strand}"
            f.write(f"{exon_id}\t{chrom}\t{start}\t{end}\t{strand}\n")
    print(f"Wrote {len(exon_entries)} unique flanking exon regions to {outpath}", file=sys.stderr)


def write_intron_metadata(intron_coords, intron_flanking, outpath):
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
    tx_exons = load_exons(EXON_GTF)
    intron_coords, intron_flanking = infer_introns_and_flanking(tx_exons)
    write_intron_saf(intron_coords, INTRON_SAF)
    write_flanking_exon_saf(intron_coords, intron_flanking, FLANKING_EXON_SAF)
    write_intron_metadata(intron_coords, intron_flanking, INTRON_META)
    print("\n=== Step 1 (hg38) complete ===", file=sys.stderr)
    print(f"  Intron SAF: {INTRON_SAF}", file=sys.stderr)
    print(f"  Flanking exon SAF: {FLANKING_EXON_SAF}", file=sys.stderr)
    print(f"  Metadata: {INTRON_META}", file=sys.stderr)


if __name__ == "__main__":
    main()
