#!/usr/bin/env python3
"""
Parse all .summary files from strandedness_verification/raw_summaries/ and
build Supplementary Table S1 as TSV + xlsx.

featureCounts .summary format (tab-delimited):
  Status                            BAM_path1   BAM_path2 ...
  Assigned                          NUM         NUM
  Unassigned_Unmapped               NUM         NUM
  Unassigned_Read_Type              NUM         NUM
  Unassigned_Singleton              NUM         NUM
  Unassigned_MappingQuality         NUM         NUM
  Unassigned_Chimera                NUM         NUM
  Unassigned_FragmentLength         NUM         NUM
  Unassigned_Duplicate              NUM         NUM
  Unassigned_MultiMapping           NUM         NUM
  Unassigned_Secondary              NUM         NUM
  Unassigned_NonSplit               NUM         NUM
  Unassigned_NoFeatures             NUM         NUM
  Unassigned_Overlapping_Length     NUM         NUM
  Unassigned_Ambiguity              NUM         NUM

We compute Assigned% = Assigned / (sum of all rows) per BAM, then average across BAMs.
"""
import os, glob, re
from pathlib import Path
import pandas as pd

RAW = Path("/mnt/HDD2/unified_results/strandedness_verification/raw_summaries")
OUT_DIR = Path("/mnt/HDD2/unified_results/strandedness_verification")
OUT_TSV = OUT_DIR / "Supplementary_Table_S1.tsv"
OUT_XLSX = OUT_DIR / "Supplementary_Table_S1.xlsx"

DATASET_META = {
    "Sakai_DHX8":    dict(label="DHX8 KD MEF (this study)", species="Mouse",
                          tissue="MEF", lib="NEBNext Ultra II Directional (dUTP)",
                          read="PE150", expected="-s 2"),
    "GSE65973":      dict(label="GSE65973 (Walker 2015)",    species="Mouse",
                          tissue="Cortex (dNLS-TDP-43)", lib="RNA ligation",
                          read="PE100", expected="-s 1"),
    "GSE124439":     dict(label="GSE124439 (Tam 2019)",      species="Human",
                          tissue="Frontal cortex (ALS)", lib="KAPA Stranded mRNA (dUTP)",
                          read="PE150", expected="-s 2"),
    "GSE116622":     dict(label="GSE116622 (Conlon 2018)",   species="Human",
                          tissue="Cerebellum (ALS)", lib="KAPA Stranded (dUTP)",
                          read="PE150", expected="-s 2"),
    "GSE126543":     dict(label="GSE126543 (Liu 2019)",      species="Human",
                          tissue="FANS-sorted nuclei (FTD-ALS)", lib="NuGEN Ovation V2 (unstranded)",
                          read="PE75",  expected="-s 0"),
    "ENCODE_PRPF8_HepG2":     dict(label="ENCODE PRPF8 KD",     species="Human",
                                   tissue="HepG2", lib="ENCODE standard (dUTP)",
                                   read="PE100", expected="-s 2"),
    "ENCODE_EFTUD2_HepG2":    dict(label="ENCODE EFTUD2 KD",    species="Human",
                                   tissue="HepG2", lib="ENCODE standard (dUTP)",
                                   read="PE100", expected="-s 2"),
    "ENCODE_SNRNP200_HepG2":  dict(label="ENCODE SNRNP200 KD",  species="Human",
                                   tissue="HepG2", lib="ENCODE standard (dUTP)",
                                   read="PE100", expected="-s 2"),
    "ENCODE_PRPF6_HepG2":     dict(label="ENCODE PRPF6 KD",     species="Human",
                                   tissue="HepG2", lib="ENCODE standard (dUTP)",
                                   read="PE100", expected="-s 2"),
    "ENCODE_RBFOX2_K562":     dict(label="ENCODE RBFOX2 KD",    species="Human",
                                   tissue="K562",  lib="ENCODE standard (dUTP)",
                                   read="PE100", expected="-s 2"),
    "ENCODE_RBM39_K562":      dict(label="ENCODE RBM39 KD",     species="Human",
                                   tissue="K562",  lib="ENCODE standard (dUTP)",
                                   read="PE100", expected="-s 2"),
}

def parse_summary(path):
    """Return Assigned% per BAM column (mean across BAMs)."""
    df = pd.read_csv(path, sep="\t")
    df = df.set_index("Status")
    cols = df.columns
    pct = []
    for c in cols:
        total = df[c].sum()
        if total <= 0:
            continue
        pct.append(100.0 * df.loc["Assigned", c] / total)
    if not pct:
        return None, 0
    return sum(pct) / len(pct), len(pct)


rows = []
for key, meta in DATASET_META.items():
    record = {"Dataset": meta["label"], "Species": meta["species"],
              "Tissue/Cell": meta["tissue"], "Library prep": meta["lib"],
              "Read": meta["read"]}
    n_bams_seen = 0
    for s in (0, 1, 2):
        # try several naming conventions
        for cand in (f"{key}_s{s}.summary",
                     f"{key}_s{s}_repBAM.summary",
                     f"{key}_s{s}_full.summary",
                     f"{key}_s{s}_fc.summary"):
            p = RAW / cand
            if p.exists():
                pct, n = parse_summary(p)
                record[f"-s {s} Assigned%"] = round(pct, 2) if pct is not None else None
                n_bams_seen = max(n_bams_seen, n)
                break
        else:
            record[f"-s {s} Assigned%"] = None
    record["n BAMs"] = n_bams_seen
    record["Setting used"] = meta["expected"]
    # Strand verification: a stranded library (-s 1 or -s 2) should keep ≥80%
    # of the unstranded (-s 0) counts in the correct direction and <20% in the
    # opposite direction. Unstranded (-s 0) libraries should give ~equal
    # counts in -s 1 and -s 2 when run with those settings (proxy for double counting).
    rates = {s: record[f"-s {s} Assigned%"] for s in (0,1,2)}
    s0, s1, s2 = rates.get(0), rates.get(1), rates.get(2)
    chosen = int(meta["expected"].replace("-s ", ""))
    verdict = "—"
    if s0 and s0 > 0:
        if chosen == 2:
            r2 = (s2 / s0) if s2 is not None else None
            r1 = (s1 / s0) if s1 is not None else None
            if r2 is not None and r2 >= 0.80 and (r1 is None or r1 < 0.20):
                verdict = f"✓ s2/s0={r2:.2f}, s1/s0={r1:.2f}" if r1 is not None else f"✓ s2/s0={r2:.2f}"
            else:
                verdict = f"✗ s2/s0={r2}, s1/s0={r1}"
        elif chosen == 1:
            r1 = (s1 / s0) if s1 is not None else None
            r2 = (s2 / s0) if s2 is not None else None
            if r1 is not None and r1 >= 0.80 and (r2 is None or r2 < 0.20):
                verdict = f"✓ s1/s0={r1:.2f}, s2/s0={r2:.2f}" if r2 is not None else f"✓ s1/s0={r1:.2f}"
            else:
                verdict = f"✗ s1/s0={r1}, s2/s0={r2}"
        elif chosen == 0:
            r1 = (s1 / s0) if s1 is not None else None
            r2 = (s2 / s0) if s2 is not None else None
            # for unstranded, s1 and s2 should be roughly equal and cover most reads
            if r1 is not None and r2 is not None and abs(r1 - r2) < 0.15 and (r1 + r2) > 0.6:
                verdict = f"✓ s1/s0={r1:.2f}, s2/s0={r2:.2f}"
            elif r1 is not None and r2 is not None:
                verdict = f"~ s1/s0={r1:.2f}, s2/s0={r2:.2f}"
    record["Verification"] = verdict
    rows.append(record)

out = pd.DataFrame(rows)
order = ["Dataset","Species","Tissue/Cell","Library prep","Read",
         "-s 0 Assigned%","-s 1 Assigned%","-s 2 Assigned%",
         "n BAMs","Setting used","Verification"]
out = out[order]

OUT_DIR.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.2f")
try:
    out.to_excel(OUT_XLSX, index=False)
except Exception as e:
    print("xlsx write failed:", e)

print(out.to_string(index=False))
print(f"\nTSV  → {OUT_TSV}")
print(f"XLSX → {OUT_XLSX}")
