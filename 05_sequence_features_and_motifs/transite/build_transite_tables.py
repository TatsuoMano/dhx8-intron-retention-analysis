#!/usr/bin/env python3
"""Transite foreground and background tables (RefSeq transcript per gene).

RECONSTRUCTED FILE (code release, 2026-10-20). The tables were written on 2026-07-08 by
/tmp/gen_transite_tables2.py on the server and copied to /mnt/HDD2/transite_fig3cd/ on
2026-07-13; the strict foreground was subset inside run_tsma_strict.R and the 3'UTR-mapped
foreground by an Rscript -e one-liner. None of these was kept as a file. The bodies were recovered
from the session logs (~/.claude-mem/logs/claude-mem-2026-07-07.log, 2026-07-08 08:38 JST, for
gen_transite_tables2.py; claude-mem-2026-07-13.log, 13:43 JST, for the mapped subset). The
selection logic below is that code unchanged; the /tmp output paths were replaced by --outdir and
a comparison with the server tables was added.

Construction:
  1. From mm10.ncbiRefSeq.gtf, per transcript: gene_name and total exon length (summed over all
     exon lines of the transcript ID).
  2. One representative coding transcript per gene: NM_ if the gene has any, otherwise XM_; among
     those the longest (total exon length; first in GTF order on ties). Genes with only NR_/XR_
     transcripts are dropped (no 3' UTR).
  3. DESeq2 v2 (deseq2_dhx8_kd_vs_ctrl.tsv, rows with padj and log2FC): up = padj < 0.05 and
     log2FC > 0.5, down = padj < 0.05 and log2FC < -0.5, background = all genes.
     -> foreground_up_only_v2.tsv, foreground_up_down_v2.tsv, foreground_up_only_v2_unversioned.tsv,
        background_all_expressed_v2.tsv
  4. Strict DEG-UP foreground (as in run_tsma_strict.R): rows of foreground_up_only_v2.tsv whose
     gene is in padj < 0.05 and log2FC > 1 -> foreground_up_strict_v2.tsv
  5. Rows of (4) whose transcript has a sequence in mm10_3utr_seqs.rds (extract_3utr.R) ->
     foreground_up_strict_mapped451.tsv (the Transite foreground, n = 451). The transcript IDs in
     the .rds are read with Rscript.

Run on the server (read-only, streamed; nothing is written unless --outdir is given):
    ssh ncnp-linux 'python3 -' < 05_sequence_features_and_motifs/transite/build_transite_tables.py
"""
import argparse
import hashlib
import os
import re
import subprocess
from collections import defaultdict

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", default="/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf")
ap.add_argument("--deseq", default="/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv")
ap.add_argument("--utr-rds", default="/mnt/HDD2/transite_fig3cd/mm10_3utr_seqs.rds")
ap.add_argument("--compare", default="/mnt/HDD2/transite_fig3cd")
ap.add_argument("--outdir", default=None)
a = ap.parse_args()

# --- 1-2: transcript per gene (verbatim logic of gen_transite_tables2.py) ---
tx_gene = {}; tx_len = {}
rgid = re.compile(r'gene_name "([^"]+)"'); rtid = re.compile(r'transcript_id "([^"]+)"')
with open(a.gtf) as f:
    for line in f:
        if "\texon\t" not in line: continue
        p = line.split("\t"); g = rgid.search(line); t = rtid.search(line)
        if not g or not t: continue
        tx_gene[t.group(1)] = g.group(1); tx_len[t.group(1)] = tx_len.get(t.group(1), 0) + (int(p[4]) - int(p[3]) + 1)
gene_tx = defaultdict(list)
for tid, gn in tx_gene.items(): gene_tx[gn].append(tid)
def best(txs):   # CODING only: NM_ preferred, then XM_; longest. Non-coding (NR_/XR_) excluded.
    nm = [t for t in txs if t.startswith("NM_")]; xm = [t for t in txs if t.startswith("XM_")]
    pool = nm if nm else xm
    return max(pool, key=lambda t: tx_len.get(t, 0)) if pool else None
gene_best = {gn: b for gn, txs in gene_tx.items() if (b := best(txs))}
print("coding genes (NM_/XM_):", len(gene_best))

# --- 3: DESeq2 groups ---
d = pd.read_csv(a.deseq, sep="\t").dropna(subset=["padj", "log2FoldChange"])
up   = sorted(d.loc[(d.padj < 0.05) & (d.log2FoldChange >  0.5), "gene"].astype(str).unique())
dn   = sorted(d.loc[(d.padj < 0.05) & (d.log2FoldChange < -0.5), "gene"].astype(str).unique())
allg = sorted(d["gene"].astype(str).unique())
def to_tx(genes): return [(g, gene_best[g]) for g in genes if g in gene_best]
up_tx = to_tx(up); dn_tx = to_tx(dn); bg_tx = to_tx(allg)
print(f"UP {len(up)} genes -> {len(up_tx)} coding tx ({100*len(up_tx)/len(up):.0f}%); dropped non-coding {len(up)-len(up_tx)}")
print(f"DOWN {len(dn)} -> {len(dn_tx)} coding; BG {len(allg)} -> {len(bg_tx)} coding")

out = {}
out["foreground_up_down_v2.tsv"] = "refseq_id\tgroup\tgene_symbol\n" + \
    "".join(f"{t}\tupregulated\t{g}\n" for g, t in up_tx) + "".join(f"{t}\tdownregulated\t{g}\n" for g, t in dn_tx)
out["foreground_up_only_v2.tsv"] = "refseq_id\tgroup\tgene_symbol\n" + "".join(f"{t}\tupregulated\t{g}\n" for g, t in up_tx)
out["background_all_expressed_v2.tsv"] = "refseq_id\tgene_symbol\n" + "".join(f"{t}\t{g}\n" for g, t in bg_tx)
out["foreground_up_only_v2_unversioned.tsv"] = "refseq_id\tgroup\tgene_symbol\n" + \
    "".join(f"{t.split('.')[0]}\tupregulated\t{g}\n" for g, t in up_tx)

# --- 4: strict DEG-UP subset (run_tsma_strict.R) ---
v2 = pd.read_csv(a.deseq, sep="\t")
strict_genes = set(v2.loc[v2.padj.notna() & (v2.padj < 0.05) & (v2.log2FoldChange > 1), "gene"].astype(str))
strict = [(g, t) for g, t in up_tx if g in strict_genes]
out["foreground_up_strict_v2.tsv"] = "refseq_id\tgroup\tgene_symbol\n" + "".join(f"{t}\tupregulated\t{g}\n" for g, t in strict)
print(f"strict DEG-UP genes (log2FC>1 & padj<0.05): {len(strict_genes)} | with a coding transcript (foreground rows): {len(strict)}")

# --- 5: transcripts with a 3'UTR sequence ---
try:
    r = subprocess.run(["Rscript", "-e", f'suppressMessages(library(Biostrings)); cat(names(readRDS("{a.utr_rds}")), sep="\\n")'],
                       capture_output=True, text=True, check=True)
    utr_ids = set(r.stdout.split())
    mapped = [(g, t) for g, t in strict if t in utr_ids]
    out["foreground_up_strict_mapped451.tsv"] = "refseq_id\tgroup\tgene_symbol\n" + "".join(f"{t}\tupregulated\t{g}\n" for g, t in mapped)
    print(f"3'UTR sequences in {os.path.basename(a.utr_rds)}: {len(utr_ids)} | strict foreground with a 3'UTR (Transite n): {len(mapped)}")
    print(f"background with a 3'UTR: {sum(t in utr_ids for _, t in bg_tx)}")
except Exception as e:
    print("Rscript step skipped:", e)

print("\nfile                                    rows  md5 (built)                       identical to server copy")
for name, text in out.items():
    md5 = hashlib.md5(text.encode()).hexdigest()
    ref = os.path.join(a.compare, name)
    same = "(no server copy)"
    if os.path.exists(ref):
        same = "yes" if open(ref, "rb").read() == text.encode() else "NO"
    print(f"{name:<40} {text.count(chr(10)) - 1:>5}  {md5}  {same}")
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        with open(os.path.join(a.outdir, name), "w") as fh:
            fh.write(text)
