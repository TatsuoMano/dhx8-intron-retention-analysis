#!/bin/bash
set -euo pipefail
ulimit -n 65536

# ============================================================
# GSE126543 — FACS-sorted neuronal nuclei RNA-seq
# TDP-43+ vs TDP-43- (7 FTD/FTD-ALS patients, frontal cortex)
# Liu/Lee 2019 Cell Rep; reused by Brown/Ma 2022 Nature
# ============================================================

PROJ=/mnt/HDD2/sorted_nuclei_GSE126543
STAR_INDEX=/mnt/HDD2/STAR_index_hg38
GTF=/mnt/NAS/database/genes/hg38.ncbiRefSeq.gtf
THREADS=48
NOTIFY=~/bin/notify_telegram.sh

mkdir -p ${PROJ}/{fastq,mapping,counts,analysis}

# Sample sheet — GSE126542 (RNA-seq SubSeries of GSE126543)
# 16 samples: 7 TDP43+, 7 TDP43-, 2 Unsorted
cat > ${PROJ}/sample_sheet.tsv << 'EOF'
SRR	Sample	Subject	Diagnosis	TDP43_status	NeuN
SRR8571937	S1_TDP43pos	S1	FTD	TDP43_positive	NeuN+
SRR8571938	S1_TDP43neg	S1	FTD	TDP43_negative	NeuN+
SRR8571939	S2_TDP43pos	S2	FTD-ALS	TDP43_positive	NeuN+
SRR8571940	S2_TDP43neg	S2	FTD-ALS	TDP43_negative	NeuN+
SRR8571941	S3_TDP43pos	S3	FTD	TDP43_positive	NeuN+
SRR8571942	S3_TDP43neg	S3	FTD	TDP43_negative	NeuN+
SRR8571943	S3_Unsorted	S3	FTD	Unsorted	Unsorted
SRR8571944	S4_TDP43pos	S4	FTD	TDP43_positive	NeuN+
SRR8571945	S4_TDP43neg	S4	FTD	TDP43_negative	NeuN+
SRR8571946	S4_Unsorted	S4	FTD	Unsorted	Unsorted
SRR8571947	S5_TDP43pos	S5	FTD	TDP43_positive	NeuN+
SRR8571948	S5_TDP43neg	S5	FTD	TDP43_negative	NeuN+
SRR8571949	S6_TDP43pos	S6	FTD-ALS	TDP43_positive	NeuN+
SRR8571950	S6_TDP43neg	S6	FTD-ALS	TDP43_negative	NeuN+
SRR8571951	S7_TDP43pos	S7	FTD-ALS	TDP43_positive	NeuN+
SRR8571952	S7_TDP43neg	S7	FTD-ALS	TDP43_negative	NeuN+
EOF

${NOTIFY} "[DHX8-Sakai] GSE126543: Pipeline started — 16 sorted nuclei RNA-seq samples"

# ============================================================
# Step 1: Download FASTQ from SRA
# ============================================================
echo "=== Step 1: Download FASTQ ==="
cd ${PROJ}/fastq

SRRS=$(tail -n+2 ${PROJ}/sample_sheet.tsv | cut -f1)
TOTAL=$(echo "$SRRS" | wc -l)
COUNT=0

SRA_DIR=${PROJ}/sra
mkdir -p ${SRA_DIR}

# Step 1a: prefetch all SRA files first
echo "--- Step 1a: prefetch SRA files ---"
for SRR in $SRRS; do
    COUNT=$((COUNT + 1))
    if [ -f "${SRR}_1.fastq.gz" ] && [ -f "${SRR}_2.fastq.gz" ]; then
        echo "  [${COUNT}/${TOTAL}] ${SRR} — FASTQ already exists, skipping"
        continue
    fi
    if [ -f "${SRA_DIR}/${SRR}/${SRR}.sra" ] || [ -f "${SRA_DIR}/${SRR}.sra" ]; then
        echo "  [${COUNT}/${TOTAL}] ${SRR} — SRA already downloaded, skipping"
        continue
    fi
    echo "  [${COUNT}/${TOTAL}] prefetch ${SRR}..."
    prefetch ${SRR} --output-directory ${SRA_DIR} --max-size 50G
    echo "  [${COUNT}/${TOTAL}] ${SRR} — prefetch done"
done

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 1a/4: prefetch complete (${TOTAL} SRA files)"

# Step 1b: convert SRA → FASTQ
echo "--- Step 1b: fasterq-dump SRA → FASTQ ---"
COUNT=0
for SRR in $SRRS; do
    COUNT=$((COUNT + 1))
    if [ -f "${SRR}_1.fastq.gz" ] && [ -f "${SRR}_2.fastq.gz" ]; then
        echo "  [${COUNT}/${TOTAL}] ${SRR} — FASTQ already exists, skipping"
        continue
    fi
    echo "  [${COUNT}/${TOTAL}] fasterq-dump ${SRR}..."
    SRA_FILE=$(find ${SRA_DIR} -name "${SRR}.sra" 2>/dev/null | head -1)
    if [ -z "${SRA_FILE}" ]; then
        echo "  ERROR: SRA file not found for ${SRR}"
        continue
    fi
    fasterq-dump ${SRA_FILE} --threads 8 --split-files --outdir .
    pigz -p 8 ${SRR}_1.fastq ${SRR}_2.fastq 2>/dev/null || gzip ${SRR}_1.fastq ${SRR}_2.fastq
    echo "  [${COUNT}/${TOTAL}] ${SRR} — FASTQ done"
done

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 1/4: FASTQ download complete (${TOTAL} samples)"

# ============================================================
# Step 2: STAR Mapping
# ============================================================
echo "=== Step 2: STAR Mapping ==="
COUNT=0

for SRR in $SRRS; do
    COUNT=$((COUNT + 1))
    OUTDIR=${PROJ}/mapping/${SRR}

    if [ -f "${OUTDIR}/${SRR}.sorted.bam" ]; then
        echo "  [${COUNT}/${TOTAL}] ${SRR} — already mapped, skipping"
        continue
    fi

    echo "  [${COUNT}/${TOTAL}] Mapping ${SRR}..."
    mkdir -p ${OUTDIR}

    STAR --runThreadN ${THREADS} \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${PROJ}/fastq/${SRR}_1.fastq.gz ${PROJ}/fastq/${SRR}_2.fastq.gz \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${OUTDIR}/ \
         --outSAMstrandField intronMotif \
         --twopassMode Basic \
         --outFilterMultimapNmax 1 \
         --outSAMattributes NH HI AS NM MD \
         --quantMode GeneCounts \
         --sjdbGTFfile ${GTF} \
         --outBAMsortingBinsN 50 \
         --limitBAMsortRAM 32000000000

    # Rename and index
    mv ${OUTDIR}/Aligned.sortedByCoord.out.bam ${OUTDIR}/${SRR}.sorted.bam
    samtools index -@ 8 ${OUTDIR}/${SRR}.sorted.bam

    echo "  [${COUNT}/${TOTAL}] ${SRR} — mapped"
done

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 2/4: STAR mapping complete (${TOTAL} samples)"

# ============================================================
# Step 3: featureCounts (gene + intron)
# ============================================================
echo "=== Step 3: featureCounts ==="

# Gene counts
BAMS=""
for SRR in $SRRS; do
    BAMS="${BAMS} ${PROJ}/mapping/${SRR}/${SRR}.sorted.bam"
done

featureCounts -T ${THREADS} -p --countReadPairs -s 0 \
    -a ${GTF} -o ${PROJ}/counts/gene_counts.txt \
    ${BAMS}

# Intron counts (for IR analysis)
INTRON_GTF=/mnt/HDD2/phase3_specificity/hg38.ncbiRefSeq.introns.gtf
if [ -f "${INTRON_GTF}" ]; then
    featureCounts -T ${THREADS} -p --countReadPairs -s 0 \
        -a ${INTRON_GTF} -o ${PROJ}/counts/intron_counts.txt \
        ${BAMS}
fi

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 3/5: featureCounts complete"

# ============================================================
# Step 4: IR analysis — TDP-43+ vs TDP-43- (paired)
# ============================================================
echo "=== Step 4: IR analysis ==="

python3 << 'IR_PYEOF'
import csv
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path

PROJ = Path('/mnt/HDD2/sorted_nuclei_GSE126543')
DHX8_IR_PATH = Path('/mnt/HDD2/ALS_patient_GSE124439/limma_ir_results_dhx8.tsv')

print('='*60)
print('GSE126543 IR Analysis — TDP-43+ vs TDP-43- (sorted nuclei)')
print('='*60)

# Load sample sheet
samples = []
with open(PROJ / 'sample_sheet.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        samples.append(row)

# Identify TDP43+ and TDP43- samples (exclude Unsorted)
tdp_pos = [s for s in samples if s['TDP43_status'] == 'TDP43_positive']
tdp_neg = [s for s in samples if s['TDP43_status'] == 'TDP43_negative']
print(f'TDP-43+ samples: {len(tdp_pos)}')
print(f'TDP-43- samples: {len(tdp_neg)}')

# Load intron counts
intron_file = PROJ / 'counts' / 'intron_counts.txt'
if not intron_file.exists():
    print('ERROR: intron_counts.txt not found — skipping IR analysis')
    exit(0)

# Parse featureCounts output
header = None
intron_ids = []
gene_names = []
counts_matrix = []

with open(intron_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if header is None:
            header = parts
            # Find sample columns (BAM paths)
            sample_cols = parts[6:]  # featureCounts: first 6 cols are annotation
            continue
        intron_id = parts[0]
        gene = parts[0].split(':')[-1] if ':' in parts[0] else parts[0]
        counts = [int(x) for x in parts[6:]]
        intron_ids.append(intron_id)
        gene_names.append(gene)
        counts_matrix.append(counts)

counts_matrix = np.array(counts_matrix, dtype=float)
print(f'Loaded {len(intron_ids)} introns, {len(sample_cols)} samples')

# Map sample columns to SRR IDs
# featureCounts uses BAM path as column name
srr_to_col = {}
for i, col in enumerate(sample_cols):
    for s in samples:
        if s['SRR'] in col:
            srr_to_col[s['SRR']] = i
            break

# Get column indices for TDP43+ and TDP43-
pos_idx = [srr_to_col[s['SRR']] for s in tdp_pos if s['SRR'] in srr_to_col]
neg_idx = [srr_to_col[s['SRR']] for s in tdp_neg if s['SRR'] in srr_to_col]
print(f'Matched columns: TDP43+ = {len(pos_idx)}, TDP43- = {len(neg_idx)}')

# Normalize: RPK (reads per kilobase)
# For simplicity, use raw counts with library size normalization
# Library sizes
lib_sizes = counts_matrix.sum(axis=0)
lib_sizes[lib_sizes == 0] = 1
# CPM-like normalization
norm_matrix = counts_matrix / lib_sizes * 1e6

# Calculate log2FC for each intron: TDP43- vs TDP43+
# TDP43- is the "disease" condition (lost TDP-43)
pos_mean = norm_matrix[:, pos_idx].mean(axis=1)
neg_mean = norm_matrix[:, neg_idx].mean(axis=1)

pseudocount = 0.1
log2fc = np.log2((neg_mean + pseudocount) / (pos_mean + pseudocount))

# Mann-Whitney U test per intron
pvals = np.ones(len(intron_ids))
for i in range(len(intron_ids)):
    pos_vals = norm_matrix[i, pos_idx]
    neg_vals = norm_matrix[i, neg_idx]
    if pos_vals.sum() + neg_vals.sum() > 0:
        try:
            _, p = stats.mannwhitneyu(neg_vals, pos_vals, alternative='two-sided')
            pvals[i] = p
        except ValueError:
            pass

# BH correction
def benjamini_hochberg(pvals):
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    adjusted = np.zeros(n)
    cummin = 1.0
    for i in range(n-1, -1, -1):
        val = min(sorted_pvals[i] * n / (i + 1), 1.0)
        cummin = min(cummin, val)
        adjusted[sorted_idx[i]] = cummin
    return adjusted

padj = benjamini_hochberg(pvals)

# Significance
sig_mask = padj < 0.05
sig_up = sig_mask & (log2fc > 0)  # IR increase in TDP43-
sig_down = sig_mask & (log2fc < 0)  # IR decrease in TDP43-

print(f'\nResults (TDP-43 negative vs positive):')
print(f'  Median log2FC: {np.median(log2fc):.4f}')
print(f'  Significant introns (FDR<0.05): {sig_mask.sum()}')
print(f'    IR increase (TDP43- > TDP43+): {sig_up.sum()}')
print(f'    IR decrease (TDP43- < TDP43+): {sig_down.sum()}')

# Gene-level summary
sig_genes = set()
sig_genes_up = set()
sig_genes_down = set()
for i in range(len(intron_ids)):
    if padj[i] < 0.05:
        g = gene_names[i].upper()
        sig_genes.add(g)
        if log2fc[i] > 0:
            sig_genes_up.add(g)
        else:
            sig_genes_down.add(g)

print(f'  Significant genes: {len(sig_genes)} (↑{len(sig_genes_up)}, ↓{len(sig_genes_down)})')

# DHX8 comparison
if DHX8_IR_PATH.exists():
    dhx8_ir_genes = set()
    with open(DHX8_IR_PATH) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('retained', '').upper() == 'TRUE':
                dhx8_ir_genes.add(row['gene_name'].upper())

    print(f'\n--- DHX8 KD IR comparison ---')
    print(f'DHX8 KD IR genes (all UP): {len(dhx8_ir_genes)}')

    # Direction-aware overlap
    overlap_all = dhx8_ir_genes & sig_genes
    overlap_up = dhx8_ir_genes & sig_genes_up
    overlap_down = dhx8_ir_genes & sig_genes_down

    print(f'Overlap (all): {len(overlap_all)}')
    print(f'Overlap (IR UP in TDP43-): {len(overlap_up)}')
    print(f'Overlap (IR DOWN in TDP43-): {len(overlap_down)}')

    # Fisher test for UP direction
    all_genes_upper = set(g.upper() for g in gene_names)
    N = len(all_genes_upper)

    for direction, overlap_set, sig_set, label in [
        ('UP', overlap_up, sig_genes_up, 'IR increase'),
        ('DOWN', overlap_down, sig_genes_down, 'IR decrease'),
        ('ALL', overlap_all, sig_genes, 'all directions'),
    ]:
        a = len(overlap_set)
        b = len(dhx8_ir_genes - sig_set)
        c = len(sig_set - dhx8_ir_genes)
        d = N - a - b - c
        if d > 0 and a > 0:
            odds, p = stats.fisher_exact([[a, b], [c, d]])
            print(f'  Fisher ({label}): OR={odds:.3f}, p={p:.2e}, overlap={a}')
        else:
            print(f'  Fisher ({label}): overlap={a} (insufficient for test)')

# Save full results
outpath = PROJ / 'analysis' / 'ir_results.tsv'
with open(outpath, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Intron_ID', 'Gene', 'TDP43pos_mean', 'TDP43neg_mean', 'log2FC', 'pval', 'padj'])
    for i in range(len(intron_ids)):
        writer.writerow([intron_ids[i], gene_names[i],
                         f'{pos_mean[i]:.4f}', f'{neg_mean[i]:.4f}',
                         f'{log2fc[i]:.4f}', f'{pvals[i]:.4e}', f'{padj[i]:.4e}'])
print(f'\nSaved: {outpath}')

# Save summary
summary_path = PROJ / 'analysis' / 'ir_summary.tsv'
with open(summary_path, 'w') as f:
    f.write(f'Comparison\tN_TDP43pos\tN_TDP43neg\tSig_introns\tSig_up\tSig_down\tSig_genes\tSig_genes_up\tSig_genes_down\tMedian_log2FC\n')
    f.write(f'TDP43neg_vs_TDP43pos\t{len(pos_idx)}\t{len(neg_idx)}\t{sig_mask.sum()}\t{sig_up.sum()}\t{sig_down.sum()}\t{len(sig_genes)}\t{len(sig_genes_up)}\t{len(sig_genes_down)}\t{np.median(log2fc):.4f}\n')
print(f'Saved: {summary_path}')

IR_PYEOF

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 4/5: IR analysis complete"

# ============================================================
# Step 5: CE v4 analysis
# ============================================================
echo "=== Step 5: CE v4 analysis ==="

# Create CE analysis script adapted for this dataset
python3 << 'PYEOF'
import csv, subprocess, re, os
import numpy as np

PROJ = '/mnt/HDD2/sorted_nuclei_GSE126543'
SAMTOOLS = '/usr/bin/samtools'
TOLERANCE = 10

TARGETS = {
    'UNC13A': {
        'chr': 'chr19', 'strand': '-',
        'upstream_exon': {'name': 'exon21', 'start': 17641393, 'end': 17641556},
        'downstream_exon': {'name': 'exon20', 'start': 17642845, 'end': 17642960},
        'ce': {'start': 17642414, 'end': 17642541},
        'normal_junction': (17641556, 17642845),
        'ce_junctions': [
            ('exon21-CE', 17641556, 17642414),
            ('CE-exon20', 17642541, 17642845),
        ],
    },
    'STMN2': {
        'chr': 'chr8', 'strand': '+',
        'upstream_exon': {'name': 'exon1', 'start': 79611117, 'end': 79611214},
        'downstream_exon': {'name': 'exon2', 'start': 79636802, 'end': 79636897},
        'ce': {'start': 79616822, 'end': 79617048},
        'normal_junction': (79611214, 79636802),
        'ce_junctions': [
            ('exon1-CE', 79611214, 79616822),
            ('CE-exon2', 79617048, 79636802),
        ],
    },
    'HDGFL2': {
        'chr': 'chr19', 'strand': '+',
        'upstream_exon': {'name': 'exon6', 'start': 4491764, 'end': 4491835},
        'downstream_exon': {'name': 'exon7', 'start': 4493703, 'end': 4493862},
        'ce': {'start': 4492015, 'end': 4492152},
        'normal_junction': (4491835, 4493703),
        'ce_junctions': [
            ('exon6-CE', 4491835, 4492015),
            ('CE-exon7', 4492152, 4493703),
        ],
    },
}

# Load sample sheet
samples = []
with open(f'{PROJ}/sample_sheet.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        samples.append(row)

print('='*60)
print('GSE126543 CE v4 Analysis — Sorted Nuclei')
print('='*60)

results = []

for s in samples:
    srr = s['SRR']
    bam = f'{PROJ}/mapping/{srr}/{srr}.sorted.bam'
    if not os.path.exists(bam):
        print(f'  SKIP {srr}: BAM not found')
        continue

    for gene, t in TARGETS.items():
        # Count CE junction reads from BAM
        chrom = t['chr']
        region_start = t['upstream_exon']['start'] - 100
        region_end = t['downstream_exon']['end'] + 100
        region = f'{chrom}:{region_start}-{region_end}'

        cmd = f'{SAMTOOLS} view -q 255 {bam} {region}'
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        ce_reads = 0
        normal_reads = 0

        for line in proc.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            pos = int(fields[3])
            cigar = fields[5]

            # Parse N junctions from CIGAR
            current_pos = pos
            for m in re.finditer(r'(\d+)([MDNISH=X])', cigar):
                length = int(m.group(1))
                op = m.group(2)
                if op == 'N':
                    intron_start = current_pos
                    intron_end = current_pos + length

                    # Check normal junction
                    norm = t['normal_junction']
                    if (abs(intron_start - norm[0]) <= TOLERANCE and
                        abs(intron_end - norm[1]) <= TOLERANCE):
                        normal_reads += 1

                    # Check CE junctions
                    for ce_name, ce_start, ce_end in t['ce_junctions']:
                        if (abs(intron_start - ce_start) <= TOLERANCE and
                            abs(intron_end - ce_end) <= TOLERANCE):
                            ce_reads += 1

                    current_pos += length
                elif op in ('M', 'D', '=', 'X'):
                    current_pos += length

        psi = 0
        if ce_reads + normal_reads > 0:
            psi = ce_reads / (ce_reads + normal_reads) * 100

        results.append({
            'SRR': srr,
            'Sample': s['Sample'],
            'Subject': s['Subject'],
            'Diagnosis': s['Diagnosis'],
            'TDP43': s['TDP43_status'],
            'Gene': gene,
            'CE_reads': ce_reads,
            'Normal_reads': normal_reads,
            'PSI': psi,
        })

        if ce_reads > 0:
            print(f'  {srr} ({s["Sample"]}) {gene}: CE={ce_reads}, Normal={normal_reads}, PSI={psi:.1f}%')

# Save results
outpath = f'{PROJ}/analysis/ce_v4_results.tsv'
with open(outpath, 'w') as f:
    if results:
        writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
print(f'\nSaved: {outpath}')

# DHX8 KD IR gene overlap with CE-positive genes
dhx8_path = '/mnt/HDD2/ALS_patient_GSE124439/limma_ir_results_dhx8.tsv'
if os.path.exists(dhx8_path):
    dhx8_ir_genes = set()
    with open(dhx8_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('retained', '').upper() == 'TRUE':
                dhx8_ir_genes.add(row['gene_name'].upper())

    print('\n--- CE genes vs DHX8 KD IR genes ---')
    for gene in TARGETS:
        is_dhx8 = gene.upper() in dhx8_ir_genes
        print(f'  {gene}: {"YES — in DHX8 KD IR gene list" if is_dhx8 else "NO"}')

# Summary by TDP43 status
print('\n' + '='*60)
print('Summary: CE detection by TDP-43 status')
print('='*60)
for gene in TARGETS:
    print(f'\n{gene}:')
    for status in ['TDP43_positive', 'TDP43_negative', 'Unsorted']:
        subset = [r for r in results if r['Gene'] == gene and r['TDP43'] == status]
        if not subset:
            continue
        n = len(subset)
        n_pos = sum(1 for r in subset if r['CE_reads'] > 0)
        mean_psi = np.mean([r['PSI'] for r in subset])
        mean_ce = np.mean([r['CE_reads'] for r in subset])
        print(f'  {status} (n={n}): {n_pos}/{n} positive ({100*n_pos/n:.0f}%), mean CE reads={mean_ce:.1f}, mean PSI={mean_psi:.1f}%')

PYEOF

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 5/6: CE v4 analysis complete"

# ============================================================
# Step 6: CE coverage + sashimi plots
# ============================================================
echo "=== Step 6: CE coverage + sashimi plots ==="

python3 << 'PLOT_PYEOF'
import csv, subprocess, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

PROJ = '/mnt/HDD2/sorted_nuclei_GSE126543'
SAMTOOLS = '/usr/bin/samtools'
FIGDIR = f'{PROJ}/analysis/figures'
os.makedirs(FIGDIR, exist_ok=True)

# Load sample sheet
samples = []
with open(f'{PROJ}/sample_sheet.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        samples.append(row)

TARGETS = {
    'UNC13A': {
        'chr': 'chr19',
        'region': 'chr19:17641200-17643100',
        'exons': [
            ('Exon 22', 17640800, 17641000),
            ('Exon 21', 17641393, 17641556),
            ('Cryptic\nExon', 17642414, 17642541),
            ('Exon 20', 17642845, 17642960),
            ('Exon 19', 17643200, 17643400),
        ],
        'ce': (17642414, 17642541),
        'strand': '-',
        'normal_junction': (17641556, 17642845),
        'ce_junctions': [(17641556, 17642414), (17642541, 17642845)],
    },
    'STMN2': {
        'chr': 'chr8',
        'region': 'chr8:79611000-79637000',
        'exons': [
            ('Exon 1', 79611117, 79611214),
            ('Cryptic\nExon', 79616822, 79617048),
            ('Exon 2', 79636802, 79636897),
        ],
        'ce': (79616822, 79617048),
        'strand': '+',
        'normal_junction': (79611214, 79636802),
        'ce_junctions': [(79611214, 79616822), (79617048, 79636802)],
    },
    'HDGFL2': {
        'chr': 'chr19',
        'region': 'chr19:4491600-4494000',
        'exons': [
            ('Exon 6', 4491764, 4491835),
            ('Cryptic\nExon', 4492015, 4492152),
            ('Exon 7', 4493703, 4493862),
        ],
        'ce': (4492015, 4492152),
        'strand': '+',
        'normal_junction': (4491835, 4493703),
        'ce_junctions': [(4491835, 4492015), (4492152, 4493703)],
    },
}

# Group samples
groups = {
    'TDP43_positive': {'samples': [], 'color': '#2ecc71', 'label': 'TDP-43+\n(NeuN+)'},
    'TDP43_negative': {'samples': [], 'color': '#e74c3c', 'label': 'TDP-43−\n(NeuN+)'},
}
for s in samples:
    status = s['TDP43_status']
    if status in groups:
        groups[status]['samples'].append(s)

# ============================================================
# Coverage plots — group mean
# ============================================================
for gene, t in TARGETS.items():
    print(f'Plotting coverage: {gene}')
    region = t['region']
    chrom, coords = region.split(':')
    rstart, rend = [int(x) for x in coords.split('-')]

    fig, axes = plt.subplots(2, 1, figsize=(14, 8),
                             gridspec_kw={'height_ratios': [1, 1]},
                             sharex=True)

    for gi, (status, ginfo) in enumerate(groups.items()):
        ax = axes[gi]
        all_depths = []

        for s in ginfo['samples']:
            srr = s['SRR']
            bam = f'{PROJ}/mapping/{srr}/{srr}.sorted.bam'
            if not os.path.exists(bam):
                continue

            cmd = f'{SAMTOOLS} depth -q 255 -r {region} {bam}'
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            positions = {}
            for line in proc.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                pos = int(parts[1])
                depth = int(parts[2])
                positions[pos] = depth

            pos_arr = np.arange(rstart, rend + 1)
            depth_arr = np.array([positions.get(p, 0) for p in pos_arr])
            all_depths.append(depth_arr)

        if all_depths:
            mean_depth = np.mean(all_depths, axis=0)
            pos_arr = np.arange(rstart, rend + 1)

            ax.fill_between(pos_arr, mean_depth, alpha=0.4, color=ginfo['color'])
            ax.plot(pos_arr, mean_depth, linewidth=0.5, color=ginfo['color'])

            # CE highlight
            ce = t['ce']
            ax.axvspan(ce[0], ce[1], alpha=0.2, color='red')

            # CE mean depth
            ce_mask = (pos_arr >= ce[0]) & (pos_arr <= ce[1])
            ce_mean = np.mean(mean_depth[ce_mask]) if ce_mask.any() else 0
            ax.text(0.98, 0.95, f'CE mean: {ce_mean:.1f}',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_ylabel(f'{ginfo["label"]}\n(n={len(ginfo["samples"])})')
        ax.set_xlim(rstart, rend)

    # Gene structure at bottom
    axes[-1].set_xlabel(f'Genomic position ({chrom})')
    fig.suptitle(f'{gene} Coverage — Sorted Nuclei (GSE126543)\n'
                 f'TDP-43+ vs TDP-43− neuronal nuclei', fontsize=13)
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/coverage_{gene.lower()}_sorted.png', dpi=150)
    fig.savefig(f'{FIGDIR}/coverage_{gene.lower()}_sorted.pdf')
    plt.close()
    print(f'  Saved: coverage_{gene.lower()}_sorted.png')

# ============================================================
# Sashimi plots
# ============================================================
import re

TOLERANCE = 10

for gene, t in TARGETS.items():
    print(f'Plotting sashimi: {gene}')
    region = t['region']
    chrom, coords = region.split(':')
    rstart, rend = [int(x) for x in coords.split('-')]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10),
                             gridspec_kw={'height_ratios': [1, 1]},
                             sharex=True)

    for gi, (status, ginfo) in enumerate(groups.items()):
        ax = axes[gi]
        all_depths = []
        junction_counts = {}

        for s in ginfo['samples']:
            srr = s['SRR']
            bam = f'{PROJ}/mapping/{srr}/{srr}.sorted.bam'
            if not os.path.exists(bam):
                continue

            # Depth
            cmd = f'{SAMTOOLS} depth -q 255 -r {region} {bam}'
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            positions = {}
            for line in proc.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                positions[int(parts[1])] = int(parts[2])
            pos_arr = np.arange(rstart, rend + 1)
            depth_arr = np.array([positions.get(p, 0) for p in pos_arr])
            all_depths.append(depth_arr)

            # Junctions from BAM
            cmd = f'{SAMTOOLS} view -q 255 {bam} {region}'
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            for line in proc.stdout.strip().split('\n'):
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 6:
                    continue
                pos = int(fields[3])
                cigar = fields[5]
                current_pos = pos
                for m in re.finditer(r'(\d+)([MDNISH=X])', cigar):
                    length = int(m.group(1))
                    op = m.group(2)
                    if op == 'N':
                        junc = (current_pos, current_pos + length)
                        junction_counts[junc] = junction_counts.get(junc, 0) + 1
                        current_pos += length
                    elif op in ('M', 'D', '=', 'X'):
                        current_pos += length

        if all_depths:
            mean_depth = np.mean(all_depths, axis=0)
            pos_arr = np.arange(rstart, rend + 1)
            ax.fill_between(pos_arr, mean_depth, alpha=0.3, color=ginfo['color'])
            ax.plot(pos_arr, mean_depth, linewidth=0.5, color=ginfo['color'])

        # CE highlight
        ce = t['ce']
        ax.axvspan(ce[0], ce[1], alpha=0.15, color='red')

        # Draw junction arcs
        for (jstart, jend), count in sorted(junction_counts.items(), key=lambda x: -x[1]):
            if count < 2:
                continue
            if jstart < rstart or jend > rend:
                continue

            # Classify junction
            is_ce = False
            for ce_s, ce_e in t['ce_junctions']:
                if abs(jstart - ce_s) <= TOLERANCE and abs(jend - ce_e) <= TOLERANCE:
                    is_ce = True
                    break

            color = '#e74c3c' if is_ce else '#7f8c8d'
            lw = 2 if is_ce else 1
            alpha = 0.8 if is_ce else 0.4

            mid = (jstart + jend) / 2
            span = jend - jstart
            height = min(span * 0.15, max(mean_depth) * 0.6) if len(all_depths) > 0 else span * 0.1
            arc_x = np.linspace(jstart, jend, 50)
            arc_y = height * np.sin(np.linspace(0, np.pi, 50)) + max(mean_depth) * 0.05
            ax.plot(arc_x, arc_y + max(mean_depth) * 0.1, color=color, lw=lw, alpha=alpha)

            mean_count = count / len(ginfo['samples'])
            ax.text(mid, max(arc_y) + max(mean_depth) * 0.15, f'{mean_count:.1f}',
                    ha='center', fontsize=7, color=color)

        n_samples = len(ginfo['samples'])
        # Count CE-positive samples
        ce_pos = 0
        for s in ginfo['samples']:
            srr = s['SRR']
            bam = f'{PROJ}/mapping/{srr}/{srr}.sorted.bam'
            if not os.path.exists(bam):
                continue
            cmd = f'{SAMTOOLS} view -q 255 {bam} {region}'
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            has_ce = False
            for line in proc.stdout.strip().split('\n'):
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) < 6:
                    continue
                pos = int(fields[3])
                cigar = fields[5]
                current_pos = pos
                for m_obj in re.finditer(r'(\d+)([MDNISH=X])', cigar):
                    length = int(m_obj.group(1))
                    op = m_obj.group(2)
                    if op == 'N':
                        for ce_s, ce_e in t['ce_junctions']:
                            if abs(current_pos - ce_s) <= TOLERANCE and abs(current_pos + length - ce_e) <= TOLERANCE:
                                has_ce = True
                        current_pos += length
                    elif op in ('M', 'D', '=', 'X'):
                        current_pos += length
                if has_ce:
                    break
            if has_ce:
                ce_pos += 1

        ax.text(0.98, 0.95, f'Cryptic junction: {ce_pos}/{n_samples} samples',
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                color='#e74c3c',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_ylabel(f'{ginfo["label"]}\n(n={n_samples})')
        ax.set_xlim(rstart, rend)

    # Add legend
    axes[0].plot([], [], color='#7f8c8d', lw=1, label='Normal junction')
    axes[0].plot([], [], color='#e74c3c', lw=2, label='Cryptic junction')
    axes[0].legend(loc='upper left', fontsize=8)

    axes[-1].set_xlabel(f'Genomic position ({chrom})')
    fig.suptitle(f'{gene} Sashimi Plot — Sorted Nuclei (GSE126543)\n'
                 f'TDP-43+ vs TDP-43− neuronal nuclei, mean junction reads per sample',
                 fontsize=13)
    plt.tight_layout()
    fig.savefig(f'{FIGDIR}/sashimi_{gene.lower()}_sorted.png', dpi=150)
    fig.savefig(f'{FIGDIR}/sashimi_{gene.lower()}_sorted.pdf')
    plt.close()
    print(f'  Saved: sashimi_{gene.lower()}_sorted.png')

print('\nAll plots saved to:', FIGDIR)

PLOT_PYEOF

${NOTIFY} "[DHX8-Sakai] GSE126543 Step 6/6: Coverage + sashimi plots complete — ALL DONE"

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "Results: ${PROJ}/analysis/"
echo "Figures: ${PROJ}/analysis/figures/"
