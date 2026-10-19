#!/bin/bash
set -euo pipefail

# DHX8 KD strandedness fix: featureCounts -s 2, then original step2 + limma
# Changes ONLY: -s 0 → -s 2 in featureCounts. Everything else is the original scripts.

PREV_DIR=/mnt/HDD2/TDP43_comparison/event_level_20260225
LIMMA_DIR=/mnt/HDD2/TDP43_comparison/validation_20260226
NOTIFY=~/bin/notify_telegram.sh

BAMS="/mnt/HDD2/Sakai_RNA-seq/output/mapping/siCont_1_Aligned.sortedByCoord.out.bam \
/mnt/HDD2/Sakai_RNA-seq/output/mapping/siCont_2_Aligned.sortedByCoord.out.bam \
/mnt/HDD2/Sakai_RNA-seq/output/mapping/siCont_3_Aligned.sortedByCoord.out.bam \
/mnt/HDD2/Sakai_RNA-seq/output/mapping/siDHX8_1_Aligned.sortedByCoord.out.bam \
/mnt/HDD2/Sakai_RNA-seq/output/mapping/siDHX8_2_Aligned.sortedByCoord.out.bam \
/mnt/HDD2/Sakai_RNA-seq/output/mapping/siDHX8_3_Aligned.sortedByCoord.out.bam"

echo "$(date) === DHX8 KD Strandfix Pipeline ==="

# Step 0: Backup old counts
echo "$(date) Backing up old counts..."
mkdir -p ${PREV_DIR}/counts_s0_backup
cp -n ${PREV_DIR}/dhx8_intron_counts.txt ${PREV_DIR}/counts_s0_backup/ 2>/dev/null || true
cp -n ${PREV_DIR}/dhx8_intron_counts.txt.summary ${PREV_DIR}/counts_s0_backup/ 2>/dev/null || true
cp -n ${PREV_DIR}/dhx8_exon_counts.txt ${PREV_DIR}/counts_s0_backup/ 2>/dev/null || true
cp -n ${PREV_DIR}/dhx8_exon_counts.txt.summary ${PREV_DIR}/counts_s0_backup/ 2>/dev/null || true

# Step 1: featureCounts -s 2 (ONLY change from original)
echo "$(date) featureCounts intron counts (-s 2)..."
featureCounts -F SAF \
    -a ${PREV_DIR}/introns.saf \
    -o ${PREV_DIR}/dhx8_intron_counts.txt \
    -p -T 8 --fracOverlap 0 -s 2 \
    ${BAMS} 2>&1 | tee ${PREV_DIR}/dhx8_intron_counts_s2.log

echo "$(date) featureCounts exon counts (-s 2)..."
featureCounts -F SAF \
    -a ${PREV_DIR}/flanking_exons.saf \
    -o ${PREV_DIR}/dhx8_exon_counts.txt \
    -p -T 8 --fracOverlap 0 -s 2 \
    ${BAMS} 2>&1 | tee ${PREV_DIR}/dhx8_exon_counts_s2.log

echo "$(date) featureCounts done"
${NOTIFY} "[DHX8-Sakai] Strandfix Step 1/3: featureCounts -s 2 done"

# Step 2: Original step2 (Welch t-test) — unchanged
echo "$(date) Running step2 (Welch t-test)..."
cd ${PREV_DIR}
python3 ${PREV_DIR}/step2_identify_retained_introns.py 2>&1 | tee ${PREV_DIR}/step2_s2.log
echo "$(date) step2 done"
${NOTIFY} "[DHX8-Sakai] Strandfix Step 2/3: Welch t-test done"

# Step 3: Original limma — unchanged
echo "$(date) Running limma eBayes..."
cd ${LIMMA_DIR}
Rscript ${LIMMA_DIR}/validation_step1_limma.R 2>&1 | tee ${LIMMA_DIR}/limma_s2.log
echo "$(date) limma done"
${NOTIFY} "[DHX8-Sakai] Strandfix Step 3/3: limma eBayes done — ALL COMPLETE"

echo "$(date) === ALL DONE ==="
