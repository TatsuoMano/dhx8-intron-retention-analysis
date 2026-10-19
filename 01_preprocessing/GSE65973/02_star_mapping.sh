#!/bin/bash
set -e
ulimit -n 4096

STAR_INDEX=/mnt/NAS/database/STAR_index/mm10
DATASET=$1  # GSE27218 or GSE65973
WORKDIR=/mnt/HDD2/TDP43_comparison/${DATASET}
mkdir -p $WORKDIR/mapping $WORKDIR/logs

echo "=== STAR mapping for ${DATASET} ==="

while IFS=$'\t' read -r srr gsm name cond; do
    BAM=$WORKDIR/mapping/${srr}_Aligned.sortedByCoord.out.bam
    if [ -f "$BAM" ] && [ -s "$BAM" ]; then
        echo "$(date) Skip ${srr} (BAM exists)"
        continue
    fi
    
    # Find FASTQ (may be .fastq.gz or _1.fastq.gz for SE)
    FQ=""
    for f in $WORKDIR/fastq/${srr}.fastq.gz $WORKDIR/fastq/${srr}_1.fastq.gz; do
        [ -f "$f" ] && FQ="$f" && break
    done
    
    if [ -z "$FQ" ]; then
        echo "$(date) WARN: No FASTQ for ${srr}, skipping"
        continue
    fi
    
    echo "$(date) Mapping ${srr} (${name})..."
    STAR --runThreadN 24 \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${FQ} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix $WORKDIR/mapping/${srr}_ \
         2>> $WORKDIR/logs/star.log
    
    # Index BAM
    samtools index $WORKDIR/mapping/${srr}_Aligned.sortedByCoord.out.bam
    echo "$(date) Done: ${srr}"
done < $WORKDIR/sample_sheet.tsv

echo "=== STAR mapping complete for ${DATASET} ==="
