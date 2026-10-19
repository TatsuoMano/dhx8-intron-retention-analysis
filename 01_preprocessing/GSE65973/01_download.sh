#!/bin/bash
set -e
WORKDIR=/mnt/HDD2/TDP43_comparison/GSE65973
cd $WORKDIR
while IFS=$'\t' read -r srr gsm name cond; do
    if [ ! -f fastq/${srr}.fastq.gz ]; then
        echo "$(date) Downloading ${srr} (${name})..."
        fasterq-dump --split-files --threads 8 -O fastq/ ${srr} 2>> logs/download.log
        gzip fastq/${srr}.fastq 2>/dev/null || gzip fastq/${srr}_1.fastq 2>/dev/null || true
        echo "$(date) Done: ${srr}"
    else
        echo "$(date) Skip ${srr} (already exists)"
    fi
done < sample_sheet.tsv
echo "=== All GSE65973 downloads complete ==="
