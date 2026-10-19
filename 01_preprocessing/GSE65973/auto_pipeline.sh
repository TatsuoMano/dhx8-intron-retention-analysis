#!/bin/bash
set -e

run_pipeline() {
    DATASET=$1
    WORKDIR=/mnt/HDD2/TDP43_comparison/${DATASET}
    EXPECTED=$2
    
    echo "$(date) === Waiting for ${DATASET} downloads (${EXPECTED} samples) ==="
    
    # Wait for all FASTQs
    while true; do
        COUNT=$(ls $WORKDIR/fastq/*.fastq.gz 2>/dev/null | wc -l)
        if [ "$COUNT" -ge "$EXPECTED" ]; then
            echo "$(date) All ${COUNT} FASTQs ready for ${DATASET}"
            break
        fi
        echo "$(date) ${DATASET}: ${COUNT}/${EXPECTED} FASTQs ready, waiting..."
        sleep 120
    done
    
    # Run STAR mapping
    echo "$(date) === Starting STAR mapping for ${DATASET} ==="
    bash /mnt/HDD2/TDP43_comparison/02_star_mapping.sh ${DATASET} 2>&1 | tee $WORKDIR/logs/star_stdout.log
    
    # Run featureCounts for intron/exon counting
    echo "$(date) === Running featureCounts for ${DATASET} ==="
    GTF=/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf
    mkdir -p $WORKDIR/counts
    
    # Get all BAM files by condition
    CTRL_BAMS=""
    KD_BAMS=""
    while IFS=$'\t' read -r srr gsm name cond; do
        BAM=$WORKDIR/mapping/${srr}_Aligned.sortedByCoord.out.bam
        if [ -f "$BAM" ]; then
            if [ "$cond" = "control" ]; then
                CTRL_BAMS="$CTRL_BAMS $BAM"
            else
                KD_BAMS="$KD_BAMS $BAM"
            fi
        fi
    done < $WORKDIR/sample_sheet.tsv
    
    # Gene-level counts
    echo "$(date) Running featureCounts (gene level)..."
    featureCounts -T 24 -a $GTF -o $WORKDIR/counts/gene_counts.txt $CTRL_BAMS $KD_BAMS 2>> $WORKDIR/logs/featurecounts.log
    
    # Intron counts (using the same approach as DHX8)
    echo "$(date) Running intron-level analysis..."
    
    echo "$(date) === ${DATASET} pipeline complete ==="
}

# Run GSE27218 (21 samples)
run_pipeline GSE27218 21 &
PID1=$!

# Run GSE65973 (8 samples)  
run_pipeline GSE65973 8 &
PID2=$!

# Wait for both
wait $PID1 $PID2
echo "$(date) === All pipelines complete ==="

# After both complete, run comparison with DHX8
echo "$(date) === Running DHX8 vs TDP-43 IR comparison ==="
python3 /mnt/HDD2/TDP43_comparison/03_ir_analysis.py /mnt/HDD2/TDP43_comparison/GSE27218
python3 /mnt/HDD2/TDP43_comparison/03_ir_analysis.py /mnt/HDD2/TDP43_comparison/GSE65973

echo "$(date) === ALL DONE ==="
