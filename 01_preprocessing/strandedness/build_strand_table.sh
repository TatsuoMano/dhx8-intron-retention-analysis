#!/bin/bash
# Build Supplementary Table S1: featureCounts strandedness verification across datasets.
# For each dataset, run featureCounts -s 0/1/2 on ONE representative BAM and report
# the Assigned% rate. Where a 3-way comparison was already done (Sakai, GSE65973),
# parse existing .summary files instead of re-running.
set -uo pipefail   # tolerate per-fc failures; we only abort on critical setup errors

OUT_DIR=/mnt/HDD2/unified_results/strandedness_verification
mkdir -p $OUT_DIR
RAW=$OUT_DIR/raw_summaries
mkdir -p $RAW
HG38_GTF=/mnt/NAS/database/genes/hg38.ncbiRefSeq.gtf
MM10_GTF=/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf

run_fc() {
    # $1=BAM  $2=GTF  $3=PE_FLAG ("-p --countReadPairs" or "")  $4=STRAND  $5=OUT_PREFIX
    local bam=$1 gtf=$2 pe=$3 s=$4 out=$5
    if [ -s "${out}.summary" ]; then
        echo "  cached: $out.summary"
        return
    fi
    echo "  running fc -s $s on $(basename $bam) ..."
    featureCounts -a "$gtf" -o "$out" -T 8 -s "$s" $pe -B "$bam" \
        > "${out}.log" 2>&1
}

# ---- 1) Sakai DHX8 KD (already has counts_s0_paired/s1/s2 — just copy summaries) ----
for s in 0 1 2; do
    src=/mnt/HDD2/Sakai_RNA-seq/counts_s${s}*/gene_counts.txt.summary
    src=$(ls $src 2>/dev/null | head -1)
    [ -n "$src" ] && cp -n "$src" "$RAW/Sakai_DHX8_s${s}.summary"
done

# ---- 2) GSE65973 (Walker mouse cortex; s0 missing → run) ----
GSE65973_DIR=/mnt/HDD2/TDP43_comparison/GSE65973
for s in 1 2; do
    cp -n $GSE65973_DIR/counts_s${s}/gene_counts.txt.summary $RAW/GSE65973_s${s}.summary 2>/dev/null || true
done
BAM=$(find $GSE65973_DIR -name '*.bam' | head -1)
# GSE65973 BAMs were aligned as single-end → no -p flag
[ -n "$BAM" ] && run_fc "$BAM" "$MM10_GTF" "" 0 "$RAW/GSE65973_s0_fc"
[ -s "$RAW/GSE65973_s0_fc.summary" ] && mv $RAW/GSE65973_s0_fc.summary $RAW/GSE65973_s0.summary

# ---- 3) GSE124439 (cortex) — only s2 was kept; verify s0 / s1 on rep BAM ----
GSE124_BAM=$(ls /mnt/HDD2/ALS_patient_GSE124439/mapping/SRR8375274/*.bam 2>/dev/null | head -1)
cp -n /mnt/HDD2/ALS_patient_GSE124439/counts/gene_counts.txt.summary $RAW/GSE124439_s2_full.summary
for s in 0 1 2; do
    run_fc "$GSE124_BAM" "$HG38_GTF" "-p --countReadPairs" "$s" "$RAW/GSE124439_s${s}_repBAM"
done

# ---- 4) GSE116622 (cerebellum) — has s0/s2; verify s1 on rep BAM ----
for s in 0 2; do
    src=/mnt/HDD2/ALS_patient_GSE116622/counts_s${s}/gene_counts.txt.summary
    [ -e "$src" ] || src=/mnt/HDD2/ALS_patient_GSE116622/counts/gene_counts.txt.summary
    cp -n "$src" "$RAW/GSE116622_s${s}.summary" 2>/dev/null || true
done
GSE116_BAM=$(find /mnt/HDD2/ALS_patient_GSE116622/mapping -name '*.bam' | head -1)
[ -n "$GSE116_BAM" ] && run_fc "$GSE116_BAM" "$HG38_GTF" "-p --countReadPairs" 1 "$RAW/GSE116622_s1_fc"
[ -s "$RAW/GSE116622_s1_fc.summary" ] && mv $RAW/GSE116622_s1_fc.summary $RAW/GSE116622_s1.summary

# ---- 5) GSE126543 (FANS nuclei) — has s0 ("counts/") + s2_wrong; run s1 ----
cp -n /mnt/HDD2/sorted_nuclei_GSE126543/counts/gene_counts.txt.summary $RAW/GSE126543_s0.summary
cp -n /mnt/HDD2/sorted_nuclei_GSE126543/counts/gene_counts_s2_wrong.txt.summary $RAW/GSE126543_s2.summary
GSE126_BAM=$(ls /mnt/HDD2/sorted_nuclei_GSE126543/mapping/SRR8571937/*.bam 2>/dev/null | head -1)
[ -n "$GSE126_BAM" ] && run_fc "$GSE126_BAM" "$HG38_GTF" "" 1 "$RAW/GSE126543_s1_fc"
[ -s "$RAW/GSE126543_s1_fc.summary" ] && mv $RAW/GSE126543_s1_fc.summary $RAW/GSE126543_s1.summary

# ---- 6) ENCODE knockdowns (PE, hg38). Run all 3 strands on one BAM per factor ----
ENCODE_DIR=/mnt/HDD2/unified_results/encode_rbp_kd/bams
declare -A ENC_BAMS=(
    ["PRPF8_HepG2"]="PRPF8_HepG2_rep1.bam"
    ["EFTUD2_HepG2"]="EFTUD2_HepG2_rep1.bam"
    ["SNRNP200_HepG2"]="SNRNP200_HepG2_rep1.bam"
    ["PRPF6_HepG2"]="PRPF6_HepG2_rep1.bam"
    ["RBFOX2_K562"]="RBFOX2_K562_rep1.bam"
    ["RBM39_K562"]="RBM39_K562_rep1.bam"
)
for label in "${!ENC_BAMS[@]}"; do
    bam_file=$ENCODE_DIR/${ENC_BAMS[$label]}
    if [ ! -e "$bam_file" ]; then
        # try fuzzy match
        prefix=${ENC_BAMS[$label]%_rep1.bam}
        bam_file=$(ls $ENCODE_DIR/${prefix}_rep1*.bam 2>/dev/null | head -1)
    fi
    if [ -z "$bam_file" ] || [ ! -e "$bam_file" ]; then
        echo "  ENCODE BAM not found for $label, skipping"
        continue
    fi
    for s in 0 1 2; do
        run_fc "$bam_file" "$HG38_GTF" "-p --countReadPairs" "$s" "$RAW/ENCODE_${label}_s${s}_fc"
        [ -s "$RAW/ENCODE_${label}_s${s}_fc.summary" ] && \
            mv "$RAW/ENCODE_${label}_s${s}_fc.summary" "$RAW/ENCODE_${label}_s${s}.summary"
    done
done

echo "All featureCounts runs complete."
ls -1 $RAW | sort
