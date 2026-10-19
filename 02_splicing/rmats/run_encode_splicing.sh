#!/bin/bash
# ENCODE RBP KD splicing analysis pipeline
# rMATS + IR quantification + CE detection
# Runs after BAM download completes
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rmats

set -euo pipefail

BAMDIR=/mnt/HDD2/unified_results/encode_rbp_kd/bams
OUTBASE=/mnt/HDD2/unified_results/encode_rbp_kd/splicing
RMATS=/home/tatsuomano-linux2/miniconda3/envs/rmats/bin/rmats.py
GTF=/mnt/NAS/database/Ensembl/GRCh38_annotation/Homo_sapiens.GRCh38.109.gtf
STAR_INDEX=/mnt/NAS/database/STAR_index/hg38
THREADS=32

notify() { ~/bin/notify_telegram.sh "$1"; }

mkdir -p $OUTBASE

# ============================================================
# Wait for BAM downloads to complete
# ============================================================
echo "Checking BAM downloads..."
EXPECTED=24
while true; do
    COMPLETED=$(ls $BAMDIR/*.bam 2>/dev/null | wc -l)
    if [ "$COMPLETED" -ge "$EXPECTED" ]; then
        echo "All $EXPECTED BAMs downloaded"
        break
    fi
    echo "$(date): $COMPLETED/$EXPECTED BAMs downloaded, waiting..."
    sleep 300
done

notify "[DHX8-Sakai] ENCODE splicing analysis started (6 KD experiments)"

# ============================================================
# Define KD-Control pairs from metadata
# ============================================================
declare -A KD_BAM1 KD_BAM2 CTRL_BAM1 CTRL_BAM2

# PRPF8
KD_BAM1[PRPF8]=$BAMDIR/PRPF8_HepG2_rep1_ENCFF673NPW.bam
KD_BAM2[PRPF8]=$BAMDIR/PRPF8_HepG2_rep2_ENCFF124KBD.bam
CTRL_BAM1[PRPF8]=$BAMDIR/Control_HepG2_ENCSR856ZRV_rep1_ENCFF528NWG.bam
CTRL_BAM2[PRPF8]=$BAMDIR/Control_HepG2_ENCSR856ZRV_rep2_ENCFF317QYL.bam

# SNRNP200
KD_BAM1[SNRNP200]=$BAMDIR/SNRNP200_HepG2_rep1_ENCFF287XNM.bam
KD_BAM2[SNRNP200]=$BAMDIR/SNRNP200_HepG2_rep2_ENCFF433UMI.bam
CTRL_BAM1[SNRNP200]=$BAMDIR/Control_HepG2_ENCSR104ABF_rep1_ENCFF710YBQ.bam
CTRL_BAM2[SNRNP200]=$BAMDIR/Control_HepG2_ENCSR104ABF_rep2_ENCFF834QZT.bam

# EFTUD2
KD_BAM1[EFTUD2]=$BAMDIR/EFTUD2_HepG2_rep1_ENCFF473RIO.bam
KD_BAM2[EFTUD2]=$BAMDIR/EFTUD2_HepG2_rep2_ENCFF009UNT.bam
CTRL_BAM1[EFTUD2]=$BAMDIR/Control_HepG2_ENCSR689PHN_rep1_ENCFF862DEF.bam
CTRL_BAM2[EFTUD2]=$BAMDIR/Control_HepG2_ENCSR689PHN_rep2_ENCFF178XWW.bam

# PRPF6
KD_BAM1[PRPF6]=$BAMDIR/PRPF6_HepG2_rep1_ENCFF283DDL.bam
KD_BAM2[PRPF6]=$BAMDIR/PRPF6_HepG2_rep2_ENCFF966GXE.bam
CTRL_BAM1[PRPF6]=$BAMDIR/Control_HepG2_ENCSR491FOC_rep1_ENCFF492EZP.bam
CTRL_BAM2[PRPF6]=$BAMDIR/Control_HepG2_ENCSR491FOC_rep2_ENCFF970TMY.bam

# RBFOX2
KD_BAM1[RBFOX2]=$BAMDIR/RBFOX2_K562_rep1_ENCFF762GIC.bam
KD_BAM2[RBFOX2]=$BAMDIR/RBFOX2_K562_rep2_ENCFF702GQW.bam
CTRL_BAM1[RBFOX2]=$BAMDIR/Control_K562_ENCSR667PLJ_rep1_ENCFF214SQV.bam
CTRL_BAM2[RBFOX2]=$BAMDIR/Control_K562_ENCSR667PLJ_rep2_ENCFF125PXK.bam

# RBM39
KD_BAM1[RBM39]=$BAMDIR/RBM39_K562_rep1_ENCFF648RHR.bam
KD_BAM2[RBM39]=$BAMDIR/RBM39_K562_rep2_ENCFF333UZU.bam
CTRL_BAM1[RBM39]=$BAMDIR/Control_K562_ENCSR164MUK_rep1_ENCFF393WNP.bam
CTRL_BAM2[RBM39]=$BAMDIR/Control_K562_ENCSR164MUK_rep2_ENCFF384QBC.bam

# ============================================================
# Run rMATS for each KD experiment
# ============================================================
EXPERIMENTS="PRPF8 SNRNP200 EFTUD2 PRPF6 RBFOX2 RBM39"
TOTAL=6
COUNT=0

for EXP in $EXPERIMENTS; do
    COUNT=$((COUNT + 1))
    OUTDIR=$OUTBASE/rmats_$EXP
    mkdir -p $OUTDIR
    TMPDIR=$OUTDIR/tmp
    mkdir -p $TMPDIR

    # Check if already completed
    if [ -f "$OUTDIR/SE.MATS.JCEC.txt" ]; then
        echo "[$COUNT/$TOTAL] $EXP: already completed, skipping"
        continue
    fi

    echo "[$COUNT/$TOTAL] Running rMATS for $EXP..."
    notify "[DHX8-Sakai] Step $COUNT/$TOTAL: rMATS $EXP started"

    # Verify BAM files exist
    for BAM in ${KD_BAM1[$EXP]} ${KD_BAM2[$EXP]} ${CTRL_BAM1[$EXP]} ${CTRL_BAM2[$EXP]}; do
        if [ ! -f "$BAM" ]; then
            echo "ERROR: Missing BAM: $BAM"
            notify "[DHX8-Sakai] ERROR: Missing BAM for $EXP: $BAM"
            continue 2
        fi
    done

    # Create BAM list files
    echo "${KD_BAM1[$EXP]},${KD_BAM2[$EXP]}" > $OUTDIR/b1.txt
    echo "${CTRL_BAM1[$EXP]},${CTRL_BAM2[$EXP]}" > $OUTDIR/b2.txt

    # Run rMATS
    python3 $RMATS \
        --b1 $OUTDIR/b1.txt \
        --b2 $OUTDIR/b2.txt \
        --gtf $GTF \
        -t paired \
        --readLength 100 \
        --nthread $THREADS \
        --od $OUTDIR \
        --tmp $TMPDIR \
        --novelSS \
        --mil 50 --mel 500 \
        2>&1 | tee $OUTDIR/rmats.log

    if [ -f "$OUTDIR/SE.MATS.JCEC.txt" ]; then
        echo "[$COUNT/$TOTAL] $EXP: rMATS completed successfully"
        notify "[DHX8-Sakai] Step $COUNT/$TOTAL: rMATS $EXP complete"
    else
        echo "[$COUNT/$TOTAL] $EXP: rMATS FAILED"
        notify "[DHX8-Sakai] ERROR: rMATS $EXP failed"
    fi
done

# ============================================================
# Summarize results
# ============================================================
notify "[DHX8-Sakai] rMATS analysis complete for all 6 experiments. Starting summary..."

echo ""
echo "============================================================"
echo "rMATS Summary: Event type distribution"
echo "============================================================"

echo -e "Experiment\tSE\tRI\tA3SS\tA5SS\tMXE\tTotal\tRI_pct\tnovelSE" > $OUTBASE/rmats_summary.tsv

for EXP in $EXPERIMENTS; do
    OUTDIR=$OUTBASE/rmats_$EXP
    if [ ! -f "$OUTDIR/SE.MATS.JCEC.txt" ]; then
        echo "$EXP: no results"
        continue
    fi

    # Count significant events (FDR < 0.05, |IncLevelDifference| > 0.1)
    SE=$(awk -F'\t' 'NR>1 && $22<0.05 && ($23>0.1 || $23<-0.1)' $OUTDIR/SE.MATS.JCEC.txt | wc -l)
    RI=$(awk -F'\t' 'NR>1 && $22<0.05 && ($23>0.1 || $23<-0.1)' $OUTDIR/RI.MATS.JCEC.txt | wc -l)
    A3SS=$(awk -F'\t' 'NR>1 && $22<0.05 && ($23>0.1 || $23<-0.1)' $OUTDIR/A3SS.MATS.JCEC.txt | wc -l)
    A5SS=$(awk -F'\t' 'NR>1 && $22<0.05 && ($23>0.1 || $23<-0.1)' $OUTDIR/A5SS.MATS.JCEC.txt | wc -l)
    MXE=$(awk -F'\t' 'NR>1 && $22<0.05 && ($23>0.1 || $23<-0.1)' $OUTDIR/MXE.MATS.JCEC.txt | wc -l)
    TOTAL=$((SE + RI + A3SS + A5SS + MXE))

    # RI percentage
    if [ "$TOTAL" -gt 0 ]; then
        RI_PCT=$(echo "scale=1; $RI * 100 / $TOTAL" | bc)
    else
        RI_PCT=0
    fi

    # Novel SE events (from novelJunction files)
    NOVEL_SE=0
    if [ -f "$OUTDIR/fromGTF.novelJunction.SE.txt" ]; then
        NOVEL_SE=$(awk -F'\t' 'NR>1' $OUTDIR/fromGTF.novelJunction.SE.txt | wc -l)
    fi

    echo -e "$EXP\t$SE\t$RI\t$A3SS\t$A5SS\t$MXE\t$TOTAL\t${RI_PCT}%\t$NOVEL_SE"
    echo -e "$EXP\t$SE\t$RI\t$A3SS\t$A5SS\t$MXE\t$TOTAL\t${RI_PCT}%\t$NOVEL_SE" >> $OUTBASE/rmats_summary.tsv
done

echo ""
echo "Summary saved to $OUTBASE/rmats_summary.tsv"

notify "[DHX8-Sakai] ENCODE splicing analysis ALL COMPLETE. Results in $OUTBASE/"
