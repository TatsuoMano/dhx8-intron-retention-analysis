#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rmats

BAMDIR=/mnt/HDD2/unified_results/encode_rbp_kd/bams
OUTBASE=/mnt/HDD2/unified_results/encode_rbp_kd/splicing
GTF=/mnt/NAS/database/Ensembl/GRCh38_annotation/Homo_sapiens.GRCh38.109.gtf

# RBFOX2
echo "[5/6] Running rMATS for RBFOX2..."
~/bin/notify_telegram.sh "[DHX8-Sakai] Step 5/6: rMATS RBFOX2 started"
mkdir -p $OUTBASE/rmats_RBFOX2/tmp
echo "$BAMDIR/RBFOX2_K562_rep1_ENCFF762GIC.bam,$BAMDIR/RBFOX2_K562_rep2_ENCFF702GQW.bam" > $OUTBASE/rmats_RBFOX2/b1.txt
echo "$BAMDIR/Control_K562_ENCSR667PLJ_rep1_ENCFF214SQV.bam,$BAMDIR/Control_K562_ENCSR667PLJ_rep2_ENCFF125PXK.bam" > $OUTBASE/rmats_RBFOX2/b2.txt

rmats.py --b1 $OUTBASE/rmats_RBFOX2/b1.txt --b2 $OUTBASE/rmats_RBFOX2/b2.txt \
  --gtf $GTF -t paired --readLength 100 --nthread 32 \
  --od $OUTBASE/rmats_RBFOX2 --tmp $OUTBASE/rmats_RBFOX2/tmp \
  --novelSS --mil 50 --mel 500 2>&1 | tee $OUTBASE/rmats_RBFOX2/rmats.log

echo "[5/6] RBFOX2: rMATS completed"
~/bin/notify_telegram.sh "[DHX8-Sakai] Step 5/6: rMATS RBFOX2 complete"

# RBM39
echo "[6/6] Running rMATS for RBM39..."
~/bin/notify_telegram.sh "[DHX8-Sakai] Step 6/6: rMATS RBM39 started"
mkdir -p $OUTBASE/rmats_RBM39/tmp
echo "$BAMDIR/RBM39_K562_rep1_ENCFF648RHR.bam,$BAMDIR/RBM39_K562_rep2_ENCFF333UZU.bam" > $OUTBASE/rmats_RBM39/b1.txt
echo "$BAMDIR/Control_K562_ENCSR164MUK_rep1_ENCFF393WNP.bam,$BAMDIR/Control_K562_ENCSR164MUK_rep2_ENCFF384QBC.bam" > $OUTBASE/rmats_RBM39/b2.txt

rmats.py --b1 $OUTBASE/rmats_RBM39/b1.txt --b2 $OUTBASE/rmats_RBM39/b2.txt \
  --gtf $GTF -t paired --readLength 100 --nthread 32 \
  --od $OUTBASE/rmats_RBM39 --tmp $OUTBASE/rmats_RBM39/tmp \
  --novelSS --mil 50 --mel 500 2>&1 | tee $OUTBASE/rmats_RBM39/rmats.log

echo "[6/6] RBM39: rMATS completed"
~/bin/notify_telegram.sh "[DHX8-Sakai] rMATS ALL 6 EXPERIMENTS COMPLETE"
