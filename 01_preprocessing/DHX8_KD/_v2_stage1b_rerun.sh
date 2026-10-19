#!/bin/bash
# Stage 1b: resume from STAR (Trimmomatic already done, full ~15M). Fix: ulimit -n.
set -uo pipefail
cd /mnt/HDD2/Sakai_RNA-seq
ulimit -n 8192
NOTIFY=~/bin/notify_telegram.sh
GTF=/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf
STAR=/usr/lib/rna-star/bin/STAR-avx2
MAP=output/mapping
BAMS_SORTED="$MAP/siCont_1_sorted.bam $MAP/siCont_2_sorted.bam $MAP/siCont_3_sorted.bam $MAP/siDHX8_1_sorted.bam $MAP/siDHX8_2_sorted.bam $MAP/siDHX8_3_sorted.bam"
L=/mnt/HDD2/Sakai_RNA-seq/_v2_stage1.log
exec >> >(tee -a "$L") 2>&1
echo "$(date) === Stage1b resume from STAR (ulimit -n $(ulimit -n)) ==="
$NOTIFY "[DHX8-Sakai] Stage1b: STAR から再開 (ulimit fix)"

# cleanup failed STAR temp + partial outputs
rm -rf $MAP/siCont_2__STARtmp 2>/dev/null || true
rm -f $MAP/siCont_2_Aligned.sortedByCoord.out.bam 2>/dev/null || true

echo "$(date) [2/4] STAR siCont_2"
$STAR --runThreadN 24 --genomeDir /mnt/NAS/database/STAR_index/mm10 \
  --readFilesIn output/trimmed/siCont_2_1_trimmed.fq.gz output/trimmed/siCont_2_2_trimmed.fq.gz \
  --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat \
  --outFileNamePrefix $MAP/siCont_2_ || { $NOTIFY "[DHX8-Sakai] Stage1b STAR FAILED"; exit 1; }
INPUT=$(grep "Number of input reads" $MAP/siCont_2_Log.final.out | awk '{print $NF}')
UNIQ=$(grep "Uniquely mapped reads %" $MAP/siCont_2_Log.final.out | awk '{print $NF}')
echo "  siCont_2 STAR input=$INPUT (was 5,291,700), uniquely-mapped=$UNIQ"
samtools sort -@ 32 $MAP/siCont_2_Aligned.sortedByCoord.out.bam -o $MAP/siCont_2_sorted.bam
samtools index -@ 24 $MAP/siCont_2_sorted.bam
samtools index -@ 24 $MAP/siCont_2_Aligned.sortedByCoord.out.bam
$NOTIFY "[DHX8-Sakai] Stage1b 2/4: STAR+sort done, siCont_2 input=$INPUT uniq=$UNIQ"

echo "$(date) [3/4] featureCounts -s2 gene"
featureCounts -T 8 -p --countReadPairs -s 2 -a $GTF -o counts_s2/gene_counts.txt $BAMS_SORTED 2>&1 | tee counts_s2/gene_counts.log | grep -iE "Assigned|Successfully" | tail -2
echo "=== SELF-VERIFY (Assigned) ==="
echo "  現行: $(grep '^Assigned' counts_s2/gene_counts.txt.summary | cut -f2-)"
echo "  v1  : 14085485  4649282  14050402  14802246  14138815  12846628"
NEW2=$(grep '^Assigned' counts_s2/gene_counts.txt.summary | cut -f3)
OTHERS_NEW=$(grep '^Assigned' counts_s2/gene_counts.txt.summary | awk -F'\t' '{print $2, $4, $5, $6, $7}')
OTHERS_V1="14085485 14050402 14802246 14138815 12846628"
if [ "$OTHERS_NEW" = "$OTHERS_V1" ]; then MATCH="YES ✅"; else MATCH="NO ⚠️ ($OTHERS_NEW)"; fi
echo "  他5サンプル v1一致: $MATCH"
$NOTIFY "[DHX8-Sakai] Stage1b 3/4: featureCounts done. siCont_2 Assigned $NEW2 (was 4,649,282). others-match=$MATCH"

echo "$(date) [4/4] DESeq2"
Rscript /mnt/HDD2/unified_results/deseq2_mouse.R 2>&1 | tail -15 || { $NOTIFY "[DHX8-Sakai] Stage1b DESeq2 FAILED"; exit 1; }
echo "$(date) === Stage1b COMPLETE ==="
$NOTIFY "[DHX8-Sakai] Stage1 COMPLETE: DESeq2 v2 done (siCont_2 full)"
