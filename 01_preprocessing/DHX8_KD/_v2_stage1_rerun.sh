#!/bin/bash
# Stage 1 v2 re-run (intact siCont_2): Trimmomatic -> STAR -> sort(both BAM names) -> featureCounts -s2 gene -> DESeq2
# Non-destructive to other 5 samples; self-verifies their counts == v1.
set -uo pipefail
cd /mnt/HDD2/Sakai_RNA-seq
NOTIFY=~/bin/notify_telegram.sh
GTF=/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf
STAR=/usr/lib/rna-star/bin/STAR-avx2
TRIM=/bin/Trimmomatic-0.39/trimmomatic-0.39.jar
RAW=/mnt/NAS/rawdata/RH23097820/rawdata
MAP=output/mapping
BAMS_SORTED="$MAP/siCont_1_sorted.bam $MAP/siCont_2_sorted.bam $MAP/siCont_3_sorted.bam $MAP/siDHX8_1_sorted.bam $MAP/siDHX8_2_sorted.bam $MAP/siDHX8_3_sorted.bam"
L=/mnt/HDD2/Sakai_RNA-seq/_v2_stage1.log
exec > >(tee -a "$L") 2>&1
echo "$(date) === Stage1 v2 re-run START ==="
$NOTIFY "[DHX8-Sakai] Stage1/3 v2 re-run START (siCont_2 full ~15.3M)"

echo "$(date) [1/4] Trimmomatic siCont_2 (intact)"
java -jar $TRIM PE -threads 24 -phred33 \
  $RAW/siCont_2_1.fq.gz $RAW/siCont_2_2.fq.gz \
  output/trimmed/siCont_2_1_trimmed.fq.gz output/trimmed/siCont_2_1_unpaired.fq.gz \
  output/trimmed/siCont_2_2_trimmed.fq.gz output/trimmed/siCont_2_2_unpaired.fq.gz \
  ILLUMINACLIP:/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:10:10 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30 2>&1 | grep -iE "Input Read Pairs|Completed" || true
$NOTIFY "[DHX8-Sakai] Stage1 1/4: Trimmomatic siCont_2 done"

echo "$(date) [2/4] STAR siCont_2 + sort both BAM names"
$STAR --runThreadN 24 --genomeDir /mnt/NAS/database/STAR_index/mm10 \
  --readFilesIn output/trimmed/siCont_2_1_trimmed.fq.gz output/trimmed/siCont_2_2_trimmed.fq.gz \
  --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat \
  --outFileNamePrefix $MAP/siCont_2_ || { $NOTIFY "[DHX8-Sakai] Stage1 STAR FAILED"; exit 1; }
INPUT=$(grep "Number of input reads" $MAP/siCont_2_Log.final.out | awk '{print $NF}')
echo "  siCont_2 STAR input reads = $INPUT (was 5,291,700)"
samtools sort -@ 32 $MAP/siCont_2_Aligned.sortedByCoord.out.bam -o $MAP/siCont_2_sorted.bam
samtools index -@ 24 $MAP/siCont_2_sorted.bam
samtools index -@ 24 $MAP/siCont_2_Aligned.sortedByCoord.out.bam
$NOTIFY "[DHX8-Sakai] Stage1 2/4: STAR+sort done, siCont_2 input reads=$INPUT (was 5.29M)"

echo "$(date) [3/4] featureCounts -s2 gene (6 sorted BAMs)"
featureCounts -T 8 -p --countReadPairs -s 2 -a $GTF -o counts_s2/gene_counts.txt $BAMS_SORTED 2>&1 | tee counts_s2/gene_counts.log | grep -iE "Process|Assigned|Successfully" | tail -3
echo "=== SELF-VERIFY (Assigned): 他5サンプルが v1 と一致すべき ==="
echo "  現行: $(grep '^Assigned' counts_s2/gene_counts.txt.summary | cut -f2-)"
echo "  v1  : 14085485  4649282  14050402  14802246  14138815  12846628"
NEW2=$(grep '^Assigned' counts_s2/gene_counts.txt.summary | cut -f3)
# check others (col 2,4,5,6,7) unchanged
OTHERS_NEW=$(grep '^Assigned' counts_s2/gene_counts.txt.summary | awk -F'\t' '{print $2, $4, $5, $6, $7}')
OTHERS_V1="14085485 14050402 14802246 14138815 12846628"
if [ "$OTHERS_NEW" = "$OTHERS_V1" ]; then echo "  ✅ 他5サンプル一致 → featureCounts コマンド正当"; else echo "  ⚠️ 他5サンプル不一致 (new=$OTHERS_NEW) → 要確認"; fi
$NOTIFY "[DHX8-Sakai] Stage1 3/4: featureCounts -s2 gene done. siCont_2 Assigned $NEW2 (was 4,649,282). others-match=$([ "$OTHERS_NEW" = "$OTHERS_V1" ] && echo YES || echo NO)"

echo "$(date) [4/4] DESeq2 (deseq2_mouse.R)"
Rscript /mnt/HDD2/unified_results/deseq2_mouse.R 2>&1 | tail -15 || { $NOTIFY "[DHX8-Sakai] Stage1 DESeq2 FAILED"; exit 1; }
echo "$(date) === Stage1 v2 re-run COMPLETE ==="
$NOTIFY "[DHX8-Sakai] Stage1 4/4 COMPLETE: DESeq2 re-run done (siCont_2 full). Stage2待ち."
