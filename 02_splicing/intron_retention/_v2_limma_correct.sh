#!/bin/bash
# Stage 2 fix: regenerate counts_s2 intron/exon (-s2 SAF, new BAMs) then run the MANUSCRIPT limma (run_limma_strandfix_v2.R).
set -uo pipefail
cd /mnt/HDD2/Sakai_RNA-seq
NOTIFY=~/bin/notify_telegram.sh
MAP=output/mapping
SAFDIR=/mnt/HDD2/TDP43_comparison/event_level_20260225
LIMMADIR=/mnt/HDD2/TDP43_comparison/validation_20260226
BAMS="$MAP/siCont_1_sorted.bam $MAP/siCont_2_sorted.bam $MAP/siCont_3_sorted.bam $MAP/siDHX8_1_sorted.bam $MAP/siDHX8_2_sorted.bam $MAP/siDHX8_3_sorted.bam"
echo "$(date) === counts_s2 intron/exon -s2 再生成 (new BAMs) ==="
$NOTIFY "[DHX8-Sakai] Stage2 limma-fix: counts_s2 intron/exon -s2 再生成"
featureCounts -F SAF -p --countReadPairs -s 2 -a $SAFDIR/introns.saf -o counts_s2/intron_counts.txt $BAMS 2>&1 | tee counts_s2/intron_counts.log | grep -iE "Assigned|Successfully" | tail -1
featureCounts -F SAF -p --countReadPairs -s 2 -a $SAFDIR/flanking_exons.saf -o counts_s2/exon_counts.txt $BAMS 2>&1 | tee counts_s2/exon_counts.log | grep -iE "Assigned|Successfully" | tail -1
echo "=== SELF-VERIFY intron Assigned (他5サンプルが v1 と一致すべき) ==="
echo "  v2: $(grep '^Assigned' counts_s2/intron_counts.txt.summary | cut -f2-)"
echo "  v1: $(grep '^Assigned' /mnt/HDD2/_v1_corrupt_snapshot_20260706/counts_s2/intron_counts.txt.summary | cut -f2-)"
NEW_OTH=$(grep '^Assigned' counts_s2/intron_counts.txt.summary | awk -F"\t" '{print $2,$4,$5,$6,$7}')
V1_OTH=$(grep '^Assigned' /mnt/HDD2/_v1_corrupt_snapshot_20260706/counts_s2/intron_counts.txt.summary | awk -F"\t" '{print $2,$4,$5,$6,$7}')
[ "$NEW_OTH" = "$V1_OTH" ] && echo "  ✅ 他5一致" || echo "  ⚠️ 他5不一致"
echo "$(date) === run_limma_strandfix_v2.R (manuscript limma) 実行 ==="
cd $LIMMADIR && Rscript run_limma_strandfix_v2.R 2>&1 | grep -iE "retained|genes|dhx8|s2|Retained" | tail -12 || { $NOTIFY "[DHX8-Sakai] limma-fix R FAILED"; exit 1; }
V2RET=$(($(wc -l < $LIMMADIR/limma_retained_dhx8_s2.tsv 2>/dev/null)-1))
echo "=== v2 _s2 retained introns = $V2RET (v1 manuscript=1233) ==="
# redeploy canonical copies (図/SG が読む)
for d in /mnt/HDD2/unified_results/sg_analysis/gsea_dhx8/data /mnt/HDD2/unified_results/sg_transcriptome_s2_reanalysis/data; do
  cp $LIMMADIR/limma_retained_dhx8_s2.tsv $d/limma_retained_dhx8.tsv 2>/dev/null && echo "  redeploy → $d"
  cp $LIMMADIR/limma_ir_results_dhx8_s2.tsv $d/limma_ir_results_dhx8.tsv 2>/dev/null
done
$NOTIFY "[DHX8-Sakai] Stage2 limma-fix DONE: v2 retained=$V2RET (v1=1233)"
echo "$(date) === limma-fix DONE ==="
