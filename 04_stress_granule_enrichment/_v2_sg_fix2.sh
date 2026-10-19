#!/bin/bash
# SG fix2: also join intron_length into the "all" IR file (limma_ir_results_dhx8.tsv), then re-run full SG.
set -uo pipefail
NOTIFY=~/bin/notify_telegram.sh
LIMMADIR=/mnt/HDD2/TDP43_comparison/validation_20260226
META=/mnt/HDD2/TDP43_comparison/event_level_20260225/intron_metadata.tsv
SGDIR=/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8
echo "$(date) === all(ir_results) ファイルに intron_length 追加 ==="
python3 /tmp/join_ir_meta.py "$LIMMADIR/limma_ir_results_dhx8_s2.tsv" "$META" "$LIMMADIR/limma_ir_results_dhx8_13col.tsv"
for d in $SGDIR/data /mnt/HDD2/unified_results/sg_transcriptome_s2_reanalysis/data; do
  cp "$LIMMADIR/limma_ir_results_dhx8_13col.tsv" "$d/limma_ir_results_dhx8.tsv" && echo "  redeploy(all,13col) → $d ($(head -1 $d/limma_ir_results_dhx8.tsv | tr '\t' '\n' | wc -l)列)"
done
echo "$(date) === SG 完全再実行 ==="
cd $SGDIR
python3 analysis_sg_transcriptome_gsea.py > _v2_sg_gsea_full.log 2>&1 && echo "  SG full OK" || { echo "  SG まだ失敗:"; tail -10 _v2_sg_gsea_full.log; $NOTIFY "[DHX8-Sakai] SG full re-run FAILED(2)"; exit 1; }
echo "=== 検証: IR overlap + DEG-UP overlap (両方 v2 か) ==="
echo "--- IR overlap ---"
awk -F"\t" 'NR>1 && ($1~/Khong 2017 SG RNA Enriched/||$1~/High Cons/||$1~/Med-High/){printf "  %-40s OR=%.2f overlap=%s p=%.1e\n",$1,$7,$6,$8}' $SGDIR/sg_rna_ir_overlap.tsv
echo "--- DEG-UP overlap (mtime=$(stat -c %y $SGDIR/sg_deg_overlap.tsv|cut -d. -f1)) ---"
awk -F"\t" 'NR>1 && tolower($2)~/up/ && ($1~/Khong 2017 SG RNA Enriched/||$1~/High Cons/||$1~/Med-High/){printf "  %-38s OR=%.2f overlap=%s p=%.1e\n",$1,$9,$8,$10}' $SGDIR/sg_deg_overlap.tsv
$NOTIFY "[DHX8-Sakai] Stage2 SG full re-run DONE (IR+DEG v2)"
echo "$(date) === SG fix2 DONE ==="
