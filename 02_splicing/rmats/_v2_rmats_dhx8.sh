#!/bin/bash
# Stage 2: rMATS DHX8 v2 (new siCont_2). Command per Obsidian 2026-02-18 pipeline note.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rmats
NOTIFY=~/bin/notify_telegram.sh
MAP=/mnt/HDD2/Sakai_RNA-seq/output/mapping
GTF=/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf
OUT=/mnt/HDD2/unified_results/rmats/dhx8_v2
TMP=/mnt/HDD2/unified_results/rmats/dhx8_v2_tmp
rm -rf "$OUT" "$TMP"; mkdir -p "$OUT" "$TMP"
# b1 = control (siCont), b2 = KD (siDHX8), comma-separated one line (per Obsidian: --b1 ctrl --b2 kd)
echo "$MAP/siCont_1_Aligned.sortedByCoord.out.bam,$MAP/siCont_2_Aligned.sortedByCoord.out.bam,$MAP/siCont_3_Aligned.sortedByCoord.out.bam" > "$OUT/rmats_ctrl.txt"
echo "$MAP/siDHX8_1_Aligned.sortedByCoord.out.bam,$MAP/siDHX8_2_Aligned.sortedByCoord.out.bam,$MAP/siDHX8_3_Aligned.sortedByCoord.out.bam" > "$OUT/rmats_kd.txt"
echo "$(date) === rMATS DHX8 v2 START ==="
$NOTIFY "[DHX8-Sakai] Stage2 rMATS v2 START"
rmats.py --b1 "$OUT/rmats_ctrl.txt" --b2 "$OUT/rmats_kd.txt" \
    --gtf "$GTF" --od "$OUT" --tmp "$TMP" \
    -t paired --readLength 150 --nthread 32 --task both \
    2>&1 | tee "$OUT/rmats_v2.log" | tail -15 || { $NOTIFY "[DHX8-Sakai] rMATS v2 FAILED"; exit 1; }
echo "=== v2 summary ==="; cat "$OUT/summary.txt" 2>/dev/null
echo "=== GTF検証: fromGTF.SE 行数 (v1=57917 なら GTF正解) ==="; wc -l "$OUT/fromGTF.SE.txt" 2>/dev/null
$NOTIFY "[DHX8-Sakai] Stage2 rMATS v2 DONE (fromGTF.SE $(wc -l < "$OUT/fromGTF.SE.txt" 2>/dev/null), v1=57917)"
echo "$(date) === rMATS v2 DONE ==="
