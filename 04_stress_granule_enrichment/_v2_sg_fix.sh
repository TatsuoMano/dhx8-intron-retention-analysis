#!/bin/bash
# Stage 2 SG fix: join intron_length/chr/start/end/strand into v2 IR file (13-col), redeploy, re-run full SG.
set -uo pipefail
NOTIFY=~/bin/notify_telegram.sh
LIMMADIR=/mnt/HDD2/TDP43_comparison/validation_20260226
META=/mnt/HDD2/TDP43_comparison/event_level_20260225/intron_metadata.tsv
SGDIR=/mnt/HDD2/unified_results/sg_analysis/gsea_dhx8
# --- 1. python join: v2 IR (8-col) + metadata -> 13-col ---
cat > /tmp/join_ir_meta.py <<'PY'
import pandas as pd, sys
ir = pd.read_csv(sys.argv[1], sep='\t')
meta = pd.read_csv(sys.argv[2], sep='\t')
addc = ['chr','start','end','strand','intron_length']
m = ir.merge(meta[['intron_id']+addc], on='intron_id', how='left')
cols = ['intron_id','logFC','avg_expr','t_stat','pvalue','adj_p','chr','start','end','strand','gene_name','intron_length','retained']
m = m[[c for c in cols if c in m.columns]]
nmiss = int(m['intron_length'].isna().sum())
m.to_csv(sys.argv[3], sep='\t', index=False)
print(f"  joined rows={len(m)} intron_length_missing={nmiss}")
PY
echo "$(date) === 13列化 (retained) ==="
python3 /tmp/join_ir_meta.py "$LIMMADIR/limma_retained_dhx8_s2.tsv" "$META" "$LIMMADIR/limma_retained_dhx8_13col.tsv"
# --- 2. redeploy 13-col to canonical locations SG reads ---
for d in $SGDIR/data /mnt/HDD2/unified_results/sg_transcriptome_s2_reanalysis/data; do
  cp "$LIMMADIR/limma_retained_dhx8_13col.tsv" "$d/limma_retained_dhx8.tsv" && echo "  redeploy(13col) → $d ($(head -1 $d/limma_retained_dhx8.tsv | tr '\t' '\n' | wc -l)列)"
done
# --- 3. re-run full SG (IR overlap + DEG overlap + length confound) ---
echo "$(date) === analysis_sg_transcriptome_gsea.py 完全再実行 ==="
cd $SGDIR
python3 analysis_sg_transcriptome_gsea.py > _v2_sg_gsea_full.log 2>&1 && echo "  SG full OK" || { echo "  SG まだ失敗:"; tail -8 _v2_sg_gsea_full.log; $NOTIFY "[DHX8-Sakai] SG full re-run FAILED"; exit 1; }
# --- 4. verify DEG overlap now v2 ---
echo "=== sg_deg_overlap.tsv 再生成確認 (DEG_UP deg_in_expressed=652なら v2) ==="
awk -F"\t" 'NR>1 && tolower($2)~/up/{print "  DEG_UP deg_in_expressed="$6; exit}' $SGDIR/sg_deg_overlap.tsv
stat -c "  mtime=%y" $SGDIR/sg_deg_overlap.tsv | cut -d. -f1
$NOTIFY "[DHX8-Sakai] Stage2 SG full re-run DONE (IR+DEG overlap v2)"
echo "$(date) === SG fix DONE ==="
