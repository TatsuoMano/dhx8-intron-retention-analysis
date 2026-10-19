#!/bin/bash
cd /mnt/HDD2/unified_results/sg_analysis/gsea_dhx8
NOTIFY=~/bin/notify_telegram.sh
$NOTIFY "[DHX8-Sakai] Stage2 SG enrichment/permutation START (IR v2=1258/1150)"
echo "=== SG transcriptome GSEA + overlap (Fig6) ==="
python3 analysis_sg_transcriptome_gsea.py > _v2_sg_gsea.log 2>&1 && echo "sg_gsea OK" || echo "sg_gsea 確認要"
echo "=== SG permutation (S10) ==="
python3 sg_permutation_test.py > _v2_sg_perm.log 2>&1 && echo "perm OK" || echo "perm 確認要"
echo "=== SG permutation lenmatched (S10) ==="
python3 sg_permutation_test_lenmatched.py > _v2_sg_perm_lm.log 2>&1 && echo "lenmatched OK" || echo "lenmatched 確認要"
$NOTIFY "[DHX8-Sakai] Stage2 SG enrichment/permutation DONE"
