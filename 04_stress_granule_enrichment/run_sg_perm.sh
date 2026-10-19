#!/bin/bash
set -e
cd /mnt/HDD2/unified_results/sg_analysis/gsea_dhx8

LOG=sg_permutation_test_lenmatched.log

~/bin/notify_telegram.sh "[DHX8-Sakai] SG perm test (length+expr matched): START on ncnp-linux"
START=$(date +%s)

/home/tatsuomano-linux2/miniconda3/bin/python sg_permutation_test_lenmatched.py > "$LOG" 2>&1
RC=$?

END=$(date +%s)
DUR=$(( (END - START) / 60 ))

if [ "$RC" -eq 0 ]; then
  ~/bin/notify_telegram.sh "[DHX8-Sakai] SG perm test (length+expr matched): DONE in ${DUR}min — outputs at /mnt/HDD2/unified_results/postfix_figures/sg_permutation_test_lenmatched.{tsv,png}"
else
  TAIL=$(tail -n 5 "$LOG" | tr '\n' ' ' | head -c 250)
  ~/bin/notify_telegram.sh "[DHX8-Sakai] SG perm test (length+expr matched): FAILED rc=${RC} after ${DUR}min. Tail: ${TAIL}"
fi
