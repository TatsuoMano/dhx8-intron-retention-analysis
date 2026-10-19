#!/bin/bash
# Stage 2b: motif S12 (clipper3) + MAJIQ Fig4B/S9 chain (majiq_env), launched as independent background jobs.
NOTIFY=~/bin/notify_telegram.sh
LIMMADIR=/mnt/HDD2/TDP43_comparison/validation_20260226

# ---------- motif S12 (repoint plain names to v2 13-col, then run) ----------
(
  cp "$LIMMADIR/limma_retained_dhx8_13col.tsv" "$LIMMADIR/limma_retained_dhx8.tsv"
  cp "$LIMMADIR/limma_ir_results_dhx8_13col.tsv" "$LIMMADIR/limma_ir_results_dhx8.tsv"
  echo "motif inputs repointed: retained=$(($(wc -l < $LIMMADIR/limma_retained_dhx8.tsv)-1)) rows (v2 1258)"
  source ~/miniconda3/etc/profile.d/conda.sh && conda activate clipper3
  cd /mnt/HDD2/phase3_specificity
  python analysis_motif_enrichment.py > _v2_motif.log 2>&1 && \
    echo "motif OK: $(grep -iE 'GC|gc_content' motif_enrichment_results.tsv | head -2)" || { echo "motif FAIL"; tail -6 _v2_motif.log; }
  $NOTIFY "[DHX8-Sakai] Fig S12 motif v2 done (1258 retained)"
) > /mnt/HDD2/_v2_motif.nohup 2>&1 &
echo "motif PID $!"

# ---------- MAJIQ Fig4B/S9 (psi-coverage -> deltapsi annotated -> integration) ----------
(
  source /mnt/HDD2/majiq_env/bin/activate
  LIC=/mnt/HDD2/majiq_env/majiq_license_academic_official.lic
  D=/mnt/HDD2/unified_results/majiq/dhx8
  SG=$D/build/splicegraph.zarr
  rm -rf $D/Control_psicov.zarr $D/DHX8_KD_psicov.zarr
  echo "=== psi-coverage Control ==="
  majiq psi-coverage --license $LIC $SG $D/Control_psicov.zarr $D/build/siCont_1_*.sj $D/build/siCont_2_*.sj $D/build/siCont_3_*.sj --nthreads 32 2>&1 | tail -3
  echo "=== psi-coverage DHX8_KD ==="
  majiq psi-coverage --license $LIC $SG $D/DHX8_KD_psicov.zarr $D/build/siDHX8_1_*.sj $D/build/siDHX8_2_*.sj $D/build/siDHX8_3_*.sj --nthreads 32 2>&1 | tail -3
  echo "=== deltapsi (annotated) ==="
  majiq deltapsi --license $LIC -grp1 $D/Control_psicov.zarr -grp2 $D/DHX8_KD_psicov.zarr -n Control DHX8_KD \
    --splicegraph $SG --annotated $SG --output-tsv $D/deltapsi_annotated.tsv --overwrite 2>&1 | tail -3
  echo "deltapsi_annotated rows=$(wc -l < $D/deltapsi_annotated.tsv 2>/dev/null)"
  echo "=== integration ==="
  cd /mnt/HDD2/unified_results/majiq && python3 analysis_deltapsi_integration.py 2>&1 | tail -5
  $NOTIFY "[DHX8-Sakai] MAJIQ v2 psicov->deltapsi->integration done (CE figure step next)"
) > /mnt/HDD2/_v2_majiq_chain.nohup 2>&1 &
echo "MAJIQ PID $!"
