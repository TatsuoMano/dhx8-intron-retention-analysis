#!/bin/bash
# RECONSTRUCTED FILE (code release, 2026-10-20).
# MAJIQ deltapsi with splice-graph annotation (columns gene_name, is_denovo, event_type, ...) for
# the comparison datasets; this wrote majiq/<dataset>/deltapsi_annotated.tsv, the input of
# analysis_deltapsi_integration.py, _bespoke_ce_fdr_arial.py (Supplementary Fig. S10) and
# 02_splicing/count_events_fig4b.py (Fig. 4B CE row; GSE126543: 160 de novo junctions,
# 24 CE candidates, 9 strong CEs).
#
# Run as one interactive ssh command on 2026-03-29 09:22 JST, after run_majiq_robust.sh had built
# the splice graphs and the psi-coverage files (majiq build, then majiq psi-coverage per group;
# the commands are in majiq/<dataset>/build.log and psicov_<group>.log). Recovered verbatim from
# ~/.claude-mem/logs/claude-mem-2026-03-29.log; it also matches the "command" field that MAJIQ
# wrote into the JSON header of each deltapsi_annotated.tsv (checked with --check).
# The same block wrote dhx8/deltapsi_annotated.tsv, which was later replaced by the v2 run in
# _v2_motif_majiq.sh (after the siCont_2 correction); the comparison datasets were not re-run.
# The deltapsi_final.tsv files next to them come from an earlier run of the same commands without
# --splicegraph (no gene annotation) and are not used.
#
# Upstream steps for GSE126543 (from build.log and psicov_*.log; run by run_majiq_robust.sh):
#   majiq build --license $LICENSE gencode.v38.hg38.gff3.gz gse126543/experiments.tsv gse126543/build --nthreads 32
#   majiq psi-coverage --license $LICENSE gse126543/build/splicegraph.zarr gse126543/TDP43_negative_psicov.zarr \
#       build/SRR8571938.sorted.sj ... SRR8571952.sorted.sj --nthreads 32      (7 TDP-43(-) samples)
#   majiq psi-coverage ... gse126543/TDP43_positive_psicov.zarr build/SRR8571937.sorted.sj ... --nthreads 32  (7 TDP-43(+))
# Group order: -psi1 TDP43_negative, -psi2 TDP43_positive.
#
# Usage:  deltapsi_annotated_2026-03-29.sh [--check]  compare the recorded command of each output
#                                                     with the command below (the default; reads
#                                                     headers only)
#         deltapsi_annotated_2026-03-29.sh --run      run the four deltapsi commands (overwrites
#                                                     the deltapsi_annotated.tsv files)
# Needs bash 4 or later (the check uses an associative array, declare -A).
if (( BASH_VERSINFO[0] < 4 )); then
    echo "bash 4 or later is required (this is $BASH_VERSION)" >&2; exit 2
fi
MODE=${1:---check}
case "$MODE" in
    --check|--run) ;;
    *) echo "usage: $0 [--check | --run]" >&2; exit 2 ;;
esac

if [ "$MODE" = "--check" ]; then
    B=/mnt/HDD2/unified_results/majiq
    L=/mnt/HDD2/majiq_env/majiq_license_academic_official.lic
    declare -A WANT=(
      [gse126543]="/mnt/HDD2/majiq_env/bin/majiq deltapsi -psi1 TDP43_negative_psicov.zarr -psi2 TDP43_positive_psicov.zarr -n TDP43_negative TDP43_positive --output-tsv deltapsi_annotated.tsv --splicegraph build/splicegraph.zarr --license $L --overwrite"
      [gse27218]="/mnt/HDD2/majiq_env/bin/majiq deltapsi -psi1 Control_psicov.zarr -psi2 TDP43_KD_psicov.zarr -n Control TDP43_KD --output-tsv deltapsi_annotated.tsv --splicegraph build/splicegraph.zarr --license $L --overwrite"
      [gse65973]="/mnt/HDD2/majiq_env/bin/majiq deltapsi -psi1 Control_psicov.zarr -psi2 dNLS_TDP_psicov.zarr -n Control dNLS_TDP --output-tsv deltapsi_annotated.tsv --splicegraph build/splicegraph.zarr --license $L --overwrite"
    )
    for d in gse126543 gse65973 gse27218; do
        have=$(grep -m1 '"command"' $B/$d/deltapsi_annotated.tsv | sed 's/^#[[:space:]]*"command": "//; s/",[[:space:]]*$//')
        [ "$have" = "${WANT[$d]}" ] && echo "OK        $d/deltapsi_annotated.tsv" || { echo "DIFFERENT $d"; echo "  have: $have"; }
    done
    exit 0
fi

# ---------------- verbatim (the ssh "..." quoting removed: \$ -> $) ----------------
source /mnt/HDD2/majiq_env/bin/activate
LICENSE=/mnt/HDD2/majiq_env/majiq_license_academic_official.lic

# Re-run all 4 with splicegraph annotation
echo '=== Starting annotated deltapsi runs ==='

(cd /mnt/HDD2/unified_results/majiq/dhx8 && \
 majiq deltapsi -psi1 Control_psicov.zarr -psi2 DHX8_KD_psicov.zarr \
   -n Control DHX8_KD --output-tsv deltapsi_annotated.tsv \
   --splicegraph build/splicegraph.zarr \
   --license $LICENSE --overwrite 2>&1 | tail -2 && \
 echo 'DHX8: DONE') &

(cd /mnt/HDD2/unified_results/majiq/gse126543 && \
 majiq deltapsi -psi1 TDP43_negative_psicov.zarr -psi2 TDP43_positive_psicov.zarr \
   -n TDP43_negative TDP43_positive --output-tsv deltapsi_annotated.tsv \
   --splicegraph build/splicegraph.zarr \
   --license $LICENSE --overwrite 2>&1 | tail -2 && \
 echo 'GSE126543: DONE') &

(cd /mnt/HDD2/unified_results/majiq/gse27218 && \
 majiq deltapsi -psi1 Control_psicov.zarr -psi2 TDP43_KD_psicov.zarr \
   -n Control TDP43_KD --output-tsv deltapsi_annotated.tsv \
   --splicegraph build/splicegraph.zarr \
   --license $LICENSE --overwrite 2>&1 | tail -2 && \
 echo 'GSE27218: DONE') &

(cd /mnt/HDD2/unified_results/majiq/gse65973 && \
 majiq deltapsi -psi1 Control_psicov.zarr -psi2 dNLS_TDP_psicov.zarr \
   -n Control dNLS_TDP --output-tsv deltapsi_annotated.tsv \
   --splicegraph build/splicegraph.zarr \
   --license $LICENSE --overwrite 2>&1 | tail -2 && \
 echo 'GSE65973: DONE') &

wait
echo ''
echo '=== Output summary ==='
for d in dhx8 gse126543 gse27218 gse65973; do
  f=/mnt/HDD2/unified_results/majiq/$d/deltapsi_annotated.tsv
  if [ -f "$f" ]; then
    lines=$(wc -l < "$f")
    echo "$d: $lines lines"
  else
    echo "$d: FAILED"
  fi
done
