#!/bin/bash
# Stage 2: MAJIQ DHX8 rebuild (whole build/, new siCont_2 BAM). ~hours.
set -uo pipefail
source /mnt/HDD2/majiq_env/bin/activate
LICENSE=/mnt/HDD2/majiq_env/majiq_license_academic_official.lic
DSDIR=/mnt/HDD2/unified_results/majiq/dhx8
NOTIFY=~/bin/notify_telegram.sh
echo "$(date) === MAJIQ dhx8 rebuild START ==="
$NOTIFY "[DHX8-Sakai] Stage2 MAJIQ dhx8 rebuild START (~hours)"
# experiments.tsv already points to the 6 _Aligned BAMs (new siCont_2). Rebuild whole build/.
rm -rf "$DSDIR/build"
echo "$(date) majiq build (splicegraph, pools all 6)..."
majiq build --license "$LICENSE" \
    /mnt/NAS/database/UCSC_GFF3/gencode.vM25.mm10.gff3.gz \
    "$DSDIR/experiments.tsv" \
    "$DSDIR/build" --nthreads 32 > "$DSDIR/build.log" 2>&1 || { $NOTIFY "[DHX8-Sakai] MAJIQ build FAILED"; exit 1; }
echo "$(date) majiq deltapsi..."
majiq deltapsi --license "$LICENSE" \
    -psi1 "$DSDIR"/build/Control_*.sj \
    -psi2 "$DSDIR"/build/DHX8_KD_*.sj \
    -n Control DHX8_KD \
    --output-tsv "$DSDIR/deltapsi.tsv" > "$DSDIR/deltapsi.log" 2>&1 || { $NOTIFY "[DHX8-Sakai] MAJIQ deltapsi FAILED"; exit 1; }
echo "$(date) === MAJIQ dhx8 DONE ==="
$NOTIFY "[DHX8-Sakai] Stage2 MAJIQ dhx8 DONE → deltapsi.tsv"
