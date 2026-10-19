#!/bin/bash
# MAJIQ Robust Pipeline — All 4 datasets
# Addresses all 8 known failure modes (see plan)
# NO set -e: each dataset runs independently

source /mnt/HDD2/majiq_env/bin/activate
LICENSE=/mnt/HDD2/majiq_env/majiq_license_academic_official.lic
OUTBASE=/mnt/HDD2/unified_results/majiq
GFF3_MM10=/mnt/NAS/database/UCSC_GFF3/gencode.vM25.mm10.gff3.gz
GFF3_HG38=/mnt/NAS/database/UCSC_GFF3/gencode.v38.hg38.gff3.gz

notify() { ~/bin/notify_telegram.sh "$1"; }

run_dataset() {
    local DS=$1
    local GFF3=$2
    local DSDIR=$OUTBASE/$DS

    echo ""
    echo "============================================================"
    echo "MAJIQ: $DS"
    echo "============================================================"

    # Ensure directories exist
    mkdir -p $DSDIR/build

    # Verify experiments.tsv
    if [ ! -f "$DSDIR/experiments.tsv" ]; then
        echo "  ERROR: $DSDIR/experiments.tsv not found"
        notify "[DHX8-Sakai] MAJIQ $DS FAILED: no experiments.tsv"
        return 1
    fi
    echo "  experiments.tsv: $(tail -n +2 $DSDIR/experiments.tsv | wc -l) samples"

    # Step 1: Build (skip if sufficient .sj dirs exist)
    local N_SJ=$(ls -d $DSDIR/build/*.sj 2>/dev/null | wc -l)
    local N_SAMPLES=$(tail -n +2 $DSDIR/experiments.tsv | wc -l)

    if [ "$N_SJ" -ge "$N_SAMPLES" ]; then
        echo "  Build: SKIPPED ($N_SJ .sj dirs already exist for $N_SAMPLES samples)"
    else
        echo "  Build: starting ($N_SJ .sj exist, need $N_SAMPLES)..."
        majiq build --license $LICENSE \
            $GFF3 \
            $DSDIR/experiments.tsv \
            $DSDIR/build \
            --nthreads 32 \
            > $DSDIR/build.log 2>&1

        N_SJ=$(ls -d $DSDIR/build/*.sj 2>/dev/null | wc -l)
        if [ "$N_SJ" -lt "$N_SAMPLES" ]; then
            echo "  Build FAILED: only $N_SJ .sj dirs (expected $N_SAMPLES)"
            echo "  Last 5 lines of build.log:"
            tail -5 $DSDIR/build.log
            notify "[DHX8-Sakai] MAJIQ $DS build FAILED ($N_SJ/$N_SAMPLES .sj)"
            return 1
        fi
        echo "  Build: OK ($N_SJ .sj dirs)"
    fi

    # Step 2: Psi-coverage (required before deltapsi in MAJIQ v3)
    local GRP1=$(awk -F'\t' 'NR>1{print $3}' $DSDIR/experiments.tsv | sort -u | head -1)
    local GRP2=$(awk -F'\t' 'NR>1{print $3}' $DSDIR/experiments.tsv | sort -u | tail -1)
    echo "  Groups: $GRP1 vs $GRP2"

    local SG=$DSDIR/build/splicegraph.zarr

    # Collect .sj dirs by group
    local SJ1="" SJ2="" N1=0 N2=0
    while IFS=$'\t' read -r sid path grp strand; do
        [ "$sid" = "sample_id" ] && continue
        local sjdir=$(ls -d $DSDIR/build/${sid}*.sj 2>/dev/null | head -1)
        [ -z "$sjdir" ] && { echo "  WARNING: no .sj for $sid"; continue; }
        if [ "$grp" = "$GRP1" ]; then SJ1="$SJ1 $sjdir"; N1=$((N1+1))
        else SJ2="$SJ2 $sjdir"; N2=$((N2+1)); fi
    done < $DSDIR/experiments.tsv

    echo "  $GRP1: $N1 samples, $GRP2: $N2 samples"

    if [ "$N1" -lt 1 ] || [ "$N2" -lt 1 ]; then
        echo "  ERROR: insufficient samples"
        notify "[DHX8-Sakai] MAJIQ $DS FAILED: $GRP1=$N1, $GRP2=$N2"
        return 1
    fi

    # Psi-coverage for group 1
    local PSICOV1=$DSDIR/${GRP1}_psicov.zarr
    if [ ! -d "$PSICOV1" ]; then
        echo "  Psi-coverage $GRP1..."
        majiq psi-coverage --license $LICENSE \
            $SG $PSICOV1 $SJ1 \
            --nthreads 32 \
            > $DSDIR/psicov_${GRP1}.log 2>&1
        if [ ! -d "$PSICOV1" ]; then
            echo "  Psi-coverage $GRP1 FAILED"
            tail -3 $DSDIR/psicov_${GRP1}.log
            notify "[DHX8-Sakai] MAJIQ $DS psicov $GRP1 FAILED"
            return 1
        fi
    fi

    # Psi-coverage for group 2
    local PSICOV2=$DSDIR/${GRP2}_psicov.zarr
    if [ ! -d "$PSICOV2" ]; then
        echo "  Psi-coverage $GRP2..."
        majiq psi-coverage --license $LICENSE \
            $SG $PSICOV2 $SJ2 \
            --nthreads 32 \
            > $DSDIR/psicov_${GRP2}.log 2>&1
        if [ ! -d "$PSICOV2" ]; then
            echo "  Psi-coverage $GRP2 FAILED"
            tail -3 $DSDIR/psicov_${GRP2}.log
            notify "[DHX8-Sakai] MAJIQ $DS psicov $GRP2 FAILED"
            return 1
        fi
    fi

    # Step 3: Deltapsi
    echo "  Deltapsi..."
    majiq heterogen --license $LICENSE \
        -psi1 $PSICOV1 \
        -psi2 $PSICOV2 \
        -n $GRP1 $GRP2 \
        --output-tsv $DSDIR/heterogen.tsv \
        > $DSDIR/heterogen.log 2>&1

    if [ $? -eq 0 ] && [ -f "$DSDIR/heterogen.tsv" ]; then
        local N_EVENTS=$(tail -n +2 $DSDIR/heterogen.tsv | wc -l)
        echo "  Heterogen: OK ($N_EVENTS events)"
        notify "[DHX8-Sakai] MAJIQ $DS DONE ($N_EVENTS events)"
    else
        echo "  Heterogen FAILED. Last 5 lines:"
        tail -5 $DSDIR/heterogen.log
        notify "[DHX8-Sakai] MAJIQ $DS heterogen FAILED"
        return 1
    fi
}

# ============================================================
# Kill any existing MAJIQ processes
# ============================================================
# Note: do NOT pkill run_majiq — it would kill this script
echo "Killing existing MAJIQ build/deltapsi processes..."
pkill -f 'majiq build' 2>/dev/null
pkill -f 'majiq deltapsi' 2>/dev/null
sleep 2

notify "[DHX8-Sakai] MAJIQ robust pipeline started (4 datasets)"

# ============================================================
# Run each dataset
# ============================================================
run_dataset dhx8 $GFF3_MM10
run_dataset gse27218 $GFF3_MM10
run_dataset gse65973 $GFF3_MM10
run_dataset gse126543 $GFF3_HG38

notify "[DHX8-Sakai] MAJIQ robust pipeline FINISHED"
echo ""
echo "============================================================"
echo "MAJIQ pipeline complete"
echo "============================================================"
