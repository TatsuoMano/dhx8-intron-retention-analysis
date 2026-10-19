#!/bin/bash
# Master pipeline: Event-level IR comparison (DHX8 vs TDP-43)
# Execute on ncnp-linux
#
# Usage:
#   # Full pipeline (from scratch):
#   bash run_pipeline.sh
#
#   # Resume from a specific step:
#   bash run_pipeline.sh --from 2
#
#   # Run only one step:
#   bash run_pipeline.sh --step 3
#
set -euo pipefail

WORKDIR="/mnt/HDD2/TDP43_comparison/event_level_20260225"
PYTHON="/home/tatsuomano-linux2/miniconda3/bin/python3"
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
LOCAL_OUTDIR="$HOME/Dropbox/Research/DHX8-Sakai/Analysis/2026-02-25_event_level_IR"

# Parse arguments
FROM_STEP=1
ONLY_STEP=0
while [[ $# -gt 0 ]]; do
    case $1 in
        --from) FROM_STEP="$2"; shift 2 ;;
        --step) ONLY_STEP="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

should_run() {
    local step=$1
    if [[ $ONLY_STEP -gt 0 ]]; then
        [[ $step -eq $ONLY_STEP ]]
    else
        [[ $step -ge $FROM_STEP ]]
    fi
}

mkdir -p "$WORKDIR" "$WORKDIR/figures"

echo "============================================================"
echo "Event-level IR comparison pipeline"
echo "Working directory: $WORKDIR"
echo "Script directory:  $SCRIPTDIR"
echo "Started at:        $(date)"
echo "============================================================"

# === Step 1: Create SAF files ===
if should_run 1; then
    echo ""
    echo ">>> STEP 1: Creating non-redundant intron SAF + flanking exon SAF"
    echo "    $(date)"
    $PYTHON "$SCRIPTDIR/step1_per_intron_count.py"
    echo "    Step 1 completed at $(date)"
fi

# === Step 1b: Run featureCounts ===
if should_run 1; then
    echo ""
    echo ">>> STEP 1b: Running featureCounts (per-intron + flanking exon)"
    echo "    $(date)"
    bash "$SCRIPTDIR/run_featurecounts.sh"
    echo "    Step 1b completed at $(date)"
fi

# === Step 2: Identify retained introns ===
if should_run 2; then
    echo ""
    echo ">>> STEP 2: Computing IR ratios + identifying retained introns"
    echo "    $(date)"
    $PYTHON "$SCRIPTDIR/step2_identify_retained_introns.py"
    echo "    Step 2 completed at $(date)"

    # Quick sanity check
    echo ""
    echo "    --- Sanity checks ---"
    for ds in dhx8 gse27218 gse65973; do
        f="$WORKDIR/retained_introns_${ds}.tsv"
        if [[ -f "$f" ]]; then
            n=$(tail -n +2 "$f" | wc -l)
            g=$(tail -n +2 "$f" | cut -f12 | sort -u | wc -l)
            echo "    $ds: $n retained introns in $g genes"
        else
            echo "    $ds: NO OUTPUT FILE"
        fi
    done
fi

# === Step 3: Cross-dataset comparison ===
if should_run 3; then
    echo ""
    echo ">>> STEP 3: Event-level cross-dataset comparison"
    echo "    $(date)"
    $PYTHON "$SCRIPTDIR/step3_compare_events.py"
    echo "    Step 3 completed at $(date)"
fi

# === Step 4: Matched-null permutation ===
if should_run 4; then
    echo ""
    echo ">>> STEP 4: Matched-null permutation (5000 iterations)"
    echo "    $(date)"
    $PYTHON "$SCRIPTDIR/step4_matched_permutation.py"
    echo "    Step 4 completed at $(date)"
fi

# === Step 5: Generate figures ===
if should_run 5; then
    echo ""
    echo ">>> STEP 5: Generating figures"
    echo "    $(date)"
    $PYTHON "$SCRIPTDIR/step5_figures.py"
    echo "    Step 5 completed at $(date)"
fi

# === Copy results to local Dropbox ===
echo ""
echo ">>> Copying results to local directory"
mkdir -p "$LOCAL_OUTDIR"
for f in retained_introns_*.tsv ir_results_all_*.tsv \
         strict_event_overlap.tsv three_way_analysis.tsv \
         strict_event_permutation_matched.tsv \
         triple_intersection_introns.tsv \
         intron_metadata.tsv \
         overlap_*.tsv; do
    src="$WORKDIR/$f"
    if [[ -f "$src" ]]; then
        cp "$src" "$LOCAL_OUTDIR/"
        echo "    Copied $f"
    fi
done
if [[ -d "$WORKDIR/figures" ]]; then
    cp -r "$WORKDIR/figures" "$LOCAL_OUTDIR/"
    echo "    Copied figures/"
fi

echo ""
echo "============================================================"
echo "Pipeline completed at $(date)"
echo "Results: $WORKDIR"
echo "Local:   $LOCAL_OUTDIR"
echo "============================================================"

# === Summary ===
echo ""
echo "=== KEY RESULTS ==="
if [[ -f "$WORKDIR/strict_event_overlap.tsv" ]]; then
    echo ""
    echo "--- Event-level overlap (strict_event_overlap.tsv) ---"
    column -t -s $'\t' "$WORKDIR/strict_event_overlap.tsv" | head -20
fi
if [[ -f "$WORKDIR/strict_event_permutation_matched.tsv" ]]; then
    echo ""
    echo "--- Permutation results (strict_event_permutation_matched.tsv) ---"
    column -t -s $'\t' "$WORKDIR/strict_event_permutation_matched.tsv"
fi
if [[ -f "$WORKDIR/three_way_analysis.tsv" ]]; then
    echo ""
    echo "--- 3-way analysis ---"
    column -t -s $'\t' "$WORKDIR/three_way_analysis.tsv"
fi
