#!/bin/bash
# RECONSTRUCTED FILE (code release, 2026-10-20).
# featureCounts recounts of the two TDP-43 comparison datasets with their verified strandedness
# (Supplementary Table S1): GSE65973 single-end, forward stranded (-s 1); GSE126543 paired-end,
# unstranded (-s 0). These counts are the inputs of deseq2_mouse.R / deseq2_unified.R
# (gene_counts.txt) and limma_ir_mouse.R / limma_ir_unified.R (intron_counts.txt, exon_counts.txt).
#
# Both recounts were run as interactive ssh commands and were not kept as files:
#   GSE126543  2026-03-27 17:53 JST (after the earlier -s 2 gene counts were renamed
#              gene_counts_s2_wrong.txt), recorded in ~/.claude-mem/logs/claude-mem-2026-03-27.log;
#   GSE65973   2026-03-28 05:41 JST (replacing the earlier -s 2 counts in counts_s2/), recorded in
#              the same log; the BAM list was read from /tmp/gse65973_bams.txt, which no longer
#              exists, and is written out below.
# The commands below are those commands with the Telegram calls removed. Each one equals, argument
# for argument, the "Command:" line that featureCounts v2.0.3 wrote into the header of the count
# file it produced; ./recount_tdp43_datasets.sh --check compares the two without running anything.
#
# Usage:  recount_tdp43_datasets.sh [--check]  compare each command with the header of its output
#                                              (the default; reads the count-file headers only)
#         recount_tdp43_datasets.sh --run      run the six featureCounts commands (overwrites the
#                                              count files)
# Needs bash 4.3 or later (the check uses name references, declare -n).
set -uo pipefail
if (( BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 3) )); then
    echo "bash 4.3 or later is required (this is $BASH_VERSION)" >&2; exit 2
fi
MODE=${1:---check}
case "$MODE" in
    --check|--run) ;;
    *) echo "usage: $0 [--check | --run]" >&2; exit 2 ;;
esac

# ---------------- GSE65973 (dNLS-TDP-43 mouse cortex), single-end, -s 1 ----------------
G65=/mnt/HDD2/TDP43_comparison/GSE65973
BAMS65=()
for s in SRR1805867 SRR1805868 SRR1805869 SRR1805870 SRR1805871 SRR1805872 SRR1805873 SRR1805874; do
    BAMS65+=("$G65/mapping/${s}_Aligned.sortedByCoord.out.bam")
done
GTF65=/mnt/HDD2/mm10.ncbiRefSeq.gtf          # same file as /mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf (md5 bf4a6d56...)
SAF=/mnt/HDD2/TDP43_comparison/event_level_20260225   # mm10 intron / flanking-exon SAF (02_splicing/intron_retention)
OUT65=$G65/counts_s1

CMD_65_gene=(featureCounts -T 32 -s 1 -a $GTF65 -o $OUT65/gene_counts.txt "${BAMS65[@]}")
CMD_65_intron=(featureCounts -F SAF -T 32 -s 1 -a $SAF/introns.saf -o $OUT65/intron_counts.txt "${BAMS65[@]}")
CMD_65_exon=(featureCounts -F SAF -T 32 -s 1 -a $SAF/flanking_exons.saf -o $OUT65/exon_counts.txt "${BAMS65[@]}")

# ---------------- GSE126543 (sorted TDP-43+/- neuronal nuclei), paired-end, -s 0 ----------------
G126=/mnt/HDD2/sorted_nuclei_GSE126543
BAMS126=()
while read -r b; do BAMS126+=("$b"); done < <(find $G126/mapping/ -name "*.sorted.bam" | sort)   # 16 BAMs
OUT126=$G126/counts
HG38SAF=/mnt/HDD2/hg38_ir_annotations   # from 02_splicing/intron_retention/step1_hg38_intron_saf.py

CMD_126_gene=(featureCounts -T 16 -p --countReadPairs -s 0 -a /mnt/NAS/database/genes/hg38.ncbiRefSeq.gtf -o $OUT126/gene_counts.txt "${BAMS126[@]}")
CMD_126_intron=(featureCounts -F SAF -T 16 -p --countReadPairs --fracOverlap 0 -s 0 -a $HG38SAF/hg38_introns.saf -o $OUT126/intron_counts.txt "${BAMS126[@]}")
CMD_126_exon=(featureCounts -F SAF -T 16 -p --countReadPairs --fracOverlap 0 -s 0 -a $HG38SAF/hg38_flanking_exons.saf -o $OUT126/exon_counts.txt "${BAMS126[@]}")

ALL="65_gene 65_intron 65_exon 126_gene 126_intron 126_exon"

if [ "$MODE" = "--check" ]; then
    fail=0
    for k in $ALL; do
        declare -n C="CMD_$k"
        out=""; for ((i = 0; i < ${#C[@]}; i++)); do [ "${C[$i]}" = "-o" ] && out=${C[$((i + 1))]}; done
        want=$(printf '"%s" ' "${C[@]}")
        have=$(head -1 "$out" | sed 's/^# Program:featureCounts v[0-9.]*; Command://')
        if [ "$want" = "$have" ]; then echo "OK        $out"; else echo "DIFFERENT $out"; fail=1; fi
        unset -n C
    done
    exit $fail
fi

mkdir -p $OUT65 $OUT126
"${CMD_65_gene[@]}"   > $OUT65/gene_counts.log 2>&1;   echo "GSE65973 gene done"
"${CMD_65_intron[@]}" > $OUT65/intron_counts.log 2>&1; echo "GSE65973 intron done"
"${CMD_65_exon[@]}"   > $OUT65/exon_counts.log 2>&1;   echo "GSE65973 exon done"
"${CMD_126_gene[@]}"   2>&1 | tail -5
"${CMD_126_intron[@]}" 2>&1 | tail -5
"${CMD_126_exon[@]}"   2>&1 | tail -5
echo "recounts done"
