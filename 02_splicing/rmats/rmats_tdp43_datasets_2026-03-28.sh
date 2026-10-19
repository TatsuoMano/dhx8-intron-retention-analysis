#!/bin/bash
# RECONSTRUCTED FILE (code release, 2026-10-19).
# rMATS for the TDP-43 comparison datasets (Fig. 4B columns GSE65973 and GSE126543) was run as an
# interactive ssh command on 2026-03-28 (about 09:54-10:00 JST) and never saved as a script.
# The BAM-list and rmats.py commands below are copied from the session log; only the Telegram
# notification calls are removed. rmats-turbo was used from a source build (/mnt/HDD2/rmats-turbo).
# Outputs: /mnt/HDD2/unified_results/rmats/{gse65973,gse126543}/ (summary.txt, *.MATS.JCEC.txt).
RMATS=/mnt/HDD2/rmats-turbo/rmats-turbo
export LD_LIBRARY_PATH=$RMATS/bamtools/lib:$LD_LIBRARY_PATH
OUTBASE=/mnt/HDD2/unified_results/rmats

# GSE65973: control vs bigenic (dNLS-TDP-43)
awk -F'\t' '$4=="control" {printf "/mnt/HDD2/TDP43_comparison/GSE65973/mapping/"$1"_Aligned.sortedByCoord.out.bam,"}' /mnt/HDD2/TDP43_comparison/GSE65973/sample_sheet.tsv | sed 's/,$//' > /tmp/rmats_65973_b1.txt
awk -F'\t' '$4!="control" {printf "/mnt/HDD2/TDP43_comparison/GSE65973/mapping/"$1"_Aligned.sortedByCoord.out.bam,"}' /mnt/HDD2/TDP43_comparison/GSE65973/sample_sheet.tsv | sed 's/,$//' > /tmp/rmats_65973_b2.txt

# GSE126543: TDP43_negative (b1) vs TDP43_positive (b2), NeuN+ nuclei only
awk -F'\t' '$5=="TDP43_negative" && $6=="NeuN+" {printf "/mnt/HDD2/sorted_nuclei_GSE126543/mapping/"$1"/"$1".sorted.bam,"}' /mnt/HDD2/sorted_nuclei_GSE126543/sample_sheet.tsv | sed 's/,$//' > /tmp/rmats_126543_b1.txt
awk -F'\t' '$5=="TDP43_positive" && $6=="NeuN+" {printf "/mnt/HDD2/sorted_nuclei_GSE126543/mapping/"$1"/"$1".sorted.bam,"}' /mnt/HDD2/sorted_nuclei_GSE126543/sample_sheet.tsv | sed 's/,$//' > /tmp/rmats_126543_b2.txt

mkdir -p $OUTBASE/gse65973 $OUTBASE/gse126543

cd $RMATS && python3 rmats.py --b1 /tmp/rmats_65973_b1.txt --b2 /tmp/rmats_65973_b2.txt \
  --gtf /mnt/HDD2/mm10.ncbiRefSeq.gtf -t single --readLength 50 --nthread 32 \
  --od $OUTBASE/gse65973 --tmp $OUTBASE/gse65973/tmp --libType fr-secondstrand \
  > $OUTBASE/gse65973/rmats.log 2>&1

cd $RMATS && python3 rmats.py --b1 /tmp/rmats_126543_b1.txt --b2 /tmp/rmats_126543_b2.txt \
  --gtf /mnt/HDD2/hg38.ncbiRefSeq.gtf -t paired --readLength 125 --nthread 32 \
  --od $OUTBASE/gse126543 --tmp $OUTBASE/gse126543/tmp --libType fr-unstranded \
  > $OUTBASE/gse126543/rmats.log 2>&1
