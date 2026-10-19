#!/bin/bash
# RECONSTRUCTED FILE (code release, 2026-10-19).
# On 2026-07-06 16:11 JST the v2 DESeq2 table was converted to the CSV layout read by
# fgsea_spliceosome_phase.R, analysis_sg_transcriptome_gsea.py and _fig6f_table_v2.py
# (data/deseq2_results.csv) with the awk one-liner below, run interactively on the server.
# The command is copied from the session log; the copy/redeploy loop that followed it is summarised.
NEW=/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv
# input columns: baseMean log2FoldChange lfcSE stat pvalue padj gene
awk -F"\t" 'NR==1{print "gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj"; next}{print $7","$1","$2","$3","$4","$5","$6}' "$NEW" > deseq2_results_v2.csv
# The CSV was then copied over every deseq2_results.csv under /mnt/HDD2/unified_results
# (including sg_analysis/gsea_dhx8/data/deseq2_results.csv and 2026-02-27_comprehensive/data/deseq2_results.csv).
