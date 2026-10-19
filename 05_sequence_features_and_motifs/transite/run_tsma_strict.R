#!/usr/bin/env Rscript
# Transite TSMA on STRICT DEG-UP foreground (log2FC>1 & padj<0.05, matching the
# manuscript Fig3A "upregulated" definition) vs both backgrounds:
#   (A) MEF-expressed 3'UTRs  (correct universe)
#   (B) all mm10 3'UTRs       (naive / web-tool-like)
# k-mer TSMA (Fig3C/D are hexamer-based) + matrix TSMA(MEF-bg) as orthogonal RBP check.
suppressMessages({library(transite); library(Biostrings)})
WD <- "/mnt/HDD2/transite_fig3cd"
KPERM <- 5000; NC <- 8

# --- build strict foreground from DESeq2 (log2FC>1 & padj<0.05), map via existing fg TSV ---
v2 <- read.delim("/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv", stringsAsFactors=FALSE)
strict_genes <- v2$gene[!is.na(v2$padj) & v2$padj<0.05 & v2$log2FoldChange>1]
fg_all <- read.delim(file.path(WD,"foreground_up_only_v2.tsv"), stringsAsFactors=FALSE)
fg <- fg_all[fg_all$gene_symbol %in% strict_genes,]
write.table(fg, file.path(WD,"foreground_up_strict_v2.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat("strict DEG-UP genes (log2FC>1 & padj<0.05):", length(strict_genes),
    " | with 3'UTR mapping (foreground):", nrow(fg), "\n")

seqs <- chartr("T","U", as.character(readRDS(file.path(WD,"mm10_3utr_seqs.rds"))))
bg   <- read.delim(file.path(WD,"background_all_expressed_v2.tsv"), stringsAsFactors=FALSE)
idc  <- if("refseq_id" %in% names(fg)) "refseq_id" else names(fg)[1]
fg_seqs  <- seqs[names(seqs) %in% unique(fg[[idc]])]
mef_bg   <- seqs[names(seqs) %in% unique(bg[[idc]])]
all_bg   <- seqs
cat("fg_seqs:", length(fg_seqs), " mef_bg:", length(mef_bg), " all_bg:", length(all_bg), "\n")

run_kmer <- function(fgs, bgs, out, lab){
  kr <- run_kmer_tsma(list(fgs), bgs, k=6, fg_permutations=KPERM,
                      produce_plot=FALSE, p_adjust_method="BH", n_cores=NC)
  kdf <- kr[[1]]$enrichment_df
  write.table(kdf, file.path(WD,out), sep="\t", row.names=FALSE, quote=FALSE)
  sig <- kdf[!is.na(kdf$adj_p_value) & kdf$adj_p_value<0.05,]
  chk <- function(k){ r<-kdf[kdf$kmer==k,]; if(nrow(r)) sprintf("%s enr=%.2f ap=%.1e",k,r$enrichment[1],r$adj_p_value[1]) else paste0(k," NA") }
  cat(sprintf("[KMER %s] -> %s | sig(adj<0.05)=%d enriched=%d | polyU:%s | UGrep:%s | GUrich:%s\n",
      lab, out, nrow(sig), sum(sig$enrichment>1), chk("UUUUUU"), chk("UGUGUG"), chk("GUUGAG")))
  if(nrow(sig)) { se<-sig[order(-sig$enrichment),]; cat("   top-enrichment:",
      paste(sprintf("%s(%.2f)",head(se$kmer,6),head(se$enrichment,6)),collapse=" "),"\n") }
}
cat("\n===== k-mer TSMA (strict fg) =====\n")
run_kmer(fg_seqs, mef_bg, "fig3c_kmer_strict_mef.tsv", "MEF-bg")
run_kmer(fg_seqs, all_bg, "fig3c_kmer_strict_all.tsv", "ALL-bg")

cat("\n===== matrix TSMA (strict fg, MEF-bg, N-filtered) =====\n")
fg_nf  <- fg_seqs[!grepl("N", fg_seqs, fixed=TRUE)]
mef_nf <- mef_bg[!grepl("N", mef_bg, fixed=TRUE)]
mr <- run_matrix_tsma(list(fg_nf), mef_nf, p_adjust_method="BH", n_cores=NC,
                      max_fg_permutations=50000, min_fg_permutations=1000)
mdf <- mr$enrichment_dfs[[1]]
write.table(mdf, file.path(WD,"fig3d_motif_strict_mef.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
ms <- mdf[!is.na(mdf$adj_p_value) & mdf$adj_p_value<0.05,]
cat(sprintf("[MATRIX strict MEF-bg] motifs=%d sig(adj<0.05)=%d enriched=%d\n", nrow(mdf), nrow(ms), sum(ms$enrichment>1)))
mdf <- mdf[order(mdf$adj_p_value),]
cat("  top8:", paste(sprintf("%s(enr=%.2f,ap=%.2f)", head(mdf$motif_rbps,8), head(mdf$enrichment,8), head(mdf$adj_p_value,8)), collapse="  "), "\n")
cat("STRICT_DONE\n")
