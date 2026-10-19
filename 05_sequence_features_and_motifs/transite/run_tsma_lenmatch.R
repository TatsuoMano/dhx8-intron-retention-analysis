#!/usr/bin/env Rscript
# Gating control: is the UG-repeat enrichment in strict DEG-UP 3'UTRs robust to
# 3'UTR LENGTH? Build a length-matched background (stratified sampling of MEF-expressed
# 3'UTRs to match the foreground 3'UTR length distribution) and re-run k-mer TSMA.
suppressMessages({library(transite)})
WD <- "/mnt/HDD2/transite_fig3cd"; NC <- 8; KPERM <- 5000
seqs <- chartr("T","U", as.character(readRDS(file.path(WD,"mm10_3utr_seqs.rds"))))
fg <- read.delim(file.path(WD,"foreground_up_strict_v2.tsv"), stringsAsFactors=FALSE)
bg <- read.delim(file.path(WD,"background_all_expressed_v2.tsv"), stringsAsFactors=FALSE)
idc <- "refseq_id"
fg_seqs <- seqs[names(seqs) %in% unique(fg[[idc]])]
mef_bg  <- seqs[names(seqs) %in% unique(bg[[idc]])]
mef_bg  <- mef_bg[!(names(mef_bg) %in% names(fg_seqs))]   # exclude fg from bg pool
Lf <- nchar(fg_seqs); Lb <- nchar(mef_bg)
cat(sprintf("fg n=%d (3'UTR len median=%d, IQR %d-%d)\n", length(fg_seqs), as.integer(median(Lf)), as.integer(quantile(Lf,.25)), as.integer(quantile(Lf,.75))))
cat(sprintf("mef_bg pool n=%d (median=%d)\n", length(mef_bg), as.integer(median(Lb))))

# stratified length matching: 20 fg-length bins, sample up to 5x fg count per bin
brks <- unique(quantile(Lf, probs=seq(0,1,length.out=21))); brks[1] <- 0; brks[length(brks)] <- Inf
fgb <- cut(Lf, brks, include.lowest=TRUE); bgb <- cut(Lb, brks, include.lowest=TRUE)
set.seed(42); pick <- character(0)
for (b in levels(fgb)) {
  nfg <- sum(fgb==b, na.rm=TRUE); pool <- names(mef_bg)[which(bgb==b)]
  take <- min(length(pool), nfg*5)
  if (take>0) pick <- c(pick, sample(pool, take))
}
mbg <- mef_bg[pick]
Lm <- nchar(mbg)
ks <- ks.test(Lf, Lm)
cat(sprintf("length-matched bg n=%d (median=%d)  KS D=%.3f p=%.3g  (matched if p>0.05)\n", length(mbg), as.integer(median(Lm)), ks$statistic, ks$p.value))

kr <- run_kmer_tsma(list(fg_seqs), mbg, k=6, fg_permutations=KPERM, produce_plot=FALSE, p_adjust_method="BH", n_cores=NC)
kdf <- kr[[1]]$enrichment_df
write.table(kdf, file.path(WD,"fig3c_kmer_strict_lenmatched.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
chk <- function(k){ r<-kdf[kdf$kmer==k,]; if(nrow(r)) sprintf("%s enr=%.2f adj_p=%.1e",k,r$enrichment[1],r$adj_p_value[1]) else paste0(k," NA") }
sig <- kdf[!is.na(kdf$adj_p_value)&kdf$adj_p_value<0.05 & kdf$enrichment>1,]; sig<-sig[order(sig$adj_p_value),]
cat("\n=== LENGTH-MATCHED k-mer result (strict fg) ===\n")
cat("  enriched&sig:", nrow(sig), "\n")
cat("  UG-repeat:", chk("UGUGUG"), "/", chk("GUGUGU"), "\n")
cat("  UC-repeat:", chk("UCUCUC"), "/", chk("CUCUCU"), "\n")
cat("  poly-U   :", chk("UUUUUU"), "\n")
cat("  CG-repeat:", chk("GCGCGC"), "/", chk("CGCGCG"), "\n")
cat("  top10 enriched-sig:\n"); print(head(sig[,c("kmer","enrichment","adj_p_value")],10), row.names=FALSE)
cat("LENMATCH_DONE\n")
