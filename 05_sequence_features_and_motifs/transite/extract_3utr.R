#!/usr/bin/env Rscript
# Extract mm10 3'UTR sequences for the RefSeq transcript sets (foreground/background)
# for Transite TSMA. No transite needed here (GenomicFeatures + Biostrings + Rsamtools).
suppressMessages({library(GenomicFeatures); library(Biostrings); library(Rsamtools)})
WD <- "/mnt/HDD2/transite_fig3cd"
GTF <- "/mnt/NAS/database/genes/mm10.ncbiRefSeq.gtf"
FA  <- "/mnt/NAS/database/UCSC_genome/mm10/mm10.fa"

if (!file.exists(paste0(FA,".fai"))) { cat("indexing genome...\n"); indexFa(FA) }
cat("building TxDb from GTF...\n")
txdb <- suppressWarnings(makeTxDbFromGFF(GTF, format="gtf"))
utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
cat("transcripts with 3'UTR:", length(utr3), "\n")
gen <- FaFile(FA)
cat("extracting 3'UTR sequences...\n")
seqs <- extractTranscriptSeqs(gen, utr3)   # DNAStringSet, names = transcript_id (with version)
cat("extracted seqs:", length(seqs), " median len:", median(width(seqs)), "\n")
saveRDS(seqs, file.path(WD,"mm10_3utr_seqs.rds"))

# match against our sets
fg <- read.delim(file.path(WD,"foreground_up_only_v2.tsv"), stringsAsFactors=FALSE)
bg <- read.delim(file.path(WD,"background_all_expressed_v2.tsv"), stringsAsFactors=FALSE)
idcol <- if("refseq_id" %in% names(fg)) "refseq_id" else names(fg)[1]
fg_ids <- unique(fg[[idcol]]); bg_ids <- unique(bg[[idcol]])
match_rate <- function(ids) c(n=length(ids), matched=sum(ids %in% names(seqs)))
cat("FOREGROUND:", match_rate(fg_ids), "\n")
cat("BACKGROUND:", match_rate(bg_ids), "\n")
# fallback: strip version and retry for unmatched
strip <- function(x) sub("\\.[0-9]+$","",x)
seqnames_nov <- strip(names(seqs))
cat("FG match (version-stripped):", sum(strip(fg_ids) %in% seqnames_nov), "/", length(fg_ids), "\n")
cat("BG match (version-stripped):", sum(strip(bg_ids) %in% seqnames_nov), "/", length(bg_ids), "\n")
cat("DONE extract_3utr\n")
