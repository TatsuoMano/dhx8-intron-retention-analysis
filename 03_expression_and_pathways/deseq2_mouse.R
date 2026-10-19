#!/usr/bin/env Rscript
# DESeq2 for mouse datasets (DHX8, GSE27218, GSE65973)
# All using correct strandedness counts

library(DESeq2)

OUTDIR <- "/mnt/HDD2/unified_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

load_featurecounts <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE)
  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta_cols)
  counts <- as.matrix(df[, count_cols])
  rownames(counts) <- df$Geneid
  clean <- gsub(".*/([^/]+)_Aligned.*", "\\1", colnames(counts))
  clean <- gsub(".*/([^/]+)_sorted\\.bam$", "\\1", clean)
  colnames(counts) <- clean
  list(counts = counts)
}

run_deseq2 <- function(counts, condition, dataset_name) {
  cat(sprintf("\n=== DESeq2: %s ===\n", dataset_name))
  cat(sprintf("  Samples: %d (%s)\n", ncol(counts),
              paste(table(condition), collapse = " vs ")))
  keep <- rowSums(counts >= 10) >= 3
  counts <- counts[keep, ]
  cat(sprintf("  Genes after filter: %d\n", nrow(counts)))

  coldata <- data.frame(condition = factor(condition),
                        row.names = colnames(counts))
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", levels(factor(condition))[2],
                                    levels(factor(condition))[1]))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[order(res_df$padj), ]

  sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("  Significant (FDR<0.05): %d (%d up, %d down)\n", sig, up, down))

  outpath <- file.path(OUTDIR, sprintf("deseq2_%s.tsv", dataset_name))
  write.table(res_df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Saved: %s\n", outpath))
  return(res_df)
}

# === DHX8 KD (MEF, mm10, -s 2) ===
cat("Loading DHX8...\n")
dhx8 <- load_featurecounts("/mnt/HDD2/Sakai_RNA-seq/counts_s2/gene_counts.txt")
# Columns: siCont_1, siCont_2, siCont_3, siDHX8_1, siDHX8_2, siDHX8_3
condition_dhx8 <- ifelse(grepl("siCont", colnames(dhx8$counts)), "Control", "DHX8_KD")
run_deseq2(dhx8$counts, condition_dhx8, "dhx8_kd_vs_ctrl")

# === GSE27218 (TDP-43 ASO KD, mm10, -s 2) ===
cat("\nLoading GSE27218...\n")
gse27218 <- load_featurecounts("/mnt/HDD2/TDP43_comparison/GSE27218/counts_s2/gene_counts.txt")
# No header in sample sheet: SRR, GSM, sample_name, condition
ss_27218 <- read.delim("/mnt/HDD2/TDP43_comparison/GSE27218/sample_sheet.tsv",
                        header = FALSE, stringsAsFactors = FALSE)
colnames(ss_27218) <- c("SRR", "GSM", "Sample", "Condition")
ctrl_27218 <- ss_27218$SRR[ss_27218$Condition == "control"]
kd_27218 <- ss_27218$SRR[ss_27218$Condition != "control"]
all_27218 <- c(ctrl_27218, kd_27218)
all_27218 <- all_27218[all_27218 %in% colnames(gse27218$counts)]
counts_27218 <- gse27218$counts[, all_27218]
condition_27218 <- ifelse(all_27218 %in% ctrl_27218, "Control", "TDP43_KD")
run_deseq2(counts_27218, condition_27218, "gse27218_tdp43kd_vs_ctrl")

# === GSE65973 (deltaNLS-hTDP-43, mm10, -s 1) ===
cat("\nLoading GSE65973...\n")
gse65973 <- load_featurecounts("/mnt/HDD2/TDP43_comparison/GSE65973/counts_s1/gene_counts.txt")
ss_65973 <- read.delim("/mnt/HDD2/TDP43_comparison/GSE65973/sample_sheet.tsv",
                        header = FALSE, stringsAsFactors = FALSE)
colnames(ss_65973) <- c("SRR", "GSM", "Sample", "Condition")
ctrl_65973 <- ss_65973$SRR[ss_65973$Condition == "control"]
kd_65973 <- ss_65973$SRR[ss_65973$Condition != "control"]
all_65973 <- c(ctrl_65973, kd_65973)
all_65973 <- all_65973[all_65973 %in% colnames(gse65973$counts)]
counts_65973 <- gse65973$counts[, all_65973]
condition_65973 <- ifelse(all_65973 %in% ctrl_65973, "Control", "dNLS_TDP43")
run_deseq2(counts_65973, condition_65973, "gse65973_dnls_vs_ctrl")

cat("\n=== Mouse DESeq2 complete ===\n")
