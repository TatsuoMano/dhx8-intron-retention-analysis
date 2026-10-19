#!/usr/bin/env Rscript
# DESeq2 Unified Analysis for ALS Datasets
# Runs DESeq2 on gene counts for GSE116622, GSE124439, GSE126543
# Each dataset has its own comparison groups defined below.

library(DESeq2)

OUTDIR <- "/mnt/HDD2/unified_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helper functions
# ============================================================
load_featurecounts <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE)
  # Extract count columns (BAM path columns)
  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta_cols)
  counts <- as.matrix(df[, count_cols])
  rownames(counts) <- df$Geneid
  # Clean column names: extract SRR from BAM path
  clean_names <- gsub(".*/mapping/([^/]+)/.*", "\\1", colnames(counts))
  # If SRR not in path, try another pattern
  if (all(grepl("^/", clean_names))) {
    clean_names <- gsub(".*/([^/]+)\\.sorted\\.bam$", "\\1", colnames(counts))
  }
  if (all(grepl("^/", clean_names))) {
    clean_names <- gsub(".*/([^/]+)_Aligned.*", "\\1", colnames(counts))
  }
  colnames(counts) <- clean_names
  list(counts = counts, gene_info = df[, intersect(meta_cols, colnames(df))])
}

run_deseq2 <- function(counts, condition, dataset_name) {
  cat(sprintf("\n=== DESeq2: %s ===\n", dataset_name))
  cat(sprintf("  Samples: %d (%s)\n", ncol(counts),
              paste(table(condition), collapse = " vs ")))
  cat(sprintf("  Genes: %d\n", nrow(counts)))

  # Filter low-count genes
  keep <- rowSums(counts >= 10) >= 3
  counts <- counts[keep, ]
  cat(sprintf("  After filtering (>=10 in >=3 samples): %d genes\n", nrow(counts)))

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

# ============================================================
# Dataset configurations
# ============================================================

# --- GSE126543: TDP43_positive vs TDP43_negative (NeuN+ sorted nuclei) ---
cat("Loading GSE126543...\n")
gse126543 <- load_featurecounts("/mnt/HDD2/sorted_nuclei_GSE126543/counts/gene_counts.txt")
sample_sheet_126 <- read.delim("/mnt/HDD2/sorted_nuclei_GSE126543/sample_sheet.tsv")
# Exclude Unsorted samples
sorted_samples <- sample_sheet_126[sample_sheet_126$NeuN == "NeuN+", ]
counts_126 <- gse126543$counts[, sorted_samples$SRR]
condition_126 <- ifelse(sorted_samples$TDP43_status == "TDP43_positive", "TDP43pos", "TDP43neg")
run_deseq2(counts_126, condition_126, "gse126543_tdp43pos_vs_neg")

# --- GSE116622: ALS vs Control (cerebellum, hnRNP H score-based) ---
cat("\nLoading GSE116622...\n")
gse116622 <- load_featurecounts("/mnt/HDD2/ALS_patient_GSE116622/counts/gene_counts.txt")
sample_sheet_116 <- read.delim("/mnt/HDD2/ALS_patient_GSE116622/sample_sheet.tsv")
# Use Control vs ALS (exclude SOD1/C9/FTD for clean comparison)
# Condition column: Control, ALS, SOD1, C9, FTD, ALS_C9high, ALS_C9low
ctrl_116 <- sample_sheet_116$SRR[sample_sheet_116$Condition == "Control"]
als_116 <- sample_sheet_116$SRR[sample_sheet_116$Condition == "ALS"]
srrs_116 <- c(ctrl_116, als_116)
counts_116 <- gse116622$counts[, srrs_116]
condition_116 <- c(rep("Control", length(ctrl_116)), rep("ALS", length(als_116)))
run_deseq2(counts_116, condition_116, "gse116622_als_vs_ctrl")

# --- GSE124439: ALS vs Control (cortex) ---
cat("\nLoading GSE124439...\n")
gse124439 <- load_featurecounts("/mnt/HDD2/ALS_patient_GSE124439/counts/gene_counts.txt")
sample_sheet_124 <- read.delim("/mnt/HDD2/ALS_patient_GSE124439/sample_sheet.tsv")
# Check sample sheet structure
cat("  Sample sheet columns: ", paste(colnames(sample_sheet_124), collapse = ", "), "\n")
# Condition column: ALS, Control; also has ALS_Subtype, Tissue, C9orf72
ctrl_124 <- sample_sheet_124$SRR[sample_sheet_124$Condition == "Control"]
als_124 <- sample_sheet_124$SRR[sample_sheet_124$Condition == "ALS"]
srrs_124 <- c(ctrl_124, als_124)
srrs_124 <- srrs_124[srrs_124 %in% colnames(gse124439$counts)]
counts_124 <- gse124439$counts[, srrs_124]
condition_124 <- ifelse(srrs_124 %in% ctrl_124, "Control", "ALS")
run_deseq2(counts_124, condition_124, "gse124439_als_vs_ctrl")

cat("\n=== DESeq2 unified analysis complete ===\n")
