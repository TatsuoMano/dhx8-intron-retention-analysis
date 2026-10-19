#!/usr/bin/env Rscript
# Limma eBayes IR analysis for mouse TDP-43 datasets
# GSE27218 (-s 2) and GSE65973 (-s 1) with correct strandedness counts

library(limma)
library(Matrix)

OUTDIR <- "/mnt/HDD2/unified_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Intron metadata (mm10)
META_PATH <- "/mnt/HDD2/TDP43_comparison/event_level_20260225/intron_metadata.tsv"

load_featurecounts <- function(path) {
  cat(sprintf("  Loading %s...\n", basename(path)))
  df <- read.delim(path, comment.char = "#", check.names = FALSE)
  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta_cols)
  counts <- as.matrix(df[, count_cols])
  rownames(counts) <- df$Geneid
  lengths <- df$Length
  names(lengths) <- df$Geneid
  clean <- gsub(".*/([^/]+)_Aligned.*", "\\1", colnames(counts))
  colnames(counts) <- clean
  list(counts = counts, lengths = lengths)
}

# Check if metadata exists
if (!file.exists(META_PATH)) {
  cat("Intron metadata not found, building from SAF...\n")
  # Build metadata from introns.saf and flanking_exons.saf
  intron_saf <- read.delim("/mnt/HDD2/TDP43_comparison/event_level_20260225/introns.saf")
  exon_saf <- read.delim("/mnt/HDD2/TDP43_comparison/event_level_20260225/flanking_exons.saf")

  # Parse intron IDs: chr:start-end:strand:gene
  parts <- strsplit(intron_saf$GeneID, ":")
  meta <- data.frame(
    intron_id = intron_saf$GeneID,
    chr = sapply(parts, `[`, 1),
    strand = sapply(parts, `[`, 3),
    gene_name = sapply(parts, `[`, 4),
    intron_length = intron_saf$End - intron_saf$Start + 1,
    stringsAsFactors = FALSE
  )
  se <- strsplit(sapply(parts, `[`, 2), "-")
  meta$start <- as.integer(sapply(se, `[`, 1))
  meta$end <- as.integer(sapply(se, `[`, 2))

  # Map flanking exons: exon_saf IDs are chr:start-end:strand
  # For each intron, find exons that flank it (end+1 = intron_start or start-1 = intron_end)
  # Simple approach: match by chr+strand, then adjacency
  exon_parts <- strsplit(exon_saf$GeneID, ":")
  exon_df <- data.frame(
    exon_id = exon_saf$GeneID,
    chr = sapply(exon_parts, `[`, 1),
    strand = sapply(exon_parts, `[`, 3),
    stringsAsFactors = FALSE
  )
  exon_se <- strsplit(sapply(exon_parts, `[`, 2), "-")
  exon_df$start <- as.integer(sapply(exon_se, `[`, 1))
  exon_df$end <- as.integer(sapply(exon_se, `[`, 2))

  # For each intron, find flanking exons (within 2bp tolerance)
  cat("  Mapping flanking exons...\n")
  flanking <- character(nrow(meta))
  for (i in seq_len(nrow(meta))) {
    chr_strand <- exon_df$chr == meta$chr[i] & exon_df$strand == meta$strand[i]
    left <- which(chr_strand & abs(exon_df$end - (meta$start[i] - 1)) <= 2)
    right <- which(chr_strand & abs(exon_df$start - (meta$end[i] + 1)) <= 2)
    ids <- unique(exon_df$exon_id[c(left, right)])
    flanking[i] <- paste(ids, collapse = ",")
    if (i %% 10000 == 0) cat(sprintf("    %d/%d\n", i, nrow(meta)))
  }
  meta$flanking_exon_ids <- flanking

  write.table(meta, META_PATH, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Saved metadata: %s\n", META_PATH))
} else {
  cat("Loading intron metadata...\n")
  meta <- read.delim(META_PATH, stringsAsFactors = FALSE)
}

meta <- meta[!is.na(meta$flanking_exon_ids) &
             meta$flanking_exon_ids != "NA" &
             meta$flanking_exon_ids != "" &
             nchar(meta$flanking_exon_ids) > 0, ]
cat(sprintf("  Introns with flanking exons: %d\n", nrow(meta)))

# Dataset configs
datasets <- list(
  gse27218 = list(
    intron_counts = "/mnt/HDD2/TDP43_comparison/GSE27218/counts_s2/intron_counts.txt",
    exon_counts   = "/mnt/HDD2/TDP43_comparison/GSE27218/counts_s2/exon_counts.txt",
    sample_sheet  = "/mnt/HDD2/TDP43_comparison/GSE27218/sample_sheet.tsv",
    ctrl_val = "control", test_val = "tdp43_kd",
    label = "GSE27218 (TDP-43 ASO KD, -s 2)"
  ),
  gse65973 = list(
    intron_counts = "/mnt/HDD2/TDP43_comparison/GSE65973/counts_s1/intron_counts.txt",
    exon_counts   = "/mnt/HDD2/TDP43_comparison/GSE65973/counts_s1/exon_counts.txt",
    sample_sheet  = "/mnt/HDD2/TDP43_comparison/GSE65973/sample_sheet.tsv",
    ctrl_val = "control", test_val = "bigenic",
    label = "GSE65973 (dNLS-TDP-43, -s 1)"
  )
)

for (ds_name in names(datasets)) {
  config <- datasets[[ds_name]]

  cat(sprintf("\n============================================================\n"))
  cat(sprintf("Processing: %s\n", config$label))
  cat(sprintf("============================================================\n"))

  if (!file.exists(config$intron_counts) || !file.exists(config$exon_counts)) {
    cat("  SKIPPED: count files not found\n")
    next
  }

  intron_data <- load_featurecounts(config$intron_counts)
  exon_data <- load_featurecounts(config$exon_counts)

  # Sample sheet (no header)
  ss <- read.delim(config$sample_sheet, header = FALSE, stringsAsFactors = FALSE)
  colnames(ss) <- c("SRR", "GSM", "Sample", "Condition")

  ctrl_srrs <- ss$SRR[ss$Condition == config$ctrl_val]
  test_srrs <- ss$SRR[ss$Condition == config$test_val]

  all_samples <- intersect(colnames(intron_data$counts), colnames(exon_data$counts))
  ctrl_srrs <- intersect(ctrl_srrs, all_samples)
  test_srrs <- intersect(test_srrs, all_samples)

  cat(sprintf("  Ctrl: %d, Test: %d\n", length(ctrl_srrs), length(test_srrs)))

  if (length(ctrl_srrs) < 2 || length(test_srrs) < 2) {
    cat("  SKIPPED: insufficient samples\n")
    next
  }

  sample_cols <- c(ctrl_srrs, test_srrs)

  common_introns <- intersect(rownames(intron_data$counts), meta$intron_id)
  cat(sprintf("  Introns in counts & metadata: %d\n", length(common_introns)))

  intron_sub <- intron_data$counts[common_introns, sample_cols]
  expressed <- rowSums(intron_sub >= 10) >= 3
  common_introns <- common_introns[expressed]
  cat(sprintf("  After expression filter: %d\n", length(common_introns)))

  meta_sub <- meta[meta$intron_id %in% common_introns, ]
  meta_sub <- meta_sub[order(match(meta_sub$intron_id, common_introns)), ]

  # Sparse matrix for exon aggregation
  exon_rownames_set <- rownames(exon_data$counts)
  exon_name_to_idx <- setNames(seq_along(exon_rownames_set), exon_rownames_set)

  all_splits <- strsplit(meta_sub$flanking_exon_ids, ",", fixed = TRUE)
  n_per_intron <- lengths(all_splits)
  all_exons <- unlist(all_splits, use.names = FALSE)
  intron_idx <- rep(seq_len(nrow(meta_sub)), n_per_intron)

  valid <- all_exons %in% exon_rownames_set
  all_exons <- all_exons[valid]
  intron_idx <- intron_idx[valid]
  exon_idx <- exon_name_to_idx[all_exons]

  mapping <- sparseMatrix(i = intron_idx, j = exon_idx,
                          dims = c(nrow(meta_sub), length(exon_rownames_set)))

  exon_sums <- as.matrix(mapping %*% exon_data$counts[exon_rownames_set, sample_cols])
  exon_len_sums <- as.numeric(mapping %*% exon_data$lengths[exon_rownames_set])

  # IR ratios
  intron_counts_sub <- intron_data$counts[meta_sub$intron_id, sample_cols]
  intron_lengths <- intron_data$lengths[meta_sub$intron_id]

  intron_density <- sweep(intron_counts_sub, 1, intron_lengths, "/")
  exon_density <- sweep(exon_sums, 1, pmax(exon_len_sums, 1), "/")
  ir_ratios <- intron_density / pmax(exon_density, 1e-10)

  valid_rows <- rowSums(exon_sums >= 10) >= 3 & rowSums(!is.na(ir_ratios) & is.finite(ir_ratios)) >= 3
  ir_ratios <- ir_ratios[valid_rows, ]
  cat(sprintf("  Valid introns for limma: %d\n", sum(valid_rows)))

  # limma
  ir_log <- log2(ir_ratios + 1e-6)
  condition <- factor(c(rep("ctrl", length(ctrl_srrs)), rep("test", length(test_srrs))))
  design <- model.matrix(~ condition)
  fit <- lmFit(ir_log, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  tt$intron_id <- rownames(tt)
  intron_info <- meta_sub[match(tt$intron_id, meta_sub$intron_id),
                          c("chr", "start", "end", "strand", "gene_name", "intron_length")]
  tt <- cbind(tt, intron_info)

  colnames(tt)[colnames(tt) == "AveExpr"] <- "avg_expr"
  colnames(tt)[colnames(tt) == "t"] <- "t_stat"
  colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "adj_p"
  tt$retained <- tt$logFC > 0 & tt$adj_p < 0.05

  sig <- sum(tt$adj_p < 0.05, na.rm = TRUE)
  retained <- sum(tt$retained, na.rm = TRUE)
  retained_genes <- length(unique(tt$gene_name[tt$retained]))
  cat(sprintf("\n  Results:\n"))
  cat(sprintf("    Introns tested: %d\n", nrow(tt)))
  cat(sprintf("    Significant (adj_p < 0.05): %d\n", sig))
  cat(sprintf("    Retained (logFC > 0, adj_p < 0.05): %d\n", retained))
  cat(sprintf("    Retained genes: %d\n", retained_genes))

  all_path <- file.path(OUTDIR, sprintf("limma_ir_results_%s.tsv", ds_name))
  ret_path <- file.path(OUTDIR, sprintf("limma_retained_%s.tsv", ds_name))
  write.table(tt, all_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(tt[tt$retained, ], ret_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("    Saved: %s\n", all_path))
}

cat("\n=== Mouse limma eBayes IR complete ===\n")
