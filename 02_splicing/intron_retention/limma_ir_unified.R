#!/usr/bin/env Rscript
# Unified limma eBayes IR analysis for human ALS datasets
# Uses unified hg38 SAF annotations (introns + flanking exons + metadata)
# Adapted from validation_step1_limma.R

library(limma)
library(Matrix)

# ============================================================
# Configuration
# ============================================================
ANNO_DIR <- "/mnt/HDD2/hg38_ir_annotations"
OUTDIR <- "/mnt/HDD2/unified_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

META_PATH <- file.path(ANNO_DIR, "hg38_intron_metadata.tsv")

datasets <- list(
  gse126543 = list(
    intron_counts = "/mnt/HDD2/sorted_nuclei_GSE126543/counts/intron_counts.txt",
    exon_counts   = "/mnt/HDD2/sorted_nuclei_GSE126543/counts/exon_counts.txt",
    sample_sheet  = "/mnt/HDD2/sorted_nuclei_GSE126543/sample_sheet.tsv",
    ctrl_col = "TDP43_status", ctrl_val = "TDP43_negative",
    kd_col   = "TDP43_status", kd_val   = "TDP43_positive",
    exclude_col = "NeuN", exclude_val = "Unsorted",
    label = "GSE126543 (sorted nuclei, TDP43pos vs neg)"
  ),
  gse116622 = list(
    intron_counts = "/mnt/HDD2/ALS_patient_GSE116622/counts/intron_counts_unified.txt",
    exon_counts   = "/mnt/HDD2/ALS_patient_GSE116622/counts/exon_counts_unified.txt",
    sample_sheet  = "/mnt/HDD2/ALS_patient_GSE116622/sample_sheet.tsv",
    ctrl_col = "Condition", ctrl_val = "Control",
    kd_col   = "Condition", kd_val   = "ALS",
    exclude_col = NULL, exclude_val = NULL,
    label = "GSE116622 (cerebellum, ALS vs Control)"
  ),
  gse124439 = list(
    intron_counts = "/mnt/HDD2/ALS_patient_GSE124439/counts/intron_counts_unified.txt",
    exon_counts   = "/mnt/HDD2/ALS_patient_GSE124439/counts/exon_counts_unified.txt",
    sample_sheet  = "/mnt/HDD2/ALS_patient_GSE124439/sample_sheet.tsv",
    ctrl_col = "Condition", ctrl_val = "Control",
    kd_col   = "Condition", kd_val   = "ALS",
    exclude_col = NULL, exclude_val = NULL,
    label = "GSE124439 (cortex, ALS vs Control)"
  )
)

# ============================================================
# Helper functions
# ============================================================
load_featurecounts <- function(path) {
  cat(sprintf("  Loading %s...\n", basename(path)))
  df <- read.delim(path, comment.char = "#", check.names = FALSE)
  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta_cols)
  counts <- as.matrix(df[, count_cols])
  rownames(counts) <- df$Geneid
  lengths <- df$Length
  names(lengths) <- df$Geneid
  # Clean column names
  clean <- gsub(".*/([^/]+)\\.sorted\\.bam$", "\\1", colnames(counts))
  clean <- gsub(".*/([^/]+)_Aligned.*", "\\1", clean)
  colnames(counts) <- clean
  list(counts = counts, lengths = lengths)
}

# ============================================================
# Load metadata
# ============================================================
cat("Loading intron metadata...\n")
meta <- read.delim(META_PATH, stringsAsFactors = FALSE)
meta <- meta[!is.na(meta$flanking_exon_ids) &
             meta$flanking_exon_ids != "NA" &
             nchar(meta$flanking_exon_ids) > 0, ]
cat(sprintf("  Introns with flanking exons: %d\n", nrow(meta)))

# ============================================================
# Process each dataset
# ============================================================
for (ds_name in names(datasets)) {
  config <- datasets[[ds_name]]

  cat(sprintf("\n============================================================\n"))
  cat(sprintf("Processing: %s\n", config$label))
  cat(sprintf("============================================================\n"))

  # Check if files exist
  if (!file.exists(config$intron_counts)) {
    cat(sprintf("  SKIPPED: %s not found\n", config$intron_counts))
    next
  }
  if (!file.exists(config$exon_counts)) {
    cat(sprintf("  SKIPPED: %s not found\n", config$exon_counts))
    next
  }

  # Load counts
  intron_data <- load_featurecounts(config$intron_counts)
  exon_data <- load_featurecounts(config$exon_counts)

  # Load sample sheet and determine groups
  ss <- read.delim(config$sample_sheet, stringsAsFactors = FALSE)

  # Exclude samples if needed
  if (!is.null(config$exclude_col)) {
    ss <- ss[ss[[config$exclude_col]] != config$exclude_val, ]
  }

  ctrl_srrs <- ss$SRR[ss[[config$ctrl_col]] == config$ctrl_val]
  kd_srrs <- ss$SRR[ss[[config$kd_col]] == config$kd_val]

  # Match to available samples
  all_samples <- intersect(colnames(intron_data$counts), colnames(exon_data$counts))
  ctrl_srrs <- intersect(ctrl_srrs, all_samples)
  kd_srrs <- intersect(kd_srrs, all_samples)

  cat(sprintf("  Ctrl: %d, Test: %d\n", length(ctrl_srrs), length(kd_srrs)))

  if (length(ctrl_srrs) < 2 || length(kd_srrs) < 2) {
    cat("  SKIPPED: insufficient samples\n")
    next
  }

  sample_cols <- c(ctrl_srrs, kd_srrs)

  # Filter to introns present in metadata AND with sufficient counts
  common_introns <- intersect(rownames(intron_data$counts), meta$intron_id)
  cat(sprintf("  Introns in counts & metadata: %d\n", length(common_introns)))

  # Pre-filter: require intron count >= 10 in at least 3 samples
  intron_sub <- intron_data$counts[common_introns, sample_cols]
  expressed <- rowSums(intron_sub >= 10) >= 3
  common_introns <- common_introns[expressed]
  cat(sprintf("  After expression filter (>=10 in >=3 samples): %d\n", length(common_introns)))

  meta_sub <- meta[meta$intron_id %in% common_introns, ]
  meta_sub <- meta_sub[order(match(meta_sub$intron_id, common_introns)), ]

  # Aggregate flanking exon counts using vectorized sparse matrix
  cat("  Aggregating flanking exon counts (vectorized)...\n")
  exon_rownames_set <- rownames(exon_data$counts)
  exon_name_to_idx <- setNames(seq_along(exon_rownames_set), exon_rownames_set)

  # Vectorized unnesting: split all at once, then build mapping
  all_exon_str <- meta_sub$flanking_exon_ids
  all_splits <- strsplit(all_exon_str, ",", fixed = TRUE)
  n_per_intron <- lengths(all_splits)
  all_exons <- unlist(all_splits, use.names = FALSE)
  intron_idx <- rep(seq_len(nrow(meta_sub)), n_per_intron)

  # Filter to exons present in count matrix
  valid <- all_exons %in% exon_rownames_set
  all_exons <- all_exons[valid]
  intron_idx <- intron_idx[valid]
  exon_idx <- exon_name_to_idx[all_exons]

  mapping <- sparseMatrix(i = intron_idx, j = exon_idx,
                          dims = c(nrow(meta_sub), length(exon_rownames_set)),
                          dimnames = list(meta_sub$intron_id, exon_rownames_set))
  cat(sprintf("    Sparse matrix: %d x %d, %d nonzero\n",
              nrow(mapping), ncol(mapping), length(intron_idx)))

  # Matrix multiply: intron x sample = mapping %*% exon_counts
  exon_sums <- as.matrix(mapping %*% exon_data$counts[exon_rownames_set, sample_cols])
  exon_len_sums <- as.numeric(mapping %*% exon_data$lengths[exon_rownames_set])
  cat("    Aggregation complete.\n")

  # Compute IR ratios
  cat("  Computing IR ratios (vectorized)...\n")
  intron_counts_sub <- intron_data$counts[meta_sub$intron_id, sample_cols]
  intron_lengths <- intron_data$lengths[meta_sub$intron_id]

  intron_density <- sweep(intron_counts_sub, 1, intron_lengths, "/")
  exon_density <- sweep(exon_sums, 1, pmax(exon_len_sums, 1), "/")

  ir_ratios <- intron_density / pmax(exon_density, 1e-10)

  # Filter: require some expression
  valid <- rowSums(exon_sums >= 10) >= 3 & rowSums(!is.na(ir_ratios) & is.finite(ir_ratios)) >= 3
  ir_ratios <- ir_ratios[valid, ]
  cat(sprintf("  Valid introns for limma: %d\n", sum(valid)))

  # limma
  cat("  Fitting limma model...\n")
  ir_log <- log2(ir_ratios + 1e-6)
  condition <- factor(c(rep("ctrl", length(ctrl_srrs)), rep("test", length(kd_srrs))))
  design <- model.matrix(~ condition)
  fit <- lmFit(ir_log, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

  # Add intron info
  tt$intron_id <- rownames(tt)
  intron_info <- meta_sub[match(tt$intron_id, meta_sub$intron_id),
                          c("chr", "start", "end", "strand", "gene_name", "intron_length")]
  tt <- cbind(tt, intron_info)

  # Rename columns for consistency
  colnames(tt)[colnames(tt) == "logFC"] <- "logFC"
  colnames(tt)[colnames(tt) == "AveExpr"] <- "avg_expr"
  colnames(tt)[colnames(tt) == "t"] <- "t_stat"
  colnames(tt)[colnames(tt) == "P.Value"] <- "pvalue"
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "adj_p"

  # Retained = logFC > 0 (IR increased in test group)
  tt$retained <- tt$logFC > 0 & tt$adj_p < 0.05

  sig <- sum(tt$adj_p < 0.05, na.rm = TRUE)
  retained <- sum(tt$retained, na.rm = TRUE)
  retained_genes <- length(unique(tt$gene_name[tt$retained]))
  cat(sprintf("\n  Results:\n"))
  cat(sprintf("    Introns tested: %d\n", nrow(tt)))
  cat(sprintf("    Significant (adj_p < 0.05): %d\n", sig))
  cat(sprintf("    Retained (logFC > 0, adj_p < 0.05): %d\n", retained))
  cat(sprintf("    Retained genes: %d\n", retained_genes))
  cat(sprintf("    Prior df (limma): %.1f\n", fit$df.prior))
  cat(sprintf("    Prior var (s2.prior): %.4f\n", fit$s2.prior))

  # Save results
  all_path <- file.path(OUTDIR, sprintf("limma_ir_results_%s.tsv", ds_name))
  ret_path <- file.path(OUTDIR, sprintf("limma_retained_%s.tsv", ds_name))
  write.table(tt, all_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(tt[tt$retained, ], ret_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("    Saved: %s\n", all_path))
}

cat("\n=== limma eBayes IR unified analysis complete ===\n")
