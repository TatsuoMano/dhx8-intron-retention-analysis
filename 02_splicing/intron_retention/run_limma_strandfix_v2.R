# Limma IR analysis with OPTIMIZED flanking exon aggregation
# Uses sparse matrix mapping instead of per-intron for loop

suppressPackageStartupMessages({ library(limma); library(edgeR); library(Matrix) })

PREV_DIR <- "/mnt/HDD2/TDP43_comparison/event_level_20260225"
WORKDIR <- "/mnt/HDD2/TDP43_comparison/validation_20260226"

load_featurecounts <- function(filepath) {
    df <- read.delim(filepath, comment.char = "#", stringsAsFactors = FALSE)
    rownames(df) <- df$Geneid
    lengths <- setNames(df$Length, df$Geneid)
    anno_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    count_cols <- setdiff(colnames(df), anno_cols)
    counts <- df[, count_cols, drop = FALSE]
    colnames(counts) <- gsub(".*\\/", "", colnames(counts))
    colnames(counts) <- gsub("_sorted\\.bam$", "", colnames(counts))
    colnames(counts) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts))
    list(counts = counts, lengths = lengths)
}

run_limma_ir <- function(intron_file, exon_file, label) {
    cat(sprintf("\n=== %s ===\n", label))
    
    intron_data <- load_featurecounts(intron_file)
    exon_data <- load_featurecounts(exon_file)
    
    # Classify samples
    snames <- colnames(intron_data$counts)
    ctrl <- grep("siCont|Cont", snames, value=TRUE)
    kd <- grep("siDHX8|DHX8", snames, value=TRUE)
    all_samples <- c(ctrl, kd)
    cat(sprintf("  Ctrl: %d, KD: %d\n", length(ctrl), length(kd)))
    
    # Load metadata + Welch-tested introns
    meta <- read.delim(file.path(PREV_DIR, "intron_metadata.tsv"), stringsAsFactors=FALSE)
    welch_res <- read.delim(file.path(PREV_DIR, "ir_results_all_dhx8.tsv"), stringsAsFactors=FALSE)
    
    meta <- meta[meta$intron_id %in% welch_res$intron_id, ]
    meta <- meta[!is.na(meta$flanking_exon_ids) & meta$flanking_exon_ids != "NA" & nchar(meta$flanking_exon_ids) > 0, ]
    valid_ids <- intersect(meta$intron_id, rownames(intron_data$counts))
    meta <- meta[meta$intron_id %in% valid_ids, ]
    
    intron_counts <- as.matrix(intron_data$counts[meta$intron_id, all_samples, drop=FALSE])
    intron_lens <- intron_data$lengths[meta$intron_id]
    cat(sprintf("  Introns to process: %d\n", nrow(intron_counts)))
    
    # === OPTIMIZED: build sparse intron->exon mapping matrix ===
    cat("  Building sparse exon mapping (vectorized)...\n")
    t0 <- proc.time()
    
    exon_splits <- strsplit(meta$flanking_exon_ids, ",")
    exon_rownames_set <- rownames(exon_data$counts)
    
    # Build triplets for sparse matrix
    ii <- integer(0)  # intron index
    jj <- integer(0)  # exon index
    
    exon_to_idx <- setNames(seq_along(exon_rownames_set), exon_rownames_set)
    
    for (i in seq_along(exon_splits)) {
        eids <- exon_splits[[i]]
        eids <- eids[eids %in% exon_rownames_set]
        if (length(eids) > 0) {
            ii <- c(ii, rep(i, length(eids)))
            jj <- c(jj, exon_to_idx[eids])
        }
    }
    
    # Sparse matrix: n_introns x n_exons, entries = 1
    mapping <- sparseMatrix(i=ii, j=jj, x=1,
                            dims=c(nrow(meta), length(exon_rownames_set)))
    
    # One matrix multiply: exon_sums = mapping %*% exon_counts_matrix
    exon_counts_mat <- as.matrix(exon_data$counts[exon_rownames_set, all_samples, drop=FALSE])
    exon_sums <- as.matrix(mapping %*% exon_counts_mat)
    rownames(exon_sums) <- meta$intron_id
    
    # Exon length sums: mapping %*% lengths
    exon_len_vec <- exon_data$lengths[exon_rownames_set]
    exon_len_sums <- as.numeric(mapping %*% exon_len_vec)
    names(exon_len_sums) <- meta$intron_id
    
    elapsed <- (proc.time() - t0)[3]
    cat(sprintf("  Exon aggregation done in %.1f seconds\n", elapsed))
    
    # === IR ratio (same as original) ===
    int_density <- sweep(intron_counts, 1, intron_lens, "/")
    ex_density <- sweep(exon_sums, 1, exon_len_sums, "/")
    ir_ratio_mat <- log2((int_density + 1e-6) / (ex_density + 1e-6))
    
    bad <- (intron_lens <= 0) | (exon_len_sums <= 0)
    ir_ratio_mat[bad, ] <- NA
    valid_rows <- apply(ir_ratio_mat, 1, function(x) all(is.finite(x)))
    ir_ratio_mat <- ir_ratio_mat[valid_rows, , drop=FALSE]
    cat(sprintf("  Valid introns for limma: %d\n", nrow(ir_ratio_mat)))
    
    # === limma (same as original) ===
    condition <- factor(ifelse(all_samples %in% ctrl, "ctrl", "kd"), levels=c("ctrl","kd"))
    design <- model.matrix(~condition)
    colnames(design) <- c("Intercept", "kd_vs_ctrl")
    
    fit <- lmFit(ir_ratio_mat, design)
    fit <- eBayes(fit)
    results <- topTable(fit, coef="kd_vs_ctrl", number=Inf, sort.by="none")
    
    final_ids <- rownames(ir_ratio_mat)
    out <- data.frame(
        intron_id = final_ids,
        logFC = results$logFC,
        avg_expr = results$AveExpr,
        t_stat = results$t,
        pvalue = results$P.Value,
        adj_p = results$adj.P.Val,
        stringsAsFactors = FALSE
    )
    
    parts <- strsplit(out$intron_id, ":")
    out$gene_name <- sapply(parts, `[`, 4)
    out$retained <- (out$logFC > 0) & (out$adj_p < 0.05)
    
    n_sig <- sum(out$adj_p < 0.05, na.rm=TRUE)
    n_retained <- sum(out$retained, na.rm=TRUE)
    n_down <- sum(out$adj_p < 0.05 & out$logFC < 0, na.rm=TRUE)
    n_genes <- length(unique(out$gene_name[out$retained]))
    
    cat(sprintf("\n  Results:\n"))
    cat(sprintf("    Introns tested: %d\n", nrow(out)))
    cat(sprintf("    FDR<0.05: %d (UP=%d, DOWN=%d)\n", n_sig, n_retained, n_down))
    cat(sprintf("    Retained genes: %d\n", n_genes))
    
    return(out)
}

# === Main ===
cat("=== Limma strandfix v2 (sparse matrix optimization) ===\n")

cat("\n--- s0 paired (unstranded) ---\n")
res_s0 <- run_limma_ir(
    "/mnt/HDD2/Sakai_RNA-seq/counts_s0_paired/intron_counts.txt",
    "/mnt/HDD2/Sakai_RNA-seq/counts_s0_paired/exon_counts.txt",
    "s0 paired (unstranded)")

cat("\n--- s2 (reverse stranded) ---\n")
res_s2 <- run_limma_ir(
    "/mnt/HDD2/Sakai_RNA-seq/counts_s2/intron_counts.txt",
    "/mnt/HDD2/Sakai_RNA-seq/counts_s2/exon_counts.txt",
    "s2 (reverse stranded)")

# Save s2 results
write.table(res_s2, file.path(WORKDIR, "limma_ir_results_dhx8_s2.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(res_s2[res_s2$retained, ], file.path(WORKDIR, "limma_retained_dhx8_s2.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Also save s0 paired for comparison
write.table(res_s0, file.path(WORKDIR, "limma_ir_results_dhx8_s0_paired.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

# === Comparison ===
cat("\n=== BEFORE/AFTER COMPARISON ===\n")
s0_retained <- res_s0$intron_id[res_s0$retained]
s2_retained <- res_s2$intron_id[res_s2$retained]
overlap <- intersect(s0_retained, s2_retained)
cat(sprintf("  s0 retained: %d introns, %d genes\n", length(s0_retained),
            length(unique(res_s0$gene_name[res_s0$retained]))))
cat(sprintf("  s2 retained: %d introns, %d genes\n", length(s2_retained),
            length(unique(res_s2$gene_name[res_s2$retained]))))
cat(sprintf("  Overlap: %d introns (%0.1f%% of s0, %0.1f%% of s2)\n",
            length(overlap),
            100*length(overlap)/max(length(s0_retained),1),
            100*length(overlap)/max(length(s2_retained),1)))

s0_genes <- unique(res_s0$gene_name[res_s0$retained])
s2_genes <- unique(res_s2$gene_name[res_s2$retained])
gene_overlap <- intersect(s0_genes, s2_genes)
cat(sprintf("\n  Gene level:\n"))
cat(sprintf("    s0: %d genes\n", length(s0_genes)))
cat(sprintf("    s2: %d genes\n", length(s2_genes)))
cat(sprintf("    Overlap: %d genes (%0.1f%%)\n", length(gene_overlap),
            100*length(gene_overlap)/max(length(s0_genes),1)))

cat("\n=== DONE ===\n")
