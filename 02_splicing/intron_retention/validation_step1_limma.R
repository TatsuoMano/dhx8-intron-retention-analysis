#!/usr/bin/env Rscript
# Validation Step 1b: limma empirical Bayes moderated t-test for IR
#
# Strategy: Load per-intron & flanking exon counts only for introns that
# were tested in the Welch analysis (ir_results_all). Compute per-sample
# IR ratios, then apply limma eBayes for variance sharing.
#
# Usage: Rscript validation_step1_limma.R
# Runs on: ncnp-linux

suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
})

# === Configuration ===
WORKDIR <- "/mnt/HDD2/TDP43_comparison/validation_20260226"
PREV_DIR <- "/mnt/HDD2/TDP43_comparison/event_level_20260225"

datasets <- list(
    dhx8 = list(
        intron_counts = file.path(PREV_DIR, "dhx8_intron_counts.txt"),
        exon_counts   = file.path(PREV_DIR, "dhx8_exon_counts.txt"),
        control_pattern = "siCont",
        kd_pattern      = "siDHX8"
    ),
    gse27218 = list(
        intron_counts = file.path(PREV_DIR, "gse27218_intron_counts.txt"),
        exon_counts   = file.path(PREV_DIR, "gse27218_exon_counts.txt"),
        control_ids = c(paste0("SRR10706", 1:9), "SRR107070", "SRR107071"),
        kd_ids      = c(paste0("SRR10707", 2:9), "SRR107080", "SRR107081")
    ),
    gse65973 = list(
        intron_counts = file.path(PREV_DIR, "gse65973_intron_counts.txt"),
        exon_counts   = file.path(PREV_DIR, "gse65973_exon_counts.txt"),
        control_ids = paste0("SRR180587", 1:4),
        kd_ids      = c(paste0("SRR180586", 7:9), "SRR1805870")
    )
)

load_featurecounts <- function(filepath) {
    df <- read.delim(filepath, comment.char = "#", stringsAsFactors = FALSE)
    rownames(df) <- df$Geneid
    lengths <- setNames(df$Length, df$Geneid)
    anno_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
    count_cols <- setdiff(colnames(df), anno_cols)
    counts <- df[, count_cols, drop = FALSE]
    colnames(counts) <- gsub(".*\\/", "", colnames(counts))
    colnames(counts) <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts))
    list(counts = counts, lengths = lengths)
}

classify_samples <- function(sample_names, config) {
    if (!is.null(config$control_pattern)) {
        ctrl <- grep(config$control_pattern, sample_names, value = TRUE)
        kd <- grep(config$kd_pattern, sample_names, value = TRUE)
    } else {
        ctrl <- sample_names[sapply(sample_names, function(s)
            any(sapply(config$control_ids, function(id) grepl(id, s))))]
        kd <- sample_names[sapply(sample_names, function(s)
            any(sapply(config$kd_ids, function(id) grepl(id, s))))]
    }
    list(ctrl = ctrl, kd = kd)
}

process_dataset <- function(ds_name, config) {
    cat(sprintf("\n%s\nProcessing: %s\n%s\n",
                paste(rep("=", 60), collapse = ""), ds_name,
                paste(rep("=", 60), collapse = "")))

    # Load Welch results to get the list of tested introns
    welch_path <- file.path(PREV_DIR, paste0("ir_results_all_", ds_name, ".tsv"))
    welch_res <- read.delim(welch_path, stringsAsFactors = FALSE)
    tested_intron_ids <- welch_res$intron_id
    cat(sprintf("  Introns tested in Welch: %d\n", length(tested_intron_ids)))

    # Load count data
    intron_data <- load_featurecounts(config$intron_counts)
    exon_data <- load_featurecounts(config$exon_counts)

    # Classify samples
    samples <- classify_samples(colnames(intron_data$counts), config)
    cat(sprintf("  Ctrl: %d, KD: %d\n", length(samples$ctrl), length(samples$kd)))
    all_samples <- c(samples$ctrl, samples$kd)

    # Load metadata - only for tested introns
    meta <- read.delim(file.path(PREV_DIR, "intron_metadata.tsv"),
                       stringsAsFactors = FALSE)
    meta <- meta[meta$intron_id %in% tested_intron_ids, ]
    meta <- meta[!is.na(meta$flanking_exon_ids) &
                 meta$flanking_exon_ids != "NA" &
                 nchar(meta$flanking_exon_ids) > 0, ]
    cat(sprintf("  Metadata rows (tested + has flanking exons): %d\n", nrow(meta)))

    # Subset count matrices to tested introns
    valid_ids <- intersect(meta$intron_id, rownames(intron_data$counts))
    meta <- meta[meta$intron_id %in% valid_ids, ]

    intron_counts <- as.matrix(intron_data$counts[meta$intron_id, all_samples, drop = FALSE])
    intron_lens <- intron_data$lengths[meta$intron_id]
    cat(sprintf("  Introns to process: %d\n", nrow(intron_counts)))

    # Compute aggregate flanking exon counts for each intron (vectorized via sparse mapping)
    cat("  Aggregating flanking exon counts...\n")
    exon_splits <- strsplit(meta$flanking_exon_ids, ",")
    exon_rownames_set <- rownames(exon_data$counts)

    # Build exon sum matrix
    exon_sums <- matrix(0, nrow = nrow(meta), ncol = length(all_samples),
                        dimnames = list(meta$intron_id, all_samples))
    exon_len_sums <- numeric(nrow(meta))
    names(exon_len_sums) <- meta$intron_id

    # Process in chunks for progress reporting
    n <- nrow(meta)
    chunk <- 5000
    for (start in seq(1, n, by = chunk)) {
        end <- min(start + chunk - 1, n)
        for (i in start:end) {
            eids <- exon_splits[[i]]
            eids <- eids[eids %in% exon_rownames_set]
            if (length(eids) > 0) {
                exon_sums[i, ] <- colSums(exon_data$counts[eids, all_samples, drop = FALSE])
                exon_len_sums[i] <- sum(exon_data$lengths[eids])
            }
        }
        cat(sprintf("    %d/%d introns...\n", min(end, n), n))
    }

    # Compute log2(IR ratio) = log2(intron_density / exon_density)
    cat("  Computing IR ratios (vectorized)...\n")
    int_density <- sweep(intron_counts, 1, intron_lens, "/")
    ex_density <- sweep(exon_sums, 1, exon_len_sums, "/")

    ir_ratio_mat <- log2((int_density + 1e-6) / (ex_density + 1e-6))

    # Remove introns with bad lengths
    bad <- (intron_lens <= 0) | (exon_len_sums <= 0)
    ir_ratio_mat[bad, ] <- NA

    # Keep only finite rows
    valid_rows <- apply(ir_ratio_mat, 1, function(x) all(is.finite(x)))
    ir_ratio_mat <- ir_ratio_mat[valid_rows, , drop = FALSE]
    final_ids <- rownames(ir_ratio_mat)
    cat(sprintf("  Valid introns for limma: %d\n", length(final_ids)))

    if (length(final_ids) < 100) {
        cat("  Too few introns, skipping.\n")
        return(NULL)
    }

    # limma: empirical Bayes moderated t-test
    condition <- factor(
        ifelse(all_samples %in% samples$ctrl, "ctrl", "kd"),
        levels = c("ctrl", "kd")
    )
    design <- model.matrix(~ condition)
    colnames(design) <- c("Intercept", "kd_vs_ctrl")

    cat("  Fitting limma model...\n")
    fit <- lmFit(ir_ratio_mat, design)
    fit <- eBayes(fit)

    results <- topTable(fit, coef = "kd_vs_ctrl", number = Inf, sort.by = "none")

    # Build output table
    out <- data.frame(
        intron_id = final_ids,
        logFC = results$logFC,
        avg_expr = results$AveExpr,
        t_stat = results$t,
        pvalue = results$P.Value,
        adj_p = results$adj.P.Val,
        stringsAsFactors = FALSE
    )

    # Parse intron_id
    parts <- strsplit(out$intron_id, ":")
    out$chr <- sapply(parts, `[`, 1)
    coord_parts <- strsplit(sapply(parts, `[`, 2), "-")
    out$start <- as.integer(sapply(coord_parts, `[`, 1))
    out$end <- as.integer(sapply(coord_parts, `[`, 2))
    out$strand <- sapply(parts, `[`, 3)
    out$gene_name <- sapply(parts, `[`, 4)
    out$intron_length <- out$end - out$start + 1
    out$retained <- (out$logFC > 0) & (out$adj_p < 0.05)

    # Summary
    cat(sprintf("\n  Results:\n"))
    cat(sprintf("    Introns tested: %d\n", nrow(out)))
    cat(sprintf("    Significant (adj_p < 0.05): %d\n", sum(out$adj_p < 0.05, na.rm = TRUE)))
    cat(sprintf("    Retained (logFC > 0, adj_p < 0.05): %d\n", sum(out$retained, na.rm = TRUE)))
    if (sum(out$retained, na.rm = TRUE) > 0) {
        cat(sprintf("    Retained genes: %d\n",
                    length(unique(out$gene_name[out$retained]))))
    }
    cat(sprintf("    Prior df (limma): %.1f\n", fit$df.prior))
    cat(sprintf("    Prior var (s2.prior): %.4f\n", fit$s2.prior))

    # Save
    all_path <- file.path(WORKDIR, paste0("limma_ir_results_", ds_name, ".tsv"))
    write.table(out, all_path, sep = "\t", row.names = FALSE, quote = FALSE)

    ret_path <- file.path(WORKDIR, paste0("limma_retained_", ds_name, ".tsv"))
    write.table(out[out$retained, ], ret_path, sep = "\t",
                row.names = FALSE, quote = FALSE)

    cat(sprintf("    Saved: %s\n", all_path))
    return(out)
}


# === Main ===
cat("=== Validation Step 1b: limma moderated IR analysis ===\n")

all_results <- list()
for (ds_name in names(datasets)) {
    result <- tryCatch(
        process_dataset(ds_name, datasets[[ds_name]]),
        error = function(e) {
            cat(sprintf("  ERROR: %s\n", e$message))
            NULL
        }
    )
    if (!is.null(result)) {
        all_results[[ds_name]] <- result
    }
}

# Welch vs limma comparison
cat(sprintf("\n%s\nWelch vs limma p-value comparison:\n", paste(rep("=", 60), collapse = "")))
for (ds_name in names(all_results)) {
    limma_res <- all_results[[ds_name]]
    welch_path <- file.path(PREV_DIR, paste0("ir_results_all_", ds_name, ".tsv"))
    welch_res <- read.delim(welch_path, stringsAsFactors = FALSE)

    common <- intersect(limma_res$intron_id, welch_res$intron_id)
    if (length(common) > 100) {
        limma_p <- -log10(limma_res$pvalue[match(common, limma_res$intron_id)])
        welch_p <- -log10(welch_res$pvalue[match(common, welch_res$intron_id)])
        valid <- is.finite(limma_p) & is.finite(welch_p)
        if (sum(valid) > 50) {
            rho <- cor(limma_p[valid], welch_p[valid], method = "spearman")
            cat(sprintf("  %s: Spearman rho(-log10p) = %.3f (n=%d)\n",
                        ds_name, rho, sum(valid)))

            # Also compare retained intron counts at various thresholds
            for (fdr_thresh in c(0.05, 0.10, 0.20)) {
                n_limma <- sum(limma_res$adj_p < fdr_thresh & limma_res$logFC > 0, na.rm = TRUE)
                n_welch <- sum(welch_res$adj_p < fdr_thresh & welch_res$delta_ir > 0, na.rm = TRUE)
                cat(sprintf("    FDR<%g: limma=%d, Welch=%d\n", fdr_thresh, n_limma, n_welch))
            }
        }
    }
}

# Cross-dataset overlap
cat("\nCross-dataset overlap (limma FDR < 0.05):\n")
if ("dhx8" %in% names(all_results)) {
    dhx8_ret <- all_results[["dhx8"]]
    dhx8_coords <- with(dhx8_ret[dhx8_ret$retained, ],
                        paste0(chr, ":", start, "-", end, ":", strand))

    for (ds_name in c("gse27218", "gse65973")) {
        if (ds_name %in% names(all_results)) {
            tdp43_ret <- all_results[[ds_name]]
            tdp43_coords <- with(tdp43_ret[tdp43_ret$retained, ],
                                 paste0(chr, ":", start, "-", end, ":", strand))
            overlap <- length(intersect(dhx8_coords, tdp43_coords))
            cat(sprintf("  DHX8 vs %s: %d vs %d, overlap=%d\n",
                        ds_name, length(dhx8_coords), length(tdp43_coords), overlap))
        }
    }
}

cat("\n=== Step 1b complete ===\n")
