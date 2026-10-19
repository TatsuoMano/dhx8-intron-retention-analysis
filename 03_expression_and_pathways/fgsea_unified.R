#!/usr/bin/env Rscript
# Unified fgsea GSEA for all datasets
# Runs Hallmark, C2:CP, GO:BP for each DESeq2 result
# Auto-detects mouse (Mm) vs human (Hs) from gene names

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
})

OUTDIR <- "/mnt/HDD2/unified_results/gsea"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# All DESeq2 result files
deseq_files <- list(
  dhx8          = "/mnt/HDD2/unified_results/deseq2_dhx8_kd_vs_ctrl.tsv",
  gse27218      = "/mnt/HDD2/unified_results/deseq2_gse27218_tdp43kd_vs_ctrl.tsv",
  gse65973      = "/mnt/HDD2/unified_results/deseq2_gse65973_dnls_vs_ctrl.tsv",
  gse126543     = "/mnt/HDD2/unified_results/deseq2_gse126543_tdp43pos_vs_neg.tsv",
  gse116622     = "/mnt/HDD2/unified_results/deseq2_gse116622_als_vs_ctrl.tsv",
  gse124439     = "/mnt/HDD2/unified_results/deseq2_gse124439_als_vs_ctrl.tsv"
)

# Species mapping
species <- list(
  dhx8 = "mouse", gse27218 = "mouse", gse65973 = "mouse",
  gse126543 = "human", gse116622 = "human", gse124439 = "human"
)

# GMT URLs
gmt_config <- list(
  mouse = list(
    base_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Mm/",
    hallmark = "mh.all.v2024.1.Mm.symbols.gmt",
    c2cp     = "m2.cp.v2024.1.Mm.symbols.gmt",
    gobp     = "m5.go.bp.v2024.1.Mm.symbols.gmt"
  ),
  human = list(
    base_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/",
    hallmark = "h.all.v2024.1.Hs.symbols.gmt",
    c2cp     = "c2.cp.v2024.1.Hs.symbols.gmt",
    gobp     = "c5.go.bp.v2024.1.Hs.symbols.gmt"
  )
)

# Download GMT files
download_gmts <- function(sp) {
  conf <- gmt_config[[sp]]
  pathways <- list()
  for (name in c("hallmark", "c2cp", "gobp")) {
    local_file <- file.path(OUTDIR, conf[[name]])
    if (!file.exists(local_file)) {
      url <- paste0(conf$base_url, conf[[name]])
      cat(sprintf("  Downloading %s %s...\n", sp, name))
      tryCatch(
        download.file(url, local_file, quiet = TRUE, method = "auto"),
        error = function(e) {
          cat(sprintf("  ERROR: %s\n", e$message))
          download.file(url, local_file, quiet = TRUE, method = "curl")
        }
      )
    }
    pathways[[name]] <- gmtPathways(local_file)
    cat(sprintf("  %s %s: %d gene sets\n", sp, name, length(pathways[[name]])))
  }
  pathways
}

# Run fgsea
run_fgsea <- function(pw, ranks, name) {
  res <- fgsea(pathways = pw, stats = ranks,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  res <- res[order(pval), ]
  res_out <- as.data.frame(res)
  res_out$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ",")
  list(res = res, res_out = res_out)
}

# Plot functions
theme_pub <- theme_classic() +
  theme(text = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold"))

make_hallmark_barplot <- function(res_df, dataset_name) {
  h <- res_df[order(res_df$NES), ]
  h$pathway_short <- gsub("^(HALLMARK_|)", "", h$pathway)
  h$pathway_short <- factor(h$pathway_short, levels = h$pathway_short)
  h$sig <- ifelse(h$padj < 0.05, "padj < 0.05",
                  ifelse(h$padj < 0.25, "padj < 0.25", "NS"))
  p <- ggplot(h, aes(x = NES, y = pathway_short, fill = sig)) +
    geom_col() +
    scale_fill_manual(values = c("padj < 0.05" = "#D62728",
                                  "padj < 0.25" = "#FF7F0E", "NS" = "grey70")) +
    labs(x = "NES", y = NULL, title = sprintf("Hallmark — %s", dataset_name)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_pub + theme(axis.text.y = element_text(size = 6))
  ggsave(file.path(OUTDIR, sprintf("hallmark_%s.pdf", dataset_name)), p, width = 8, height = 10)
  ggsave(file.path(OUTDIR, sprintf("hallmark_%s.png", dataset_name)), p, width = 8, height = 10, dpi = 150)
}

make_top_dotplot <- function(res_df, title, filename, n = 30) {
  sig <- res_df[res_df$padj < 0.25, ]
  if (nrow(sig) < n) {
    top <- res_df[order(res_df$padj), ][1:min(n, nrow(res_df)), ]
  } else {
    pos <- sig[sig$NES > 0, ][order(sig[sig$NES > 0, ]$padj), ]
    neg <- sig[sig$NES < 0, ][order(sig[sig$NES < 0, ]$padj), ]
    n_each <- floor(n / 2)
    top <- rbind(head(pos, n_each), head(neg, n - n_each))
  }
  top <- top[order(top$NES), ]
  top$pathway_short <- gsub("^(REACTOME_|KEGG_|WP_|BIOCARTA_|PID_|GOBP_|GO_)", "", top$pathway)
  top$pathway_short <- ifelse(nchar(top$pathway_short) > 60,
                               paste0(substr(top$pathway_short, 1, 57), "..."),
                               top$pathway_short)
  top$pathway_short <- factor(top$pathway_short, levels = top$pathway_short)
  p <- ggplot(top, aes(x = NES, y = pathway_short)) +
    geom_point(aes(size = size, color = -log10(padj))) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(2, 6)) +
    labs(x = "NES", y = NULL, title = title) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_pub + theme(axis.text.y = element_text(size = 5.5))
  ggsave(file.path(OUTDIR, filename), p, width = 10, height = 8)
  ggsave(file.path(OUTDIR, gsub("\\.pdf$", ".png", filename)), p, width = 10, height = 8, dpi = 150)
}

# Download GMT files for both species
cat("=== Downloading GMT files ===\n")
gmt_mouse <- download_gmts("mouse")
gmt_human <- download_gmts("human")

# Process each dataset
for (ds_name in names(deseq_files)) {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("GSEA: %s (%s)\n", ds_name, species[[ds_name]]))
  cat(sprintf("============================================================\n"))

  fpath <- deseq_files[[ds_name]]
  if (!file.exists(fpath)) {
    cat("  SKIPPED: file not found\n")
    next
  }

  deseq <- read.delim(fpath, stringsAsFactors = FALSE)

  # Determine stat column
  if ("stat" %in% colnames(deseq)) {
    stat_col <- "stat"
  } else if ("log2FoldChange" %in% colnames(deseq) && "lfcSE" %in% colnames(deseq)) {
    deseq$stat <- deseq$log2FoldChange / deseq$lfcSE
    stat_col <- "stat"
  } else {
    cat("  SKIPPED: no stat column\n")
    next
  }

  # Gene name column
  gene_col <- if ("gene" %in% colnames(deseq)) "gene" else "Geneid"
  deseq_filt <- deseq[!is.na(deseq[[stat_col]]), ]
  ranks <- setNames(deseq_filt[[stat_col]], deseq_filt[[gene_col]])
  ranks <- sort(ranks, decreasing = TRUE)
  # Remove duplicates
  ranks <- ranks[!duplicated(names(ranks))]
  cat(sprintf("  Genes ranked: %d (range: %.1f to %.1f)\n",
              length(ranks), min(ranks), max(ranks)))

  pw <- if (species[[ds_name]] == "mouse") gmt_mouse else gmt_human

  results <- list()
  for (collection in names(pw)) {
    cat(sprintf("  Running %s...\n", collection))
    results[[collection]] <- run_fgsea(pw[[collection]], ranks, collection)

    res <- as.data.frame(results[[collection]]$res)
    n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
    cat(sprintf("    Significant (padj<0.05): %d\n", n_sig))

    # Save table
    outfile <- file.path(OUTDIR, sprintf("fgsea_%s_%s.tsv", collection, ds_name))
    write.table(results[[collection]]$res_out, outfile,
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Plots
  cat("  Generating plots...\n")
  make_hallmark_barplot(as.data.frame(results$hallmark$res), ds_name)
  make_top_dotplot(as.data.frame(results$c2cp$res),
                   sprintf("C2:CP Top 30 — %s", ds_name),
                   sprintf("c2cp_top30_%s.pdf", ds_name))
  make_top_dotplot(as.data.frame(results$gobp$res),
                   sprintf("GO:BP Top 30 — %s", ds_name),
                   sprintf("gobp_top30_%s.pdf", ds_name))
}

cat("\n=== Unified fgsea GSEA complete ===\n")
