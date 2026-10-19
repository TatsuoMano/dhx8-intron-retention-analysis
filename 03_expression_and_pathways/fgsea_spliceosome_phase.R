#!/usr/bin/env Rscript
# Phase-based fGSEA: Spliceosome cycle phase gene sets vs DHX8 KD ranking
# Tests whether early-stage spliceosome factors are preferentially upregulated
#
# Input:  deseq2_results.csv + data/spliceosome_phase_genesets.gmt
# Output: fgsea results table + enrichment plots

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
})

cat("=== Phase-based Spliceosome fGSEA ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- 1. Load DESeq2 results ---
deseq_path <- "../2026-02-27_comprehensive/data/deseq2_results.csv"
cat("Loading DESeq2 results from:", deseq_path, "\n")
deseq <- read.csv(deseq_path, row.names = 1)
deseq_filt <- deseq[!is.na(deseq$stat), ]
ranks <- setNames(deseq_filt$stat, rownames(deseq_filt))
ranks <- sort(ranks, decreasing = TRUE)
cat("  Genes with rank:", length(ranks), "\n")
cat("  Rank range:", round(min(ranks), 2), "to", round(max(ranks), 2), "\n\n")

# --- 2. Load custom GMT ---
gmt_path <- "data/spliceosome_phase_genesets.gmt"
cat("Loading custom GMT:", gmt_path, "\n")
pathways <- gmtPathways(gmt_path)
for (name in names(pathways)) {
  # Check overlap with ranking
  overlap <- sum(pathways[[name]] %in% names(ranks))
  cat(sprintf("  %s: %d genes (%d in ranking)\n", name, length(pathways[[name]]), overlap))
}
cat("\n")

# --- 3. Run fGSEA ---
cat("Running fGSEA...\n")
res <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 5,       # smaller min because some phase sets are small
  maxSize = 500,
  nPermSimple = 50000  # more permutations for stability
)
res <- res[order(pval), ]

cat("\n=== Results ===\n")
for (i in 1:nrow(res)) {
  r <- res[i, ]
  sig_mark <- ifelse(r$padj < 0.05, "***", ifelse(r$padj < 0.1, "**", ifelse(r$padj < 0.25, "*", "")))
  cat(sprintf("  %-45s NES=%6.2f  padj=%.3e  size=%3d  %s\n",
              r$pathway, r$NES, r$padj, r$size, sig_mark))
}

# --- 4. Save results ---
res_out <- as.data.frame(res)
res_out$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ",")
write.table(res_out, "fgsea_spliceosome_phase.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved: fgsea_spliceosome_phase.tsv\n")

# --- 5. Visualization ---
cat("\nCreating visualizations...\n")

theme_pub <- theme_classic() +
  theme(
    text = element_text(family = "sans", size = 9),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 8)
  )

# 5a. NES bar chart colored by phase
phase_colors <- c(
  "SPLICEOSOME_PHASE1_EA_COMPLEX" = "#4DAF4A",
  "SPLICEOSOME_PHASE2_B_COMPLEX" = "#377EB8",
  "SPLICEOSOME_PHASE3_BACT_NTC" = "#FF7F00",
  "SPLICEOSOME_PHASE4_CATALYTIC" = "#E41A1C",
  "SPLICEOSOME_PHASE5_POST_CATALYTIC" = "#984EA3",
  "SPLICEOSOME_PHASE6_REGULATORY" = "#A65628",
  "SPLICEOSOME_ALL_CORE" = "#666666"
)

res_plot <- as.data.frame(res)
res_plot$pathway_short <- gsub("^SPLICEOSOME_", "", res_plot$pathway)
res_plot$pathway_short <- gsub("_", " ", res_plot$pathway_short)
res_plot <- res_plot[order(res_plot$NES), ]
res_plot$pathway_short <- factor(res_plot$pathway_short, levels = res_plot$pathway_short)
res_plot$sig <- ifelse(res_plot$padj < 0.05, "padj < 0.05",
                       ifelse(res_plot$padj < 0.25, "padj < 0.25", "NS"))

p_nes <- ggplot(res_plot, aes(x = NES, y = pathway_short)) +
  geom_col(aes(fill = pathway), show.legend = FALSE) +
  scale_fill_manual(values = phase_colors) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = "Spliceosome Phase Gene Sets - DHX8 KD",
       subtitle = "fGSEA: early assembly factors show strongest positive enrichment") +
  theme_pub +
  theme(axis.text.y = element_text(size = 8))

# Add significance stars
for (i in 1:nrow(res_plot)) {
  if (res_plot$padj[i] < 0.05) {
    p_nes <- p_nes + annotate("text", x = res_plot$NES[i] + 0.05 * sign(res_plot$NES[i]),
                               y = i, label = "***", size = 4, hjust = ifelse(res_plot$NES[i] > 0, 0, 1))
  } else if (res_plot$padj[i] < 0.25) {
    p_nes <- p_nes + annotate("text", x = res_plot$NES[i] + 0.05 * sign(res_plot$NES[i]),
                               y = i, label = "*", size = 4, hjust = ifelse(res_plot$NES[i] > 0, 0, 1))
  }
}

ggsave("figures/splicing_factor_phase_gsea_nes.pdf", p_nes, width = 9, height = 5)
ggsave("figures/splicing_factor_phase_gsea_nes.png", p_nes, width = 9, height = 5, dpi = 300)
cat("  Saved: figures/splicing_factor_phase_gsea_nes.pdf/png\n")

# 5b. Running enrichment plots for each phase
pdf("figures/splicing_factor_phase_gsea_running.pdf", width = 7, height = 5)
for (pw_name in names(pathways)) {
  if (pw_name %in% res$pathway) {
    tryCatch({
      p <- plotEnrichment(pathways[[pw_name]], ranks) +
        labs(title = gsub("^SPLICEOSOME_", "", pw_name)) +
        theme_pub
      print(p)
    }, error = function(e) {
      cat("  Warning: could not plot", pw_name, "\n")
    })
  }
}
dev.off()
cat("  Saved: figures/splicing_factor_phase_gsea_running.pdf\n")

cat("\nDone!\n")
