# =============================================================================
# AfriOmics | scripts/deseq2_analysis.R
# =============================================================================
# Differential expression of microbial genes between disease groups
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

cat("[AfriOmics] DESeq2 differential expression analysis\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
counts_raw <- read.delim(snakemake@input[["counts"]], comment.char="#",
                          row.names = 1, check.names = FALSE)
metadata   <- read.delim(snakemake@input[["metadata"]], row.names = 1)

group_col  <- snakemake@params[["group_col"]]
ref_level  <- snakemake@params[["reference"]]
padj_cut   <- snakemake@params[["padj_cutoff"]]
lfc_cut    <- snakemake@params[["lfc_cutoff"]]

# Remove Geneid length columns if present (featureCounts output)
count_cols <- setdiff(colnames(counts_raw), c("Chr","Start","End","Strand","Length"))
counts     <- counts_raw[, count_cols, drop = FALSE]
counts     <- counts[, colnames(counts) %in% rownames(metadata), drop = FALSE]

# Align metadata to count columns
metadata   <- metadata[colnames(counts), , drop = FALSE]

cat(sprintf("[AfriOmics] %d genes × %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("[AfriOmics] Group column: %s | Reference: %s\n", group_col, ref_level))

# ---------------------------------------------------------------------------
# DESeq2 object
# ---------------------------------------------------------------------------
coldata           <- data.frame(row.names = colnames(counts),
                                 group     = factor(metadata[[group_col]]))
coldata$group     <- relevel(coldata$group, ref = ref_level)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = ~ group
)

# Filter low-count genes: keep if ≥10 reads in ≥20% of samples
keep <- rowSums(counts(dds) >= 10) >= max(2, ncol(dds) * 0.2)
dds  <- dds[keep, ]
cat(sprintf("[AfriOmics] After low-count filtering: %d genes retained\n", nrow(dds)))

# ---------------------------------------------------------------------------
# Run DESeq2
# ---------------------------------------------------------------------------
cat("[AfriOmics] Running DESeq2...\n")
dds    <- DESeq(dds, fitType = "parametric", quiet = TRUE)
saveRDS(dds, snakemake@output[["rds"]])

# Collect results for all group comparisons vs reference
all_results <- list()
groups      <- levels(coldata$group)
comparisons <- groups[groups != ref_level]

for (grp in comparisons) {
  cat(sprintf("[AfriOmics] Computing: %s vs %s\n", grp, ref_level))
  res <- results(dds,
                  contrast  = c("group", grp, ref_level),
                  alpha     = padj_cut,
                  lfcThreshold = 0,
                  independentFiltering = TRUE)
  res <- lfcShrink(dds,
                    contrast = c("group", grp, ref_level),
                    res = res,
                    type = "ashr",
                    quiet = TRUE)
  res_df           <- as.data.frame(res)
  res_df$gene      <- rownames(res_df)
  res_df$comparison <- paste0(grp, "_vs_", ref_level)
  all_results[[paste0(grp, "_vs_", ref_level)]] <- res_df
}

results_merged <- do.call(rbind, all_results)
results_merged <- results_merged %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  mutate(
    significant = padj < padj_cut & abs(log2FoldChange) > lfc_cut,
    direction   = case_when(
      significant & log2FoldChange > 0 ~ "UP",
      significant & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

write.table(results_merged, snakemake@output[["results"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("[AfriOmics] Significant genes (padj<%.2f, |LFC|>%.1f): %d\n",
            padj_cut, lfc_cut, sum(results_merged$significant, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------
top_genes <- results_merged %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(20)

volcano_colours <- c("UP" = "#ff6b6b", "DOWN" = "#5bc4ff", "NS" = "#3a4a3a")

p_volcano <- ggplot(results_merged %>% filter(comparison == comparisons[1]),
                    aes(x = log2FoldChange, y = -log10(padj),
                        colour = direction, size = significant)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed",
             colour = "#f5c842", linewidth = 0.6) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
             colour = "#f5c842", linewidth = 0.6) +
  geom_text_repel(data = top_genes %>% filter(comparison == comparisons[1]),
                  aes(label = gene), size = 3, max.overlaps = 15,
                  colour = "white") +
  scale_colour_manual(values = volcano_colours) +
  scale_size_manual(values = c("TRUE" = 2, "FALSE" = 0.8)) +
  labs(
    title    = paste("Differential Gene Expression:", comparisons[1], "vs", ref_level),
    subtitle = "AfriOmics metatranscriptomics module",
    x        = "Log2 Fold Change",
    y        = "-Log10 (Adjusted p-value)",
    colour   = "Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background  = element_rect(fill = "#0a0e0c", colour = NA),
    panel.background = element_rect(fill = "#111710", colour = NA),
    panel.grid.major = element_line(colour = "#2a3828"),
    text             = element_text(colour = "#e8f0e9"),
    axis.text        = element_text(colour = "#7a9b7e"),
    legend.position  = "right"
  )

ggsave(snakemake@output[["volcano"]], p_volcano, width = 10, height = 7, dpi = 300)

# ---------------------------------------------------------------------------
# Heatmap — top 50 significant genes
# ---------------------------------------------------------------------------
top50_genes <- results_merged %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)

if (length(top50_genes) >= 5) {
  vsd       <- vst(dds, blind = FALSE)
  mat       <- assay(vsd)[top50_genes[top50_genes %in% rownames(vsd)], ]
  mat_scaled <- t(scale(t(mat)))

  ann_col <- data.frame(Group = coldata$group, row.names = colnames(mat))
  pal     <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

  pdf(snakemake@output[["heatmap"]], width = 12, height = 10)
  pheatmap(mat_scaled,
           annotation_col = ann_col,
           color          = pal,
           show_rownames  = TRUE,
           show_colnames  = TRUE,
           cluster_rows   = TRUE,
           cluster_cols   = TRUE,
           fontsize_row   = 7,
           main           = "Top 50 Differentially Expressed Genes\nAfriOmics Metatranscriptomics")
  dev.off()
}

cat("[AfriOmics] DESeq2 analysis complete.\n")
