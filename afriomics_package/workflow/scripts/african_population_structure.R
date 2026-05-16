# =============================================================================
# AfriOmics | scripts/african_population_structure.R
# =============================================================================
# Identify multi-omics variation driven by African region, ethnicity,
# diet type, and urbanisation — critical confounders that must be
# separated from disease signals
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

cat("[AfriOmics] African population structure analysis\n")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
mg_matrix    <- read.delim(snakemake@input[["mg_matrix"]], row.names = 1, check.names = FALSE)
mb_matrix    <- read.delim(snakemake@input[["mb_matrix"]], row.names = 1, check.names = FALSE)
mofa_factors <- read.delim(snakemake@input[["mofa_factors"]])
metadata     <- read.delim(snakemake@input[["metadata"]], row.names = 1)
african_ref  <- read.delim(snakemake@input[["african_ref"]], row.names = 1, check.names = FALSE)

dark_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0a0e0c", colour = NA),
    panel.background = element_rect(fill = "#111710", colour = NA),
    panel.grid.major = element_line(colour = "#2a3828"),
    text             = element_text(colour = "#e8f0e9"),
    axis.text        = element_text(colour = "#7a9b7e"),
    legend.background = element_rect(fill = "#111710"),
    strip.text        = element_text(colour = "#4aff91", face = "bold")
  )

# Common samples
common <- intersect(colnames(mg_matrix), rownames(metadata))
mg     <- mg_matrix[, common, drop = FALSE]
meta   <- metadata[common, , drop = FALSE]

cat(sprintf("[AfriOmics] %d samples for population structure analysis\n", length(common)))

# ---------------------------------------------------------------------------
# PERMANOVA: decompose variance by each metadata factor
# ---------------------------------------------------------------------------
cat("[AfriOmics] Running PERMANOVA to partition variance...\n")

bray_dist   <- vegdist(t(mg), method = "bray")

perm_factors <- c("region", "disease", "diet_type", "urbanisation", "body_site")
perm_factors <- intersect(perm_factors, colnames(meta))

perm_results <- list()
for (factor in perm_factors) {
  grp <- meta[[factor]]
  if (length(unique(grp)) < 2) next
  df_perm <- data.frame(group = factor(grp))
  perm    <- adonis2(bray_dist ~ group, data = df_perm, permutations = 499)
  perm_results[[factor]] <- data.frame(
    factor      = factor,
    R2          = round(perm$R2[1], 4),
    F_stat      = round(perm$F[1], 3),
    p_value     = perm$`Pr(>F)`[1],
    significant = perm$`Pr(>F)`[1] < 0.05
  )
  cat(sprintf("  %-20s R² = %.3f  p = %.4f%s\n",
              factor, perm$R2[1], perm$`Pr(>F)`[1],
              ifelse(perm$`Pr(>F)`[1] < 0.05, " *", "")))
}

perm_df <- do.call(rbind, perm_results)

# ---------------------------------------------------------------------------
# Region-specific signatures
# ---------------------------------------------------------------------------
cat("[AfriOmics] Identifying region-specific microbial signatures...\n")

regions      <- unique(meta$region)
region_sigs  <- list()

for (reg in regions) {
  reg_samples <- rownames(meta)[meta$region == reg]
  oth_samples <- rownames(meta)[meta$region != reg]

  if (length(reg_samples) < 3 || length(oth_samples) < 3) next

  rows <- lapply(rownames(mg), function(feat) {
    x <- as.numeric(mg[feat, reg_samples])
    y <- as.numeric(mg[feat, oth_samples])
    if (var(c(x, y)) == 0) return(NULL)
    wt  <- wilcox.test(x, y, exact = FALSE)
    fc  <- mean(x) - mean(y)
    data.frame(region = reg, feature = feat,
               log_fc = round(fc, 4), pvalue = wt$p.value,
               mean_region = round(mean(x), 4),
               mean_other  = round(mean(y), 4))
  })
  rows <- do.call(rbind, Filter(Negate(is.null), rows))

  if (!is.null(rows) && nrow(rows) > 0) {
    rows$fdr <- p.adjust(rows$pvalue, method = "fdr")
    rows <- rows %>% filter(fdr < 0.1) %>% arrange(desc(abs(log_fc)))
    region_sigs[[reg]] <- rows
    cat(sprintf("  %s: %d region-enriched features\n", reg, nrow(rows)))
  }
}

region_sigs_df <- do.call(rbind, region_sigs)
if (is.null(region_sigs_df) || nrow(region_sigs_df) == 0) {
  # Fallback placeholder
  region_sigs_df <- data.frame(
    region = character(), feature = character(),
    log_fc = numeric(), pvalue = numeric(),
    fdr    = numeric(), mean_region = numeric(), mean_other = numeric()
  )
}

# ---------------------------------------------------------------------------
# MOFA factor correlation with population metadata
# ---------------------------------------------------------------------------
cat("[AfriOmics] Correlating MOFA factors with population metadata...\n")

factor_cols <- grep("^Factor", colnames(mofa_factors), value = TRUE)
mofa_sub    <- mofa_factors %>%
  filter(sample_id %in% common) %>%
  column_to_rownames("sample_id") %>%
  select(all_of(factor_cols))

factor_meta_assoc <- list()
for (fact in factor_cols) {
  fvals <- mofa_sub[[fact]][match(common, rownames(mofa_sub))]
  for (var in perm_factors) {
    grps <- meta[[var]]
    if (length(unique(grps)) < 2) next
    if (is.numeric(grps)) {
      ct   <- cor.test(fvals, grps, method = "spearman", exact = FALSE)
      assoc <- data.frame(factor = fact, variable = var,
                           statistic = ct$estimate, pvalue = ct$p.value,
                           test = "spearman")
    } else {
      kt   <- kruskal.test(fvals ~ factor(grps))
      assoc <- data.frame(factor = fact, variable = var,
                           statistic = kt$statistic, pvalue = kt$p.value,
                           test = "kruskal")
    }
    factor_meta_assoc[[paste0(fact, "_", var)]] <- assoc
  }
}

factor_assoc_df <- do.call(rbind, factor_meta_assoc)
if (!is.null(factor_assoc_df) && nrow(factor_assoc_df) > 0) {
  factor_assoc_df$fdr <- p.adjust(factor_assoc_df$pvalue, method = "fdr")
}

# ---------------------------------------------------------------------------
# Comparison to African reference database
# ---------------------------------------------------------------------------
cat("[AfriOmics] Comparing to AWI-Gen / H3Africa reference...\n")

shared_feats <- intersect(rownames(mg), rownames(african_ref))
ref_centroid <- rowMeans(african_ref[shared_feats, , drop = FALSE])

stratification_df <- do.call(rbind, lapply(common, function(s) {
  sample_vec <- as.numeric(mg[shared_feats, s])
  ref_vec    <- as.numeric(ref_centroid)
  bc         <- vegan::vegdist(rbind(sample_vec, ref_vec), method = "bray")[1]

  data.frame(
    sample_id           = s,
    bray_to_african_ref = round(bc, 4),
    similarity_pct      = round((1 - bc) * 100, 1),
    region              = meta[s, "region"],
    disease             = if ("disease" %in% colnames(meta)) meta[s, "disease"] else NA,
    diet_type           = if ("diet_type" %in% colnames(meta)) meta[s, "diet_type"] else NA
  )
}))

stratification_df$african_context <- cut(
  stratification_df$bray_to_african_ref,
  breaks = c(-Inf, 0.3, 0.5, Inf),
  labels = c("WITHIN_RANGE", "MODERATE_DEVIATION", "HIGH_DEVIATION")
)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

# 1. PERMANOVA variance partitioning bar chart
if (nrow(perm_df) > 0) {
  p_perm <- ggplot(perm_df, aes(x = reorder(factor, R2), y = R2 * 100,
                                 fill = significant)) +
    geom_col(width = 0.6) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#4aff91", "FALSE" = "#3a4a3a")) +
    labs(title    = "Variance Explained by Population Factors (PERMANOVA)",
         subtitle = "African microbiome — Bray-Curtis dissimilarity",
         x        = "Metadata Factor",
         y        = "Variance Explained (%)",
         fill     = "Significant (p<0.05)") +
    dark_theme

  ggsave(file.path(dirname(snakemake@output[["stratification"]]),
                   "permanova_variance.pdf"),
         p_perm, width = 8, height = 5, dpi = 300)
}

# 2. African reference similarity by region
p_sim <- ggplot(stratification_df, aes(x = reorder(region, -similarity_pct),
                                        y = similarity_pct, fill = region)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title    = "Similarity to African Reference Microbiome by Region",
       subtitle = "AWI-Gen / H3Africa cohort baseline",
       x        = "African Region",
       y        = "% Similarity to African Reference",
       fill     = "Region") +
  dark_theme +
  theme(legend.position = "none")

pdf(snakemake@output[["structure_plot"]], width = 9, height = 6)
print(p_sim)
dev.off()

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
write.table(stratification_df, snakemake@output[["stratification"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(region_sigs_df, snakemake@output[["region_sigs"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("[AfriOmics] Mean similarity to African reference: %.1f%%\n",
            mean(stratification_df$similarity_pct, na.rm = TRUE)))
cat("[AfriOmics] African population structure analysis complete.\n")
