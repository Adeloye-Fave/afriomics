# =============================================================================
# AfriOmics | scripts/mofa_analysis.R
# =============================================================================
# MOFA+ Multi-Omics Factor Analysis
# Decomposes shared variation across metagenomics, metatranscriptomics,
# and metabolomics into interpretable latent factors
# =============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(corrplot)
})

cat("[AfriOmics] MOFA+ Multi-Omics Factor Analysis\n")

# ---------------------------------------------------------------------------
# Load harmonised data matrices
# ---------------------------------------------------------------------------
mg_matrix <- read.delim(snakemake@input[["mg_matrix"]], row.names = 1, check.names = FALSE)
mt_matrix <- read.delim(snakemake@input[["mt_matrix"]], row.names = 1, check.names = FALSE)
mb_matrix <- read.delim(snakemake@input[["mb_matrix"]], row.names = 1, check.names = FALSE)
metadata  <- read.delim(snakemake@input[["metadata"]], row.names = 1)

n_factors <- snakemake@params[["n_factors"]]
seed      <- snakemake@params[["seed"]]

cat(sprintf("[AfriOmics] Metagenomics  : %d features × %d samples\n", nrow(mg_matrix), ncol(mg_matrix)))
cat(sprintf("[AfriOmics] Metatranscript: %d features × %d samples\n", nrow(mt_matrix), ncol(mt_matrix)))
cat(sprintf("[AfriOmics] Metabolomics  : %d features × %d samples\n", nrow(mb_matrix), ncol(mb_matrix)))

# Common samples across all three layers
common_samples <- Reduce(intersect, list(
  colnames(mg_matrix),
  colnames(mt_matrix),
  colnames(mb_matrix),
  rownames(metadata)
))
cat(sprintf("[AfriOmics] Common samples across all omics: %d\n", length(common_samples)))

# Subset to common samples
mg <- mg_matrix[, common_samples]
mt <- mt_matrix[, common_samples]
mb <- mb_matrix[, common_samples]

# ---------------------------------------------------------------------------
# Create MOFA object
# ---------------------------------------------------------------------------
# MOFA expects a named list: view name → features × samples matrix
mofa_data <- list(
  "Metagenomics"        = as.matrix(mg),
  "Metatranscriptomics" = as.matrix(mt),
  "Metabolomics"        = as.matrix(mb)
)

cat("[AfriOmics] Creating MOFA object...\n")
mofa_obj <- create_mofa(mofa_data)

# ---------------------------------------------------------------------------
# Set MOFA options
# ---------------------------------------------------------------------------
data_opts <- get_default_data_options(mofa_obj)
data_opts$scale_views <- TRUE   # Scale views to comparable variance

model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- n_factors
model_opts$spikeslab_factors <- FALSE
model_opts$spikeslab_weights <- TRUE   # Sparse weights per view

train_opts <- get_default_training_options(mofa_obj)
train_opts$convergence_mode <- snakemake@params[["convergence_mode"]]
train_opts$seed             <- seed
train_opts$verbose          <- FALSE
train_opts$maxiter          <- 1000

mofa_obj <- prepare_mofa(mofa_obj,
  data_options  = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# ---------------------------------------------------------------------------
# Train MOFA+ model
# ---------------------------------------------------------------------------
cat("[AfriOmics] Training MOFA+ model...\n")
mofa_obj <- run_mofa(mofa_obj, outfile = snakemake@output[["model"]], use_basilisk = FALSE)
cat("[AfriOmics] MOFA+ training complete.\n")

# ---------------------------------------------------------------------------
# Add metadata to MOFA object
# ---------------------------------------------------------------------------
meta_for_mofa <- metadata[common_samples, , drop = FALSE]
meta_for_mofa$sample <- rownames(meta_for_mofa)
samples_metadata(mofa_obj) <- meta_for_mofa

# ---------------------------------------------------------------------------
# Extract results
# ---------------------------------------------------------------------------

# Variance explained per factor per view
var_expl <- get_variance_explained(mofa_obj)
var_r2   <- var_expl$r2_per_factor[[1]]
var_df   <- as.data.frame(var_r2)
var_df$Factor <- rownames(var_df)

write.table(var_df, snakemake@output[["var_expl"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

# Factor values per sample
factor_vals <- get_factors(mofa_obj, factors = "all")[[1]]
factor_df   <- as.data.frame(factor_vals)
factor_df$sample_id <- rownames(factor_df)
factor_df   <- merge(factor_df, meta_for_mofa, by.x = "sample_id", by.y = "sample", all.x = TRUE)

write.table(factor_df, snakemake@output[["factors"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

# Feature weights per factor per view
weight_list <- get_weights(mofa_obj, factors = "all", as.data.frame = TRUE)
write.table(weight_list, snakemake@output[["weights"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("[AfriOmics] Extracted %d factors, weights, and variance explained.\n", n_factors))

# ---------------------------------------------------------------------------
# GENERATE PLOTS
# ---------------------------------------------------------------------------
dir.create(snakemake@output[["plots_dir"]], showWarnings = FALSE, recursive = TRUE)
plots_dir <- snakemake@output[["plots_dir"]]

dark_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0a0e0c", colour = NA),
    panel.background = element_rect(fill = "#111710", colour = NA),
    panel.grid.major = element_line(colour = "#2a3828"),
    text             = element_text(colour = "#e8f0e9"),
    axis.text        = element_text(colour = "#7a9b7e"),
    legend.background = element_rect(fill = "#111710"),
    strip.text       = element_text(colour = "#4aff91")
  )

# Plot 1: Variance explained heatmap
p_var <- plot_variance_explained(mofa_obj, x = "view", y = "factor") +
  dark_theme +
  scale_fill_gradient(low = "#111710", high = "#4aff91") +
  labs(title = "MOFA+ Variance Explained per Factor per Omics Layer",
       subtitle = "AfriOmics Multi-Omics Integration")

ggsave(file.path(plots_dir, "variance_explained_heatmap.pdf"),
       p_var, width = 8, height = 5, dpi = 300)

# Plot 2: Factor scatter — coloured by disease
if ("disease" %in% colnames(meta_for_mofa)) {
  p_f1f2 <- plot_factors(mofa_obj, factors = c(1, 2),
                          color_by = "disease",
                          shape_by = if ("body_site" %in% colnames(meta_for_mofa)) "body_site" else NULL) +
    dark_theme +
    scale_colour_brewer(palette = "Set2") +
    labs(title = "MOFA+ Factor 1 vs Factor 2",
         subtitle = "Coloured by disease group")

  ggsave(file.path(plots_dir, "factor1_vs_factor2_disease.pdf"),
         p_f1f2, width = 7, height = 6, dpi = 300)
}

# Plot 3: Factor scatter — coloured by region
if ("region" %in% colnames(meta_for_mofa)) {
  p_f1f2_reg <- plot_factors(mofa_obj, factors = c(1, 2), color_by = "region") +
    dark_theme +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = "MOFA+ Factor 1 vs Factor 2",
         subtitle = "Coloured by African region")

  ggsave(file.path(plots_dir, "factor1_vs_factor2_region.pdf"),
         p_f1f2_reg, width = 7, height = 6, dpi = 300)
}

# Plot 4: Top weights per factor — metagenomics view
for (f in 1:min(5, n_factors)) {
  p_weights <- plot_top_weights(mofa_obj,
                                 view   = "Metagenomics",
                                 factor = f,
                                 nfeatures = 15) +
    dark_theme +
    labs(title  = paste0("Factor ", f, " — Top Microbial Taxa Weights"),
         subtitle = "Metagenomics view")

  ggsave(file.path(plots_dir, paste0("factor", f, "_metagenomics_weights.pdf")),
         p_weights, width = 7, height = 5, dpi = 300)
}

# Plot 5: Top weights — metabolomics view
for (f in 1:min(3, n_factors)) {
  p_met_weights <- plot_top_weights(mofa_obj,
                                     view   = "Metabolomics",
                                     factor = f,
                                     nfeatures = 15) +
    dark_theme +
    labs(title  = paste0("Factor ", f, " — Top Metabolite Weights"),
         subtitle = "Metabolomics view")

  ggsave(file.path(plots_dir, paste0("factor", f, "_metabolomics_weights.pdf")),
         p_met_weights, width = 7, height = 5, dpi = 300)
}

# Plot 6: Factor correlation with metadata
if ("disease" %in% colnames(meta_for_mofa)) {
  p_factor_meta <- plot_factor(mofa_obj,
                                factors = 1:min(5, n_factors),
                                color_by = "disease",
                                dodge = TRUE,
                                add_boxplot = TRUE) +
    dark_theme +
    scale_colour_brewer(palette = "Set2") +
    labs(title = "MOFA+ Factor Values by Disease Group",
         subtitle = "Factors correlated with disease phenotype are key disease signatures")

  ggsave(file.path(plots_dir, "factor_values_by_disease.pdf"),
         p_factor_meta, width = 10, height = 6, dpi = 300)
}

# Plot 7: Correlation between factors and metadata (heatmap)
factor_meta_cor <- cor(
  factor_df[, paste0("Factor", 1:n_factors), drop = FALSE],
  model.matrix(~ 0 + ., data = meta_for_mofa %>%
    select(where(is.numeric)) %>%
    select(where(~ var(., na.rm = TRUE) > 0))),
  use = "pairwise.complete.obs"
)

if (nrow(factor_meta_cor) > 0 && ncol(factor_meta_cor) > 0) {
  pdf(file.path(plots_dir, "factor_metadata_correlation.pdf"), width = 10, height = 6)
  corrplot(factor_meta_cor, method = "color", tl.cex = 0.7,
           col = colorRampPalette(c("#5bc4ff", "#111710", "#4aff91"))(200),
           bg = "#111710", tl.col = "#e8f0e9",
           title = "Factor–Metadata Correlations",
           mar = c(0, 0, 2, 0))
  dev.off()
}

cat(sprintf("[AfriOmics] MOFA+ plots saved to: %s\n", plots_dir))
cat("[AfriOmics] MOFA+ analysis complete.\n")
