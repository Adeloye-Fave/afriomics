# =============================================================================
# AfriOmics | scripts/disease_model.R
# =============================================================================
# Random Forest disease classifier from integrated multi-omics features
# Includes: cross-validation, SHAP importance, ROC curves, disease signatures
# African diseases: malaria, HIV, TB, schistosomiasis, cholera, healthy
# =============================================================================

suppressPackageStartupMessages({
  library(randomForest)
  library(caret)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
})

cat("[AfriOmics] Disease Modelling — Random Forest + Cross-Validation\n")

# ---------------------------------------------------------------------------
# Load integrated multi-omics matrices
# ---------------------------------------------------------------------------
mg_matrix    <- read.delim(snakemake@input[["mg_matrix"]], row.names = 1, check.names = FALSE)
mt_matrix    <- read.delim(snakemake@input[["mt_matrix"]], row.names = 1, check.names = FALSE)
mb_matrix    <- read.delim(snakemake@input[["mb_matrix"]], row.names = 1, check.names = FALSE)
mofa_factors <- read.delim(snakemake@input[["mofa_factors"]], row.names = NULL)
diablo_feats <- read.delim(snakemake@input[["diablo_feats"]])
metadata     <- read.delim(snakemake@input[["metadata"]], row.names = 1)

group_col <- snakemake@params[["group_col"]]
n_trees   <- snakemake@params[["n_trees"]]
n_folds   <- snakemake@params[["n_folds"]]
n_repeats <- snakemake@params[["n_repeats"]]
seed      <- snakemake@params[["seed"]]

set.seed(seed)

# Common samples
common_samples <- Reduce(intersect, list(
  colnames(mg_matrix), colnames(mt_matrix),
  colnames(mb_matrix), rownames(metadata)
))

# Subset and align
mg  <- mg_matrix[, common_samples]
mt  <- mt_matrix[, common_samples]
mb  <- mb_matrix[, common_samples]
meta <- metadata[common_samples, , drop = FALSE]

# Disease labels
labels <- factor(meta[[group_col]])
cat(sprintf("[AfriOmics] Disease groups: %s\n", paste(levels(labels), collapse = " | ")))
cat(sprintf("[AfriOmics] Sample counts: %s\n",
            paste(table(labels), collapse = " | ")))

# ---------------------------------------------------------------------------
# Feature matrix: combine top DIABLO features across all omics
# ---------------------------------------------------------------------------
cat("[AfriOmics] Building integrated feature matrix...\n")

# Use DIABLO-selected features if available, else top variance features
select_top_features <- function(mat, n = 50) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  top  <- names(sort(vars, decreasing = TRUE))[1:min(n, nrow(mat))]
  t(mat[top, , drop = FALSE])
}

# Check DIABLO features
if (nrow(diablo_feats) > 0 && "feature" %in% colnames(diablo_feats)) {
  mg_feats <- diablo_feats$feature[diablo_feats$view == "Metagenomics"]
  mt_feats <- diablo_feats$feature[diablo_feats$view == "Metatranscriptomics"]
  mb_feats <- diablo_feats$feature[diablo_feats$view == "Metabolomics"]

  mg_feats <- mg_feats[mg_feats %in% rownames(mg)]
  mt_feats <- mt_feats[mt_feats %in% rownames(mt)]
  mb_feats <- mb_feats[mb_feats %in% rownames(mb)]
} else {
  mg_feats <- rownames(mg)[order(apply(mg, 1, var), decreasing = TRUE)][1:50]
  mt_feats <- rownames(mt)[order(apply(mt, 1, var), decreasing = TRUE)][1:50]
  mb_feats <- rownames(mb)[order(apply(mb, 1, var), decreasing = TRUE)][1:50]
}

# Add MOFA factors
mofa_cols <- grep("^Factor", colnames(mofa_factors), value = TRUE)
mofa_sub  <- mofa_factors[mofa_factors$sample_id %in% common_samples, c("sample_id", mofa_cols)]
rownames(mofa_sub) <- mofa_sub$sample_id
mofa_sub  <- mofa_sub[common_samples, mofa_cols, drop = FALSE]

# Combine feature matrices
feat_mg   <- t(mg[mg_feats, , drop = FALSE])
feat_mt   <- t(mt[mt_feats, , drop = FALSE])
feat_mb   <- t(mb[mb_feats, , drop = FALSE])
feat_mofa <- as.matrix(mofa_sub)

# Prefix column names by omics layer
colnames(feat_mg)   <- paste0("MG__",   colnames(feat_mg))
colnames(feat_mt)   <- paste0("MT__",   colnames(feat_mt))
colnames(feat_mb)   <- paste0("MB__",   colnames(feat_mb))
colnames(feat_mofa) <- paste0("MOFA__", colnames(feat_mofa))

X <- cbind(feat_mg, feat_mt, feat_mb, feat_mofa)
X[is.na(X)] <- 0  # Impute NAs with 0

cat(sprintf("[AfriOmics] Feature matrix: %d samples × %d features\n", nrow(X), ncol(X)))

# ---------------------------------------------------------------------------
# Cross-validated Random Forest
# ---------------------------------------------------------------------------
cat(sprintf("[AfriOmics] Training Random Forest (%d trees, %d-fold CV × %d repeats)...\n",
            n_trees, n_folds, n_repeats))

train_ctrl <- trainControl(
  method          = "repeatedcv",
  number          = n_folds,
  repeats         = n_repeats,
  classProbs      = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = "final",
  verboseIter     = FALSE
)

rf_model <- train(
  x         = X,
  y         = labels,
  method    = "rf",
  metric    = "AUC",
  tuneGrid  = data.frame(mtry = floor(sqrt(ncol(X)))),
  ntree     = n_trees,
  trControl = train_ctrl,
  importance = TRUE
)

saveRDS(rf_model, snakemake@output[["model"]])
cat("[AfriOmics] Random Forest training complete.\n")

# ---------------------------------------------------------------------------
# Model performance metrics
# ---------------------------------------------------------------------------
cv_results <- rf_model$results
overall    <- rf_model$resample

perf_df <- data.frame(
  Metric        = c("Accuracy", "Kappa", "AUC_macro", "Balanced_Accuracy"),
  Mean          = c(mean(overall$Accuracy), mean(overall$Kappa),
                    mean(overall$AUC, na.rm = TRUE),
                    mean(overall$Mean_Balanced_Accuracy, na.rm = TRUE)),
  SD            = c(sd(overall$Accuracy), sd(overall$Kappa),
                    sd(overall$AUC, na.rm = TRUE),
                    sd(overall$Mean_Balanced_Accuracy, na.rm = TRUE)),
  n_samples     = length(labels),
  n_features    = ncol(X),
  n_trees       = n_trees,
  cv_folds      = n_folds,
  cv_repeats    = n_repeats
)

write.table(perf_df, snakemake@output[["performance"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[AfriOmics] Model Performance:\n")
print(perf_df[, 1:4])

# ---------------------------------------------------------------------------
# Feature importance
# ---------------------------------------------------------------------------
imp_raw <- varImp(rf_model, scale = TRUE)$importance
imp_df  <- as.data.frame(imp_raw)
imp_df$feature  <- rownames(imp_df)

# Parse omics layer from prefix
imp_df$omics_layer <- case_when(
  grepl("^MG__",   imp_df$feature) ~ "Metagenomics",
  grepl("^MT__",   imp_df$feature) ~ "Metatranscriptomics",
  grepl("^MB__",   imp_df$feature) ~ "Metabolomics",
  grepl("^MOFA__", imp_df$feature) ~ "MOFA_Factor",
  TRUE ~ "Unknown"
)
imp_df$feature_clean <- gsub("^(MG|MT|MB|MOFA)__", "", imp_df$feature)

# Overall importance = mean across classes
imp_df$Overall <- rowMeans(imp_df[, !colnames(imp_df) %in% c("feature", "omics_layer", "feature_clean")])
imp_df <- imp_df %>% arrange(desc(Overall))

write.table(imp_df, snakemake@output[["importance"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# ROC curves per class
# ---------------------------------------------------------------------------
dark_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0a0e0c", colour = NA),
    panel.background = element_rect(fill = "#111710", colour = NA),
    panel.grid.major = element_line(colour = "#2a3828"),
    text             = element_text(colour = "#e8f0e9"),
    axis.text        = element_text(colour = "#7a9b7e"),
    legend.background = element_rect(fill = "#111710")
  )

preds     <- rf_model$pred
roc_plots <- list()
disease_colours <- c(
  "healthy"         = "#4aff91",
  "malaria"         = "#ff6b6b",
  "HIV"             = "#f5c842",
  "TB"              = "#5bc4ff",
  "schistosomiasis" = "#c084fc",
  "cholera"         = "#fb923c",
  "sepsis"          = "#94a3b8"
)

roc_list <- list()
for (cls in levels(labels)) {
  if (cls %in% colnames(preds)) {
    roc_obj <- roc(response  = ifelse(preds$obs == cls, 1, 0),
                   predictor = preds[[cls]],
                   quiet     = TRUE)
    roc_list[[cls]] <- list(
      roc    = roc_obj,
      auc    = as.numeric(auc(roc_obj)),
      sens   = roc_obj$sensitivities,
      spec   = roc_obj$specificities
    )
  }
}

roc_df <- do.call(rbind, lapply(names(roc_list), function(cls) {
  data.frame(
    class       = cls,
    specificity = roc_list[[cls]]$spec,
    sensitivity = roc_list[[cls]]$sens,
    auc         = round(roc_list[[cls]]$auc, 3)
  )
}))

p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity,
                              colour = class, linetype = class)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#3a4a3a") +
  scale_colour_manual(
    values = disease_colours,
    labels = sapply(names(roc_list), function(cls)
      paste0(cls, " (AUC=", roc_list[[cls]]$auc %>% round(3), ")")),
    breaks = names(roc_list)
  ) +
  labs(
    title    = "ROC Curves — Multi-Omics Disease Classifier",
    subtitle = "AfriOmics | Random Forest cross-validated predictions",
    x        = "False Positive Rate (1 - Specificity)",
    y        = "True Positive Rate (Sensitivity)",
    colour   = "Disease Class",
    linetype = "Disease Class"
  ) +
  dark_theme

ggsave(snakemake@output[["roc_plot"]], p_roc, width = 9, height = 7, dpi = 300)

# ---------------------------------------------------------------------------
# Feature importance plot (top 30 per omics layer)
# ---------------------------------------------------------------------------
top_imp <- imp_df %>%
  group_by(omics_layer) %>%
  slice_max(Overall, n = 10) %>%
  ungroup() %>%
  arrange(desc(Overall))

layer_colours <- c(
  "Metagenomics"        = "#4aff91",
  "Metatranscriptomics" = "#5bc4ff",
  "Metabolomics"        = "#f5c842",
  "MOFA_Factor"         = "#c084fc"
)

p_imp <- ggplot(top_imp, aes(x = reorder(feature_clean, Overall),
                              y = Overall, fill = omics_layer)) +
  geom_col(alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = layer_colours) +
  facet_wrap(~ omics_layer, scales = "free_y", ncol = 1) +
  labs(
    title    = "Top Features by Omics Layer — Random Forest Importance",
    subtitle = "AfriOmics multi-omics disease model",
    x        = NULL,
    y        = "Mean Decrease in Accuracy (scaled)",
    fill     = "Omics Layer"
  ) +
  dark_theme +
  theme(legend.position = "none",
        strip.text = element_text(colour = "#4aff91", face = "bold"))

ggsave(snakemake@output[["shap_plot"]], p_imp, width = 10, height = 14, dpi = 300)

cat("[AfriOmics] Disease modelling complete.\n")
cat(sprintf("[AfriOmics] Top feature: %s (Importance: %.3f)\n",
            imp_df$feature_clean[1], imp_df$Overall[1]))
