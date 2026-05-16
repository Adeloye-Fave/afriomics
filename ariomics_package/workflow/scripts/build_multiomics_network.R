# =============================================================================
# AfriOmics | scripts/build_multiomics_network.R
# =============================================================================
# Build a unified multi-omics interaction network where:
# Nodes = taxa | active pathways | metabolites | disease phenotypes
# Edges = SparCC correlations | RNA-metabolite links | MOFA loadings
# Export = GraphML (Cytoscape) + interactive HTML (vis.js via htmlwidgets)
# =============================================================================

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(htmlwidgets)
  library(visNetwork)
})

cat("[AfriOmics] Building multi-omics interaction network\n")

# ---------------------------------------------------------------------------
# Load correlation data
# ---------------------------------------------------------------------------
sparcc_corr  <- read.delim(snakemake@input[["sparcc_corr"]])
sparcc_pval  <- read.delim(snakemake@input[["sparcc_pval"]])
rna_met_corr <- read.delim(snakemake@input[["rna_met_corr"]])
mofa_weights <- read.delim(snakemake@input[["mofa_weights"]])
metadata     <- read.delim(snakemake@input[["metadata"]], row.names = 1)

corr_thresh <- snakemake@params[["corr_threshold"]]
pval_thresh <- snakemake@params[["pval_threshold"]]
min_mod_sz  <- snakemake@params[["min_module_size"]]

# ---------------------------------------------------------------------------
# Build edge list from SparCC: taxa ↔ metabolites
# ---------------------------------------------------------------------------
cat("[AfriOmics] Processing SparCC edges (taxa ↔ metabolites)...\n")

# Melt correlation and p-value matrices
sparcc_long <- sparcc_corr %>%
  pivot_longer(-feature, names_to = "target", values_to = "correlation") %>%
  filter(abs(correlation) >= corr_thresh)

pval_long <- sparcc_pval %>%
  pivot_longer(-feature, names_to = "target", values_to = "pvalue") %>%
  filter(pvalue < pval_thresh)

# Merge
sparcc_edges <- inner_join(sparcc_long, pval_long, by = c("feature", "target")) %>%
  filter(feature != target) %>%
  mutate(
    edge_type  = "taxa_metabolite",
    weight     = abs(correlation),
    direction  = ifelse(correlation > 0, "positive", "negative")
  ) %>%
  select(from = feature, to = target, weight, correlation, direction, edge_type)

cat(sprintf("[AfriOmics] SparCC edges: %d\n", nrow(sparcc_edges)))

# ---------------------------------------------------------------------------
# Build edges from RNA-metabolite correlations
# ---------------------------------------------------------------------------
cat("[AfriOmics] Processing RNA–metabolite edges...\n")

rna_met_edges <- rna_met_corr %>%
  filter(!is.na(pvalue), pvalue < pval_thresh, abs(correlation) >= corr_thresh) %>%
  mutate(
    edge_type = "rna_metabolite",
    weight    = abs(correlation),
    direction = ifelse(correlation > 0, "positive", "negative")
  ) %>%
  select(from = pathway, to = metabolite, weight, correlation, direction, edge_type)

cat(sprintf("[AfriOmics] RNA-metabolite edges: %d\n", nrow(rna_met_edges)))

# ---------------------------------------------------------------------------
# Build edges from MOFA+ high-weight features (co-loading = shared factor)
# ---------------------------------------------------------------------------
cat("[AfriOmics] Processing MOFA+ co-loading edges...\n")

mofa_thresh <- 0.4   # Minimum absolute weight to consider a strong loading

# Find features with high loading on same factor across different views
mofa_high <- mofa_weights %>%
  filter(abs(value) >= mofa_thresh) %>%
  group_by(factor) %>%
  filter(n_distinct(view) >= 2) %>%   # Must span ≥2 omics views
  ungroup()

# Create edges between co-loaded features from different views
mofa_edges <- mofa_high %>%
  inner_join(mofa_high, by = "factor", suffix = c("_a", "_b")) %>%
  filter(view_a != view_b, feature_a != feature_b) %>%
  mutate(
    from       = feature_a,
    to         = feature_b,
    weight     = abs(value_a * value_b),   # Product of absolute weights
    correlation = value_a * value_b,
    direction  = ifelse(sign(value_a) == sign(value_b), "positive", "negative"),
    edge_type  = paste0("mofa_factor_", factor)
  ) %>%
  select(from, to, weight, correlation, direction, edge_type) %>%
  distinct(from, to, .keep_all = TRUE)

cat(sprintf("[AfriOmics] MOFA+ co-loading edges: %d\n", nrow(mofa_edges)))

# ---------------------------------------------------------------------------
# Combine all edges
# ---------------------------------------------------------------------------
all_edges <- bind_rows(sparcc_edges, rna_met_edges, mofa_edges) %>%
  distinct(from, to, edge_type, .keep_all = TRUE)

cat(sprintf("[AfriOmics] Total edges in network: %d\n", nrow(all_edges)))

# ---------------------------------------------------------------------------
# Build node table
# ---------------------------------------------------------------------------
all_nodes <- data.frame(
  id = unique(c(all_edges$from, all_edges$to)),
  stringsAsFactors = FALSE
)

# Classify each node by omics layer
classify_node <- function(node, mg_features, mt_features, mb_features) {
  if (node %in% mg_features) return("Microbial Taxon")
  if (node %in% mt_features) return("Active Pathway/Gene")
  if (node %in% mb_features) return("Metabolite")
  return("Unknown")
}

# Get feature lists (from column names of harmonised matrices or edge source)
mg_nodes <- unique(sparcc_edges$from)
mt_nodes <- unique(rna_met_edges$from)
mb_nodes <- unique(c(sparcc_edges$to, rna_met_edges$to))

all_nodes <- all_nodes %>%
  mutate(
    node_type = case_when(
      id %in% mg_nodes ~ "Microbial Taxon",
      id %in% mt_nodes ~ "Active Pathway",
      id %in% mb_nodes ~ "Metabolite",
      TRUE ~ "Multi-omics Feature"
    ),
    degree = sapply(id, function(n)
      sum(all_edges$from == n | all_edges$to == n)
    )
  )

cat(sprintf("[AfriOmics] Network nodes: %d total\n", nrow(all_nodes)))
print(table(all_nodes$node_type))

# ---------------------------------------------------------------------------
# Build igraph network
# ---------------------------------------------------------------------------
cat("[AfriOmics] Constructing igraph object...\n")

g <- graph_from_data_frame(
  d        = all_edges[, c("from", "to", "weight", "direction", "edge_type")],
  vertices = all_nodes[, c("id", "node_type", "degree")],
  directed = FALSE
)

# Remove self-loops and duplicate edges
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# Community detection (Louvain algorithm)
cat("[AfriOmics] Detecting network modules (Louvain)...\n")
set.seed(42)
communities <- cluster_louvain(g, weights = E(g)$weight)
V(g)$module <- membership(communities)

cat(sprintf("[AfriOmics] Detected %d modules\n", max(membership(communities))))

# Filter small modules
module_sizes <- table(V(g)$module)
valid_modules <- as.integer(names(module_sizes[module_sizes >= min_mod_sz]))
V(g)$module_filtered <- ifelse(V(g)$module %in% valid_modules, V(g)$module, 0)

# Module table
module_df <- data.frame(
  feature     = V(g)$name,
  node_type   = V(g)$node_type,
  module      = V(g)$module,
  degree      = degree(g),
  betweenness = betweenness(g, normalized = TRUE),
  closeness   = closeness(g, normalized = TRUE)
) %>% arrange(module, desc(degree))

write.table(module_df, snakemake@output[["modules"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# Export as GraphML (Cytoscape-compatible)
# ---------------------------------------------------------------------------
write_graph(g, snakemake@output[["network_graphml"]], format = "graphml")

# Also export edge table
edge_table <- as_data_frame(g, what = "edges")
write.table(edge_table, snakemake@output[["network_tsv"]],
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[AfriOmics] Network exported as GraphML and TSV.\n")

# ---------------------------------------------------------------------------
# Interactive HTML network (visNetwork)
# ---------------------------------------------------------------------------
cat("[AfriOmics] Building interactive HTML network visualisation...\n")

node_colours <- c(
  "Microbial Taxon"      = "#4aff91",
  "Active Pathway"       = "#5bc4ff",
  "Metabolite"           = "#f5c842",
  "Multi-omics Feature"  = "#c084fc"
)
edge_colours <- c(
  "positive" = "#4aff91",
  "negative" = "#ff6b6b"
)

# visNetwork node dataframe
vis_nodes <- data.frame(
  id    = V(g)$name,
  label = ifelse(degree(g) >= 5, V(g)$name, ""),   # Only label hub nodes
  title = paste0("<b>", V(g)$name, "</b><br>",
                 "Type: ", V(g)$node_type, "<br>",
                 "Module: ", V(g)$module, "<br>",
                 "Degree: ", degree(g)),
  color = node_colours[V(g)$node_type],
  size  = pmax(8, pmin(30, degree(g) * 2)),
  group = V(g)$node_type,
  font  = "14px Arial #e8f0e9"
)

# visNetwork edge dataframe
edge_df     <- as_data_frame(g, what = "edges")
vis_edges   <- data.frame(
  from   = edge_df$from,
  to     = edge_df$to,
  width  = pmax(1, pmin(5, edge_df$weight * 3)),
  color  = ifelse(edge_df$direction == "positive", "#4aff91", "#ff6b6b"),
  title  = paste0("Type: ", edge_df$edge_type, "<br>",
                  "Weight: ", round(edge_df$weight, 3), "<br>",
                  "Direction: ", edge_df$direction),
  smooth = TRUE
)

net_viz <- visNetwork(vis_nodes, vis_edges,
  main = list(text = "AfriOmics | Multi-Omics Interaction Network",
              style = "font-family:Sora; color:#4aff91; font-size:18px; font-weight:800;"),
  background = "#0a0e0c"
) %>%
  visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    selectedBy       = list(variable = "group", highlight = TRUE),
    nodesIdSelection = TRUE
  ) %>%
  visPhysics(
    solver           = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -50, springLength = 100),
    stabilization    = list(iterations = 200)
  ) %>%
  visLegend(
    useGroups  = TRUE,
    position   = "right",
    main       = "Node Type"
  ) %>%
  visInteraction(
    navigationButtons = TRUE,
    tooltipDelay      = 100
  )

saveWidget(net_viz, snakemake@output[["viz_html"]],
           selfcontained = TRUE,
           title = "AfriOmics Multi-Omics Network")

cat("[AfriOmics] Interactive network saved.\n")
cat("[AfriOmics] Network construction complete.\n")
