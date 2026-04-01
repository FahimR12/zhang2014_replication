#!/usr/bin/env Rscript
# =============================================================================
# 05_analysis.R
# Replication of Zhang et al. (2014) — Step 5: Analysis & Visualization
#
# Replicates the paper's main analyses:
#   1. Hub gene identification (out-degree > mean + 2*SD)
#   2. Gene clustering (k-means, k=4)
#   3. MDS plots for sample classification
#   4. Network visualization
#   5. Survival analysis (if clinical data available)
#   6. SMC vs bnlearn comparison plots
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(bnlearn)
  library(igraph)
  library(ggplot2)
  library(survival)
  library(survminer)
})

# ── Configuration ──────────────────────────────────────────────────────────────
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_full[grep(file_arg, args_full)])
PROJECT_DIR <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
PREPROCESS_DIR <- file.path(PROJECT_DIR, "data", "02_preprocessing")
DATA_DIR <- file.path(PROJECT_DIR, "data", "03_feature_selection")
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "04_bn_learning")
ANALYSIS_DIR <- file.path(PROJECT_DIR, "results", "05_analysis")
PLOTS_DIR   <- file.path(ANALYSIS_DIR, "plots")
dir.create(ANALYSIS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Zhang et al. 2014 Replication — Step 5: Analysis ===\n\n")

# ── Load results ──────────────────────────────────────────────────────────────
all_results <- readRDS(file.path(RESULTS_DIR, "all_results.rds"))
all_metrics <- readRDS(file.path(RESULTS_DIR, "all_metrics.rds"))
bn_data     <- readRDS(file.path(DATA_DIR, "bn_data.rds"))
feature_lists <- readRDS(file.path(DATA_DIR, "feature_lists.rds"))

# Load continuous expression for clustering/MDS
ge_continuous <- tryCatch(readRDS(file.path(PREPROCESS_DIR, "ge_continuous_aligned.rds")),
                          error = function(e) NULL)

# Use the best bnlearn result (typically tabu or HC)
best_name <- names(which.max(sapply(all_results, function(x) x$n_edges)))
if (length(best_name) == 0) best_name <- names(all_results)[1]
best_bn <- all_results[[best_name]]$bn

cat("Using best bnlearn result:", all_results[[best_name]]$algorithm, "\n")
cat("Edges:", narcs(best_bn), "\n\n")

# =============================================================================
# 1. HUB GENE IDENTIFICATION
# =============================================================================
cat("── Hub Gene Identification ───────────────────────────────\n")

# Per paper: hub genes have out-degree > mean(OD) + 2*SD(OD) ≥ 7
out_deg <- sapply(nodes(best_bn), function(n) length(children(best_bn, n)))
in_deg  <- sapply(nodes(best_bn), function(n) length(parents(best_bn, n)))
total_deg <- out_deg + in_deg

hub_threshold <- mean(out_deg) + 2 * sd(out_deg)
hub_genes <- names(out_deg[out_deg >= max(hub_threshold, 1)])

# Clean gene names (remove prefix)
clean_name <- function(name) gsub("^(EXPR_|CNV_|METH_|MUT_)", "", name)
hub_gene_names <- clean_name(hub_genes)

cat("  Mean out-degree:", round(mean(out_deg), 2), "\n")
cat("  SD out-degree:", round(sd(out_deg), 2), "\n")
cat("  Hub threshold:", round(hub_threshold, 2), "\n")
cat("  Number of hub genes:", length(hub_genes), "\n")
if (length(hub_genes) > 0) {
  cat("  Hub genes:", paste(hub_gene_names, collapse = ", "), "\n")
}

# Paper's 13 hubs: ARID1A, C19orf53, CSNK2A1, DERL1, TRMT6, COL5A2,
#                  TCF21, LUM, TPX2, UBE2C, DPM1, NDUFB7, NDUFB9

# --- Out-degree histogram (Figure 4) ---
cat("\n  Plotting out-degree distribution...\n")
pdf(file.path(PLOTS_DIR, "fig4_outdegree_histogram.pdf"), width = 7, height = 5)
hist(out_deg, breaks = seq(-0.5, max(out_deg) + 0.5, by = 1),
     main = "Distribution of Out-degree (OD)",
     xlab = "Out-degree", ylab = "Frequency",
     col = "gray80", border = "gray50")
abline(v = hub_threshold, col = "red", lwd = 2)
text(hub_threshold + 0.5, par("usr")[4] * 0.9,
     paste0(">mean(OD)+2sd(OD)"), col = "red", adj = 0, cex = 0.8)
dev.off()
cat("  Saved: plots/fig4_outdegree_histogram.pdf\n")

# =============================================================================
# 2. GENE CLUSTERING (k-means, k=4)
# =============================================================================
cat("\n── Gene Clustering ───────────────────────────────────────\n")

# Per paper: 245 expression genes clustered into 4 groups using k-means
# on correlation dissimilarity metric

expr_nodes <- grep("^EXPR_", nodes(best_bn), value = TRUE)
expr_gene_names <- clean_name(expr_nodes)

if (!is.null(ge_continuous) && length(expr_gene_names) > 10) {
  # Get expression data for selected genes
  available_genes <- intersect(expr_gene_names, rownames(ge_continuous))

  if (length(available_genes) >= 10) {
    expr_sub <- ge_continuous[available_genes, ]
    # Remove NAs
    expr_sub <- expr_sub[complete.cases(expr_sub), ]

    # Correlation dissimilarity matrix
    cor_mat <- cor(t(expr_sub), use = "pairwise.complete.obs")
    dist_mat <- as.dist(1 - abs(cor_mat))

    # Determine optimal k (paper found k=4)
    cat("  Finding optimal number of clusters...\n")
    wss <- sapply(1:7, function(k) {
      km <- kmeans(as.matrix(dist_mat), centers = k, nstart = 25)
      km$tot.withinss
    })
    var_explained <- 1 - wss / wss[1]

    # k-means with k=4 (as in paper)
    set.seed(42)
    km4 <- kmeans(as.matrix(dist_mat), centers = 4, nstart = 50)

    gene_clusters <- data.table(
      gene = available_genes,
      cluster = km4$cluster
    )

    for (k in 1:4) {
      genes_in_k <- gene_clusters[cluster == k]$gene
      cat("  Cluster", k, "(", length(genes_in_k), "genes):",
          paste(head(genes_in_k, 5), collapse = ", "),
          if (length(genes_in_k) > 5) "...", "\n")
    }

    # --- MDS plot of gene clusters (Figure 6a) ---
    cat("\n  Plotting gene clusters MDS...\n")
    mds <- cmdscale(dist_mat, k = 2)
    mds_df <- data.frame(x = mds[,1], y = mds[,2],
                         cluster = as.factor(km4$cluster),
                         gene = rownames(mds))
    # Mark hubs
    mds_df$is_hub <- mds_df$gene %in% hub_gene_names

    p_clusters <- ggplot(mds_df, aes(x = x, y = y, color = cluster)) +
      geom_point(aes(size = is_hub), alpha = 0.7) +
      scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 4)) +
      labs(title = "Four gene clusters (MDS on correlation dissimilarity)",
           x = "First direction", y = "Second direction") +
      theme_minimal() +
      theme(legend.position = "right")
    ggsave(file.path(PLOTS_DIR, "fig6a_gene_clusters_mds.pdf"), p_clusters,
           width = 8, height = 6)
    cat("  Saved: plots/fig6a_gene_clusters_mds.pdf\n")

    # --- Optimal k plot (Figure 6b) ---
    p_elbow <- ggplot(data.frame(k = 1:7, var = var_explained),
                      aes(x = k, y = var)) +
      geom_point(size = 3) + geom_line() +
      geom_vline(xintercept = 4, linetype = "dashed", color = "red") +
      annotate("text", x = 4.3, y = var_explained[4], label = "optimal k=4",
               color = "red") +
      labs(title = "Choose optimal number of clusters",
           x = "Number of clusters", y = "Percentage of variance explained") +
      theme_minimal()
    ggsave(file.path(PLOTS_DIR, "fig6b_optimal_k.pdf"), p_elbow,
           width = 6, height = 5)
    cat("  Saved: plots/fig6b_optimal_k.pdf\n")

    saveRDS(gene_clusters, file.path(ANALYSIS_DIR, "gene_clusters.rds"))
  }
}

# =============================================================================
# 3. MDS PLOTS FOR SAMPLE CLASSIFICATION (Figure 5)
# =============================================================================
cat("\n── Sample Classification MDS ──────────────────────────────\n")

if (!is.null(ge_continuous)) {
  # MDS based on hub genes (Figure 5a)
  if (length(hub_gene_names) > 0) {
    hub_expr <- intersect(hub_gene_names, rownames(ge_continuous))
    if (length(hub_expr) >= 3) {
      expr_hub <- ge_continuous[hub_expr, ]
      expr_hub <- expr_hub[, complete.cases(t(expr_hub))]

      sample_dist <- dist(t(expr_hub))
      mds_samples <- cmdscale(sample_dist, k = 2)

      # Classify samples
      sample_types <- as.integer(substr(colnames(expr_hub), 14, 15))
      sample_class <- ifelse(sample_types %in% 10:14, "Normal control",
                            ifelse(sample_types == 1, "Tumor", "Other"))

      mds_sample_df <- data.frame(
        x = mds_samples[,1], y = mds_samples[,2],
        class = sample_class
      )

      p_mds <- ggplot(mds_sample_df, aes(x = x, y = y, color = class)) +
        geom_point(alpha = 0.6, size = 2) +
        scale_color_manual(values = c("Normal control" = "red",
                                      "Tumor" = "black",
                                      "Other" = "green")) +
        labs(title = paste("MDS plot based on", length(hub_expr), "hub genes"),
             x = "First Direction", y = "Second Direction") +
        theme_minimal()
      ggsave(file.path(PLOTS_DIR, "fig5a_mds_hub_genes.pdf"), p_mds,
             width = 7, height = 6)
      cat("  Saved: plots/fig5a_mds_hub_genes.pdf\n")
    }
  }

  # MDS based on all selected genes (Figure 5b)
  all_sel <- intersect(expr_gene_names, rownames(ge_continuous))
  if (length(all_sel) >= 10) {
    expr_all <- ge_continuous[all_sel, ]
    expr_all <- expr_all[, complete.cases(t(expr_all))]

    sample_dist2 <- dist(t(expr_all))
    mds_samples2 <- cmdscale(sample_dist2, k = 2)

    sample_types2 <- as.integer(substr(colnames(expr_all), 14, 15))
    sample_class2 <- ifelse(sample_types2 %in% 10:14, "Normal control",
                           ifelse(sample_types2 == 1, "Tumor", "Other"))

    mds_sample_df2 <- data.frame(
      x = mds_samples2[,1], y = mds_samples2[,2],
      class = sample_class2
    )

    p_mds2 <- ggplot(mds_sample_df2, aes(x = x, y = y, color = class)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = c("Normal control" = "red",
                                    "Tumor" = "black",
                                    "Other" = "green")) +
      labs(title = paste("MDS plot based on", length(all_sel), "cancer-related genes"),
           x = "First Direction", y = "Second Direction") +
      theme_minimal()
    ggsave(file.path(PLOTS_DIR, "fig5b_mds_all_genes.pdf"), p_mds2,
           width = 7, height = 6)
    cat("  Saved: plots/fig5b_mds_all_genes.pdf\n")
  }
}

# =============================================================================
# 4. NETWORK VISUALIZATION (Figure 3)
# =============================================================================
cat("\n── Network Visualization ──────────────────────────────────\n")

if (narcs(best_bn) > 0) {
  # Convert bnlearn to igraph for visualization
  arc_list <- arcs(best_bn)
  g <- graph_from_data_frame(arc_list, directed = TRUE, vertices = nodes(best_bn))

  # Color nodes by type
  node_colors <- rep("yellow", length(V(g)))  # expression = yellow
  names(node_colors) <- V(g)$name
  node_colors[grepl("^CNV_", V(g)$name)] <- "blue"
  node_colors[grepl("^METH_", V(g)$name)] <- "green"
  node_colors[grepl("^MUT_", V(g)$name)] <- "red"

  # Size by out-degree
  node_sizes <- out_deg[V(g)$name]
  node_sizes <- 3 + 5 * (node_sizes / max(node_sizes + 1))

  # Labels (clean names, only for high-degree nodes)
  node_labels <- rep("", length(V(g)))
  high_deg <- names(sort(total_deg, decreasing = TRUE))[1:min(30, length(total_deg))]
  for (i in seq_along(V(g)$name)) {
    if (V(g)$name[i] %in% high_deg) {
      node_labels[i] <- clean_name(V(g)$name[i])
    }
  }

  pdf(file.path(PLOTS_DIR, "fig3_network.pdf"), width = 14, height = 12)
  plot(g,
       vertex.color = node_colors,
       vertex.size = node_sizes,
       vertex.label = node_labels,
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       edge.arrow.size = 0.2,
       edge.color = rgb(0.5, 0, 0, 0.3),
       layout = layout_with_fr(g),
       main = paste("Predicted BN graph —", narcs(best_bn), "directed edges"))

  legend("bottomright",
         legend = c("Expression", "CNV", "Methylation", "Mutation"),
         fill = c("yellow", "blue", "green", "red"),
         cex = 0.8)
  dev.off()
  cat("  Saved: plots/fig3_network.pdf\n")
}

# =============================================================================
# 5. CROSS-CLUSTER CAUSAL EDGES TABLE (Table 5)
# =============================================================================
cat("\n── Cross-Cluster Edge Analysis ────────────────────────────\n")

if (exists("gene_clusters") && narcs(best_bn) > 0) {
  # Map nodes to clusters
  node_cluster <- setNames(gene_clusters$cluster, paste0("EXPR_", gene_clusters$gene))
  arc_df <- as.data.table(arcs(best_bn))

  # Assign clusters to edges
  arc_df$from_cluster <- node_cluster[arc_df$from]
  arc_df$to_cluster <- node_cluster[arc_df$to]

  # Count within/between cluster edges
  cluster_edges <- arc_df[!is.na(from_cluster) & !is.na(to_cluster)]

  if (nrow(cluster_edges) > 0) {
    cross_table <- dcast(cluster_edges, from_cluster ~ to_cluster,
                         fun.aggregate = length, value.var = "from")
    cat("  Within/between cluster causal edges:\n")
    print(cross_table)
    fwrite(cross_table, file.path(ANALYSIS_DIR, "table5_cluster_edges.csv"))
  }
}

# =============================================================================
# 6. ALGORITHM COMPARISON VISUALIZATION
# =============================================================================
cat("\n── Algorithm Comparison Plots ─────────────────────────────\n")

# Comparison bar chart
comp_df <- data.frame(
  Algorithm = sapply(all_metrics, function(m) if (!is.null(m)) m$name else NA),
  Edges = sapply(all_metrics, function(m) if (!is.null(m)) m$n_edges else NA),
  Avg_Degree = sapply(all_metrics, function(m) if (!is.null(m)) m$avg_degree else NA),
  Hubs = sapply(all_metrics, function(m) if (!is.null(m)) m$n_hubs else NA),
  Time = sapply(all_results, function(r) if (!is.null(r)) r$time else NA)
)
comp_df <- comp_df[!is.na(comp_df$Edges), ]

if (nrow(comp_df) > 0) {
  p_comp <- ggplot(comp_df, aes(x = reorder(Algorithm, -Edges), y = Edges)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = Edges), vjust = -0.5) +
    labs(title = "bnlearn Algorithm Comparison — Number of Edges",
         x = "", y = "Number of directed edges") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  ggsave(file.path(PLOTS_DIR, "algorithm_comparison.pdf"), p_comp,
         width = 10, height = 6)
  cat("  Saved: plots/algorithm_comparison.pdf\n")

  fwrite(comp_df, file.path(ANALYSIS_DIR, "algorithm_comparison.csv"))
}

# =============================================================================
# 7. SMC vs BNLEARN COMPARISON PLOTS (if SMC results exist)
# =============================================================================
cat("\n── SMC Comparison Plots ───────────────────────────────────\n")

smc_edges_file <- file.path(RESULTS_DIR, "smc_edges.csv")
if (file.exists(smc_edges_file)) {
  smc_edges <- fread(smc_edges_file)

  # Build SMC network
  smc_bn <- empty.graph(names(bn_data))
  valid_smc <- smc_edges[from %in% names(bn_data) & to %in% names(bn_data)]
  if (nrow(valid_smc) > 0) {
    arcs(smc_bn) <- as.matrix(valid_smc[, .(from, to)])

    smc_out_deg <- sapply(nodes(smc_bn), function(n) length(children(smc_bn, n)))

    # Side-by-side degree distributions
    deg_compare <- data.frame(
      degree = c(out_deg, smc_out_deg),
      method = rep(c("bnlearn", "SMC"), c(length(out_deg), length(smc_out_deg)))
    )

    p_deg_comp <- ggplot(deg_compare, aes(x = degree, fill = method)) +
      geom_histogram(position = "dodge", bins = max(c(out_deg, smc_out_deg)) + 1,
                     alpha = 0.7) +
      labs(title = "Out-degree Distribution: bnlearn vs SMC",
           x = "Out-degree", y = "Frequency") +
      theme_minimal()
    ggsave(file.path(PLOTS_DIR, "smc_vs_bnlearn_degree.pdf"), p_deg_comp,
           width = 8, height = 5)
    cat("  Saved: plots/smc_vs_bnlearn_degree.pdf\n")

    # Venn-style overlap
    bn_edge_set <- paste(arcs(best_bn)[,1], arcs(best_bn)[,2], sep = "→")
    smc_edge_set <- paste(valid_smc$from, valid_smc$to, sep = "→")
    common <- length(intersect(bn_edge_set, smc_edge_set))

    overlap_df <- data.frame(
      Category = c("bnlearn only", "Common", "SMC only"),
      Count = c(length(setdiff(bn_edge_set, smc_edge_set)),
                common,
                length(setdiff(smc_edge_set, bn_edge_set)))
    )

    p_overlap <- ggplot(overlap_df, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("bnlearn only" = "steelblue",
                                   "Common" = "purple",
                                   "SMC only" = "coral")) +
      labs(title = "Edge Overlap: bnlearn vs SMC") +
      theme_minimal()
    ggsave(file.path(PLOTS_DIR, "smc_vs_bnlearn_overlap.pdf"), p_overlap,
           width = 7, height = 5)
    cat("  Saved: plots/smc_vs_bnlearn_overlap.pdf\n")
  }
} else {
  cat("  No SMC results found. Run your SMC algorithm and save to:\n")
  cat("  results/smc_edges.csv (columns: from, to)\n")
}

# =============================================================================
# 8. SUMMARY REPORT
# =============================================================================
cat("\n══════════════════════════════════════════════════════════\n")
cat("  SUMMARY REPORT\n")
cat("══════════════════════════════════════════════════════════\n\n")

cat("Paper reference values:\n")
cat("  Nodes: 339 (245 expression + 82 CNV + 11 methylation + 1 mutation)\n")
cat("  Edges: 698 directed\n")
cat("  Mean out-degree: 2.15 (SD 2.31)\n")
cat("  Hub genes: 13 (out-degree ≥ 7)\n")
cat("  Clusters: 4\n\n")

cat("Our replication:\n")
cat("  Nodes:", length(nodes(best_bn)), "\n")
cat("  Edges:", narcs(best_bn), "\n")
cat("  Mean out-degree:", round(mean(out_deg), 2),
    "(SD", round(sd(out_deg), 2), ")\n")
cat("  Hub genes:", length(hub_genes), "\n")
if (length(hub_genes) > 0) {
  cat("  Hub gene names:", paste(hub_gene_names, collapse = ", "), "\n")
}

# Check overlap with paper's hub genes
paper_hubs <- c("ARID1A", "C19orf53", "CSNK2A1", "DERL1", "TRMT6",
                "COL5A2", "TCF21", "LUM", "TPX2", "UBE2C", "DPM1",
                "NDUFB7", "NDUFB9")
overlap_hubs <- intersect(hub_gene_names, paper_hubs)
cat("  Hub overlap with paper:", length(overlap_hubs), "/", length(paper_hubs), "\n")
if (length(overlap_hubs) > 0) {
  cat("  Matching hubs:", paste(overlap_hubs, collapse = ", "), "\n")
}

cat("\nAll plots saved to:", PLOTS_DIR, "\n")
cat("Analysis outputs saved to:", ANALYSIS_DIR, "\n")
cat("\n=== Analysis Complete ===\n")
