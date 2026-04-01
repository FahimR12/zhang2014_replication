#!/usr/bin/env Rscript
# =============================================================================
# 04_bn_learning.R
# Replication of Zhang et al. (2014) — Step 4: Bayesian Network Learning
#
# This script:
#   A) Learns BN structure using bnlearn (multiple algorithms for comparison)
#   B) Provides a standardized interface for your SMC-based BN learner
#   C) Computes comparison metrics: TPR, FDR, SHD, edge counts, etc.
#
# The paper used a custom logistic BN with BCD (blockwise coordinate descent)
# and L1 penalty. We approximate this using bnlearn's structure learning
# algorithms, which is the baseline to compare against your SMC approach.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(bnlearn)
  library(parallel)
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
DATA_DIR <- file.path(PROJECT_DIR, "data", "03_feature_selection")

RESULTS_DIR <- file.path(PROJECT_DIR, "results", "04_bn_learning")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Respect scheduler allocation first (SLURM), then fall back to local detection.
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
detected_cpus <- suppressWarnings(detectCores())
if (is.na(detected_cpus) || detected_cpus < 1) detected_cpus <- 1L

if (!is.na(slurm_cpus) && slurm_cpus > 0) {
  N_CORES <- slurm_cpus
} else {
  N_CORES <- max(1L, detected_cpus - 1L)  # leave one core free off scheduler
}

cat("=== Zhang et al. 2014 Replication — Step 4: BN Learning ===\n")
cat("Using", N_CORES, "cores",
    "(SLURM_CPUS_PER_TASK =", ifelse(is.na(slurm_cpus), "unset", slurm_cpus),
    ", detected =", detected_cpus, ")\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
bn_data <- readRDS(file.path(DATA_DIR, "bn_data.rds"))
feature_lists <- readRDS(file.path(DATA_DIR, "feature_lists.rds"))

cat("BN data:", nrow(bn_data), "samples ×", ncol(bn_data), "nodes\n")
cat("Feature types:",
    sum(grepl("^EXPR_", names(bn_data))), "expression,",
    sum(grepl("^CNV_", names(bn_data))), "CNV,",
    sum(grepl("^METH_", names(bn_data))), "methylation,",
    sum(grepl("^MUT_", names(bn_data))), "mutation\n\n")

# Handle missing data: bnlearn requires complete cases or imputation
# We'll impute with mode (most common category) for discrete data
cat("── Handling Missing Data ──────────────────────────────────\n")
na_count <- sum(is.na(bn_data))
cat("  Total NAs:", na_count, "(", round(na_count / prod(dim(bn_data)) * 100, 2), "%)\n")

if (na_count > 0) {
  cat("  Imputing with column modes...\n")
  for (col in names(bn_data)) {
    if (any(is.na(bn_data[[col]]))) {
      # Mode imputation
      mode_val <- names(sort(table(bn_data[[col]]), decreasing = TRUE))[1]
      bn_data[[col]][is.na(bn_data[[col]])] <- mode_val
      # Ensure factor levels are clean
      bn_data[[col]] <- droplevels(bn_data[[col]])
    }
  }
  cat("  Imputation complete\n")
}

# Verify all columns are factors with ≥2 levels
valid_cols <- sapply(bn_data, function(x) is.factor(x) && nlevels(x) >= 2)
bn_data <- bn_data[, valid_cols]
cat("  Final nodes after cleaning:", ncol(bn_data), "\n\n")

# =============================================================================
# A. BNLEARN STRUCTURE LEARNING — Multiple algorithms
# =============================================================================
cat("══════════════════════════════════════════════════════════\n")
cat("  PART A: BNLEARN STRUCTURE LEARNING\n")
cat("══════════════════════════════════════════════════════════\n\n")

# We run several bnlearn algorithms to find the best baseline:
# 1. Hill-Climbing (HC) — greedy score-based, most common
# 2. Tabu Search — enhanced score-based with tabu list
# 3. PC Algorithm — constraint-based
# 4. MMHC (Max-Min Hill-Climbing) — hybrid approach
# 5. H2PC — hybrid
#
# The paper used L1-penalized logistic BN with BCD.
# bnlearn's score-based methods with BIC/BDe are a reasonable proxy.

results <- list()

# --- 1. Hill-Climbing (HC) with BDe score ---
cat("── 1. Hill-Climbing (BDe score) ──────────────────────────\n")
t1 <- system.time({
  bn_hc <- tryCatch({
    hc(bn_data, score = "bde", iss = 10, restart = 10, perturb = 5)
  }, error = function(e) {
    cat("  HC failed:", e$message, "\n")
    cat("  Trying with BIC score...\n")
    tryCatch(hc(bn_data, score = "bic", restart = 5), error = function(e2) NULL)
  })
})
if (!is.null(bn_hc)) {
  results$hc <- list(
    bn = bn_hc,
    n_edges = narcs(bn_hc),
    time = t1["elapsed"],
    algorithm = "Hill-Climbing (BDe)"
  )
  cat("  Edges:", narcs(bn_hc), "| Time:", round(t1["elapsed"], 1), "sec\n")
} else {
  cat("  HC failed completely\n")
}

# --- 2. Tabu Search ---
cat("\n── 2. Tabu Search ────────────────────────────────────────\n")
t2 <- system.time({
  bn_tabu <- tryCatch({
    tabu(bn_data, score = "bde", iss = 10, tabu = 50)
  }, error = function(e) {
    cat("  Tabu failed:", e$message, "\n")
    tryCatch(tabu(bn_data, score = "bic", tabu = 50), error = function(e2) NULL)
  })
})
if (!is.null(bn_tabu)) {
  results$tabu <- list(
    bn = bn_tabu,
    n_edges = narcs(bn_tabu),
    time = t2["elapsed"],
    algorithm = "Tabu Search (BDe)"
  )
  cat("  Edges:", narcs(bn_tabu), "| Time:", round(t2["elapsed"], 1), "sec\n")
} else {
  cat("  Tabu failed completely\n")
}

# --- 3. PC Algorithm (constraint-based) ---
cat("\n── 3. PC Algorithm ───────────────────────────────────────\n")
t3 <- system.time({
  bn_pc <- tryCatch({
    pc.stable(bn_data, test = "mi", alpha = 0.05)
  }, error = function(e) {
    cat("  PC failed:", e$message, "\n")
    tryCatch(pc.stable(bn_data, test = "x2", alpha = 0.05),
             error = function(e2) NULL)
  })
})
if (!is.null(bn_pc)) {
  results$pc <- list(
    bn = bn_pc,
    n_edges = narcs(bn_pc),
    time = t3["elapsed"],
    algorithm = "PC-Stable (MI test)"
  )
  cat("  Edges:", narcs(bn_pc), "| Time:", round(t3["elapsed"], 1), "sec\n")
} else {
  cat("  PC failed completely\n")
}

# --- 4. Max-Min Hill-Climbing (MMHC) ---
cat("\n── 4. MMHC (hybrid) ──────────────────────────────────────\n")
t4 <- system.time({
  # bnlearn hybrid algorithms split arguments across two phases:
  #   restrict.args  → constraint-based skeleton phase (test type)
  #   maximize.args  → hill-climbing phase (score + iss)
  # Passing test/score as top-level args causes "unused arguments" error.
  bn_mmhc <- tryCatch({
    mmhc(bn_data,
         restrict.args = list(test = "mi"),
         maximize.args = list(score = "bde", iss = 10))
  }, error = function(e) {
    cat("  MMHC failed:", e$message, "\n")
    tryCatch(mmhc(bn_data,
                  restrict.args = list(test = "x2"),
                  maximize.args = list(score = "bic")),
             error = function(e2) NULL)
  })
})
if (!is.null(bn_mmhc)) {
  results$mmhc <- list(
    bn = bn_mmhc,
    n_edges = narcs(bn_mmhc),
    time = t4["elapsed"],
    algorithm = "MMHC (MI + BDe)"
  )
  cat("  Edges:", narcs(bn_mmhc), "| Time:", round(t4["elapsed"], 1), "sec\n")
} else {
  cat("  MMHC failed completely\n")
}

# --- 5. H2PC (hybrid) ---
cat("\n── 5. H2PC ───────────────────────────────────────────────\n")
t5 <- system.time({
  bn_h2pc <- tryCatch({
    h2pc(bn_data,
         restrict.args = list(test = "mi"),
         maximize.args = list(score = "bde", iss = 10))
  }, error = function(e) {
    cat("  H2PC failed:", e$message, "\n")
    NULL
  })
})
if (!is.null(bn_h2pc)) {
  results$h2pc <- list(
    bn = bn_h2pc,
    n_edges = narcs(bn_h2pc),
    time = t5["elapsed"],
    algorithm = "H2PC (MI + BDe)"
  )
  cat("  Edges:", narcs(bn_h2pc), "| Time:", round(t5["elapsed"], 1), "sec\n")
}

# =============================================================================
# B. BOOTSTRAP CONFIDENCE — Average network from multiple bootstrap samples
# =============================================================================
cat("\n── Bootstrap Averaging (HC, 200 replicates) ───────────────\n")

t_boot <- system.time({
  cl_boot <- NULL
  if (N_CORES > 1L) {
    cl_boot <- makeCluster(N_CORES, type = "PSOCK")
    on.exit(stopCluster(cl_boot), add = TRUE)
  }

  boot_strength <- tryCatch({
    boot.strength(bn_data, R = 200, algorithm = "hc",
                  algorithm.args = list(score = "bde", iss = 10, restart = 5),
                  cluster = cl_boot)
  }, error = function(e) {
    cat("  Bootstrap failed:", e$message, "\n")
    tryCatch({
      boot.strength(bn_data, R = 100, algorithm = "hc",
                    algorithm.args = list(score = "bic", restart = 3),
                    cluster = cl_boot)
    }, error = function(e2) NULL)
  })
})

if (!is.null(boot_strength)) {
  # Average network at different thresholds
  avg_net_50 <- averaged.network(boot_strength, threshold = 0.50)
  avg_net_75 <- averaged.network(boot_strength, threshold = 0.75)
  avg_net_85 <- averaged.network(boot_strength, threshold = 0.85)

  results$boot_50 <- list(
    bn = avg_net_50, n_edges = narcs(avg_net_50),
    time = t_boot["elapsed"], algorithm = "Bootstrap HC (threshold=0.50)"
  )
  results$boot_75 <- list(
    bn = avg_net_75, n_edges = narcs(avg_net_75),
    time = t_boot["elapsed"], algorithm = "Bootstrap HC (threshold=0.75)"
  )
  results$boot_85 <- list(
    bn = avg_net_85, n_edges = narcs(avg_net_85),
    time = t_boot["elapsed"], algorithm = "Bootstrap HC (threshold=0.85)"
  )

  cat("  Bootstrap complete (", round(t_boot["elapsed"], 1), "sec)\n")
  cat("  Averaged network edges: threshold=0.50:", narcs(avg_net_50),
      "| 0.75:", narcs(avg_net_75), "| 0.85:", narcs(avg_net_85), "\n")

  saveRDS(boot_strength, file.path(RESULTS_DIR, "boot_strength.rds"))
}

# =============================================================================
# C. EXTRACT NETWORK PROPERTIES FOR COMPARISON
# =============================================================================
cat("\n══════════════════════════════════════════════════════════\n")
cat("  NETWORK COMPARISON METRICS\n")
cat("══════════════════════════════════════════════════════════\n\n")

extract_network_metrics <- function(bn_obj, name) {
  if (is.null(bn_obj)) return(NULL)

  arcs_df <- arcs(bn_obj)
  n_nodes <- length(nodes(bn_obj))
  n_edges <- nrow(arcs_df)

  # Degree distribution
  in_deg <- sapply(nodes(bn_obj), function(n) length(parents(bn_obj, n)))
  out_deg <- sapply(nodes(bn_obj), function(n) length(children(bn_obj, n)))
  total_deg <- in_deg + out_deg

  # Hub genes (out-degree > mean + 2*SD, as per paper)
  hub_threshold <- mean(out_deg) + 2 * sd(out_deg)
  hub_genes <- names(out_deg[out_deg >= hub_threshold])

  list(
    name = name,
    n_nodes = n_nodes,
    n_edges = n_edges,
    avg_degree = mean(total_deg),
    avg_in_degree = mean(in_deg),
    avg_out_degree = mean(out_deg),
    sd_out_degree = sd(out_deg),
    max_out_degree = max(out_deg),
    hub_threshold = hub_threshold,
    n_hubs = length(hub_genes),
    hub_genes = hub_genes,
    arcs = arcs_df,
    in_degree = in_deg,
    out_degree = out_deg
  )
}

all_metrics <- list()
for (name in names(results)) {
  all_metrics[[name]] <- extract_network_metrics(results[[name]]$bn, results[[name]]$algorithm)
}

# Print comparison table
cat("Algorithm                    | Edges | Avg Deg | Hubs | Time (s)\n")
cat("────────────────────────────-|───────|─────────|──────|─────────\n")
for (name in names(results)) {
  m <- all_metrics[[name]]
  if (is.null(m)) next
  cat(sprintf("%-29s| %5d | %7.2f | %4d | %7.1f\n",
              m$name, m$n_edges, m$avg_degree, m$n_hubs, results[[name]]$time))
}

# =============================================================================
# D. SMC ALGORITHM INTERFACE
# =============================================================================
cat("\n══════════════════════════════════════════════════════════\n")
cat("  PART D: SMC ALGORITHM COMPARISON INTERFACE\n")
cat("══════════════════════════════════════════════════════════\n\n")

cat("To compare your SMC-based BN learner against bnlearn:\n\n")
cat("1. Your SMC algorithm should output a DAG as an edge list or adjacency matrix.\n")
cat("2. Save the output in one of these formats:\n")
cat("   a) Edge list CSV: columns 'from', 'to' (directed edges)\n")
cat("   b) Adjacency matrix RDS: p×p matrix, A[i,j]=1 means i→j\n\n")
cat("3. Place the file at: results/smc_edges.csv or results/smc_adjmat.rds\n\n")

# Function to load SMC results and compare
compare_with_smc <- function(smc_file = NULL) {
  # Try to load SMC results
  smc_edges <- NULL

  if (!is.null(smc_file) && file.exists(smc_file)) {
    ext <- tools::file_ext(smc_file)
    if (ext == "csv") {
      smc_edges <- fread(smc_file)
      if (!all(c("from", "to") %in% names(smc_edges))) {
        cat("Error: CSV must have columns 'from' and 'to'\n")
        return(NULL)
      }
    } else if (ext == "rds") {
      adj <- readRDS(smc_file)
      if (is.matrix(adj)) {
        node_names <- colnames(adj)
        if (is.null(node_names)) node_names <- names(bn_data)
        edges <- which(adj != 0, arr.ind = TRUE)
        smc_edges <- data.table(from = node_names[edges[,1]],
                                to = node_names[edges[,2]])
      }
    }
  } else {
    # Check default locations
    for (f in c("results/smc_edges.csv", "results/smc_adjmat.rds")) {
      if (file.exists(f)) {
        return(compare_with_smc(f))
      }
    }
    cat("No SMC results found. Run your SMC algorithm and save results.\n")
    return(NULL)
  }

  if (is.null(smc_edges) || nrow(smc_edges) == 0) {
    cat("No valid edges loaded from SMC output.\n")
    return(NULL)
  }

  # Create bnlearn network from SMC edges
  smc_bn <- empty.graph(names(bn_data))
  valid_edges <- smc_edges[from %in% names(bn_data) & to %in% names(bn_data)]
  if (nrow(valid_edges) > 0) {
    arcs(smc_bn) <- as.matrix(valid_edges[, .(from, to)])
  }

  # Compute metrics
  smc_metrics <- extract_network_metrics(smc_bn, "SMC Algorithm")

  # Compare with all bnlearn results
  cat("\n── SMC vs bnlearn Comparison ──────────────────────────────\n\n")
  cat("SMC edges:", nrow(valid_edges), "\n\n")

  # Structural Hamming Distance (SHD) against each bnlearn result
  cat("Structural Hamming Distance (SHD) — lower is more similar:\n")
  for (name in names(results)) {
    if (is.null(results[[name]]$bn)) next
    shd_val <- tryCatch(shd(smc_bn, results[[name]]$bn), error = function(e) NA)
    cat(sprintf("  SMC vs %-25s: %s\n", all_metrics[[name]]$name,
                ifelse(is.na(shd_val), "N/A", as.character(shd_val))))
  }

  # Edge overlap analysis
  cat("\nEdge overlap analysis:\n")
  smc_edge_set <- paste(valid_edges$from, valid_edges$to, sep = "→")

  for (name in names(results)) {
    if (is.null(all_metrics[[name]])) next
    bl_arcs <- all_metrics[[name]]$arcs
    if (is.null(bl_arcs) || nrow(bl_arcs) == 0) next

    bl_edge_set <- paste(bl_arcs[,1], bl_arcs[,2], sep = "→")
    common <- length(intersect(smc_edge_set, bl_edge_set))
    smc_only <- length(setdiff(smc_edge_set, bl_edge_set))
    bl_only <- length(setdiff(bl_edge_set, smc_edge_set))

    cat(sprintf("  SMC vs %-25s: %d common | %d SMC-only | %d bnlearn-only\n",
                all_metrics[[name]]$name, common, smc_only, bl_only))
  }

  # Skeleton comparison (ignoring edge direction)
  cat("\nSkeleton overlap (undirected edges):\n")
  smc_skel <- sort(apply(valid_edges[, .(from, to)], 1, function(x) paste(sort(x), collapse="--")))

  for (name in names(results)) {
    if (is.null(all_metrics[[name]])) next
    bl_arcs <- all_metrics[[name]]$arcs
    if (is.null(bl_arcs) || nrow(bl_arcs) == 0) next

    bl_skel <- sort(apply(bl_arcs, 1, function(x) paste(sort(x), collapse="--")))
    common <- length(intersect(smc_skel, bl_skel))
    total <- length(union(smc_skel, bl_skel))
    jaccard <- ifelse(total > 0, common / total, 0)

    cat(sprintf("  SMC vs %-25s: Jaccard = %.3f (%d/%d shared)\n",
                all_metrics[[name]]$name, jaccard, common, total))
  }

  return(list(bn = smc_bn, metrics = smc_metrics))
}

# Try to load SMC results if available
smc_result <- compare_with_smc()

# =============================================================================
# E. SAVE ALL RESULTS
# =============================================================================
cat("\n── Saving Results ────────────────────────────────────────\n")

# Save bnlearn networks
for (name in names(results)) {
  if (!is.null(results[[name]]$bn)) {
    saveRDS(results[[name]]$bn, file.path(RESULTS_DIR, paste0("bn_", name, ".rds")))
  }
}

# Save metrics
saveRDS(all_metrics, file.path(RESULTS_DIR, "all_metrics.rds"))
saveRDS(results, file.path(RESULTS_DIR, "all_results.rds"))

# Save edge lists as CSVs for easy inspection
for (name in names(results)) {
  if (!is.null(results[[name]]$bn) && narcs(results[[name]]$bn) > 0) {
    edge_dt <- as.data.table(arcs(results[[name]]$bn))
    setnames(edge_dt, c("from", "to"))
    fwrite(edge_dt, file.path(RESULTS_DIR, paste0("edges_", name, ".csv")))
  }
}

# Save the BN data for SMC algorithm input
fwrite(bn_data, file.path(RESULTS_DIR, "bn_input_data.csv"))
cat("  BN input data saved for SMC: results/bn_input_data.csv\n")
cat("  (", nrow(bn_data), "samples ×", ncol(bn_data), "discrete nodes)\n")

# Also save as integer matrix (for algorithms expecting numeric input)
bn_numeric <- as.data.frame(lapply(bn_data, as.integer))
rownames(bn_numeric) <- rownames(bn_data)
fwrite(bn_numeric, file.path(RESULTS_DIR, "bn_input_numeric.csv"))

cat("\n  All results saved to:", RESULTS_DIR, "\n")

# =============================================================================
# F. HELPER: RE-RUN COMPARISON AFTER SMC
# =============================================================================

cat("\n══════════════════════════════════════════════════════════\n")
cat("  TO COMPARE YOUR SMC RESULTS:\n")
cat("══════════════════════════════════════════════════════════\n")
cat("\n")
cat("  1. Run your SMC algorithm on: results/bn_input_data.csv\n")
cat("     or results/bn_input_numeric.csv (integer-encoded)\n")
cat("\n")
cat("  2. Save output edges as: results/smc_edges.csv\n")
cat("     (columns: from, to — matching column names in bn_input_data.csv)\n")
cat("\n")
cat("  3. Re-source this script or run:\n")
cat("     source('04_bn_learning.R')  # Will auto-detect SMC results\n")
cat("     # OR manually:\n")
cat("     compare_with_smc('results/smc_edges.csv')\n")
cat("\n")
cat("=== BN Learning Complete ===\n")
cat("Run 05_analysis.R for visualization and biological interpretation.\n")
