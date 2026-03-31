#!/usr/bin/env Rscript
# =============================================================================
# 03_feature_selection.R
# Replication of Zhang et al. (2014) — Step 3: Feature Selection
#
# Per the paper:
#   1. Define seed genes:
#      a) 48 tumor suppressors/oncogenes from TCGA data (Table 2)
#      b) 36 tumor suppressors/oncogenes from literature (Table 3)
#      c) Union → 68 unique seed genes
#   2. Apply Stepwise Correlation-Based Selection (SCBS):
#      - Start from seed genes (phenotype node)
#      - Iteratively select k most correlated features with current node set
#      - Filter by BH-adjusted p-value ≤ 0.05
#      - Continue until desired number of features reached
#   3. Final feature set: 177 genes + 82 CNV sites + 11 methylation sites + 1 mutation
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
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
DATA_DIR <- file.path(PROJECT_DIR, "data")
N_CORES  <- min(8L, max(1L, detectCores() - 1L))

cat("=== Zhang et al. 2014 Replication — Step 3: Feature Selection ===\n")
cat("Parallel workers:", N_CORES, "\n\n")

# ── Load preprocessed data ────────────────────────────────────────────────────
ge_discrete    <- readRDS(file.path(DATA_DIR, "ge_aligned.rds"))
meth_discrete  <- readRDS(file.path(DATA_DIR, "meth_aligned.rds"))
mut_matrix     <- readRDS(file.path(DATA_DIR, "mut_aligned.rds"))
ge_continuous  <- readRDS(file.path(DATA_DIR, "ge_continuous_aligned.rds"))
meth_continuous <- readRDS(file.path(DATA_DIR, "meth_continuous_aligned.rds"))
sample_info    <- readRDS(file.path(DATA_DIR, "sample_info.rds"))

# Try to load CNV
cnv_discrete <- tryCatch(readRDS(file.path(DATA_DIR, "cnv_aligned.rds")),
                         error = function(e) NULL)

cat("Loaded aligned matrices:\n")
cat("  Expression:", nrow(ge_discrete), "genes ×", ncol(ge_discrete), "samples\n")
cat("  Methylation:", nrow(meth_discrete), "genes ×", ncol(meth_discrete), "samples\n")
cat("  Mutation:", nrow(mut_matrix), "genes ×", ncol(mut_matrix), "samples\n")
if (!is.null(cnv_discrete)) {
  cat("  CNV:", nrow(cnv_discrete), "genes ×", ncol(cnv_discrete), "samples\n")
}

# =============================================================================
# 1. DEFINE SEED GENES (Tables 2 and 3 from the paper)
# =============================================================================
cat("\n── Defining Seed Genes ────────────────────────────────────\n")

# Table 2: 48 tumor suppressors and oncogenes from TCGA data
table2_suppressors <- c("CDKN2A", "MAP2K4", "MAGEC1", "RIMBP2", "DIRAS3", "PEG3",
                        "DAB2", "NF1", "ARID1A", "OPCML", "PLAGL1", "CASP9",
                        "WWOX", "RPS6KA2", "SPARC", "DLEC1")
table2_oncogenes <- c("THY1", "ALG3", "ATP5E", "ATP6V1C1", "C19orf53", "CSNK2A1",
                      "CTSF1", "DERL1", "HSF1", "ITPA", "MRPL34", "NCBP2",
                      "NDUFA13", "NDUFB7", "NDUFB9", "OSBPL2", "POLR2H",
                      "PIK3R1", "AKT2", "ERG", "PTK2", "RAE1", "RIOK1",
                      "SNRPB2", "SNX5", "SRXN1", "STX10", "TRMT1", "TRMT6",
                      "WDR53", "YWHAZ", "RAB25")

# Table 3: 36 tumor suppressors and oncogenes from literature
table3_suppressors <- c("RB1", "PTEN", "DAB2", "DLEC1", "TP53", "NF1",
                        "SPARC", "TMPRSS2", "CASP9", "PLAGL1", "WWOX",
                        "RPS6KA2", "BRCA1", "BRCA2", "DIRAS3", "PEG3",
                        "ARID1A", "OPCML")
table3_oncogenes <- c("MYC", "CDC25A", "PIK3CA", "NOTCH3", "EIF5A2",
                      "STAT3", "ETV6", "EGFR", "FGF1", "AKT2", "KRAS",
                      "RAB25", "AURKA", "PIK3R1", "ERG", "ATAD2",
                      "PDGFRA", "ERBB2")

# Union of all seed genes
all_seed_genes <- unique(c(table2_suppressors, table2_oncogenes,
                           table3_suppressors, table3_oncogenes))
cat("  Total unique seed genes:", length(all_seed_genes), "\n")

# Check which seed genes are in our expression data
seed_in_expr <- all_seed_genes[all_seed_genes %in% rownames(ge_discrete)]
seed_missing <- all_seed_genes[!all_seed_genes %in% rownames(ge_discrete)]
cat("  Seed genes found in expression data:", length(seed_in_expr), "\n")
if (length(seed_missing) > 0) {
  cat("  Missing seed genes:", paste(seed_missing, collapse = ", "), "\n")
}

# =============================================================================
# 2. DIFFERENTIAL EXPRESSION ANALYSIS (Wilcoxon rank-sum test)
# =============================================================================
cat("\n── Differential Expression Analysis ──────────────────────\n")

# Per the paper: use Wilcoxon rank-sum test to identify differentially
# expressed genes between cancer and control groups
# with BH correction, FDR ≤ 0.05

# Determine tumor vs normal from sample barcodes
sample_types <- as.integer(substr(colnames(ge_continuous), 14, 15))
is_tumor <- sample_types %in% 1:9
is_normal <- sample_types %in% 10:14

n_tumor <- sum(is_tumor, na.rm = TRUE)
n_normal <- sum(is_normal, na.rm = TRUE)
cat("  Tumor samples:", n_tumor, "| Normal samples:", n_normal, "\n")

if (n_normal >= 3) {
  cat("  Running Wilcoxon rank-sum tests using", N_CORES, "cores...\n")

  # ── Parallel Wilcoxon: one job per gene ───────────────────────────────────
  cl <- makeCluster(N_CORES, type = "PSOCK")
  clusterExport(cl, c("ge_continuous", "is_tumor", "is_normal"))

  wilcox_raw <- parLapply(cl, seq_len(nrow(ge_continuous)), function(i) {
    tv <- ge_continuous[i, is_tumor];  tv <- tv[!is.na(tv)]
    nv <- ge_continuous[i, is_normal]; nv <- nv[!is.na(nv)]
    if (length(tv) < 3 || length(nv) < 3) return(c(NA_real_, NA_real_))
    p  <- tryCatch(wilcox.test(tv, nv)$p.value, error = function(e) NA_real_)
    c(p, abs(mean(tv) - mean(nv)))
  })

  stopCluster(cl)

  wilcox_results <- data.table(
    gene        = rownames(ge_continuous),
    p_value     = sapply(wilcox_raw, `[`, 1),
    effect_size = sapply(wilcox_raw, `[`, 2)
  )
  rm(wilcox_raw); gc()

  wilcox_results$fdr <- p.adjust(wilcox_results$p_value, method = "BH")
  sig_genes_expr <- wilcox_results[fdr <= 0.05]$gene
  cat("  Significantly DE genes (FDR ≤ 0.05):", length(sig_genes_expr), "\n")
} else {
  cat("  Insufficient normal samples for Wilcoxon test\n")
  cat("  Using variance-based gene ranking as alternative\n")
  gene_vars <- apply(ge_continuous, 1, var, na.rm = TRUE)
  sig_genes_expr <- names(sort(gene_vars, decreasing = TRUE))[1:min(5000, length(gene_vars))]
}

# =============================================================================
# 3. DIFFERENTIAL ANALYSIS FOR OTHER DATA TYPES
# =============================================================================
cat("\n── Differential Analysis: Mutation, Methylation, CNV ─────\n")

# --- Mutation: genes mutated in at least 1 sample ---
# Per paper: Wilcoxon test on somatic mutation / methylation / CNV
# For mutation, genes with significant mutation frequency differences

mut_genes_all <- rownames(mut_matrix)
mut_freq <- rowSums(mut_matrix, na.rm = TRUE) / ncol(mut_matrix)
# Keep genes mutated in at least 2% of samples
sig_genes_mut <- names(mut_freq[mut_freq >= 0.02])
cat("  Genes with mutation freq ≥ 2%:", length(sig_genes_mut), "\n")

# --- Methylation: differentially methylated genes ---
if (n_normal >= 3 && nrow(meth_continuous) > 0) {
  meth_sample_types <- as.integer(substr(colnames(meth_continuous), 14, 15))
  is_tumor_meth <- meth_sample_types %in% 1:9
  is_normal_meth <- meth_sample_types %in% 10:14

  if (sum(is_normal_meth) >= 3) {
    cl <- makeCluster(N_CORES, type = "PSOCK")
    clusterExport(cl, c("meth_continuous", "is_tumor_meth", "is_normal_meth"))
    meth_p_raw <- parLapply(cl, seq_len(nrow(meth_continuous)), function(i) {
      tv <- meth_continuous[i, is_tumor_meth];  tv <- tv[!is.na(tv)]
      nv <- meth_continuous[i, is_normal_meth]; nv <- nv[!is.na(nv)]
      if (length(tv) < 3 || length(nv) < 3) return(NA_real_)
      tryCatch(wilcox.test(tv, nv)$p.value, error = function(e) NA_real_)
    })
    stopCluster(cl)
    wilcox_meth <- data.table(probe   = rownames(meth_continuous),
                              p_value = unlist(meth_p_raw))
    rm(meth_p_raw); gc()
    wilcox_meth$fdr <- p.adjust(wilcox_meth$p_value, method = "BH")
    sig_genes_meth <- wilcox_meth[fdr <= 0.05]$probe
    cat("  Significantly DM genes (FDR ≤ 0.05):", length(sig_genes_meth), "\n")
  } else {
    meth_vars <- apply(meth_continuous, 1, var, na.rm = TRUE)
    sig_genes_meth <- names(sort(meth_vars, decreasing = TRUE))[1:min(1000, length(meth_vars))]
  }
} else {
  meth_vars <- apply(meth_continuous, 1, var, na.rm = TRUE)
  sig_genes_meth <- names(sort(meth_vars, decreasing = TRUE))[1:min(1000, length(meth_vars))]
  cat("  Using top variable methylation genes:", length(sig_genes_meth), "\n")
}

# --- CNV: genes with significant copy number changes ---
if (!is.null(cnv_discrete)) {
  cnv_freq_gain <- rowSums(cnv_discrete == 1, na.rm = TRUE) / ncol(cnv_discrete)
  cnv_freq_loss <- rowSums(cnv_discrete == 0, na.rm = TRUE) / ncol(cnv_discrete)
  # Keep genes with gain or loss in ≥10% of samples
  sig_genes_cnv <- names(which(cnv_freq_gain >= 0.1 | cnv_freq_loss >= 0.1))
  cat("  Genes with CNV freq ≥ 10%:", length(sig_genes_cnv), "\n")
} else {
  sig_genes_cnv <- character(0)
}

# =============================================================================
# 4. CANDIDATE TUMOR SUPPRESSORS / ONCOGENES (Paper's criteria)
# =============================================================================
cat("\n── Identifying Candidate Tumor Suppressors/Oncogenes ─────\n")

# Per paper: a gene is a candidate if ALL THREE conditions hold:
# (1) Gene expression p-value significant (BH FDR ≤ 0.05)
# (2) Somatic mutation/methylation/CNV p-value significant (BH FDR ≤ 0.05)
# (3) |correlation| between expression and mutation/methylation/CNV > 0.4

# Genes meeting condition 1
cond1_genes <- sig_genes_expr

# Genes meeting condition 2 (significant in any epigenetic modality)
cond2_genes <- unique(c(sig_genes_mut, sig_genes_meth, sig_genes_cnv))

# Genes meeting both 1 and 2
cond12_genes <- intersect(cond1_genes, cond2_genes)
cat("  Genes meeting conditions 1 & 2:", length(cond12_genes), "\n")

# Check condition 3: correlation between expression and epigenetic features
cat("  Checking cross-modal correlations (condition 3)...\n")
candidate_genes <- character(0)

for (gene in cond12_genes) {
  if (!(gene %in% rownames(ge_continuous))) next

  expr_vals <- ge_continuous[gene, ]
  pass_cor <- FALSE

  # Check correlation with mutation
  if (gene %in% rownames(mut_matrix)) {
    mut_vals <- mut_matrix[gene, ]
    # Align samples
    common <- intersect(names(expr_vals[!is.na(expr_vals)]),
                        names(mut_vals[!is.na(mut_vals)]))
    if (length(common) >= 20) {
      cor_val <- tryCatch(cor(expr_vals[common], mut_vals[common], method = "pearson"),
                          error = function(e) 0)
      if (!is.na(cor_val) && abs(cor_val) > 0.4) pass_cor <- TRUE
    }
  }

  # Check correlation with methylation
  if (!pass_cor && gene %in% rownames(meth_continuous)) {
    meth_vals <- meth_continuous[gene, ]
    common <- intersect(names(expr_vals[!is.na(expr_vals)]),
                        names(meth_vals[!is.na(meth_vals)]))
    if (length(common) >= 20) {
      cor_val <- tryCatch(cor(expr_vals[common], meth_vals[common], method = "pearson"),
                          error = function(e) 0)
      if (!is.na(cor_val) && abs(cor_val) > 0.4) pass_cor <- TRUE
    }
  }

  # Check correlation with CNV
  if (!pass_cor && !is.null(cnv_discrete) && gene %in% rownames(cnv_discrete)) {
    cnv_vals <- cnv_discrete[gene, ]
    common <- intersect(names(expr_vals[!is.na(expr_vals)]),
                        names(cnv_vals[!is.na(cnv_vals)]))
    if (length(common) >= 20) {
      cor_val <- tryCatch(cor(as.numeric(expr_vals[common]),
                              as.numeric(cnv_vals[common]), method = "pearson"),
                          error = function(e) 0)
      if (!is.na(cor_val) && abs(cor_val) > 0.4) pass_cor <- TRUE
    }
  }

  if (pass_cor) candidate_genes <- c(candidate_genes, gene)
}

cat("  Candidate tumor suppressors/oncogenes (all 3 conditions):",
    length(candidate_genes), "\n")

# If stringent criteria yield few candidates, relax the correlation threshold
if (length(candidate_genes) < 20) {
  cat("  Relaxing correlation threshold to 0.2...\n")
  candidate_genes_relaxed <- character(0)

  for (gene in cond12_genes) {
    if (!(gene %in% rownames(ge_continuous))) next
    expr_vals <- ge_continuous[gene, ]
    pass_cor <- FALSE

    if (gene %in% rownames(mut_matrix)) {
      common <- intersect(names(expr_vals[!is.na(expr_vals)]),
                          colnames(mut_matrix))
      if (length(common) >= 20) {
        cor_val <- tryCatch(cor(expr_vals[common], mut_matrix[gene, common]),
                            error = function(e) 0)
        if (!is.na(cor_val) && abs(cor_val) > 0.2) pass_cor <- TRUE
      }
    }
    if (!pass_cor && gene %in% rownames(meth_continuous)) {
      common <- intersect(names(expr_vals[!is.na(expr_vals)]),
                          colnames(meth_continuous))
      if (length(common) >= 20) {
        cor_val <- tryCatch(cor(expr_vals[common], meth_continuous[gene, common]),
                            error = function(e) 0)
        if (!is.na(cor_val) && abs(cor_val) > 0.2) pass_cor <- TRUE
      }
    }
    if (pass_cor) candidate_genes_relaxed <- c(candidate_genes_relaxed, gene)
  }

  candidate_genes <- unique(c(candidate_genes, candidate_genes_relaxed))
  cat("  Candidates after relaxation:", length(candidate_genes), "\n")
}

# =============================================================================
# 5. FORM SEED GENE SET (union of identified candidates + literature genes)
# =============================================================================
cat("\n── Forming Final Seed Gene Set ────────────────────────────\n")

seed_genes_final <- unique(c(candidate_genes, seed_in_expr))
cat("  Data-derived candidates:", length(candidate_genes), "\n")
cat("  Literature seed genes (in data):", length(seed_in_expr), "\n")
cat("  Union (final seed set):", length(seed_genes_final), "\n")

# =============================================================================
# 6. STEPWISE CORRELATION-BASED SELECTION (SCBS)
# =============================================================================
cat("\n── Stepwise Correlation-Based Selection (SCBS) ───────────\n")

# Per the paper:
# Step 1: For current node Xi, calculate correlation with all other nodes.
#         Keep k most correlated for further filtering.
# Step 2: For each of the k nodes, test if p-value is significant under BH (FDR ≤ 0.05).
#         Select if significant.
# Step 3: Repeat until p nodes are selected.
# Recommended k = 4, 5, or 6.

run_scbs <- function(seed_genes, feature_matrix, k = 5, max_features = 300,
                     fdr_threshold = 0.05, n_cores = 1L) {
  # seed_genes    : character vector of starting gene names
  # feature_matrix: continuous gene × sample matrix (all candidate features)
  # k             : top-k most correlated neighbours considered per node
  # max_features  : stop when this many features selected
  # fdr_threshold : BH FDR cutoff
  # n_cores       : PSOCK workers for the correlation sweep (Windows-safe)

  cat("    SCBS: starting with", length(seed_genes), "seeds, k =", k,
      ", cores =", n_cores, "\n")

  selected  <- seed_genes[seed_genes %in% rownames(feature_matrix)]
  remaining <- setdiff(rownames(feature_matrix), selected)

  # ── One-time parallel cluster for correlation sweeps ─────────────────────
  cl <- makeCluster(n_cores, type = "PSOCK")
  clusterEvalQ(cl, {})   # warm up
  on.exit(stopCluster(cl), add = TRUE)

  iteration <- 0
  while (length(selected) < max_features && length(remaining) > 0) {
    iteration  <- iteration + 1
    candidates <- character(0)

    for (node in selected) {
      if (length(remaining) == 0) break
      node_vals <- feature_matrix[node, ]

      # ── Parallel correlation of node vs all remaining features ──────────
      clusterExport(cl, c("node_vals", "feature_matrix", "remaining"),
                    envir = environment())

      cors_raw <- parLapply(cl, remaining, function(r) {
        rv    <- feature_matrix[r, ]
        valid <- !is.na(node_vals) & !is.na(rv)
        if (sum(valid) < 10L) return(0)
        tryCatch(cor(node_vals[valid], rv[valid]), error = function(e) 0)
      })

      cors      <- setNames(unlist(cors_raw), remaining)
      top_k     <- names(sort(abs(cors), decreasing = TRUE))[1:min(k, length(cors))]

      # Quick pre-filter (p < 0.05 unadjusted)
      for (cand in top_k) {
        cv    <- feature_matrix[cand, ]
        valid <- !is.na(node_vals) & !is.na(cv)
        if (sum(valid) < 10L) next
        p <- tryCatch(cor.test(node_vals[valid], cv[valid])$p.value,
                      error = function(e) 1)
        if (!is.na(p) && p < 0.05) candidates <- c(candidates, cand)
      }
    }

    candidates <- unique(candidates)
    if (length(candidates) == 0) break

    # ── BH correction across all candidates this iteration ────────────────
    # Parallel: compute min-p over all selected nodes for each candidate
    clusterExport(cl, c("candidates", "selected", "feature_matrix"),
                  envir = environment())

    minp_raw <- parLapply(cl, candidates, function(cand) {
      cv   <- feature_matrix[cand, ]
      minp <- 1
      for (nd in selected) {
        nv    <- feature_matrix[nd, ]
        valid <- !is.na(nv) & !is.na(cv)
        if (sum(valid) < 10L) next
        p <- tryCatch(cor.test(nv[valid], cv[valid])$p.value,
                      error = function(e) 1)
        if (!is.na(p)) minp <- min(minp, p)
      }
      minp
    })

    fdr     <- p.adjust(unlist(minp_raw), method = "BH")
    sig_new <- candidates[fdr <= fdr_threshold]
    if (length(sig_new) == 0) break

    selected  <- c(selected, sig_new)
    remaining <- setdiff(remaining, sig_new)

      if (iteration %% 5 == 0 || iteration <= 3) {
        cat("    Iteration", iteration, ": selected", length(selected), "features\n")
      }
    } else {
      break  # No new candidates found
    }

    if (length(selected) >= max_features) break
  }

  cat("    SCBS complete:", length(selected), "features selected in",
      iteration, "iterations\n")
  return(selected)
}

# Build combined feature matrix for SCBS
# Include genes that have expression data AND at least one epigenetic feature
cat("  Building combined feature pool...\n")

# Identify genes with data in expression + at least one other modality
genes_with_epi <- unique(c(
  intersect(rownames(ge_continuous), rownames(meth_continuous)),
  intersect(rownames(ge_continuous), rownames(mut_matrix)),
  if (!is.null(cnv_discrete)) intersect(rownames(ge_continuous), rownames(cnv_discrete))
))

# Per paper: "12,000 genes that have records of expression level and at least
# one of the three (epi)genetic factors"
cat("  Genes with expression + ≥1 epigenetic feature:", length(genes_with_epi), "\n")

# Use continuous expression for SCBS correlation analysis
scbs_feature_matrix <- ge_continuous[genes_with_epi[genes_with_epi %in% rownames(ge_continuous)], ]

# Remove genes with too many NAs
na_frac <- rowMeans(is.na(scbs_feature_matrix))
scbs_feature_matrix <- scbs_feature_matrix[na_frac < 0.3, ]
cat("  Features available for SCBS:", nrow(scbs_feature_matrix), "\n")

# Run SCBS  (parallelised correlation sweeps)
# Paper used k = 5 and selected ~271 additional features beyond the 68 seeds
selected_genes <- run_scbs(
  seed_genes     = seed_genes_final,
  feature_matrix = scbs_feature_matrix,
  k              = 5,
  max_features   = 350,   # Paper got 339 total nodes
  fdr_threshold  = 0.05,
  n_cores        = N_CORES
)

cat("\n  SCBS selected", length(selected_genes), "gene expression features\n")

# =============================================================================
# 7. SELECT CORRESPONDING EPIGENETIC FEATURES
# =============================================================================
cat("\n── Selecting Epigenetic Features ──────────────────────────\n")

# For the selected genes, include their CNV, methylation, and mutation features

# CNV features
if (!is.null(cnv_discrete)) {
  cnv_features <- intersect(selected_genes, rownames(cnv_discrete))
  # Filter: keep CNV features with sufficient variation
  if (length(cnv_features) > 0) {
    cnv_var <- apply(cnv_discrete[cnv_features, , drop = FALSE], 1, function(x) {
      valid <- x[!is.na(x)]
      if (length(valid) < 10) return(0)
      var(valid)
    })
    cnv_features <- cnv_features[cnv_var > 0.01]
  }
  cat("  CNV features:", length(cnv_features), "(paper had 82)\n")
} else {
  cnv_features <- character(0)
  cat("  CNV features: 0 (CNV data not mapped to genes)\n")
}

# Methylation features
meth_features <- intersect(selected_genes, rownames(meth_discrete))
if (length(meth_features) > 0) {
  meth_var <- apply(meth_discrete[meth_features, , drop = FALSE], 1, function(x) {
    valid <- x[!is.na(x)]
    if (length(valid) < 10) return(0)
    var(valid)
  })
  meth_features <- meth_features[meth_var > 0.01]
}
cat("  Methylation features:", length(meth_features), "(paper had 11)\n")

# Mutation features — per paper: only TP53 somatic mutation was included
# (the one gene with sufficient mutation frequency in the selected set)
mut_features <- intersect(selected_genes, rownames(mut_matrix))
if (length(mut_features) > 0) {
  mut_freq_sel <- rowSums(mut_matrix[mut_features, , drop = FALSE]) / ncol(mut_matrix)
  # Keep mutations present in ≥5% of samples
  mut_features <- mut_features[mut_freq_sel >= 0.05]
}
cat("  Mutation features:", length(mut_features), "(paper had 1: TP53)\n")

# Expression features (the selected genes themselves)
expr_features <- intersect(selected_genes, rownames(ge_discrete))
cat("  Expression features:", length(expr_features), "(paper had 177→245)\n")

# =============================================================================
# 8. BUILD FINAL FEATURE MATRIX FOR BAYESIAN NETWORK
# =============================================================================
cat("\n── Building Final Feature Matrix for BN ──────────────────\n")

# Combine all features into a single sample × node data frame
# Using discretized values as per the paper

# Align all features to same sample set
all_barcodes <- colnames(ge_discrete)

# Expression nodes (3-level: low=1, medium=2, high=3)
bn_data <- as.data.frame(t(ge_discrete[expr_features, ]))
colnames(bn_data) <- paste0("EXPR_", expr_features)

# CNV nodes (binary: gain=1, loss=0)
if (length(cnv_features) > 0) {
  cnv_sub <- t(cnv_discrete[cnv_features, , drop = FALSE])
  # Align samples
  common_samples <- intersect(rownames(bn_data), rownames(cnv_sub))
  if (length(common_samples) > 0) {
    cnv_df <- as.data.frame(cnv_sub[common_samples, , drop = FALSE])
    colnames(cnv_df) <- paste0("CNV_", cnv_features)
    bn_data <- cbind(bn_data[common_samples, , drop = FALSE], cnv_df)
  }
}

# Methylation nodes (binary: hypo=1, hyper=2)
if (length(meth_features) > 0) {
  meth_sub <- t(meth_discrete[meth_features, , drop = FALSE])
  common_samples <- intersect(rownames(bn_data), rownames(meth_sub))
  if (length(common_samples) > 0) {
    meth_df <- as.data.frame(meth_sub[common_samples, , drop = FALSE])
    colnames(meth_df) <- paste0("METH_", meth_features)
    bn_data <- cbind(bn_data[common_samples, , drop = FALSE], meth_df)
  }
}

# Mutation nodes (binary: 0/1)
if (length(mut_features) > 0) {
  mut_sub <- t(mut_matrix[mut_features, , drop = FALSE])
  common_samples <- intersect(rownames(bn_data), rownames(mut_sub))
  if (length(common_samples) > 0) {
    mut_df <- as.data.frame(mut_sub[common_samples, , drop = FALSE])
    colnames(mut_df) <- paste0("MUT_", mut_features)
    bn_data <- cbind(bn_data[common_samples, , drop = FALSE], mut_df)
  }
}

# Convert all columns to factors (required for discrete BN in bnlearn)
for (col in names(bn_data)) {
  bn_data[[col]] <- as.factor(bn_data[[col]])
}

# Remove columns with only one level (no variation)
single_level <- sapply(bn_data, function(x) length(levels(x)) < 2)
bn_data <- bn_data[, !single_level]

# Remove rows with too many NAs
row_na_frac <- rowMeans(is.na(bn_data))
bn_data <- bn_data[row_na_frac < 0.3, ]

cat("  Final BN data: ", nrow(bn_data), "samples ×", ncol(bn_data), "nodes\n")
cat("  (Paper had ~580 samples × 339 nodes)\n")

# Feature type breakdown
n_expr <- sum(grepl("^EXPR_", names(bn_data)))
n_cnv  <- sum(grepl("^CNV_", names(bn_data)))
n_meth <- sum(grepl("^METH_", names(bn_data)))
n_mut  <- sum(grepl("^MUT_", names(bn_data)))
cat("  Breakdown: ", n_expr, "expression |", n_cnv, "CNV |",
    n_meth, "methylation |", n_mut, "mutation\n")

# =============================================================================
# 9. SAVE FEATURE SELECTION RESULTS
# =============================================================================
cat("\n── Saving Feature Selection Results ──────────────────────\n")

saveRDS(bn_data, file.path(DATA_DIR, "bn_data.rds"))
saveRDS(selected_genes, file.path(DATA_DIR, "selected_genes.rds"))
saveRDS(seed_genes_final, file.path(DATA_DIR, "seed_genes.rds"))

# Save feature lists
feature_lists <- list(
  seed_genes = seed_genes_final,
  candidate_genes = candidate_genes,
  scbs_selected = selected_genes,
  expr_features = expr_features,
  cnv_features = cnv_features,
  meth_features = meth_features,
  mut_features = mut_features
)
saveRDS(feature_lists, file.path(DATA_DIR, "feature_lists.rds"))

cat("  Saved to:", DATA_DIR, "\n")
cat("\n=== Feature Selection Complete ===\n")
cat("  Total nodes for BN:", ncol(bn_data), "\n")
cat("  Total samples:", nrow(bn_data), "\n")
cat("Run 04_bn_learning.R next.\n")
