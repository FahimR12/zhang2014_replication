#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R — Master runner for Zhang et al. (2014) replication
#
# Usage:
#   Rscript run_pipeline.R              # Run all steps
#   Rscript run_pipeline.R 1            # Run only step 1
#   Rscript run_pipeline.R 3 5          # Run steps 3 through 5
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Determine which steps to run
if (length(args) == 0) {
  steps <- 1:5
} else if (length(args) == 1) {
  steps <- as.integer(args[1])
} else {
  steps <- as.integer(args[1]):as.integer(args[2])
}

cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  Zhang et al. (2014) TCGA-OV Replication Pipeline      ║\n")
cat("║  BMC Systems Biology 8:1338                             ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat("Running steps:", paste(steps, collapse = ", "), "\n\n")

# Set working directory to script location
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_full[grep(file_arg, args_full)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
setwd(script_dir)

total_time <- system.time({
  if (1 %in% steps) {
    cat("\n▶ Step 1/5: Data Loading\n")
    cat(strrep("─", 60), "\n")
    source("01_data_loading.R", local = TRUE)
  }

  if (2 %in% steps) {
    cat("\n▶ Step 2/5: Preprocessing\n")
    cat(strrep("─", 60), "\n")
    source("02_preprocessing.R", local = TRUE)
  }

  if (3 %in% steps) {
    cat("\n▶ Step 3/5: Feature Selection\n")
    cat(strrep("─", 60), "\n")
    source("03_feature_selection.R", local = TRUE)
  }

  if (4 %in% steps) {
    cat("\n▶ Step 4/5: Bayesian Network Learning\n")
    cat(strrep("─", 60), "\n")
    source("04_bn_learning.R", local = TRUE)
  }

  if (5 %in% steps) {
    cat("\n▶ Step 5/5: Analysis & Visualization\n")
    cat(strrep("─", 60), "\n")
    source("05_analysis.R", local = TRUE)
  }
})

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  PIPELINE COMPLETE                                      ║\n")
cat("║  Total time:", sprintf("%-40s", paste(round(total_time["elapsed"]/60, 1), "minutes")), "║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")