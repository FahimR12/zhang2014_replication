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

status_file <- Sys.getenv("PIPELINE_STATUS_FILE", unset = file.path("logs", "pipeline_status.txt"))
dir.create(dirname(status_file), recursive = TRUE, showWarnings = FALSE)

status_log <- function(msg) {
  line <- sprintf("%s | %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = status_file, append = TRUE)
  flush.console()
}

run_step <- function(step_num, title, script_name) {
  cat(sprintf("\n▶ Step %d/5: %s\n", step_num, title))
  cat(strrep("─", 60), "\n")
  status_log(sprintf("START step %d: %s", step_num, title))
  # Use global scope so existing scripts that rely on clusterExport(..., envir=.GlobalEnv)
  # keep working exactly as before.
  t <- system.time(source(script_name, local = .GlobalEnv))
  status_log(sprintf("DONE  step %d: %s (elapsed %.1f min)",
                     step_num, title, unname(t["elapsed"]) / 60))
}

status_log(sprintf("PIPELINE START (steps=%s)", paste(steps, collapse = ",")))

total_time <- system.time({
  if (1 %in% steps) {
    run_step(1, "Data Loading", "01_data_loading.R")
  }

  if (2 %in% steps) {
    run_step(2, "Preprocessing", "02_preprocessing.R")
  }

  if (3 %in% steps) {
    run_step(3, "Feature Selection", "03_feature_selection.R")
  }

  if (4 %in% steps) {
    run_step(4, "Bayesian Network Learning", "04_bn_learning.R")
  }

  if (5 %in% steps) {
    run_step(5, "Analysis & Visualization", "05_analysis.R")
  }
})

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  PIPELINE COMPLETE                                      ║\n")
cat("║  Total time:", sprintf("%-40s", paste(round(total_time["elapsed"]/60, 1), "minutes")), "║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")
status_log(sprintf("PIPELINE COMPLETE (elapsed %.1f min)", unname(total_time["elapsed"]) / 60))
