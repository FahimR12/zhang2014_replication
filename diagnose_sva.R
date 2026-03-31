#!/usr/bin/env Rscript
# =============================================================================
# diagnose_sva.R — Pinpoints exactly why sva won't load in this environment
# Run: Rscript diagnose_sva.R
# =============================================================================

cat("=== SVA / mgcv Diagnostic ===\n\n")
cat("R version:", R.version$version.string, "\n")
cat("R home:", R.home(), "\n")
cat(".libPaths():\n")
cat(paste(" ", .libPaths(), collapse = "\n"), "\n\n")

# Check every package sva depends on
sva_deps <- c("mgcv", "nlme", "MASS", "lattice", "Matrix",
              "BiocGenerics", "genefilter", "limma", "edgeR")

cat("── Dependency load test ────────────────────────────────────\n")
for (pkg in sva_deps) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  loadable  <- tryCatch({ library(pkg, character.only = TRUE); TRUE },
                         error = function(e) e$message)
  status <- if (isTRUE(loadable)) "OK" else paste("FAIL:", loadable)
  cat(sprintf("  %-30s installed=%-5s  load=%s\n",
              pkg, installed, status))
}

cat("\n── Attempting to load sva directly ───────────────────────\n")
err <- tryCatch({
  library(sva)
  cat("  sva loaded successfully!\n")
  NULL
}, error = function(e) {
  cat("  FAILED:", e$message, "\n")
  e$message
})

if (!is.null(err)) {
  cat("\n── Recommended fix ────────────────────────────────────────\n")

  if (grepl("mgcv", err)) {
    cat("  mgcv is the blocking dependency.\n\n")
    cat("  FASTEST FIX (run in your conda PowerShell, NOT inside R):\n")
    cat("    conda install -c conda-forge r-mgcv\n\n")
    cat("  ALTERNATIVE (run inside R if conda is not available):\n")
    cat("    install.packages('mgcv', type='binary', repos='https://cran.r-project.org')\n")
    cat("    # Note: type='binary' avoids Fortran compilation on Windows\n\n")
  } else {
    cat("  Run: install.packages(c('", paste(sva_deps, collapse="','"), "'),\n")
    cat("         type='binary', repos='https://cran.r-project.org')\n")
  }
}

cat("\n── Quick mgcv install attempt ─────────────────────────────\n")
cat("  Trying install.packages('mgcv', type='binary')...\n")
result <- tryCatch({
  install.packages("mgcv", type = "binary", repos = "https://cran.r-project.org",
                   quiet = FALSE)
  library(mgcv)
  cat("  mgcv installed and loaded via binary!\n")
  TRUE
}, warning = function(w) {
  cat("  Warning:", w$message, "\n")
  FALSE
}, error = function(e) {
  cat("  Binary install also failed:", e$message, "\n")
  cat("  You must install via conda: conda install -c conda-forge r-mgcv\n")
  FALSE
})

if (result) {
  cat("\n  Now re-trying sva...\n")
  tryCatch({
    library(sva)
    cat("  sva now loads! Run install_dependencies.R again to finalise.\n")
  }, error = function(e) cat("  Still failing:", e$message, "\n"))
}
