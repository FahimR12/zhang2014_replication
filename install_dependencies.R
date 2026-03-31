#!/usr/bin/env Rscript
# =============================================================================
# install_dependencies.R — Install all required R packages
#
# NOTE FOR CONDA USERS: conda R environments often omit "recommended" base-R
# packages (mgcv, nlme, MASS, lattice, etc.) that are normally pre-installed.
# This script installs them explicitly so downstream packages like sva work.
# =============================================================================

cat("Installing required packages for Zhang et al. (2014) replication...\n")
cat("R version:", R.version$version.string, "\n\n")

CRAN <- "https://cran.r-project.org"

# ── Helper: try to load (not just find) a package ────────────────────────────
can_load <- function(pkg) {
  isTRUE(tryCatch({
    library(pkg, character.only = TRUE)
    TRUE
  }, error = function(e) FALSE))
}

install_if_needed <- function(pkgs, installer = "cran") {
  for (pkg in pkgs) {
    needs_install <- !requireNamespace(pkg, quietly = TRUE) || !can_load(pkg)
    if (needs_install) {
      cat("  Installing:", pkg, "\n")
      if (installer == "cran") {
        install.packages(pkg, repos = CRAN, quiet = TRUE)
      } else {
        BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      }
      # Verify
      if (!can_load(pkg)) {
        cat("  WARNING: Could not load", pkg, "after installation.\n")
      } else {
        cat("  OK:", pkg, "\n")
      }
    } else {
      cat("  OK (already loaded):", pkg, "\n")
    }
  }
}

# =============================================================================
# STEP 1 — Base-R recommended packages that conda environments often omit
# These are normally bundled with R but are absent in many conda setups.
# sva depends on: mgcv → nlme; also uses MASS, lattice, Matrix
# =============================================================================
cat("── Recommended base-R packages (often missing in conda) ──\n")
base_recommended <- c(
  "mgcv",      # Required by sva (the ComBat batch correction package)
  "nlme",      # Required by mgcv
  "MASS",      # Required by sva & many others
  "lattice",   # Required by nlme / sva
  "Matrix",    # Required by sva
  "codetools", # Required by BiocParallel (dependency of sva)
  "nnet",      # Occasionally needed
  "survival"   # Also in our direct deps; ensure it's loadable
)
install_if_needed(base_recommended, "cran")

# =============================================================================
# STEP 2 — CRAN packages
# =============================================================================
cat("\n── CRAN Packages ──────────────────────────────────────────\n")
cran_pkgs <- c(
  "data.table",
  "jsonlite",
  "httr",
  "bnlearn",
  "igraph",
  "ggplot2",
  "survival",
  "survminer",
  "parallel"
)
install_if_needed(cran_pkgs, "cran")

# =============================================================================
# STEP 3 — Bioconductor packages
# =============================================================================
cat("\n── Bioconductor Packages ──────────────────────────────────\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = CRAN)

bioc_required <- c(
  "sva",      # ComBat batch correction
  "biomaRt"   # Gene coordinate annotation
)

bioc_optional <- c(
  "IlluminaHumanMethylation27kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

install_if_needed(bioc_required, "bioc")

cat("\n── Optional Annotation Packages ──────────────────────────\n")
for (pkg in bioc_optional) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "(may take a minute)\n")
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE),
      error = function(e) cat("  Warning: could not install", pkg, "\n")
    )
  } else {
    cat("  OK:", pkg, "\n")
  }
}

# =============================================================================
# STEP 4 — Full load-based verification
# =============================================================================
cat("\n── Load Verification ──────────────────────────────────────\n")

must_load <- c(cran_pkgs, bioc_required)
failed <- character(0)

for (pkg in must_load) {
  ok <- can_load(pkg)
  cat(sprintf("  %-45s %s\n", pkg, if (ok) "✓" else "✗  FAILED"))
  if (!ok) failed <- c(failed, pkg)
}

cat("\n")
if (length(failed) == 0) {
  cat("All required packages load successfully. Run: Rscript run_pipeline.R\n")
} else {
  cat("FAILED to load:", paste(failed, collapse = ", "), "\n\n")
  cat("Troubleshooting for conda environments:\n")
  cat("  Option A — Install missing packages manually inside R:\n")
  cat("    install.packages(c('", paste(failed, collapse = "', '"), "'))\n\n")
  cat("  Option B — Install system libraries then retry:\n")
  cat("    conda install -c conda-forge r-mgcv r-nlme r-mass\n\n")
  cat("  Option C — Use the plain system R instead of the conda R:\n")
  cat("    deactivate your conda env, then: Rscript install_dependencies.R\n")
}
