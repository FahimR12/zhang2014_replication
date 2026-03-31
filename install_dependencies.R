#!/usr/bin/env Rscript
# =============================================================================
# install_dependencies.R — Install all required R packages
# =============================================================================

cat("Installing required packages for Zhang et al. (2014) replication...\n\n")

# CRAN packages
cran_pkgs <- c(
  "data.table",    # Fast data manipulation
  "jsonlite",      # JSON parsing for GDC API
  "httr",          # HTTP requests for GDC API
  "bnlearn",       # Bayesian network learning
  "igraph",        # Network analysis & visualization
  "ggplot2",       # Plotting
  "survival",      # Survival analysis
  "survminer",     # Survival plot helpers
  "parallel"       # Multi-core processing
)

cat("── CRAN Packages ──\n")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

# Bioconductor packages
cat("\n── Bioconductor Packages ──\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

bioc_pkgs <- c(
  "sva",           # ComBat for batch correction
  "biomaRt"        # Gene coordinate annotation
)

# Optional Illumina annotation packages (try both 27K and 450K)
bioc_optional <- c(
  "IlluminaHumanMethylation27kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

cat("\n── Optional Annotation Packages ──\n")
for (pkg in bioc_optional) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, " (this may take a while)\n")
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE),
      error = function(e) cat("  Warning: Could not install", pkg, ":", e$message, "\n")
    )
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

cat("\n── Verification ──\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) == 0) {
  cat("  All required packages installed successfully!\n")
} else {
  cat("  WARNING: Missing packages:", paste(missing, collapse = ", "), "\n")
  cat("  The pipeline may fail without these.\n")
}

cat("\nDone. You can now run: Rscript run_pipeline.R\n")
