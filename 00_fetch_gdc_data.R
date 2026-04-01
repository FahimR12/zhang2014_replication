#!/usr/bin/env Rscript
# =============================================================================
# 00_fetch_gdc_data.R
# Download TCGA-OV raw files from GDC into the folder layout expected by
# 01_data_loading.R.
#
# Usage:
#   Rscript 00_fetch_gdc_data.R
#   TCGA_OV_DATA_ROOT=/path/to/TCGA-OV Rscript 00_fetch_gdc_data.R
#   GDC_FETCH_LIMIT=20 Rscript 00_fetch_gdc_data.R         # test mode
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(httr)
})

args_full <- commandArgs(trailingOnly = FALSE)
file_arg  <- "--file="
script_path <- sub(file_arg, "", args_full[grep(file_arg, args_full)])
PROJECT_DIR <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

DATA_ROOT_DEFAULT <- file.path(PROJECT_DIR, "data", "GDCdata", "TCGA-OV")
DATA_ROOT <- Sys.getenv("TCGA_OV_DATA_ROOT", unset = DATA_ROOT_DEFAULT)
FETCH_LIMIT_RAW <- Sys.getenv("GDC_FETCH_LIMIT", unset = "")
FETCH_LIMIT <- if (nzchar(FETCH_LIMIT_RAW)) as.integer(FETCH_LIMIT_RAW) else NA_integer_

if (!is.na(FETCH_LIMIT) && FETCH_LIMIT <= 0) {
  stop("GDC_FETCH_LIMIT must be a positive integer if set.")
}

PROFILES <- list(
  list(
    label = "Gene Expression",
    out_subdir = file.path("Transcriptome_Profiling", "Gene_Expression_Quantification"),
    filters = list(
      data_category = "Transcriptome Profiling",
      data_type = "Gene Expression Quantification"
    )
  ),
  list(
    label = "Copy Number Variation",
    out_subdir = file.path("Copy_Number_Variation", "Copy_Number_Segment"),
    filters = list(
      data_category = "Copy Number Variation",
      data_type = "Copy Number Segment"
    )
  ),
  list(
    label = "DNA Methylation",
    out_subdir = file.path("DNA_Methylation", "Methylation_Beta_Value"),
    filters = list(
      data_category = "DNA Methylation",
      data_type = "Methylation Beta Value"
    )
  ),
  list(
    label = "Somatic Mutation",
    out_subdir = file.path("Simple_Nucleotide_Variation", "Masked_Somatic_Mutation"),
    filters = list(
      data_category = "Simple Nucleotide Variation",
      data_type = "Masked Somatic Mutation"
    )
  )
)

gdc_files_endpoint <- "https://api.gdc.cancer.gov/files"
gdc_data_endpoint  <- "https://api.gdc.cancer.gov/data"
dir.create(DATA_ROOT, recursive = TRUE, showWarnings = FALSE)

cat("=== TCGA-OV GDC Fetch ===\n")
cat("Project dir:", PROJECT_DIR, "\n")
cat("Data root:", DATA_ROOT, "\n")
if (!is.na(FETCH_LIMIT)) cat("Per-profile fetch limit:", FETCH_LIMIT, "\n")
cat("\n")

build_filter <- function(extra_filters) {
  content <- list(
    list(op = "=", content = list(field = "cases.project.project_id", value = "TCGA-OV")),
    list(op = "=", content = list(field = "access", value = "open"))
  )
  for (nm in names(extra_filters)) {
    content[[length(content) + 1L]] <- list(
      op = "=",
      content = list(field = nm, value = extra_filters[[nm]])
    )
  }
  list(op = "and", content = content)
}

query_files <- function(extra_filters, page_size = 200L) {
  all_rows <- list()
  from <- 1L

  repeat {
    body <- list(
      filters = build_filter(extra_filters),
      format = "JSON",
      fields = paste(
        c("file_id", "file_name", "data_category", "data_type", "cases.submitter_id"),
        collapse = ","
      ),
      size = as.character(page_size),
      from = as.character(from)
    )

    resp <- POST(
      gdc_files_endpoint,
      body = toJSON(body, auto_unbox = TRUE),
      content_type_json(),
      timeout(120)
    )
    stop_for_status(resp)

    parsed <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyDataFrame = FALSE)
    hits <- parsed$data$hits
    if (is.null(hits) || length(hits) == 0) break

    dt <- rbindlist(lapply(hits, function(hit) {
      data.table(
        file_id = hit$file_id %||% NA_character_,
        file_name = hit$file_name %||% NA_character_,
        data_category = hit$data_category %||% NA_character_,
        data_type = hit$data_type %||% NA_character_
      )
    }), fill = TRUE)

    all_rows[[length(all_rows) + 1L]] <- dt
    if (nrow(dt) < page_size) break
    from <- from + page_size
  }

  if (length(all_rows) == 0) return(data.table())
  unique(rbindlist(all_rows, fill = TRUE), by = "file_id")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

download_gdc_file <- function(file_id, out_file, attempts = 3L) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  url <- paste0(gdc_data_endpoint, "/", file_id)

  for (try_i in seq_len(attempts)) {
    ok <- tryCatch({
      resp <- GET(url, write_disk(out_file, overwrite = TRUE), timeout(600))
      status_code(resp) == 200L
    }, error = function(e) FALSE)

    if (ok && file.exists(out_file) && file.info(out_file)$size > 0) return(TRUE)
    if (try_i < attempts) Sys.sleep(1 + try_i)
  }
  FALSE
}

manifest <- data.table(
  profile = character(),
  file_id = character(),
  file_name = character(),
  local_file = character(),
  status = character()
)

for (profile in PROFILES) {
  cat("──", profile$label, "──\n")
  out_base <- file.path(DATA_ROOT, profile$out_subdir)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  files_dt <- query_files(profile$filters)
  if (!is.na(FETCH_LIMIT)) files_dt <- head(files_dt, FETCH_LIMIT)
  cat("Files discovered:", nrow(files_dt), "\n")

  if (nrow(files_dt) == 0) {
    cat("No files found for this profile.\n\n")
    next
  }

  for (i in seq_len(nrow(files_dt))) {
    fid <- files_dt$file_id[i]
    fname <- files_dt$file_name[i]
    out_file <- file.path(out_base, fid, fname)

    if (file.exists(out_file) && file.info(out_file)$size > 0) {
      status <- "skipped_existing"
    } else {
      status <- if (download_gdc_file(fid, out_file)) "downloaded" else "failed"
    }

    manifest <- rbind(
      manifest,
      data.table(
        profile = profile$label,
        file_id = fid,
        file_name = fname,
        local_file = out_file,
        status = status
      )
    )

    if (i %% 50 == 0 || i == nrow(files_dt)) {
      cat("  Progress:", i, "/", nrow(files_dt), "\n")
    }
  }
  cat("\n")
}

manifest_path <- file.path(DATA_ROOT, "gdc_fetch_manifest.csv")
fwrite(manifest, manifest_path)

cat("=== Fetch Summary ===\n")
if (nrow(manifest) > 0) {
  print(manifest[, .N, by = .(profile, status)][order(profile, status)])
} else {
  cat("No records in manifest.\n")
}
cat("Manifest:", manifest_path, "\n")
cat("Done.\n")
