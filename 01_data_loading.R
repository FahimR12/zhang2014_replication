#!/usr/bin/env Rscript
# =============================================================================
# 01_data_loading.R
# Replication of Zhang et al. (2014) "Integrative network analysis of TCGA data
# for ovarian cancer" — BMC Systems Biology 8:1338
#
# Step 1: Load raw GDC-downloaded TCGA-OV data and assemble sample × feature
#         matrices for gene expression, CNV, DNA methylation, somatic mutation.
#         Map file-level UUIDs to TCGA patient barcodes via the GDC API.
#
# Parallelization: all four data types use Windows-safe PSOCK cluster
# (parallel::parLapply), avoiding the fork-based mclapply that fails on Windows.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(httr)
  library(parallel)
})

# ── Configuration ──────────────────────────────────────────────────────────────
DATA_ROOT <- "C:/Users/fahim/Desktop/scripts/BNPipeline/GDCdata/TCGA-OV"
args_full <- commandArgs(trailingOnly = FALSE)
file_arg  <- "--file="
script_path <- sub(file_arg, "", args_full[grep(file_arg, args_full)])
PROJECT_DIR <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
OUT_DIR  <- file.path(PROJECT_DIR, "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GE_DIR   <- file.path(DATA_ROOT, "Transcriptome_Profiling",    "Gene_Expression_Quantification")
CNV_DIR  <- file.path(DATA_ROOT, "Copy_Number_Variation",      "Copy_Number_Segment")
METH_DIR <- file.path(DATA_ROOT, "DNA_Methylation",            "Methylation_Beta_Value")
MUT_DIR  <- file.path(DATA_ROOT, "Simple_Nucleotide_Variation","Masked_Somatic_Mutation")

# Number of parallel workers.  On Windows, PSOCK clusters can saturate disk I/O
# past ~8 workers; tune down if reads become slower rather than faster.
N_CORES <- min(8L, max(1L, detectCores() - 1L))

cat("=== Zhang et al. 2014 Replication — Step 1: Data Loading ===\n")
cat("Data root:", DATA_ROOT, "\n")
cat("Output dir:", OUT_DIR, "\n")
cat("Parallel workers:", N_CORES, "\n\n")

# =============================================================================
# HELPERS
# =============================================================================

# ── TCGA barcode utilities ────────────────────────────────────────────────────
patient_from_barcode    <- function(bc) substr(bc, 1, 12)
sample_type_from_barcode <- function(bc) as.integer(substr(bc, 14, 15))

# ── Platform-safe .gz reader ─────────────────────────────────────────────────
# Replaces fread(cmd = "zcat …") which is unavailable on Windows.
# Strategy (in order of preference):
#   1. data.table native gz support  – fread(path)  works on v1.14.0+ (2021)
#   2. R gzfile() connection fallback – slower but universal
read_gz_maf <- function(maf_file, select_cols) {
  if (!file.exists(maf_file)) return(NULL)

  # --- Attempt 1: data.table native (recommended, fast) ---
  dat <- tryCatch({
    d <- fread(maf_file, sep = "\t", skip = "Hugo_Symbol",
               select = select_cols, showProgress = FALSE)
    if (nrow(d) > 0) d else NULL
  }, warning = function(w) NULL,
     error   = function(e) NULL)

  if (!is.null(dat)) return(dat)

  # --- Attempt 2: R gzfile() connection (universal fallback) ---
  tryCatch({
    con <- if (grepl("\\.gz$", maf_file)) gzfile(maf_file, "rt") else file(maf_file, "rt")
    on.exit(close(con), add = TRUE)

    lines <- readLines(con)
    header_idx <- which(startsWith(lines, "Hugo_Symbol"))[1]
    if (is.na(header_idx)) return(NULL)

    fread(text = paste(lines[header_idx:length(lines)], collapse = "\n"),
          sep = "\t", select = select_cols, showProgress = FALSE)
  }, error = function(e) NULL)
}

# ── UUID → TCGA barcode mapping via GDC API ───────────────────────────────────
map_uuids_to_barcodes <- function(uuid_dirs) {
  file_uuids <- basename(uuid_dirs)
  cat("  Querying GDC API for", length(file_uuids), "file UUIDs...\n")

  batch_size <- 100
  all_mappings <- data.table(file_id = character(), barcode = character(),
                             sample_type = character())

  for (i in seq(1, length(file_uuids), by = batch_size)) {
    batch <- file_uuids[i:min(i + batch_size - 1, length(file_uuids))]

    body <- list(
      filters = list(op = "in", content = list(field = "files.file_id", value = batch)),
      format  = "JSON",
      fields  = paste0("file_id,cases.samples.sample_type,",
                       "cases.samples.portions.analytes.aliquots.submitter_id"),
      size    = as.character(length(batch))
    )

    resp <- tryCatch(
      POST("https://api.gdc.cancer.gov/files",
           body = toJSON(body, auto_unbox = TRUE),
           content_type_json(), timeout(60)),
      error = function(e) NULL
    )

    if (is.null(resp) || status_code(resp) != 200) next

    result <- tryCatch(
      fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE),
      error = function(e) NULL
    )
    if (is.null(result)) next

    hits <- result$data$hits
    if (is.null(hits) || length(hits) == 0) next

    for (j in seq_len(nrow(hits))) {
      fid   <- hits$file_id[j]
      cases <- hits$cases[[j]]
      if (is.null(cases)) next
      for (k in seq_len(nrow(cases))) {
        stype    <- cases$sample_type[k]
        portions <- cases$portions[[k]]
        if (is.null(portions)) next
        analytes <- portions$analytes[[1]]
        if (is.null(analytes)) next
        aliquots <- analytes$aliquots[[1]]
        if (is.null(aliquots)) next
        bc <- aliquots$submitter_id[1]
        all_mappings <- rbind(all_mappings,
          data.table(file_id = fid, barcode = bc, sample_type = stype))
      }
    }

    if (i %% 500 == 1 && i > 1)
      cat("    Progress:", min(i + batch_size - 1, length(file_uuids)),
          "/", length(file_uuids), "\n")
    Sys.sleep(0.3)
  }

  cat("  Mapped", nrow(all_mappings), "UUID → barcode entries\n")
  all_mappings
}

# =============================================================================
# 1. GENE EXPRESSION  (426 samples, parallelised)
# =============================================================================
cat("── Loading Gene Expression Data ──────────────────────────\n")

ge_uuid_dirs <- list.dirs(GE_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(ge_uuid_dirs), "expression samples\n")

ge_mapping <- map_uuids_to_barcodes(ge_uuid_dirs)

# Reference gene list from first file
first_tsv <- list.files(ge_uuid_dirs[1], full.names = TRUE, pattern = "\\.tsv$")[1]
template   <- fread(first_tsv, skip = "gene_id", sep = "\t",
                    select = c("gene_id", "gene_name", "gene_type"))
gene_info  <- template[grepl("^ENSG", gene_id)]
cat("  Genes in annotation:", nrow(gene_info), "\n")

# ── Parallel read ─────────────────────────────────────────────────────────────
cat("  Reading", length(ge_uuid_dirs), "files using", N_CORES, "cores...\n")

cl <- makeCluster(N_CORES, type = "PSOCK")
clusterEvalQ(cl, suppressPackageStartupMessages(library(data.table)))
clusterExport(cl, c("ge_uuid_dirs", "gene_info", "ge_mapping"))

ge_results <- parLapply(cl, seq_along(ge_uuid_dirs), function(i) {
  uuid     <- basename(ge_uuid_dirs[i])
  tsv_file <- list.files(ge_uuid_dirs[i], full.names = TRUE, pattern = "\\.tsv$")[1]
  if (is.na(tsv_file)) return(list(vals = rep(NA_real_, nrow(gene_info)), bc = uuid))

  dat <- tryCatch(
    fread(tsv_file, skip = "gene_id", sep = "\t",
          select = c("gene_name", "fpkm_unstranded"), showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat)) return(list(vals = rep(NA_real_, nrow(gene_info)), bc = uuid))

  dat <- dat[grepl("^(TSPAN|ENSG)", gene_name) | gene_name %in% gene_info$gene_name]
  idx  <- match(gene_info$gene_name, dat$gene_name)
  vals <- dat$fpkm_unstranded[idx]

  bc_row <- ge_mapping[file_id == uuid]
  bc     <- if (nrow(bc_row) > 0 && !is.na(bc_row$barcode[1])) bc_row$barcode[1] else uuid

  list(vals = vals, bc = bc)
})

stopCluster(cl)

# Assemble matrix
ge_matrix            <- matrix(NA_real_, nrow = nrow(gene_info), ncol = length(ge_uuid_dirs))
rownames(ge_matrix)  <- gene_info$gene_name
sample_barcodes_ge   <- sapply(ge_results, `[[`, "bc")
for (i in seq_along(ge_results)) ge_matrix[, i] <- ge_results[[i]]$vals
rm(ge_results); gc()

colnames(ge_matrix) <- sample_barcodes_ge
ge_matrix           <- log2(ge_matrix + 1)                        # log2(FPKM+1)

gene_vars           <- apply(ge_matrix, 1, var, na.rm = TRUE)
ge_matrix           <- ge_matrix[!is.na(gene_vars) & gene_vars > 0, ]
gene_info_filtered  <- gene_info[gene_name %in% rownames(ge_matrix)]

sample_types_ge  <- sample_type_from_barcode(colnames(ge_matrix))
is_tumor_ge      <- sample_types_ge %in% 1:9
is_normal_ge     <- sample_types_ge %in% 10:14

cat("  Expression matrix:", nrow(ge_matrix), "genes ×", ncol(ge_matrix), "samples\n")
cat("  Tumor:", sum(is_tumor_ge, na.rm=TRUE), "| Normal:", sum(is_normal_ge, na.rm=TRUE), "\n")

# =============================================================================
# 2. COPY NUMBER VARIATION  (1564 samples, parallelised)
# =============================================================================
cat("\n── Loading Copy Number Variation Data ─────────────────────\n")

cnv_uuid_dirs <- list.dirs(CNV_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(cnv_uuid_dirs), "CNV samples\n")

cnv_mapping <- map_uuids_to_barcodes(cnv_uuid_dirs)

cat("  Reading", length(cnv_uuid_dirs), "segment files using", N_CORES, "cores...\n")

cl <- makeCluster(N_CORES, type = "PSOCK")
clusterEvalQ(cl, suppressPackageStartupMessages(library(data.table)))
clusterExport(cl, c("cnv_uuid_dirs", "cnv_mapping"))

cnv_results <- parLapply(cl, seq_along(cnv_uuid_dirs), function(i) {
  uuid     <- basename(cnv_uuid_dirs[i])
  seg_file <- list.files(cnv_uuid_dirs[i], full.names = TRUE, pattern = "\\.txt$")[1]
  if (is.na(seg_file)) return(NULL)

  seg <- tryCatch(fread(seg_file, sep = "\t", showProgress = FALSE), error = function(e) NULL)
  if (is.null(seg) || nrow(seg) == 0) return(NULL)

  bc_row <- cnv_mapping[file_id == uuid]
  bc     <- if (nrow(bc_row) > 0 && !is.na(bc_row$barcode[1])) bc_row$barcode[1] else uuid

  seg[, sample_barcode := bc]
  seg[, file_uuid      := uuid]
  seg
})

stopCluster(cl)

all_segments <- rbindlist(Filter(Negate(is.null), cnv_results), fill = TRUE)
rm(cnv_results); gc()

cat("  Segments:", nrow(all_segments), "| Unique samples:", length(unique(all_segments$sample_barcode)), "\n")

# =============================================================================
# 3. DNA METHYLATION  (592 samples, parallelised)
# =============================================================================
cat("\n── Loading DNA Methylation Data ───────────────────────────\n")

meth_uuid_dirs <- list.dirs(METH_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(meth_uuid_dirs), "methylation samples\n")

meth_mapping <- map_uuids_to_barcodes(meth_uuid_dirs)

# Probe list from first file
first_meth_file <- list.files(meth_uuid_dirs[1], full.names = TRUE, pattern = "\\.txt$")[1]
first_meth      <- fread(first_meth_file, header = FALSE, col.names = c("probe_id", "beta"))
probe_ids       <- first_meth$probe_id
cat("  Probes per sample:", length(probe_ids), "\n")

cat("  Reading", length(meth_uuid_dirs), "files using", N_CORES, "cores...\n")

cl <- makeCluster(N_CORES, type = "PSOCK")
clusterEvalQ(cl, suppressPackageStartupMessages(library(data.table)))
clusterExport(cl, c("meth_uuid_dirs", "probe_ids", "meth_mapping"))

meth_results <- parLapply(cl, seq_along(meth_uuid_dirs), function(i) {
  uuid      <- basename(meth_uuid_dirs[i])
  meth_file <- list.files(meth_uuid_dirs[i], full.names = TRUE, pattern = "\\.txt$")[1]
  if (is.na(meth_file)) return(list(vals = rep(NA_real_, length(probe_ids)), bc = uuid))

  dat <- tryCatch(
    fread(meth_file, header = FALSE, col.names = c("probe_id", "beta"),
          showProgress = FALSE, colClasses = c("character", "numeric")),
    error = function(e) NULL
  )
  if (is.null(dat)) return(list(vals = rep(NA_real_, length(probe_ids)), bc = uuid))

  idx  <- match(probe_ids, dat$probe_id)
  vals <- dat$beta[idx]

  bc_row <- meth_mapping[file_id == uuid]
  bc     <- if (nrow(bc_row) > 0 && !is.na(bc_row$barcode[1])) bc_row$barcode[1] else uuid

  list(vals = vals, bc = bc)
})

stopCluster(cl)

# Assemble matrix
meth_matrix           <- matrix(NA_real_, nrow = length(probe_ids), ncol = length(meth_uuid_dirs))
rownames(meth_matrix) <- probe_ids
sample_barcodes_meth  <- sapply(meth_results, `[[`, "bc")
for (i in seq_along(meth_results)) meth_matrix[, i] <- meth_results[[i]]$vals
rm(meth_results); gc()

colnames(meth_matrix) <- sample_barcodes_meth
cat("  Methylation matrix:", nrow(meth_matrix), "probes ×", ncol(meth_matrix), "samples\n")

# =============================================================================
# 4. SOMATIC MUTATION  (482 files, parallelised, platform-safe .gz reading)
# =============================================================================
cat("\n── Loading Somatic Mutation Data ──────────────────────────\n")

mut_uuid_dirs <- list.dirs(MUT_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(mut_uuid_dirs), "mutation files\n")
cat("  Using platform-safe .gz reader (no zcat required)\n")

MAF_COLS <- c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "Variant_Type")

cat("  Reading", length(mut_uuid_dirs), "MAF files using", N_CORES, "cores...\n")

cl <- makeCluster(N_CORES, type = "PSOCK")
clusterEvalQ(cl, suppressPackageStartupMessages(library(data.table)))
clusterExport(cl, c("mut_uuid_dirs", "MAF_COLS", "read_gz_maf"))  # export helper

mut_results <- parLapply(cl, seq_along(mut_uuid_dirs), function(i) {
  uuid     <- basename(mut_uuid_dirs[i])
  maf_file <- list.files(mut_uuid_dirs[i], full.names = TRUE,
                          pattern = "\\.maf\\.gz$|\\.maf$")[1]
  if (is.na(maf_file)) return(NULL)

  dat <- read_gz_maf(maf_file, MAF_COLS)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)

  dat[, file_uuid := uuid]
  dat
})

stopCluster(cl)

all_mutations <- rbindlist(Filter(Negate(is.null), mut_results), fill = TRUE)
rm(mut_results); gc()

if (nrow(all_mutations) == 0) {
  stop("FATAL: No mutations loaded. Check your MAF files and data.table version (need >= 1.14).\n",
       "  Run: packageVersion('data.table')  — update if < 1.14.0")
}

cat("  Total mutations:", nrow(all_mutations), "\n")
cat("  Unique genes:", length(unique(all_mutations$Hugo_Symbol)), "\n")
cat("  Unique samples:", length(unique(all_mutations$Tumor_Sample_Barcode)), "\n")

# Binary variable: 1 = non-silent, 0 = silent / not mutated  (per paper)
SILENT_CLASSES <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
                    "IGR", "RNA", "lincRNA")
all_mutations[, is_nonsilent := !(Variant_Classification %in% SILENT_CLASSES)]

nonsilent         <- all_mutations[is_nonsilent == TRUE]
unique_genes_mut  <- unique(nonsilent$Hugo_Symbol)
unique_samples_mut <- unique(nonsilent$Tumor_Sample_Barcode)

cat("  Non-silent mutations — genes:", length(unique_genes_mut),
    "| samples:", length(unique_samples_mut), "\n")

# Build binary matrix (gene × sample)
mut_matrix            <- matrix(0L, nrow = length(unique_genes_mut),
                                    ncol = length(unique_samples_mut))
rownames(mut_matrix)  <- unique_genes_mut
colnames(mut_matrix)  <- unique_samples_mut

# Fast fill via data.table cross-join
nons_dt <- nonsilent[, .(Hugo_Symbol, Tumor_Sample_Barcode)]
nons_dt <- unique(nons_dt)
row_idx  <- match(nons_dt$Hugo_Symbol,         unique_genes_mut)
col_idx  <- match(nons_dt$Tumor_Sample_Barcode, unique_samples_mut)
valid    <- !is.na(row_idx) & !is.na(col_idx)
mut_matrix[cbind(row_idx[valid], col_idx[valid])] <- 1L

cat("  Mutation matrix:", nrow(mut_matrix), "genes ×", ncol(mut_matrix), "samples\n")
cat("  Sparsity:", round(1 - mean(mut_matrix), 4), "\n")

# =============================================================================
# 5. SAVE ALL MATRICES
# =============================================================================
cat("\n── Saving Matrices ───────────────────────────────────────\n")

saveRDS(ge_matrix,          file.path(OUT_DIR, "ge_matrix.rds"))
saveRDS(gene_info_filtered, file.path(OUT_DIR, "gene_info.rds"))
saveRDS(all_segments,       file.path(OUT_DIR, "cnv_segments.rds"))
saveRDS(meth_matrix,        file.path(OUT_DIR, "meth_matrix.rds"))
saveRDS(mut_matrix,         file.path(OUT_DIR, "mut_matrix.rds"))

# UUID mappings
saveRDS(ge_mapping,   file.path(OUT_DIR, "ge_uuid_mapping.rds"))
saveRDS(cnv_mapping,  file.path(OUT_DIR, "cnv_uuid_mapping.rds"))
saveRDS(meth_mapping, file.path(OUT_DIR, "meth_uuid_mapping.rds"))

# Sample metadata
sample_info <- data.table(
  barcode          = colnames(ge_matrix),
  patient          = patient_from_barcode(colnames(ge_matrix)),
  sample_type_code = sample_type_from_barcode(colnames(ge_matrix)),
  is_tumor         = is_tumor_ge,
  is_normal        = is_normal_ge
)
saveRDS(sample_info, file.path(OUT_DIR, "sample_info.rds"))

cat("  All matrices saved to:", OUT_DIR, "\n")

cat("\n=== Summary ===\n")
cat("Gene Expression:  ", nrow(ge_matrix), "genes ×", ncol(ge_matrix), "samples\n")
cat("CNV Segments:     ", nrow(all_segments), "segments,",
    length(unique(all_segments$sample_barcode)), "samples\n")
cat("DNA Methylation:  ", nrow(meth_matrix), "probes ×", ncol(meth_matrix), "samples\n")
cat("Somatic Mutation: ", nrow(mut_matrix), "genes ×", ncol(mut_matrix), "samples\n")
cat("\nStep 1 complete. Run 02_preprocessing.R next.\n")
