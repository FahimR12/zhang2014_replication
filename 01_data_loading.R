#!/usr/bin/env Rscript
# =============================================================================
# 01_data_loading.R
# Replication of Zhang et al. (2014) "Integrative network analysis of TCGA data
# for ovarian cancer" — BMC Systems Biology 8:1338
#
# Step 1: Load raw GDC-downloaded TCGA-OV data and assemble sample × feature
#         matrices for gene expression, CNV, DNA methylation, somatic mutation.
#         Map file-level UUIDs to TCGA patient barcodes via the GDC API.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(httr)
})

# ── Configuration ──────────────────────────────────────────────────────────────
# CHANGE THIS to your local TCGA-OV download path
DATA_ROOT <- "C:/Users/fahim/Desktop/scripts/BNPipeline/GDCdata/TCGA-OV"
OUT_DIR   <- file.path(dirname(DATA_ROOT), "..", "zhang2014_replication", "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GE_DIR   <- file.path(DATA_ROOT, "Transcriptome_Profiling", "Gene_Expression_Quantification")
CNV_DIR  <- file.path(DATA_ROOT, "Copy_Number_Variation", "Copy_Number_Segment")
METH_DIR <- file.path(DATA_ROOT, "DNA_Methylation", "Methylation_Beta_Value")
MUT_DIR  <- file.path(DATA_ROOT, "Simple_Nucleotide_Variation", "Masked_Somatic_Mutation")

cat("=== Zhang et al. 2014 Replication — Step 1: Data Loading ===\n")
cat("Data root:", DATA_ROOT, "\n")
cat("Output dir:", OUT_DIR, "\n\n")

# =============================================================================
# HELPER: Map GDC file UUIDs → TCGA barcodes via the GDC API
# =============================================================================
map_uuids_to_barcodes <- function(uuid_dirs, data_dir) {
  # Each UUID folder name is the file UUID in GDC
  file_uuids <- basename(uuid_dirs)

  cat("  Querying GDC API for", length(file_uuids), "file UUIDs...\n")

  # Query in batches of 100
  batch_size <- 100
  all_mappings <- data.table(file_id = character(), case_id = character(),
                             barcode = character(), sample_type = character())

  for (i in seq(1, length(file_uuids), by = batch_size)) {
    batch <- file_uuids[i:min(i + batch_size - 1, length(file_uuids))]

    body <- list(
      filters = list(
        op = "in",
        content = list(
          field = "files.file_id",
          value = batch
        )
      ),
      format = "JSON",
      fields = paste0(
        "file_id,cases.case_id,",
        "cases.samples.sample_type,",
        "cases.samples.portions.analytes.aliquots.submitter_id"
      ),
      size = as.character(length(batch))
    )

    resp <- tryCatch({
      POST("https://api.gdc.cancer.gov/files",
           body = toJSON(body, auto_unbox = TRUE),
           content_type_json(),
           timeout(60))
    }, error = function(e) {
      cat("    GDC API error (batch", ceiling(i/batch_size), "):", e$message, "\n")
      return(NULL)
    })

    if (is.null(resp) || status_code(resp) != 200) {
      cat("    Warning: GDC API returned status",
          ifelse(is.null(resp), "NULL", status_code(resp)),
          "for batch", ceiling(i/batch_size), "\n")
      next
    }

    result <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
    hits <- result$data$hits

    if (length(hits) == 0 || nrow(hits) == 0) next

    for (j in seq_len(nrow(hits))) {
      fid <- hits$file_id[j]
      cases <- hits$cases[[j]]
      if (is.null(cases) || length(cases) == 0) next

      for (k in seq_len(nrow(cases))) {
        case_id <- cases$case_id[k]
        samples <- cases$samples[[k]]
        if (is.null(samples) || length(samples) == 0) next

        sample_type <- samples$sample_type[1]
        # Extract TCGA barcode from aliquot submitter_id
        barcode <- NA_character_
        portions <- samples$portions[[1]]
        if (!is.null(portions) && length(portions) > 0) {
          analytes <- portions$analytes[[1]]
          if (!is.null(analytes) && length(analytes) > 0) {
            aliquots <- analytes$aliquots[[1]]
            if (!is.null(aliquots) && length(aliquots) > 0) {
              barcode <- aliquots$submitter_id[1]
            }
          }
        }
        all_mappings <- rbind(all_mappings, data.table(
          file_id = fid, case_id = case_id,
          barcode = barcode, sample_type = sample_type
        ))
      }
    }

    if (i %% 500 == 1) cat("    Processed", min(i + batch_size - 1, length(file_uuids)),
                           "/", length(file_uuids), "\n")
    Sys.sleep(0.3)  # Be polite to the API
  }

  cat("  Mapped", nrow(all_mappings), "file-to-barcode entries\n")
  return(all_mappings)
}

# Fallback: use filename patterns to extract TCGA barcodes from MAF files
extract_barcode_from_maf <- function(maf_file) {
  # MAF files contain TCGA barcodes in the Tumor_Sample_Barcode column
  lines <- fread(maf_file, nrows = 1, skip = "Hugo_Symbol", sep = "\t",
                 select = "Tumor_Sample_Barcode")
  if (nrow(lines) > 0) return(lines$Tumor_Sample_Barcode[1])
  return(NA_character_)
}

# Helper: extract patient ID from TCGA barcode (first 12 chars: TCGA-XX-XXXX)
patient_from_barcode <- function(barcode) {
  substr(barcode, 1, 12)
}

# Helper: extract sample type code from barcode (chars 14-15)
# 01 = Primary Tumor, 10 = Blood Normal, 11 = Solid Tissue Normal
sample_type_from_barcode <- function(barcode) {
  as.integer(substr(barcode, 14, 15))
}

# =============================================================================
# 1. GENE EXPRESSION — Build sample × gene matrix (FPKM values)
# =============================================================================
cat("── Loading Gene Expression Data ──────────────────────────\n")

ge_uuid_dirs <- list.dirs(GE_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(ge_uuid_dirs), "gene expression samples\n")

# Map UUIDs to barcodes
ge_mapping <- map_uuids_to_barcodes(ge_uuid_dirs, GE_DIR)

# Read first file to get gene list
first_file <- list.files(ge_uuid_dirs[1], full.names = TRUE, pattern = "\\.tsv$")[1]
template <- fread(first_file, skip = "gene_id", sep = "\t")
# Keep only ENSG entries (protein coding genes etc.), skip summary rows
template <- template[grepl("^ENSG", gene_id)]
gene_info <- template[, .(gene_id, gene_name, gene_type)]

cat("  Total genes in annotation:", nrow(gene_info), "\n")

# Build expression matrix — use FPKM (unstranded) as in the paper (Agilent microarray equivalent)
# The paper used microarray; we use RNA-seq FPKM as the modern equivalent
ge_matrix <- matrix(NA_real_, nrow = nrow(gene_info), ncol = length(ge_uuid_dirs))
rownames(ge_matrix) <- gene_info$gene_name
sample_barcodes_ge <- character(length(ge_uuid_dirs))

cat("  Reading expression files...\n")
pb <- txtProgressBar(min = 0, max = length(ge_uuid_dirs), style = 3)

for (i in seq_along(ge_uuid_dirs)) {
  uuid <- basename(ge_uuid_dirs[i])
  tsv_file <- list.files(ge_uuid_dirs[i], full.names = TRUE, pattern = "\\.tsv$")[1]

  if (is.na(tsv_file) || !file.exists(tsv_file)) {
    setTxtProgressBar(pb, i)
    next
  }

  dat <- fread(tsv_file, skip = "gene_id", sep = "\t",
               select = c("gene_id", "gene_name", "fpkm_unstranded"))
  dat <- dat[grepl("^ENSG", gene_id)]

  # Match to template gene order
  idx <- match(gene_info$gene_name, dat$gene_name)
  ge_matrix[, i] <- dat$fpkm_unstranded[idx]

  # Get barcode
  bc_row <- ge_mapping[file_id == uuid]
  if (nrow(bc_row) > 0 && !is.na(bc_row$barcode[1])) {
    sample_barcodes_ge[i] <- bc_row$barcode[1]
  } else {
    sample_barcodes_ge[i] <- uuid  # fallback
  }

  setTxtProgressBar(pb, i)
}
close(pb)

colnames(ge_matrix) <- sample_barcodes_ge

# Log2 transform (add small constant to avoid log(0))
ge_matrix <- log2(ge_matrix + 1)

# Remove genes with zero variance
gene_vars <- apply(ge_matrix, 1, var, na.rm = TRUE)
ge_matrix <- ge_matrix[!is.na(gene_vars) & gene_vars > 0, ]
gene_info_filtered <- gene_info[gene_name %in% rownames(ge_matrix)]

cat("  Expression matrix:", nrow(ge_matrix), "genes ×", ncol(ge_matrix), "samples\n")

# Separate tumor vs normal
sample_types_ge <- sample_type_from_barcode(colnames(ge_matrix))
is_tumor_ge <- sample_types_ge %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9)
is_normal_ge <- sample_types_ge %in% c(10, 11, 12, 13, 14)
cat("  Tumor samples:", sum(is_tumor_ge, na.rm = TRUE),
    "| Normal samples:", sum(is_normal_ge, na.rm = TRUE), "\n")

# =============================================================================
# 2. COPY NUMBER VARIATION — Build sample × gene matrix (segment means)
# =============================================================================
cat("\n── Loading Copy Number Variation Data ─────────────────────\n")

cnv_uuid_dirs <- list.dirs(CNV_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(cnv_uuid_dirs), "CNV samples\n")

# Map UUIDs to barcodes
cnv_mapping <- map_uuids_to_barcodes(cnv_uuid_dirs, CNV_DIR)

# We need gene coordinates to map segments to genes
# Download gene coordinates from UCSC or use a built-in reference
# For reproducibility, we'll create a minimal gene location table from Ensembl
cat("  Downloading gene coordinates from Ensembl BioMart...\n")

get_gene_coordinates <- function() {
  # Try to download from Ensembl BioMart REST API
  url <- paste0(
    "https://rest.ensembl.org/biomart/query?query=",
    URLencode(paste0(
      '<?xml version="1.0" encoding="UTF-8"?>',
      '<!DOCTYPE Query>',
      '<Query virtualSchemaName="default" formatter="TSV" header="1">',
      '<Dataset name="hsapiens_gene_ensembl" interface="default">',
      '<Attribute name="hgnc_symbol"/>',
      '<Attribute name="chromosome_name"/>',
      '<Attribute name="start_position"/>',
      '<Attribute name="end_position"/>',
      '<Attribute name="strand"/>',
      '</Dataset>',
      '</Query>'
    ))
  )

  # Alternative: use a simpler approach with gencode info from our gene expression files
  # We already have gene_info with gene_id (ENSG) and gene_name
  # We'll use the GDC API to get gene coordinates

  # Simplest approach: read coordinates from the GDC gene model annotation
  # For now, use a pre-built mapping
  cache_file <- file.path(OUT_DIR, "gene_coordinates.rds")
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  # Fetch from Ensembl REST API (works reliably)
  cat("    Fetching gene coordinates from Ensembl REST API...\n")
  all_genes <- unique(gene_info_filtered$gene_name)

  # Query in batches
  gene_coords <- data.table(
    gene_name = character(), chr = character(),
    start = integer(), end = integer()
  )

  batch_size <- 200
  for (i in seq(1, length(all_genes), by = batch_size)) {
    batch_genes <- all_genes[i:min(i + batch_size - 1, length(all_genes))]

    for (g in batch_genes) {
      resp <- tryCatch({
        GET(paste0("https://rest.ensembl.org/lookup/symbol/homo_sapiens/", g),
            content_type("application/json"),
            timeout(10))
      }, error = function(e) NULL)

      if (!is.null(resp) && status_code(resp) == 200) {
        info <- fromJSON(content(resp, "text", encoding = "UTF-8"))
        gene_coords <- rbind(gene_coords, data.table(
          gene_name = g,
          chr = as.character(info$seq_region_name),
          start = info$start,
          end = info$end
        ))
      }
      Sys.sleep(0.05)
    }

    if (i %% 1000 == 1) {
      cat("    Gene coordinates:", nrow(gene_coords), "/", length(all_genes), "\n")
    }
  }

  saveRDS(gene_coords, cache_file)
  return(gene_coords)
}

# Alternative faster approach: parse gene locations from GENCODE GTF
# or use the simpler segment-level summary
cat("  Building gene-level CNV from segment data...\n")
cat("  (This maps each gene to its overlapping CNV segment mean)\n")

# For efficiency, we first collect all segments, then map to genes
# Step 1: Read all segment files
all_segments <- list()
sample_barcodes_cnv <- character(length(cnv_uuid_dirs))

pb <- txtProgressBar(min = 0, max = length(cnv_uuid_dirs), style = 3)
for (i in seq_along(cnv_uuid_dirs)) {
  uuid <- basename(cnv_uuid_dirs[i])
  seg_file <- list.files(cnv_uuid_dirs[i], full.names = TRUE, pattern = "\\.txt$")[1]

  if (is.na(seg_file) || !file.exists(seg_file)) {
    setTxtProgressBar(pb, i)
    next
  }

  seg <- fread(seg_file, sep = "\t")
  seg$file_uuid <- uuid

  # Get barcode
  bc_row <- cnv_mapping[file_id == uuid]
  if (nrow(bc_row) > 0 && !is.na(bc_row$barcode[1])) {
    sample_barcodes_cnv[i] <- bc_row$barcode[1]
  } else {
    # Try to get from GDC_Aliquot column
    sample_barcodes_cnv[i] <- uuid
  }
  seg$sample_barcode <- sample_barcodes_cnv[i]

  all_segments[[i]] <- seg
  setTxtProgressBar(pb, i)
}
close(pb)

all_segments <- rbindlist(all_segments, fill = TRUE)
cat("  Total segments:", nrow(all_segments), "\n")
cat("  Unique samples:", length(unique(all_segments$sample_barcode)), "\n")

# Save raw segments for gene mapping in preprocessing step
# (Gene mapping requires coordinates which we'll handle in script 02)

# =============================================================================
# 3. DNA METHYLATION — Build sample × probe matrix (beta values)
# =============================================================================
cat("\n── Loading DNA Methylation Data ───────────────────────────\n")

meth_uuid_dirs <- list.dirs(METH_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(meth_uuid_dirs), "methylation samples\n")

# Map UUIDs to barcodes
meth_mapping <- map_uuids_to_barcodes(meth_uuid_dirs, METH_DIR)

# Read first file to get probe list
first_meth_file <- list.files(meth_uuid_dirs[1], full.names = TRUE, pattern = "\\.txt$")[1]
first_meth <- fread(first_meth_file, header = FALSE, col.names = c("probe_id", "beta"))
probe_ids <- first_meth$probe_id

cat("  Total probes:", length(probe_ids), "\n")

# Build methylation matrix
meth_matrix <- matrix(NA_real_, nrow = length(probe_ids), ncol = length(meth_uuid_dirs))
rownames(meth_matrix) <- probe_ids
sample_barcodes_meth <- character(length(meth_uuid_dirs))

cat("  Reading methylation files...\n")
pb <- txtProgressBar(min = 0, max = length(meth_uuid_dirs), style = 3)

for (i in seq_along(meth_uuid_dirs)) {
  uuid <- basename(meth_uuid_dirs[i])
  meth_file <- list.files(meth_uuid_dirs[i], full.names = TRUE, pattern = "\\.txt$")[1]

  if (is.na(meth_file) || !file.exists(meth_file)) {
    setTxtProgressBar(pb, i)
    next
  }

  dat <- fread(meth_file, header = FALSE, col.names = c("probe_id", "beta"))
  idx <- match(probe_ids, dat$probe_id)
  meth_matrix[, i] <- dat$beta[idx]

  # Get barcode
  bc_row <- meth_mapping[file_id == uuid]
  if (nrow(bc_row) > 0 && !is.na(bc_row$barcode[1])) {
    sample_barcodes_meth[i] <- bc_row$barcode[1]
  } else {
    sample_barcodes_meth[i] <- uuid
  }

  setTxtProgressBar(pb, i)
}
close(pb)

colnames(meth_matrix) <- sample_barcodes_meth
cat("  Methylation matrix:", nrow(meth_matrix), "probes ×", ncol(meth_matrix), "samples\n")

# =============================================================================
# 4. SOMATIC MUTATION — Build sample × gene binary matrix
# =============================================================================
cat("\n── Loading Somatic Mutation Data ──────────────────────────\n")

mut_uuid_dirs <- list.dirs(MUT_DIR, full.names = TRUE, recursive = FALSE)
cat("  Found", length(mut_uuid_dirs), "mutation files\n")

# Read all MAF files and combine
# Each MAF file may contain mutations for multiple samples (aliquot-merged)
all_mutations <- list()

cat("  Reading MAF files...\n")
pb <- txtProgressBar(min = 0, max = length(mut_uuid_dirs), style = 3)

for (i in seq_along(mut_uuid_dirs)) {
  uuid <- basename(mut_uuid_dirs[i])
  maf_file <- list.files(mut_uuid_dirs[i], full.names = TRUE, pattern = "\\.maf\\.gz$|\\.maf$")[1]

  if (is.na(maf_file) || !file.exists(maf_file)) {
    setTxtProgressBar(pb, i)
    next
  }

  # Read MAF — select only key columns for memory efficiency
  dat <- tryCatch({
    fread(cmd = paste("zcat", shQuote(maf_file)),
          sep = "\t", skip = "Hugo_Symbol",
          select = c("Hugo_Symbol", "Variant_Classification",
                     "Tumor_Sample_Barcode", "Variant_Type"))
  }, error = function(e) {
    tryCatch({
      fread(maf_file, sep = "\t", skip = "Hugo_Symbol",
            select = c("Hugo_Symbol", "Variant_Classification",
                       "Tumor_Sample_Barcode", "Variant_Type"))
    }, error = function(e2) NULL)
  })

  if (!is.null(dat) && nrow(dat) > 0) {
    dat$file_uuid <- uuid
    all_mutations[[i]] <- dat
  }

  setTxtProgressBar(pb, i)
}
close(pb)

all_mutations <- rbindlist(all_mutations, fill = TRUE)
cat("  Total mutations:", nrow(all_mutations), "\n")
cat("  Unique genes with mutations:", length(unique(all_mutations$Hugo_Symbol)), "\n")
cat("  Unique samples:", length(unique(all_mutations$Tumor_Sample_Barcode)), "\n")

# Per the paper: binary variable — "1" for non-silent mutations, "0" for silent/none
# Non-silent = Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, Frame_Shift_Ins,
#              In_Frame_Del, In_Frame_Ins, Splice_Site, Translation_Start_Site, Nonstop_Mutation
silent_classes <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
                    "IGR", "RNA", "lincRNA")
all_mutations$is_nonsilent <- !(all_mutations$Variant_Classification %in% silent_classes)

# Create binary matrix: genes × samples
nonsilent <- all_mutations[is_nonsilent == TRUE]
unique_genes_mut <- unique(nonsilent$Hugo_Symbol)
unique_samples_mut <- unique(nonsilent$Tumor_Sample_Barcode)

cat("  Genes with non-silent mutations:", length(unique_genes_mut), "\n")
cat("  Samples with non-silent mutations:", length(unique_samples_mut), "\n")

# Build sparse binary matrix
mut_matrix <- matrix(0L, nrow = length(unique_genes_mut), ncol = length(unique_samples_mut))
rownames(mut_matrix) <- unique_genes_mut
colnames(mut_matrix) <- unique_samples_mut

# Fill in mutations
for (g in unique_genes_mut) {
  samples_with_mut <- unique(nonsilent[Hugo_Symbol == g]$Tumor_Sample_Barcode)
  mut_matrix[g, samples_with_mut] <- 1L
}

cat("  Mutation matrix:", nrow(mut_matrix), "genes ×", ncol(mut_matrix), "samples\n")
cat("  Sparsity:", round(1 - sum(mut_matrix) / length(mut_matrix), 4), "\n")

# =============================================================================
# 5. SAVE ALL MATRICES
# =============================================================================
cat("\n── Saving Matrices ───────────────────────────────────────\n")

saveRDS(ge_matrix, file.path(OUT_DIR, "ge_matrix.rds"))
saveRDS(gene_info_filtered, file.path(OUT_DIR, "gene_info.rds"))
saveRDS(all_segments, file.path(OUT_DIR, "cnv_segments.rds"))
saveRDS(meth_matrix, file.path(OUT_DIR, "meth_matrix.rds"))
saveRDS(mut_matrix, file.path(OUT_DIR, "mut_matrix.rds"))

# Save barcode mappings
saveRDS(ge_mapping, file.path(OUT_DIR, "ge_uuid_mapping.rds"))
saveRDS(cnv_mapping, file.path(OUT_DIR, "cnv_uuid_mapping.rds"))
saveRDS(meth_mapping, file.path(OUT_DIR, "meth_uuid_mapping.rds"))

# Save sample type info
sample_info <- data.table(
  barcode = colnames(ge_matrix),
  patient = patient_from_barcode(colnames(ge_matrix)),
  sample_type_code = sample_type_from_barcode(colnames(ge_matrix)),
  is_tumor = is_tumor_ge,
  is_normal = is_normal_ge
)
saveRDS(sample_info, file.path(OUT_DIR, "sample_info.rds"))

cat("  All matrices saved to:", OUT_DIR, "\n")

cat("\n=== Summary ===\n")
cat("Gene Expression:   ", nrow(ge_matrix), "genes ×", ncol(ge_matrix), "samples\n")
cat("CNV Segments:      ", nrow(all_segments), "segments across",
    length(unique(all_segments$sample_barcode)), "samples\n")
cat("DNA Methylation:   ", nrow(meth_matrix), "probes ×", ncol(meth_matrix), "samples\n")
cat("Somatic Mutation:  ", nrow(mut_matrix), "genes ×", ncol(mut_matrix), "samples\n")
cat("\nStep 1 complete. Run 02_preprocessing.R next.\n")
