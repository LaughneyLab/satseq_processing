rm(list = ls())

# ============================================================
#  run_MAST_celltype_variant_vs_REF.R
#
#  Usage (called by SLURM array job via submit_MAST_celltype.sh):
#    Rscript --vanilla run_MAST_celltype_variant_vs_REF.R <task_id>
#
#  task_id encodes both the cell type AND the comparison:
#    task_id = (variant_index - 1) * n_celltypes + celltype_index
#  where:
#    variant_index  in 1..3  (P13, K337, STOP  – all vs REF)
#    celltype_index in 1..11 (see CELL_TYPES below)
#  Total tasks = 3 * 11 = 33
#
#  Outputs:
#    MAST_output/raw/<celltype>_<variant>_vs_REF.csv
#
#  After all 33 jobs finish, run the companion collation script:
#    Rscript collate_MAST_to_excel.R
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

# ---- Libraries ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(MAST)
  library(rhdf5)
  library(Matrix)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(parallel)
})

options(mc.cores = 4)
set.seed(123)

# ---- Parameters --------------------------------------------------------------

H5FILE <- "../mt269/satseq/matt_cbsulaughney/10x/lung_mets/adata_save_2026_07_05.h5ad"

# Variants to compare against REF (order defines variant_index 1,2,3)
VARIANTS <- c("P13", "K337", "STOP")

# Cell types (order defines celltype_index 1..11)
# These must match values in adata.obs["detailed_celltype"] exactly.
# Adjust spellings if your h5ad uses different labels.
CELL_TYPES <- c(
"Neutrophil",
  "Interstitial macrophages",
  "Alveolar macrophages",
  "Naive CD4+ T cells",
  "B cells",
  "Fibroblasts",
  "cDCs",
  "Naive CD8+ T cells",
  "Tregs",
  "Inflammatory monocytes",
  "Activated CD4+ T cells",
  "pDCs",
  "Exhausted CD8+ T cells",
  "γδ T cells",
  "NK cells")


N_CT  <- length(CELL_TYPES)
N_VAR <- length(VARIANTS)

# Decode task_id -> (variant_index, celltype_index)
# task_id runs 1..(N_VAR * N_CT)
# variant_index  = ceiling(task_id / N_CT)
# celltype_index = ((task_id - 1) %% N_CT) + 1
variant_index  <- ceiling(task_id / N_CT)
celltype_index <- ((task_id - 1) %% N_CT) + 1

variant_name  <- VARIANTS[variant_index]
celltype_name <- CELL_TYPES[celltype_index]

cat("=== Task", task_id, "===\n")
cat("  Variant  :", variant_name, "\n")
cat("  Cell type:", celltype_name, "\n\n")

# ---- Output directory --------------------------------------------------------
out_dir <- "MAST_output/raw"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_ct <- gsub("[^A-Za-z0-9_]", "_", celltype_name)
out_csv <- file.path(out_dir, paste0(safe_ct, "_", variant_name, "_vs_REF.csv"))

if (file.exists(out_csv)) {
  cat("Output already exists, skipping:", out_csv, "\n")
  quit(status = 0)
}

# ---- Helper functions --------------------------------------------------------

# Read a CSR sparse matrix from an h5ad /X slot
csr_list_to_sparse_anndata <- function(file, path) {
  lst   <- h5read(file, path)
  attrs <- h5readAttributes(file, path)

  x      <- as.numeric(lst$data)
  j      <- as.integer(lst$indices)
  indptr <- as.integer(lst$indptr)

  if (!is.null(attrs$shape)) {
    shape <- as.integer(attrs$shape)
  } else {
    X_attrs <- h5readAttributes(file, "/X")
    shape   <- as.integer(X_attrs$shape)
  }

  nrow <- shape[1L]; ncol <- shape[2L]
  row_counts <- diff(indptr)
  i0  <- rep(seq_len(nrow) - 1L, row_counts)
  i   <- i0 + 1L
  j   <- j  + 1L

  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(nrow, ncol))
}

# Read a column from /obs (handles plain dataset, factor, or encoded category)
readObsFromH5 <- function(infile, colname) {
  tmp <- h5ls(infile)
  w   <- which(tmp$name == colname & tmp$otype == "H5I_DATASET")
  if (length(w) == 1)
    return(h5read(infile, sprintf("/obs/%s", colname)))
  w <- which(tmp$group == sprintf("/obs/%s", colname))
  if (length(w) == 1)
    return(h5read(infile, sprintf("/obs/%s", colname)))
  if (length(w) == 2) {
    w2       <- which(tmp[w, ]$name == "categories")
    catNames <- h5read(infile, sprintf("%s/%s",
                                       tmp[w, ][w2, "group"],
                                       tmp[w, ][w2, "name"]))
    w2 <- which(tmp[w, ]$name == "codes")
    return(catNames[as.numeric(h5read(infile,
                                       sprintf("%s/%s",
                                               tmp[w, ][w2, "group"],
                                               tmp[w, ][w2, "name"]))) + 1])
  }
  stop("Cannot read ", colname, " from ", infile)
}

# Fix NA coefficients for genes with 0%/100% expression (from Melissa)
fixNAs <- function(fcHurdle, cluster, clusters, geneNames, counts) {
  pcol <- grep("Chisq", names(fcHurdle), value = TRUE)
  if (length(pcol) != 1) {
    cat("Warning: could not locate p-value column in fixNAs\n")
    return(fcHurdle)
  }
  fixNA    <- is.na(fcHurdle$coef) & fcHurdle[, pcol] < 1
  coefRange <- c(-10000, 10000)
  geneOrder <- sapply(fcHurdle$primerid, function(x) which(geneNames == x))
  fracIn    <- apply(counts[geneOrder, clusters == cluster, drop = FALSE], 1, sum) /
               sum(clusters == cluster)
  fracOut   <- apply(counts[geneOrder, clusters != cluster, drop = FALSE], 1, sum) /
               sum(clusters != cluster)
  fcHurdle[fixNA & fracIn <= fracOut, "coef"] <- coefRange[1]
  fcHurdle[fixNA & fracIn >  fracOut, "coef"] <- coefRange[2]
  fcHurdle
}

# ---- Read h5ad ---------------------------------------------------------------
cat("Reading h5ad file...\n")
counts    <- csr_list_to_sparse_anndata(H5FILE, "X")
geneNames <- as.character(h5read(H5FILE, "/var/_index"))
cellNames <- as.character(h5read(H5FILE, "/obs/_index"))

# --- Derive variant from 'sample', exactly as the Python script does:
#     adata.obs["variant"] = adata.obs["sample"].str.split("_").str.get(0)
sample_raw  <- readObsFromH5(H5FILE, "sample")
variant_raw <- sub("_.*", "", sample_raw)   # keep everything before the first "_"

# --- Derive cell type, mirroring the Python logic exactly:
#     Start from celltype; overwrite with T_cell_subtype where non-null;
#     then apply the case-insensitive mapping.
celltype_raw <- trimws(readObsFromH5(H5FILE, "celltype"))

has_t_subtype <- tryCatch({
  tmp <- h5ls(H5FILE)
  any(tmp$name == "T_cell_subtype" | tmp$group == "/obs/T_cell_subtype")
}, error = function(e) FALSE)

if (has_t_subtype) {
  t_sub  <- readObsFromH5(H5FILE, "T_cell_subtype")
  t_mask <- !is.na(t_sub) & trimws(as.character(t_sub)) != "nan"
  celltype_raw[t_mask] <- trimws(as.character(t_sub))[t_mask]
}

celltype_map <- c(
  "am"                = "Alveolar macrophages",
  "endothelial cells" = "Endothelial",
  "endothelial"       = "Endothelial",
  "im"                = "Interstitial macrophages",
  "ilc"               = "Innate lymphoid cells",
  "imon"              = "Inflammatory monocytes",
  "scmf"              = "Secondary crest myofibroblasts",
  "4t1"               = "Tumor cells"
)
mapped <- celltype_map[tolower(celltype_raw)]
keep   <- !is.na(mapped)
celltype_raw[keep] <- mapped[keep]

# ---- Summary: variants and cell types in the data ---------------------------
var_counts <- sort(table(variant_raw), decreasing = TRUE)
cat("Variants found in data (n cells):\n")
cat(paste0("  ", names(var_counts), ": ", var_counts, collapse = "\n"), "\n\n")

ct_counts <- sort(table(celltype_raw), decreasing = TRUE)
cat("Cell types found in data (n cells):\n")
cat(paste0("  ", names(ct_counts), ": ", ct_counts, collapse = "\n"), "\n\n")

# ---- Filter to REF + focal variant, focal cell type -------------------------
filter_variant <- which(variant_raw %in% c("REF", variant_name))
filter_ct      <- which(trimws(celltype_raw) == celltype_name)
filter_cells   <- intersect(filter_variant, filter_ct)

n_focal <- sum(variant_raw[filter_cells] == variant_name)
n_ref   <- sum(variant_raw[filter_cells] == "REF")
cat(sprintf("Cells after filtering: %d (%s=%d, REF=%d)\n",
            length(filter_cells), variant_name, n_focal, n_ref))

if (n_focal < 10 || n_ref < 10) {
  cat("WARNING: Too few cells for reliable MAST analysis. Skipping.\n")
  empty <- data.frame(primerid = character(0), `Pr(>Chisq)` = numeric(0),
                      coef = numeric(0), ci.hi = numeric(0), ci.lo = numeric(0),
                      fdr = numeric(0), check.names = FALSE)
  write.csv(empty, out_csv, row.names = FALSE)
  quit(status = 0)
}

counts_sub     <- as.matrix(counts[filter_cells, ])
cellNames_sub  <- cellNames[filter_cells]
variant_sub    <- variant_raw[filter_cells]

rm(counts); gc()

# ---- Gene filtering: expressed in > 20% of selected cells -------------------
expr_frac     <- colSums(counts_sub != 0) / nrow(counts_sub)
gene_keep     <- expr_frac > 0.20
counts_sub    <- counts_sub[, gene_keep]
geneNames_sub <- geneNames[gene_keep]
cat(sprintf("Genes after filtering (>20%% expressed): %d\n", sum(gene_keep)))

# ---- Build MAST SingleCellAssay ---------------------------------------------
sca <- FromMatrix(
  t(counts_sub),
  cData = data.frame(
    condition = variant_sub,
    wellKey   = cellNames_sub,
    barcode   = cellNames_sub
  ),
  fData = data.frame(
    symbolid  = geneNames_sub,
    primerid  = geneNames_sub,
    geneName  = geneNames_sub
  )
)

# Cellular detection rate covariate
colData(sca)$cngeneson <- scale(colSums(assay(sca) > 0))

# Set REF as reference level
cond <- factor(colData(sca)$condition)
cond <- relevel(cond, "REF")
colData(sca)$condition <- cond

# ---- Run MAST hurdle model --------------------------------------------------
cat("Fitting zlm...\n")
zlmRes <- zlm(~ condition + cngeneson, sca)

cat("Summarising LRT...\n")
contrast_name <- paste0("condition", variant_name)
zt <- summary(zlmRes, doLRT = contrast_name)$datatable

# ---- Extract hurdle test results --------------------------------------------
fcHurdle <- merge(
  zt[contrast == contrast_name & component == "H",
     .(primerid, `Pr(>Chisq)`)],
  zt[contrast == contrast_name & component == "logFC",
     .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, "fdr")]
setorder(fcHurdle, fdr)

fcHurdle <- fixNAs(
  fcHurdle  = as.data.frame(fcHurdle),
  cluster   = variant_name,
  clusters  = variant_sub,
  geneNames = geneNames_sub,
  counts    = t(counts_sub)
)

# ---- Save raw output ---------------------------------------------------------
fcHurdle$celltype <- celltype_name
fcHurdle$variant  <- variant_name

write.csv(fcHurdle, out_csv, row.names = FALSE)
cat("Saved:", out_csv, "\n")
cat("Done.\n")
