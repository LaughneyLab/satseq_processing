rm(list=ls()) #Clear environment

#################
## libraries ####
#################
library(MAST)
library(rhdf5)
library(Matrix)
library(data.table)
library(stringr)
library(dplyr)
library(parallel)

options(mc.cores = 4)
set.seed(123)

#######################
## USER CONFIG ########
#######################

h5file <- "../mt269/satseq/matt_cbsulaughney/10x/lung_mets/adata_save_2026_07_05.h5ad"

# Raw adata.obs$celltype category strings for the 3 cell types of interest.
# "4T1" and "AM" are taken directly from sc_analysis_lung_mets.py (AM is the
# RAW label there — the script only renames it to "Alveolar macrophages" for
# plotting, after loading). "T cells" is inferred from your message, not
# confirmed against the raw data — the script prints every celltype it
# actually finds below, so check that table on first run and fix this line
# if the literal string doesn't match.
celltypes_of_interest <- c("4T1", "T cells", "AM")

variants_of_interest  <- c("P13", "K337", "STOP")   # each compared vs "REF"

out_root <- "MAST_output_lung_mets"

# Combinations with fewer cells than this in either group are skipped (with
# a warning) instead of being allowed to crash zlm() on too little data
min_cells_per_group <- 10

#######Functions##########
csr_list_to_sparse_anndata <- function(file, path) {
  lst   <- h5read(file, path)
  attrs <- h5readAttributes(file, path)

  x      <- as.numeric(lst$data)
  j      <- as.integer(lst$indices)
  indptr <- as.integer(lst$indptr)

  # shape might be an attribute of the layer; if not, fall back to /X
  if (!is.null(attrs$shape)) {
    shape <- as.integer(attrs$shape)
  } else {
    X_attrs <- h5readAttributes(file, "/X")
    shape   <- as.integer(X_attrs$shape)
  }

  nrow <- shape[1L]
  ncol <- shape[2L]

  row_counts <- diff(indptr)
  i0 <- rep(seq_len(nrow) - 1L, row_counts)

  i <- i0 + 1L       # 1-based rows
  j <- j    + 1L     # 0-based -> 1-based cols

  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(nrow, ncol))
}

## given h5ad anndata file, return given column from adata.obs
readObsFromh5File <- function(infile, colname) {
  tmp <- h5ls(infile)
  w <- which(tmp$name == colname & tmp$otype=="H5I_DATASET")
  if (length(w) == 1) {
    return(h5read(infile, sprintf("/obs/%s", colname)))
  }
  w <- which(tmp$group == sprintf("/obs/%s", colname))
  if (length(w) == 1) {
    return(h5read(infile, sprintf("/obs/%s", colname)))
  }
  if (length(w) == 2) {
    w2 <- which(tmp[w,]$name=="categories")
    catNames <- h5read(infile, sprintf("%s/%s", tmp[w,][w2,"group"], tmp[w,][w2,"name"]))
    w2 <- which(tmp[w,]$name == "codes")
    return(catNames[as.numeric(h5read(infile, sprintf("%s/%s", tmp[w,][w2, "group"], tmp[w,][w2, "name"])))+1])
  }
  stop("Error reading ", colname, " from ", infile)
}

# For zlm output file (Melissa provided)
fixNAs <- function(fcHurdle, cluster, clusters, geneNames, counts) {

  # 1. Identify the p-value column (usually 'Pr(>Chisq)') using a keyword search
  pcol <- grep("Chisq", names(fcHurdle), value=TRUE)
  if (length(pcol) != 1) {
    cat("Error finding p-value column in fixNAs")
    return(fcHurdle)
  }

  # 2. Find rows where the coefficient is missing (NA) but the test actually ran (p-value < 1)
  # This targets genes where the Hurdle model failed to converge due to 0/100% expression
  fixNA <- (is.na(fcHurdle$coef) & fcHurdle[,pcol] < 1)

  # 3. Define arbitrary "infinite" bounds for the missing coefficients
  coefRange <- c(-10000, 10000)

  # 4. Map the genes in the results (primerid) to their row indices in the raw count matrix
  geneOrder <- sapply(fcHurdle$primerid, function(x) {which(geneNames == x)})

  # 5. Calculate 'fracIn': Fraction of cells expressing the gene INSIDE the cluster of interest
  # sum() is used on the counts; if counts are binarized or counts > 0, this identifies frequency
  fracIn <- apply(counts[geneOrder, clusters == cluster, drop=FALSE], 1, sum) / sum(clusters == cluster)

  # 6. Calculate 'fracOut': Fraction of cells expressing the gene in ALL OTHER clusters (background)
  fracOut <- apply(counts[geneOrder, clusters != cluster, drop=FALSE], 1, sum) / sum(clusters != cluster)

  # 7. Logic Check: If the gene is more frequent in the background, assign the large negative value
  fcHurdle[fixNA & fracIn <= fracOut, "coef"] <- coefRange[1]

  # 8. Logic Check: If the gene is more frequent in the target cluster, assign the large positive value
  fcHurdle[fixNA & fracIn > fracOut, "coef"] <- coefRange[2]

  # Return the data frame with the formerly NA coefficients now filled with +/- 10000
  fcHurdle
}


## Read Data ######

counts <- csr_list_to_sparse_anndata(h5file, "X")

# Extract cell and gene names
geneNames <- as.character(h5read(h5file, "/var/_index"))
cellNames <- as.character(h5read(h5file, "/obs/_index"))

# Variant: first "_"-delimited token of `sample`, mirroring
# sc_analysis_lung_mets.py's
#   adata.obs["variant"] = adata.obs["sample"].astype(str).str.split("_").str.get(0)
sample_col <- as.character(readObsFromh5File(h5file, "sample"))
pvariant   <- sub("_.*$", "", sample_col)

celltype <- as.character(readObsFromh5File(h5file, "celltype"))

# Keep only cells that are (a) one of the 4 known variant labels and
# (b) one of the 3 cell types we're actually going to test — trims the
# matrix before the expensive per-combination steps below
known_variants <- c("P13", "K337", "STOP", "REF")
keep <- !is.na(pvariant) & pvariant %in% known_variants &
        !is.na(celltype) & celltype %in% celltypes_of_interest

counts    <- counts[keep, ]
cellNames <- cellNames[keep]
pvariant  <- pvariant[keep]
celltype  <- celltype[keep]

# Drop genes never expressed in any retained cell (cheap global prefilter —
# the real >20%-of-selected-cells filter is reapplied per celltype/variant
# combination below, on the much smaller subset)
filter_unexpressed_genes <- which(Matrix::colSums(counts) != 0)
geneNames <- geneNames[filter_unexpressed_genes]
counts <- counts[, filter_unexpressed_genes]

cat("Cells retained after variant/celltype filtering:", length(cellNames), "\n")
cat("Celltype counts (verify these match `celltypes_of_interest` above):\n")
print(table(celltype))
cat("Variant counts:\n")
print(table(pvariant))

### run MAST, per cell type x variant ######

dir.create(out_root, showWarnings = FALSE)

for (ct in celltypes_of_interest) {

  ct_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  out_dir <- file.path(out_root, paste0(ct_safe, "_variant_vs_REF_lung_mets"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!(ct %in% celltype)) {
    warning("Celltype '", ct, "' not found in data — check spelling against the table printed above. Skipping.")
    next
  }

  for (v in variants_of_interest) {

    filter_MAST   <- which(celltype == ct & pvariant %in% c("REF", v))
    pvariant_MAST <- pvariant[filter_MAST]

    n_ref <- sum(pvariant_MAST == "REF")
    n_var <- sum(pvariant_MAST == v)
    cat(sprintf("\n[%s vs REF | %s] REF=%d cells, %s=%d cells\n", v, ct, n_ref, v, n_var))

    if (n_ref < min_cells_per_group || n_var < min_cells_per_group) {
      warning(sprintf("Skipping %s vs REF for %s: too few cells (REF=%d, %s=%d, min=%d)",
                       v, ct, n_ref, v, n_var, min_cells_per_group))
      next
    }

    counts_MAST    <- as.matrix(counts[filter_MAST, ])
    cellNames_MAST <- cellNames[filter_MAST]

    ## Subset genes expressed in >20% of selected cells
    filter_MAST_genes <- colSums(counts_MAST != 0) / length(cellNames_MAST) > 0.2
    counts_MAST_f  <- counts_MAST[, filter_MAST_genes]
    geneNames_MAST <- geneNames[filter_MAST_genes]

    ## Construct sca for MAST
    sca_MAST <- FromMatrix(t(counts_MAST_f),
                         cData = data.frame(condition=pvariant_MAST,
                                    wellKey=cellNames_MAST,
                                    barcode=cellNames_MAST),
                         fData = data.frame(symbolid=geneNames_MAST,
                                    primerid=geneNames_MAST,
                                    geneName=geneNames_MAST)
    )
    ## Calculate cngeneson
    colData(sca_MAST)$cngeneson <- scale(colSums(assay(sca_MAST)>0))

    # Set reference level
    cond <- factor(colData(sca_MAST)$condition)
    cond <- relevel(cond, "REF")
    colData(sca_MAST)$condition <- cond

    # Run MAST (https://rglab.github.io/MAST/articles/MAITAnalysis.html)
    zlmResults <- zlm(~condition + cngeneson, sca_MAST)

    # Summarize and print results
    zt <- summary(zlmResults, doLRT=paste0('condition', v))$datatable

    fcHurdle <- merge(zt[contrast==paste0('condition', v) & component=='H',.(primerid, `Pr(>Chisq)`)],
                      zt[contrast==paste0('condition', v) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')

    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    setorder(fcHurdle, fdr)

    fcHurdle <- fixNAs(
      fcHurdle = as.data.frame(fcHurdle),
      cluster = v,
      clusters = pvariant_MAST,
      geneNames = geneNames_MAST,
      counts = t(counts_MAST_f)
    )

    # Add to full data and save
    fcHurdle$variant  <- v
    fcHurdle$celltype <- ct

    out_file <- file.path(out_dir, paste0("fcHurdle_", ct_safe, "_", v, ".csv"))
    write.csv(fcHurdle, out_file, row.names = FALSE)
    cat("Wrote", out_file, "\n")
  }
}
