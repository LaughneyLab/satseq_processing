rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
pvariant_index <- as.integer(args[1]) 

library(MAST)
library(rhdf5)
library(Matrix)
library(data.table)
library(stringr)
library(dplyr)
library(parallel)

options(mc.cores = 4)
set.seed(123)

####### Functions (Keep existing csr_list_to_sparse_anndata, readObsFromh5File, and fixNAs) ##########
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

# Process variant data
load_processed_data <- function(h5file, label) {
  counts <- csr_list_to_sparse_anndata(h5file, "X")
  genes  <- as.character(h5read(h5file, "/var/_index"))
  cells  <- as.character(h5read(h5file, "/obs/_index"))
  var_sub <- readObsFromh5File(h5file, 'variant_substitution')
  
  # Pool per position
  parts <- str_match(var_sub, '^([a-zA-Z]+)([0-9.]+)([a-zA-Z]+)$')
  pvar <- str_c(parts[, 2], parts[, 3])
  pvar[var_sub == 'mSTING_WT'] <- "REF"
  pvar[var_sub == 'STOP'] <- "STOP"
  
  # Filter "None"
  valid_idx <- which(var_sub != "None")
  return(list(counts = counts[valid_idx, ], genes = genes,
              cells = cells[valid_idx],
              pvariant = pvar[valid_idx], cond = label))
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

# 1. Load Data and Metadata for files
## Paths
file_cgamp <- "../mt269/satseq/10x/adata_EO771_cGAMP.h5ad"
file_h2o   <- "../mt269/satseq/10x/adata_EO771_H2O.h5ad"

data_cgamp <- load_processed_data(file_cgamp, "cGAMP")
data_h2o   <- load_processed_data(file_h2o, "H2O")

# 2. Find common genes and subset
common_genes <- intersect(data_cgamp$genes, data_h2o$genes)
idx_cgamp_genes <- match(common_genes, data_cgamp$genes)
idx_h2o_genes   <- match(common_genes, data_h2o$genes)

# 3. Identify the variant to analyze
# Note: We assume the variant index refers to the set of variants present in both or the cGAMP file
pvariant_list <- unique(data_cgamp$pvariant)
target_variant <- pvariant_list[pvariant_index]
print(paste("Analyzing Variant:", target_variant))

# 4. SUBSET MODIFICATION: Now include BOTH the variant and the REF (reference)
# Assuming data_cgamp and data_h2o have a 'pvariant' column containing "REF" or the target
sub_cgamp_idx <- which(data_cgamp$pvariant %in% c(target_variant, "REF"))
sub_h2o_idx   <- which(data_h2o$pvariant %in% c(target_variant, "REF"))

# Create merged matrix and metadata
counts_MAST <- rbind(
  data_cgamp$counts[sub_cgamp_idx, idx_cgamp_genes],
  data_h2o$counts[sub_h2o_idx, idx_h2o_genes]
)

meta_MAST <- data.frame(
  wellKey   = c(data_cgamp$cells[sub_cgamp_idx], data_h2o$cells[sub_h2o_idx]),
  condition = c(rep("cGAMP", length(sub_cgamp_idx)), rep("H2O", length(sub_h2o_idx))),
  mutation  = c(data_cgamp$pvariant[sub_cgamp_idx], data_h2o$pvariant[sub_h2o_idx])
)

geneNames_MAST <- common_genes
counts_MAST <- as.matrix(counts_MAST)

# 5. Filter genes expressed in >20% of selected cells
filter_genes <- colSums(counts_MAST != 0) / nrow(counts_MAST) > 0.2
counts_MAST <- counts_MAST[, filter_genes]
geneNames_MAST <- geneNames_MAST[filter_genes]

## Construct sca for MAST
sca_MAST <- FromMatrix(t(counts_MAST),
                       cData = meta_MAST,
                       fData = data.frame(symbolid=geneNames_MAST, primerid=geneNames_MAST)
)

colData(sca_MAST)$cngeneson <- scale(colSums(assay(sca_MAST) > 0))

## Prepare factors for interaction
# Set H2O as reference for condition
colData(sca_MAST)$condition <- relevel(factor(colData(sca_MAST)$condition), "H2O")
# Set REF as reference for mutation
colData(sca_MAST)$mutation  <- relevel(factor(colData(sca_MAST)$mutation), "REF")

# 6. Run MAST with Interaction
# Formula: condition * mutation includes Main Effects + Interaction
zlmResults <- zlm(~condition * mutation + cngeneson, sca_MAST)

# Define the Interaction Term Name 
# MAST combines factor names/levels: 'conditioncGAMP:mutation[VariantName]'
intTerm <- paste0("conditioncGAMP:mutation", target_variant)

# Test the interaction effect specifically
zt <- summary(zlmResults, doLRT = intTerm)$datatable

fcHurdle <- merge(zt[contrast == intTerm & component == 'H', .(primerid, `Pr(>Chisq)`)],
                  zt[contrast == intTerm & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], 
                  by = 'primerid')

fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
setorder(fcHurdle, fdr)

# Use your fixNAs function
fcHurdle <- fixNAs(
  fcHurdle = as.data.frame(fcHurdle), 
  cluster = "cGAMP", 
  clusters = meta_MAST$condition, 
  geneNames = geneNames_MAST, 
  counts = t(counts_MAST)
)

fcHurdle$variant <- target_variant

# 7. Save
write.csv(fcHurdle, paste0('MAST_output_cGAMP_vs_H2O/EO771_interaction/fcHurdle_', target_variant, '.csv'))

  




