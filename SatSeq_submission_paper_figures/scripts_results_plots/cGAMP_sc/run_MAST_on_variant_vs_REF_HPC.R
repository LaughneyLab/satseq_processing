rm(list=ls())#Clear environment
args = commandArgs(trailingOnly=TRUE)

# Inputs
pvariant_index <- as.integer(args[1]) #input variant ndex

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


# Read gene expression matrix
h5file <- "../mt269/satseq/10x/adata_EO771_cGAMP.h5ad"

counts <- csr_list_to_sparse_anndata(h5file,"X")

# Extract cell, gene, variant names
geneNames <- as.character(h5read(h5file, "/var/_index"))
cellNames <- as.character(h5read(h5file, "/obs/_index"))

variant <- readObsFromh5File(h5file,'variant_substitution')

# Pool per position
parts <- str_match(variant, '^([a-zA-Z]+)([0-9.]+)([a-zA-Z]+)$')
pvariant <- str_c(parts[, 2], parts[, 3])
pvariant[variant=='mSTING_WT'] <- "REF"
pvariant[variant=='STOP'] <- "STOP"
pvariant[variant=='dCTT'] <- "dCTT"
pvariant[variant=='RRDD'] <- "RRDD"

# Filterings
# Variants which are None
filter_none_variant <- which(variant != "None")
counts <- counts[filter_none_variant,]
cellNames <- cellNames[filter_none_variant]
variant <- variant[filter_none_variant]
pvariant <- pvariant[filter_none_variant]

# Filter genes which are never expressed in any cell
filter_unexpressed_genes <- which(colSums(counts)!=0)
geneNames <- geneNames[filter_unexpressed_genes]
counts <- counts[,filter_unexpressed_genes]

### run MAST ######

## subset to REF and a variant
pvariant_list <- unique(pvariant)
fcHurdle_df <- data.frame()
filter_MAST <- which(pvariant %in% c("REF",pvariant_list[pvariant_index]))

counts_MAST <- as.matrix(counts[filter_MAST,])
cellNames_MAST <- cellNames[filter_MAST]
pvariant_MAST <- pvariant[filter_MAST]

print(c(pvariant_index,pvariant_list[pvariant_index],length(pvariant_MAST)))

rm(counts)

## Subset genes expressed in >20% of selected cells
filter_MAST_genes <- colSums(counts_MAST!=0)/(length(cellNames_MAST))>0.2
counts_MAST <- counts_MAST[,filter_MAST_genes]
geneNames_MAST <- geneNames[filter_MAST_genes]

## Construct sca for MAST
sca_MAST <- FromMatrix(t(counts_MAST),
                     cData = data.frame(condition=pvariant_MAST,
                                wellKey=cellNames_MAST,
                                barcode=cellNames_MAST),
                     fData = data.frame(symbolid=geneNames_MAST,
                                primerid=geneNames_MAST,
                                geneName=geneNames_MAST)
)
## Calculate cngeneson
colData(sca_MAST)$cngeneson <- scale(colSums(assay(sca_MAST)>0))

# Set refrence level
cond<-factor(colData(sca_MAST)$condition)
cond<-relevel(cond,"REF")
colData(sca_MAST)$condition<-cond

# Run MAST (https://rglab.github.io/MAST/articles/MAITAnalysis.html)
zlmResults <- zlm(~condition + cngeneson, sca_MAST)

# Summarize and print results
zt <- summary(zlmResults, doLRT=paste0('condition',pvariant_list[pvariant_index]))$datatable

fcHurdle <- merge(zt[contrast==paste0('condition',pvariant_list[pvariant_index]) & component=='H',.(primerid, `Pr(>Chisq)`)],
                  zt[contrast==paste0('condition',pvariant_list[pvariant_index]) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
setorder(fcHurdle, fdr)

fcHurdle <- fixNAs(
  fcHurdle = as.data.frame(fcHurdle), 
  cluster = pvariant_list[pvariant_index], 
  clusters = pvariant_MAST, 
  geneNames = geneNames_MAST, 
  counts = t(counts_MAST)
)

# Add to full data and save
fcHurdle$variant <- pvariant_list[pvariant_index]

fcHurdle_df <- rbind(fcHurdle_df,fcHurdle)

write.csv(fcHurdle_df,paste0('MAST_output_cGAMP_vs_H2O/EO771_variant_vs_REF_cGAMP/fcHurdle_EO771_',
                             pvariant_list[pvariant_index],'.csv'))






