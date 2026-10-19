rm(list=ls())#Clear environment

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
# 1. Load contact matrices
bonds_AA  <- read.csv('bonds_monomer.csv', row.names = 1)
bonds_AB  <- read.csv('bonds_dimer.csv', row.names = 1)
bonds_ACD <- read.csv('bonds_tetramer.csv', row.names = 1)

# 2. Load fitness data
fitness_data <- read.csv("FACS_mets_weights_all.csv")

# Mets cutoff for hi/lo
mets_cutoff <- sd(fitness_data$mets_add)

# 3. Build the wide dataframe using fast column-summing
mets_bonds <- data.frame(
  position = 1:378,
  monomer  = colSums(bonds_AA[, 1:378, drop = FALSE], na.rm = TRUE),
  dimer    = colSums(bonds_AB[, 1:378, drop = FALSE], na.rm = TRUE),
  tetramer = colSums(bonds_ACD[, 1:378, drop = FALSE], na.rm = TRUE)
) %>% 
  # 4. Bring in the fitness column directly
  left_join(select(fitness_data, position, variant, mets_add,
                   STING_add, IFN_add, NFkB_add, ATF3_add), by = "position") %>%
# 1. Calculate fractions and categorize groups
mutate(
    total_bonds = monomer + dimer + tetramer,
    monomer_fraction = ifelse(total_bonds > 0, monomer / total_bonds, 0),
    dimer_fraction   = ifelse(total_bonds > 0, dimer / total_bonds, 0),
    
    # Classify based on your fractional rules
    bond_group = case_when(
      monomer_fraction > dimer_fraction & mets_add > mets_cutoff ~ "mets_hi_monomer",
      dimer_fraction > monomer_fraction & mets_add > mets_cutoff ~ "mets_hi_dimer",
      mets_add <= mets_cutoff ~ "mets_lo")) %>%
  filter(total_bonds > 0, !is.na(bond_group)) 

# Read gene expression matrix
h5file <- "../mt269/satseq/10x/adata_4T1_tumor.h5ad"
counts <- csr_list_to_sparse_anndata(h5file,"X")

# Extract cell, gene, variant names
geneNames <- as.character(h5read(h5file, "/var/_index"))
cellNames <- as.character(h5read(h5file, "/obs/_index"))
variant <- readObsFromh5File(h5file,'variant_substitution')

# 3. Process variant labels into a pooled condition map
# Extract the numeric position from strings like "P13A" or "E71G"
parts <- str_match(variant, '^([a-zA-Z]+)([0-9.]+)([a-zA-Z]+)$')
res_positions <- as.integer(parts[, 3]) # Pulls out the numeric residue sequence index

# Build condition classifications
pvariant <- rep("Other", length(variant))
pvariant[variant == 'mSTING_WT'] <- "REF"

# Map variants to their specific bond_group
for (i in seq_along(variant)) {
  if (variant[i] != "mSTING_WT" && variant[i] != "None") {
    pos <- res_positions[i]
    # Match the numeric position in mets_bonds to assign the bond_group
    match_row <- mets_bonds[mets_bonds$position == pos, ]
    if (nrow(match_row) > 0) {
      pvariant[i] <- match_row$bond_group[1]
    }
  }
}

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

##################
## Run MAST ######
##################

# --- MAST setup function ---
run_mast_test <- function(counts_matrix, cells, variants, group1, group2, suffix) {
  
  # Filter down specifically to the two groups being compared
  filter_MAST <- which(variants %in% c(group1, group2))
  if (length(filter_MAST) == 0) {
    message(paste("No cells found for comparison:", group1, "vs", group2))
    return(NULL)
  }
  
  counts_MAST <- as.matrix(counts_matrix[filter_MAST, ])
  cellNames_MAST <- cells[filter_MAST]
  pvariant_MAST <- variants[filter_MAST]
  
  print(paste("Running MAST for", group1, "vs", group2, ". Total Cells:", length(cellNames_MAST)))
  
  # Subset genes expressed in >20% of selected cells
  filter_MAST_genes <- colSums(counts_MAST != 0) / length(cellNames_MAST) > 0.2
  counts_MAST <- counts_MAST[, filter_MAST_genes]
  geneNames_MAST <- geneNames[filter_MAST_genes]
  
  if (ncol(counts_MAST) == 0) {
    message("No genes passed the 20% expression threshold.")
    return(NULL)
  }
  
  # Construct SingleCellAssay Object with explicit factor ordering
  sca_MAST <- FromMatrix(
    t(counts_MAST), 
    cData = data.frame(
      # Force the condition to have group1 (REF) as the explicit first level
      condition = factor(pvariant_MAST, levels = c(group1, group2)), 
      wellKey = cellNames_MAST, 
      barcode = cellNames_MAST
    ), 
    fData = data.frame(
      symbolid = geneNames_MAST, 
      primerid = geneNames_MAST, 
      geneName = geneNames_MAST
    )
  )
  
  # Calculate Cellular Detection Rate
  colData(sca_MAST)$cngeneson <- scale(colSums(assay(sca_MAST) > 0))
  
  # Fit the hurdle model
  zlmResults <- zlm(~condition + cngeneson, sca_MAST)
  
  # Compare conditiongroup2 (since group2 is the test)
  lrt_contrast <- paste0("condition", group2) 
  zt <- summary(zlmResults, doLRT = lrt_contrast)$datatable
  
  fcHurdle <- merge(
    zt[contrast == lrt_contrast & component == 'H', .(primerid, `Pr(>Chisq)`)],
    zt[contrast == lrt_contrast & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
    by = 'primerid'
  )
  
  fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
  setorder(fcHurdle, fdr)
  
  # Adjust missing model values using your fixNAs rules
  fcHurdle <- fixNAs(
    fcHurdle = as.data.frame(fcHurdle),
    cluster = group2,
    clusters = pvariant_MAST,
    geneNames = geneNames_MAST,
    counts = t(counts_MAST)
  )
  
  # Append tracking descriptor and export
  fcHurdle$variant <- paste0(group2, "_vs_", group1)
  write.csv(fcHurdle, paste0('dimer_monomer_MAST_results/fcHurdle_4T1_', suffix, '.csv'), row.names = FALSE)
}

# --- 1. monomer_hi vs REF ---
run_mast_test(counts, cellNames, pvariant, "REF", "mets_hi_monomer", "mets_hi_monomer_vs_REF")

# --- 2. dimer_hi vs REF ---
run_mast_test(counts, cellNames, pvariant, "REF", "mets_hi_dimer", "mets_hi_dimer_vs_REF")

# --- 3. monomer_hi vs REF ---
run_mast_test(counts, cellNames, pvariant, "mets_lo", "mets_hi_monomer", "mets_hi_monomer_vs_mets_lo")

# --- 4. dimer_hi vs REF ---
run_mast_test(counts, cellNames, pvariant, "mets_lo", "mets_hi_dimer", "mets_hi_dimer_vs_mets_lo")

# --- 5. dimer_hi vs monomer_hi ---
run_mast_test(counts, cellNames, pvariant, "mets_hi_monomer", "mets_hi_dimer", "mets_hi_dimer_vs_mets_hi_monomer")







