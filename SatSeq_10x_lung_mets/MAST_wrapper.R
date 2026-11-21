require("MAST")
require("SingleCellExperiment", quietly=TRUE)
require("rhdf5")
require("data.table")
require("argparse")
require("parallel")
require("Matrix")
options(mc.cores=detectCores())


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

## Get infile from command line
parser <- ArgumentParser()
parser$add_argument("-i", "--infile", type="character",
                    default="/workdir/ContactTracing/tutorial/test.h5ad",
                    help="input file (h5ad).  Should have fields: .layers[\"logX\"].h5ad with log-normalized counts for each cell/gene")
parser$add_argument("-I", "--infile2", type="character", default=NULL, help="optional second input h5ad file. Will concat with infile if provided")
parser$add_argument("-o", "--outdir", type="character",
                    default="/workdir/ContactTracing/tutorial/outdir",
                    help="output directory")
parser$add_argument("-g", "--groups", type="character",
                    default="condition",
                    help="column in .obs that will be used to group cells")
parser$add_argument("-c", "--comp-groups", nargs=2)
parser$add_argument("-C", "--numcore", type="integer", default=1, help="number of cores to use")
parser$add_argument("-f", "--force", action="store_true", default=FALSE, help="Overwrite output files (default: do not re-run MAST if output already exists)")
parser$add_argument("-V", "--covariates", type="character", default="cngeneson", help="Other columns in .obs to be used as covariate (comma-delimited)")

parser$add_argument("-F", "--filter", type="character", nargs=2, help="if provided, the first argument is column name of .obs matrix, and the second is the value, only use cells with this value")

args <- parser$parse_args()

infile <- args$infile
cat("Input file = ", infile, "\n")

infile2 <- args$infile2
if (!is.null(infile2)) {
   cat("Input file 2 = ", infile2, "\n")
} else {
   cat("No second infile\n")
}
outdir <- args$outdir
cat("Output directory = ", outdir, "\n")
system(sprintf("mkdir -p %s", outdir))

groupName <- args$groups
cat("grouping by ", groupName, "\n")

compGroups <- args$comp_groups
if (is.null(compGroups)) {
    cat("Comparing all values of", groupName, "to each other\n")
} else {
    cat("Comparing", compGroups[1], "to", compGroups[2], "\n")
}

filterCells <- args$filter
if (!is.null(filterCells)) {
  cat("Only using cells with", filterCells[1], "value of ", filterCells[2], "\n")
}

cat("Using ", args$numcore, "cores\n")
options(mc.cores=args$numcore)

counts <- h5read(infile, "/layers/logX")
if (is.null(dim(counts))) {
   if (length(counts) == 3) {
      counts = t(csr_list_to_sparse_anndata(infile, "/layers/logX"))
      counts <- as.matrix(counts)
   }
}
if (is.null(counts)) {
    stop("Error: no data in /layers/logX")
}
if (!is.null(infile2)) {
   counts2 <- h5read(infile2, "/layers/logX")
   if (is.null(counts2)) {
      stop("Error: no data in /layers/logX for infile2")
   }
   if (is.null(dim(counts2))) {
     if (length(counts2) == 3) {
	counts2 <- t(csr_list_to_sparse_anndata(infile2, "/layers/logX"))
        counts2 <- as.matrix(counts2)
     }
   }
   counts <- cbind(counts, counts2)
}


## given h5ad anndata file, return given column from adata.obs
readObsFromh5File <- function(infile, colname, errorOK=FALSE) {
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
    if (errorOK) return(NULL)
    stop("Error reading ", colname, " from ", infile)
}


cellNames <- as.character(h5read(infile, "/obs/_index"))
geneNames <- as.character(h5read(infile, "/var/_index"))
clusters <- readObsFromh5File(infile, groupName)
if (!is.null(infile2)) {
   cellNames <- c(cellNames, as.character(h5read(infile2, "/obs/_index")))
   geneNames2 <- as.character(h5read(infile2, "/var/_index"))
   if ((length(geneNames) != length(geneNames2)) || (sum(geneNames != geneNames2) != 0)) {
      stop("Error: infile2 does not have same genes as infile")
   }
   clusters <- c(clusters, readObsFromh5File(infile2, groupName))
}
names(clusters) <- cellNames

keep <- rep(TRUE, length(cellNames))
if (!is.null(filterCells)) {
  filterCol <- readObsFromh5File(infile, filterCells[1])
  keep = keep & (filterCol == filterCells[2])
  cat("Discarding ", sum(!keep), " cells with ", filterCells[1], "!=", filterCells[2], "total kept=", sum(keep),"\n")
}

if (length(compGroups) == 2) {
    if (compGroups[2] == 'others') {
#        keep <- rep(TRUE, length(cellNames))
    } else {
        keep <-as.vector((keep) & (clusters == compGroups[1] | clusters == compGroups[2]))
        clusters <- clusters[keep]
        counts <- counts[,keep]
        cellNames <- cellNames[keep]
    }
    clusterList <- compGroups[1]
} else {
    clusterList <- unique(as.character(clusters))
    #keep <- rep(TRUE, length(cellNames))
}



covariates <- list()
if (!is.null(args$covariates)) {
    covariateList <- strsplit(args$covariates, ",")[[1]]
    for (v in covariateList) {
        covariates[[v]] <- readObsFromh5File(infile, v, errorOK=(v=="cngeneson"))
        if (v == "cngeneson" && is.null(covariates[[v]])) {
            cat("cngeneson not found... computing from counts\n")
            covariates[[v]] <- scale(colSums(counts > 0))
        } else if (!is.null(infile2)) {
	    covariates[[v]] <- c(covariates[[v]], readObsFromh5File(infile2, v))[keep]
	} else {
	    covariates[[v]] <- covariates[[v]][keep]
	}
        cat("Using", v, "as a covariate\n")
    }
}


fixNAs <- function(fcHurdle, cluster, clusters, geneNames, counts) {
    pcol <- grep("Chisq", names(fcHurdle), value=TRUE)
    if (length(pcol) != 1) {
        cat("Error finding p-value column in fixNAs")
        return(fcHurdle)
    }
    fixNA <- (is.na(fcHurdle$coef) & fcHurdle[,pcol] < 1)
                                        #    coefRange <- range(fcHurdle$coef, na.rm=TRUE)
    coefRange <- c(-10000, 10000)
    geneOrder <- sapply(fcHurdle$primerid, function(x) {which(geneNames ==x)})
    fracIn <- apply(counts[geneOrder,clusters==cluster, drop=FALSE], 1, sum)/sum(clusters==cluster)
    fracOut <- apply(counts[geneOrder,clusters != cluster, drop=FALSE], 1, sum)/sum(clusters!=cluster)
    fcHurdle[fixNA & fracIn <= fracOut,"coef"] <- coefRange[1]
    fcHurdle[fixNA & fracIn > fracOut, "coef"] <- coefRange[2]
    fcHurdle
}


writeResults <- function(outdir, testname, zt, zlr=NULL) {
  cons <- unique(zt$contrast)
  if (!is.null(zlr)) {
    zlr0 <- zlr[,"hurdle","Pr(>Chisq)"]
    df <- data.frame(primerid=names(zlr0), pval=as.numeric(zlr0))
  } else {
    pcontrast <- unique(zt[!is.na(zt$`Pr(>Chisq)`),'contrast']$contrast)
    if (length(pcontrast) != 1) {
      stop("Multiple contrasts have p-values, need to pass zlr to writeResults")
    }
    f <- zt$contrast == pcontrast & zt$component=='H'
    df <- data.frame(primerid=zt[f,]$primerid, pval=zt[f,]$`Pr(>Chisq)`)
  }
  for (con in cons) {
    if (con == '(Intercept)') next
    currdf <- zt[contrast == con & component=='logFC', .(primerid, coef, ci.hi, ci.lo)]
    names(currdf)[2:4] <- sprintf("%s_%s", names(currdf)[2:4], con)
    df <- merge(df, currdf, by='primerid')
  }
  df <- df[order(df$pval),]
  if (!is.null(outdir)) {
    write.table(df, file=sprintf("%s/%s.txt", outdir, testname), quote=FALSE, sep="\t", row.names=FALSE)
  }
  invisible(df)
}



for (i in 1:length(clusterList)) {
    if (class(clusterList) == "list") {
        clusters <- clusterList[[i]]
        cluster <- names(clusterList)[i]
    } else {
        cluster <- clusterList[i]
    }
    clusterName <- gsub("/","_", cluster, fixed=TRUE)
    clusterName2 <- sprintf("cluster%s", cluster)
    clusterAssign <- ifelse(clusters == cluster, cluster, "background")
                                        # get outfile name
    if (length(compGroups) == 2) {
        outbn <- sprintf("%s_vs_%s", compGroups[1], compGroups[2])
    } else {
	outbn <- clusterName
    }
    if (!is.null(filterCells)) {
	    outbn <- sprintf("%s_%s", outbn, filterCells[2])
    }
    outbn <- gsub(" ", "_", outbn, fixed=TRUE)
    outfile <- sprintf("%s/%s.csv", outdir, outbn)
    cat("outfile = \"", outfile, "\"\n", sep="")
    
    if ((!args$force) && file.exists(outfile)) {
        cat("outfile already exists, skipping\n")
        next
    }
    
    cat("Running MAST on ", clusterName, "\n")
    sca <- FromMatrix(counts, 
                      cData=data.frame(wellKey=cellNames, barcode=cellNames, cluster=clusterAssign), 
                      fData=data.frame(primerid=geneNames, geneName=geneNames))
    tmpc <- factor(colData(sca)$cluster)
    tmpc <- relevel(tmpc, "background")
    colData(sca)$cluster <- tmpc
                                        #  colData(sca)$cngeneson <- scale(colSums(assay(sca) > 0))
    for (v in names(covariates)) {
        colData(sca)[[v]] <- covariates[[v]]
    }


    zlmstr <- ("~cluster")
    for (v in names(covariates)) {
        zlmstr <- sprintf("%s + %s", zlmstr, v)
    }
    cat("zlmstr = ", zlmstr, "\n")
    system.time(zlmResults <- zlm(formula(zlmstr), sca))
    system.time(zt <- summary(zlmResults, doLRT=clusterName2)$datatable)
    
    fcHurdle <- merge(zt[contrast==clusterName2 & component=='H',.(primerid, `Pr(>Chisq)`)],
                      zt[contrast==clusterName2 & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
    
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    setorder(fcHurdle, fdr)
    fcHurdle <- fixNAs(as.data.frame(fcHurdle), cluster, clusters, geneNames, counts)
    writeResults(outdir, outbn, zt)
    write.csv(fcHurdle, outfile, quote=FALSE)
}



