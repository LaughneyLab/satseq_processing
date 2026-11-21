options(stringsAsFactors=FALSE)
require("SingleCellExperiment", quietly=TRUE)
require("cellassign")
require("rhdf5")
args <- commandArgs(trailingOnly=TRUE)
cat("args=", args, "\n")

## given count matrix and markerMat, call cellAssign
## - remove genes in markerMat not in count matrix
## - remove genes from both with no counts
## - compute size factors before removin genes
### ... arguments get passed to cell assign
### returns cellAssign fit object
runCellAssign <- function(mat, markerMat, outroot, ...) {
    ## get cell sizes before removing any genes
    sizes <- apply(mat, 2, sum)
    normSizes <- sizes/sum(sizes)
    rownames(mat) <- tolower(rownames(mat))
    rownames(markerMat) <- tolower(rownames(markerMat))
    geneNames <- rownames(mat)
    f <- is.element(rownames(markerMat), geneNames)
    if (sum(!f) > 0) {
        cat(sum(!f), " genes out of ", length(f), "not found in data matrix.. removing\n")
        print(rownames(markerMat)[!f])
        cat("\n")
        cat(rownames(markerMat)[!f], sep="\n", file=sprintf("%s.gene_not_found.txt", outroot))
        markerMat <- markerMat[f,]
    }
    
    # reduce data count matrix to marker genes (and reorder to match markerMat)
    smallMat <- mat[rownames(markerMat),]
    
    ## also get rid of genes that have no counts in the data
    f <- apply(smallMat, 1, sum) > 0
    cat("Removing ", length(f) - sum(f), " genes with no counts\n")
    smallMat <- smallMat[f,]
    markerMat <- markerMat[f,]

    data <- SingleCellExperiment(assays = List(counts = smallMat))
    sizeFactors(data) <- normSizes

    # This command takes awhile (10 minutes for 3000 cells/40 genes/10 cell types... using ~40 CPUs, very little memory)
    cat("Calling Cell Assign with ", nrow(smallMat),
        " genes and ", ncol(smallMat), " cells\n")
    system.time(fit <- cellassign(exprs_obj=data, marker_gene_info=markerMat, s=sizeFactors(data), ...))
    write.table(data.frame(barcode=colnames(data), celltype=celltypes(fit)),
                row.names=FALSE, quote=FALSE, sep="\t", file=sprintf("%s_assignments.txt", outroot))
    write.table(cellprobs(fit), row.names=FALSE, quote=FALSE, sep="\t", file=sprintf("%s_probs.txt", outroot))
    save(fit, smallMat, markerMat, file=sprintf("%s.Rdata", outroot))
    return(fit)
}


############ Parse arguments #############

## -i or --infile
w <- which(args == "--infile" | args == "-i")
if (length(w) == 1 && w < length(args)) {
  infile <- args[w+1]
} else {
  infile <- "/workdir/CIN/4T1_LUC911_pipeline_out/tmp_results/adata_cell_filtered_doublets.h5ad"
}
cat("infile =", infile, "\n")

## -o or --outroot <outroot>
w <- which(args == "--outroot" | args == "-o")
if (length(w) == 1 && w < length(args)) {
  outroot <- args[w+1]
} else {
  outroot <- sprintf("/home/ruser/%s_cellassign", "4T1_LUC911")
  outroot <- sprintf("/home/ruser/%s_cellassign", "4T1_CGAMP_PROCESSING")
}
system(sprintf("mkdir -p %s", dirname(outroot)))
cat("infile =", infile, "\n")
cat("outroot = ", outroot, "\n")

w <- which(args == "--seed")
if (length(w) == 1 && w < length(args)) {
    cat("Setting seed to ", as.numeric(args[w+1]), "\n")
    set.seed(as.numeric(args[w+1]))
} else {
  set.seed(100)
}

seed2 <- 100
w <- which(args == "--seed2")
if (length(w) == 1 && w < length(args)) {
   seed2 <- as.numeric(args[w+1])
} 
w <- which(args == "--sample-cells")
if (length(w) ==1 && w < length(args)-1) {
   random.sample = as.numeric(args[w+2])
   random.sample.num = as.numeric(args[w+1])
   cat("Dividing cells into group of ", random.sample, ", keeping sample ", random.sample.num, "\n")
} else {
   random.sample = 0
}

maxCellTypes <- NA
w <- which(args == "--maxCellTypes")
if (length(w) == 1 && w < length(args)) {
  maxCellTypes <- as.numeric(args[w+1])
  cat("Using at most ", maxCellTypes, " cell types", "\n")
}


## -c or --clusters: combine all cells in the same cluster
## (infile should have /obs/cluster defined)
groupCluster <-  length(which(args == "--cluster" | args == "-c") == 1)

## read in data
h5ls(infile)
geneNames <- trimws(tolower(h5read(infile, "/var/_index")))
ngene <- length(geneNames)
geneIds <- h5read(infile, "/var/gene_ids")
mat <- h5read(infile, "/layers/X")
lsize <- as.numeric(h5read(infile, "/obs/library_size"))
cellNames <- h5read(infile,"/obs/_index")
dimnames(mat) <- list(geneNames, cellNames)

if (random.sample > 1) {
     sampleNum <- sample(1:random.sample, length(cellNames), replace=TRUE)
     keepCells <- (sampleNum == random.sample.num)
     mat <- mat[,keepCells]
     cellNames <- cellNames[keepCells]
     lsize <- lsize[keepCells]
     cat("Keeping ", length(cellNames), " out of ", length(keepCells), " cells\n")
}

mat <- t(t(mat)*lsize)/10000.0

if (groupCluster) {
    categories <- NULL
    try (
        categories <- as.numeric(h5read(infile, "/obs/cluster"))
    )
    if (is.null(categories)) {
        stop()
    }
    cats <- sort(unique(categories))
    cats <- cats[cats!=-1]
    tinyMat <- matrix(nrow=nrow(mat), ncol=length(cats))
    rownames(tinyMat) <- rownames(mat)
    for (cat in cats) {
        tinyMat[,cat+1] <- rowSums(mat[,categories==cat])
    }
    mat <- tinyMat
    colnames(mat) <- cats
}

## --markers or -m to specify marker set
## depending on set, different file format expected
## for some sets, --minScore <column> <score> or --maxScore <column> <score>
## supported
w <- which(args=="--markers" | args == "-m")
if (length(w) == 1 && w < length(args)) {
     markerFile = args[w+1]
     cat("markerFile = ", markerFile, "\n")
     if (grepl("mcaMarkers_.*_DI.*", basename(markerFile))) {
         tmp <- read.table(markerFile, header=FALSE, sep=",")
         names(tmp) <- c("CellType", "genes")
         markers <- list()
         for (i in 1:nrow(tmp)) {
             markers[[tmp$CellType[i]]] <- strsplit(gsub(' ', '', gsub('\'', '', gsub('\\}', '', gsub('\\{', '', tmp$genes[i], perl=FALSE), perl=FALSE), perl=FALSE), perl=FALSE), split=',', fixed=TRUE)[[1]]
         }
     } else if (grepl("JK.*.csv", basename(markerFile)) || basename(markerFile)=="markers_filtered_CLEANED_ete4001_2021_06_24.csv" || basename(markerFile)=="markers_2021_10_14.csv"){
         m0 <- read.table(markerFile, header=TRUE, sep=",", check.names=FALSE)
         markers <- list()
         for (i in 1:ncol(m0)) {
             if (class(m0[,i]) == "integer") next
             f <- !(is.na(m0[,i]) | m0[,i] == "")
             if (sum(f) == 0) next
             celltype <- names(m0)[i]
             markers[[celltype]] <- unique(m0[f,i])
             cat("Added ", length(markers[[celltype]]),
                 " markers for cell type ", celltype, "\n")
         }
     }  else if (grepl("markers_filtered_CLEANED.*.csv", basename(markerFile)) || basename(markerFile) == "markers_2021_11_18.csv") {
       m0 <- read.table(markerFile, header=TRUE, sep=",", check.names=FALSE)
       genes <- tolower(m0[,1])
       markers <- list()
       for (i in 2:ncol(m0)) {
         markers[[colnames(m0)[i]]] <- genes[m0[,i] > 0]
       }
     } else {
         markers <- read.table(markerFile, header=TRUE, sep="\t")
         markers$gene <- tolower(markers$gene)
         cat("Read ", markerFile, " with ", nrow(markers), " markers\n")
         w <- which(args == "--minScore")
         if (length(w) == 1 && w < length(args)-1) {
             cat("Subsetting by minScore\n")
             markers <- markers[markers[,args[w+1]] >= as.numeric(args[w+2]),]
         }
         w <- which(args == "--maxScore")
         if (length(w) == 1 && w < length(args)-1) {
             cat("Subsetting by max Score\n")
             markers <- markers[markers[,args[w+1]] <= as.numeric(args[w+2]),]
         }
	 if (is.element("avg_logFC", names(markers))) {
           markers <- markers[order(markers$avg_logFC, decreasing=TRUE),]
	 }
         print(names(markers))
         cat("\n")
         origMarkers <- markers
         cellTypes <- unique(origMarkers[,"Cell.Type"])
         markers <- list()
         for (ct in cellTypes) {
             f <- origMarkers[,"Cell.Type"] == ct
             if (sum(f) > 0) {
                 markers[[ct]] <- unique(origMarkers[f,"gene"])
             }
             cat(sum(f), " markers for ", ct, "\n")
         }
     }
} else if (FALSE) {
    ## This is an attempt to use CellMarker, it is way too much though, and messy. Try something simpler for now
    markers <- read.table("/home/mt269/Mouse_cell_markers.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
## first column is always "mouse"
    markerGenesLower <- unique(sort(tolower(sub("\\]", "", sub("\\[", "", trimws(unlist(strsplit(markers$geneSymbol, split=','))))))))
    markerGenes <- unique(sort(sub("\\]", "", sub("\\[", "", trimws(unlist(strsplit(markers$geneSymbol, split=',')))))))
} else {
     cat("Using default mmLung markers from Garnett\n")
     markers <- list()
     markers[["Endothelial"]] <- c("Lyve1", "Tek", "Kdr", "Ramp2", "Flt1")
     markers[["B cells"]] <- c("Cd19", "Ms4a1", "Cd79a")
     markers[["T cells"]] <- c("Cd3d", "Cd3e", "Cd3g")
     markers[["Fibroblasts"]] <- c("Des", "Pdgfra", "Pdgfrb", "Fap", "Thy1")
     markers[["Alveolar"]] <- c("Sftpa1", "Sftpb", "Ager", "Aqp4")
     markers[["NK"]] <- c("Ncam1", "Cd160", "Itga2", "Klrb1a")
     markers[["Granulocytes"]] <- c("Psg26", "Fut4", "Siglece", "Itgam", "Ccr3", "Il5ra")
     markers[["Ciliated cells"]] <- c("Cyp2f2", "Ccdc153")
     markers[["Macs/Monos/DCs"]] <- c("Cd14", "Cd68", "Cd33", "Tlr7", "Tlr9", "Itgax", "Flt3")
}

w <- which(args == "--max-markers-per-ct")
if (length(w) == 1 && w < length(args)) {
    maxMarkers = as.numeric(args[w+1])
    for (i in 1:length(markers)) {
        markers[[i]] <- markers[[i]][is.element(tolower(markers[[i]]), tolower(rownames(mat)))]
        if (length(markers[[i]]) > maxMarkers) {
            cat("Reducing to ", maxMarkers, " from ", length(markers[[i]]), " for cell type ", names(markers)[i], "\n")
            markers[[i]] <- markers[[i]][1:maxMarkers]
        }
    }
}

w <- which(args == "--subset-markers")
if (length(w) == 1 && w < length(args)) {
    subsetMarkers = as.numeric(args[w+1])
    cat("Subsetting to ", subsetMarkers, " random markers\n")
    for (i in 1:length(markers)) {
        if (length(markers[[i]]) > subsetMarkers) {
            w <- sample(length(markers[[i]]), subsetMarkers, replace=FALSE)
            markers[[i]] <- markers[[i]][w]
         }
     }
    cat("marker length: ", sapply(markers, length), "\n")
}

if (!is.na(maxCellTypes)) {
  if (length(markers) > maxCellTypes) {
     markers <- markers[sample(1:length(markers), maxCellTypes, replace=FALSE)]
  }
cat("Number of cell Types: ", length(markers), "\n")
}

markerMat <- marker_list_to_mat(markers, include_other=FALSE)
rownames(markerMat) <- trimws(rownames(markerMat))

if (!is.null(seed2)) {
   cat("Setting seed2 to ", seed2, "\n")
   set.seed(seed2)
}
system.time(fit <- runCellAssign(mat, markerMat, outroot=outroot, shrinkage=TRUE, verbose=TRUE))

## for each cell type, compute fraction of cells assigned to type that have
## (TODO)


