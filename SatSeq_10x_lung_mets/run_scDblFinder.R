
############ Load Libraries ############
cat('Loading scDblFinder')
library("scDblFinder")
require("rhdf5")

############ Parse arguments #############
cat('Parsing Arguments')

args <- commandArgs(trailingOnly=TRUE)
cat("args=", args, "\n")

## -i or --infile
w <- which(args == "--infile" | args == "-i")
if (length(w) == 1 && w < length(args)) {
  infile <- args[w+1]
} else {
  infile <- "ERROR. Input file required."
}
cat("infile =", infile, "\n")

## -o or --outfile
w <- which(args == "--outfile" | args == "-o")
if (length(w) == 1 && w < length(args)) {
  outfile <- args[w+1]
} else {
  outfile <- "ERROR. Output file required."
}
cat("outfile =", outfile, "\n")

## -l or --layer
layer <- "layers/rawX"
w <- which(args == "--layer" | args == "-l")
if (length(w) == 1 && w < length(args)) {
  layer <- args[w+1]
} 
cat("layer =", layer, "\n")

clusters <- FALSE
w <- which(args == "--clusters" | args == "-c")
if (length(w) ==1) {
   clusters <- TRUE
}
cat("clusters=", clusters, "\n")


## -d or --dbr
w <- which(args == "--dbr" | args == "-d")
if (length(w) == 1 && w < length(args)) {
  dbr <- as.numeric(args[w+1])
} else {
  dbr <- NULL
}
cat("dbr =", dbr, "\n")

samples <- NULL
w <- which(args == "--samples" | args == "-s")
if (length(w) == 1 && w < length(args)) {
  tmp <- h5read(infile, sprintf("/obs/%s", args[w+1]))
  samples <- tmp$categories[tmp$codes+1]
}
print(head(samples))
#cat("samples =", samples, "\n")




############ Run scDblFinder ############
cat('Running scDblFinder')

# Read in count data:
counts <- h5read(infile, layer)

# Run scDblFinder:
results <- scDblFinder(counts, dbr=dbr, clusters=clusters, samples=samples) # The default doublet rate (dbr) value is NULL and the algorithm infers this parameter

# Save results:
write.csv(colData(results), outfile)

cat('scDblFinder Complete.')
				
##### End of Script #####
