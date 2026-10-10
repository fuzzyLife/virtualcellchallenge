# Inspect serialized Seurat structure without modifying inputs or guessing assay axes.
suppressPackageStartupMessages(library(Matrix))
Sys.setlocale("LC_ALL", "C")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
x <- readRDS(args[1])
sink(args[2])
cat("Object class:",class(x),"\nAttributes:",names(attributes(x)),"\n")
meta <- attr(x,"meta.data")
cat("Metadata dimensions:",dim(meta),"\nColumns:",names(meta),"\n")
for (nm in names(meta)) {
 cat("\n",nm,"class",class(meta[[nm]]),"unique",length(unique(meta[[nm]])),"values:")
 print(utils::head(unique(meta[[nm]]),12))
}
assays <- attr(x,"assays")
for(nm in names(assays)) {
 a <- assays[[nm]];cat("\nAssay",nm,"class",class(a),"attributes",names(attributes(a)),"\n")
 if (!is.null(attr(a,"counts"))) cat("counts dimensions:",dim(attr(a,"counts")),"\n")
 layers <- attr(a,"layers")
 if (!is.null(layers)) for (layer in names(layers)) cat("Layer",layer,"dimensions",dim(layers[[layer]]),"rownames?",!is.null(rownames(layers[[layer]])),"colnames?",!is.null(colnames(layers[[layer]])),"\n")
}
sink()
