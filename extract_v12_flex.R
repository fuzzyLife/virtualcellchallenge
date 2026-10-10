# Extract one verified VIPerturb bin with base R/Matrix, without Seurat dependencies.
suppressPackageStartupMessages(library(Matrix))
Sys.setlocale("LC_ALL","C")
a <- commandArgs(trailingOnly=TRUE);stopifnot(length(a)==5)
source <- a[1];out <- a[2];panel <- read.csv(a[3],stringsAsFactors=FALSE)[[1]];bin <- a[4]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readRDS(source);meta <- attr(x,"meta.data");rna <- attr(x,"assays")[["RNA"]]
counts <- attr(rna,"layers")[["counts"]]
genes <- attr(attr(rna,"features"),"dimnames")[[1]];cells <- attr(attr(rna,"cells"),"dimnames")[[1]]
stopifnot(inherits(counts,"dgCMatrix"),length(genes)==nrow(counts),length(cells)==ncol(counts),identical(cells,rownames(meta)),!anyDuplicated(genes),!anyDuplicated(cells))
reference <- read.csv(a[5],header=FALSE,stringsAsFactors=FALSE,col.names=c("guide_id","protospacer","gene"))
control <- meta$gene == "NO-TARGET"
expected <- sub("NO-TARGET_","NO-TARGETg",reference$guide_id[reference$gene=="NO-TARGET"],fixed=TRUE)
stopifnot(length(expected)==496,all(meta$guide[control] %in% expected))
# Some shards may contain only perturbation genes; controls are shared globally
# after combining verified bins. Never substitute a human gene called CTRL.
selected <- which(meta$gene %in% panel | control)
sel <- counts[,selected,drop=FALSE];library <- Matrix::colSums(sel)
stopifnot(all(sel@x>=0),all(is.finite(sel@x)),all(sel@x==floor(sel@x)),all(library>0),length(sel@x)<2147483647)
meta <- meta[selected,,drop=FALSE];meta$cell_id <- cells[selected];meta$source_row <- selected-1L;meta$bin <- bin
meta$gene_target <- ifelse(control[selected],"Non-Targeting",meta$gene);meta$guide_id <- meta$guide
meta$total_counts <- library;meta$original_nCount_RNA <- meta$nCount_RNA
meta$ntc <- control[selected]
# RNA is a fixed 19,068-feature assay; preserve the actual complete assay totals
# before projecting to the challenge gene axis.
write.csv(data.frame(gene_name=genes),file.path(out,"gene_names.csv"),row.names=FALSE)
write.csv(meta,file.path(out,"selected_cells.csv"),row.names=FALSE)
writeLines(c(as.character(ncol(sel)),as.character(nrow(sel)),as.character(length(sel@x))),file.path(out,"shape.txt"))
con <- file(file.path(out,"counts_data.bin"),"wb");writeBin(as.integer(sel@x),con,size=4,endian="little");close(con)
con <- file(file.path(out,"counts_indices.bin"),"wb");writeBin(as.integer(sel@i),con,size=4,endian="little");close(con)
con <- file(file.path(out,"counts_indptr.bin"),"wb");writeBin(as.integer(sel@p),con,size=4,endian="little");close(con)
write.csv(as.data.frame(table(meta$gene_target)),file.path(out,"target_cell_counts.csv"),row.names=FALSE)
writeLines(c(paste("source_cells",ncol(counts)),paste("selected_cells",ncol(sel)),paste("source_genes",nrow(counts)),paste("controls",sum(meta$ntc)),paste("metadata_library_max_error",max(abs(library-meta$nCount_RNA)))),file.path(out,"extraction_summary.txt"))
cat("Extracted",bin,"cells",ncol(sel),"features",nrow(sel),"controls",sum(meta$ntc),"\n")
