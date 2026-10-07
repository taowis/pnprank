# Extract raw counts and cell metadata from a Seurat object saved as RDS.
# Usage: Rscript seurat_to_mtx.R <in.rds[.gz]> <outdir> [assay]
# Writes matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz, meta.tsv.gz.
# Fails if the counts are not integer-valued (i.e. not raw counts).
suppressPackageStartupMessages({library(Matrix); library(SeuratObject)})
args <- commandArgs(trailingOnly = TRUE)
assay <- if (length(args) >= 3) args[3] else "RNA"
obj <- tryCatch(readRDS(args[1]), error = function(e) {
  con <- gzcon(gzfile(args[1], "rb")); on.exit(close(con)); readRDS(con)
})
cat("class:", class(obj), "\n")
a <- obj[[assay]]
counts <- if (inherits(a, "Assay5")) {
  layers <- grep("^counts", Layers(a), value = TRUE)
  cat("counts layers:", layers, "\n")
  if (length(layers) > 1) do.call(cbind, lapply(layers, function(l) LayerData(a, layer = l))) else LayerData(a, layer = layers)
} else {
  slot(a, "counts")
}
counts <- as(counts, "CsparseMatrix")
cat("dim:", dim(counts), " nnz:", length(counts@x), "\n")
if (length(counts@x) == 0) stop("counts slot is empty")
if (!all(counts@x == round(counts@x))) stop("counts are not integers")
md <- obj@meta.data
md$active_ident <- as.character(Idents(obj))
md <- md[colnames(counts), , drop = FALSE]
rm(obj, a); invisible(gc())  # free the full object before writing (large atlases)
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
writeMM(counts, file.path(args[2], "matrix.mtx"))
system2("gzip", c("-f", shQuote(file.path(args[2], "matrix.mtx"))))
writeLines(rownames(counts), gzfile(file.path(args[2], "features.tsv.gz")))
writeLines(colnames(counts), gzfile(file.path(args[2], "barcodes.tsv.gz")))
write.table(cbind(barcode = rownames(md), md), gzfile(file.path(args[2], "meta.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("meta columns:", colnames(md), "\n")
