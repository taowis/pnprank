# Convert an RDS holding a sparse count matrix (genes x cells) to Matrix Market.
# Usage: Rscript rds_to_mtx.R <in.RDS[.gz]> <outdir>
# Writes matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz and prints the object class.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
# Some GEO uploads gzip an already-compressed RDS; fall back to stripping one layer.
obj <- tryCatch(readRDS(args[1]), error = function(e) {
  message("readRDS failed (", conditionMessage(e), "); retrying with one gzip layer removed")
  con <- gzcon(gzfile(args[1], "rb"))
  on.exit(close(con))
  readRDS(con)
})
cat("class:", class(obj), "\n")
if (!inherits(obj, "dgCMatrix")) obj <- as(obj, "CsparseMatrix")
cat("dim:", dim(obj), "\n")
cat("integer-valued:", all(obj@x == round(obj@x)), " min:", min(obj@x), " max:", max(obj@x), "\n")
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
writeMM(obj, file.path(args[2], "matrix.mtx"))
system2("gzip", c("-f", shQuote(file.path(args[2], "matrix.mtx"))))
writeLines(rownames(obj), gzfile(file.path(args[2], "features.tsv.gz")))
writeLines(colnames(obj), gzfile(file.path(args[2], "barcodes.tsv.gz")))
