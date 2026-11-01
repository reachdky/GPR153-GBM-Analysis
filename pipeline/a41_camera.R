# CAMERA (limma 3.68.5), correlation-aware competitive gene-set test (PREREGISTRATION.md, Section 10).
# Design ~ GPR153 + batch on log2(TPM + 1) of eligible genes (GPR153 removed); test of the GPR153 coefficient.
# inter.gene.cor = NA (estimated per set) as primary; inter.gene.cor = 0.01 (limma default) as sensitivity.
# MSigDB 2026.2.Hs Hallmark and Reactome; sets of 15-500 genes present in the data; FDR per library.
# Usage: Rscript pipeline/a41_camera.R <camera_expression.csv.gz> <camera_covariates.csv> <genesets dir> <outdir>
suppressPackageStartupMessages(library(limma))
a <- commandArgs(trailingOnly = TRUE)
y <- as.matrix(read.csv(gzfile(a[1]), row.names = 1, check.names = FALSE))
cv <- read.csv(a[2])
stopifnot(identical(colnames(y), cv$aliquot))
design <- model.matrix(~ gpr153 + batch_gdc, data = cv)
read_gmt <- function(p) {
  l <- strsplit(readLines(p), "\t")
  setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1))
}
libs <- c("h.all.v2026.2.Hs.symbols", "c2.cp.reactome.v2026.2.Hs.symbols")
targeted <- list()
for (lib in libs) {
  sets <- read_gmt(file.path(a[3], paste0(lib, ".gmt")))
  idx <- ids2indices(sets, rownames(y))
  idx <- idx[lengths(idx) >= 15 & lengths(idx) <= 500]
  for (igc in c("estimated", "0.01")) {
    r <- camera(y, idx, design, contrast = 2, inter.gene.cor = if (igc == "estimated") NA else 0.01)
    r$FDR <- p.adjust(r$PValue, "BH")
    if (!"Correlation" %in% names(r)) r <- cbind(Correlation = 0.01, r)
    r <- r[, c("NGenes", "Correlation", "Direction", "PValue", "FDR")]
    r <- data.frame(term = rownames(r), r, library = lib, inter_gene_cor = igc, check.names = FALSE)
    write.csv(r, file.path(a[4], paste0("camera_", lib, "_igc_", igc, ".csv")), row.names = FALSE)
    targeted[[paste(lib, igc)]] <- r[r$term %in% c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_DNA_REPAIR"), ]
  }
}
write.csv(do.call(rbind, targeted), file.path(a[4], "camera_targeted_sets.csv"), row.names = FALSE)
