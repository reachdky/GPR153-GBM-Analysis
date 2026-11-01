# T7: primary screen re-run in R: rank, lm residuals on batch, cor (PREREGISTRATION.md, Section 13.1).
# Usage: Rscript tests/t7_screen.R <expression.csv.gz> <batch.csv> <output.csv>
a <- commandArgs(trailingOnly = TRUE)
x <- as.matrix(read.csv(gzfile(a[1]), row.names = 1, check.names = FALSE))
b <- read.csv(a[2])
stopifnot(identical(colnames(x), b$aliquot))
r <- t(apply(x, 1, rank))
res <- t(apply(r, 1, function(v) residuals(lm(v ~ b$batch_gdc))))
m <- res["MGMT_ROW", ]
g <- res[rownames(res) != "MGMT_ROW", , drop = FALSE]
rho <- as.vector(cor(t(g), m))
write.csv(data.frame(gene_id = rownames(g), rho_r = rho), a[3], row.names = FALSE)
