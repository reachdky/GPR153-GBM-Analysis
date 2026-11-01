# DESeq2 1.52.0 differential expression, tumor vs normal, design ~ condition, apeglm 1.34.0 log2FC shrinkage
# (PREREGISTRATION.md, Section 12). Usage: Rscript pipeline/a60_deseq2.R <counts.csv.gz> <samples.csv> <output.csv>
# counts: genes x samples (first column gene_id); samples: columns sample, condition (tumor / normal).
suppressPackageStartupMessages({library(DESeq2); library(apeglm)})
a <- commandArgs(trailingOnly = TRUE)
cts <- as.matrix(read.csv(gzfile(a[1]), row.names = 1, check.names = FALSE))
cd <- read.csv(a[2])
stopifnot(identical(colnames(cts), cd$sample))
cd$condition <- factor(cd$condition, levels = c("normal", "tumor"))
dds <- DESeqDataSetFromMatrix(cts, cd, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, name = "condition_tumor_vs_normal")
shr <- lfcShrink(dds, coef = "condition_tumor_vs_normal", type = "apeglm", quiet = TRUE)
out <- data.frame(gene_id = rownames(res), baseMean = res$baseMean, log2FoldChange = res$log2FoldChange,
                  lfcSE = res$lfcSE, log2FoldChange_apeglm = shr$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
write.csv(out, a[3], row.names = FALSE)
