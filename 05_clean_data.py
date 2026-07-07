import os
import pandas as pd

# ==== PATHS ====
TCGA_EXPR_PATH = "data/raw/tcga/expression/tcga_gbm_expression_tpm.csv"
GTEX_EXPR_PATH = "data/raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"

OUTPUT_DIR = "data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==== LOAD TCGA EXPRESSION ====
print("Loading TCGA expression data...")
tcga_expr = pd.read_csv(TCGA_EXPR_PATH, index_col=0)
print(f"TCGA shape: {tcga_expr.shape}")

# ==== LOAD GTEX EXPRESSION ====
print("Loading GTEx expression data...")
# The GTEx gct.gz file has 2 header rows; skiprows=2
gtex_expr = pd.read_csv(GTEX_EXPR_PATH, sep="\t", skiprows=2, index_col=0)
print(f"GTEx shape: {gtex_expr.shape}")

# ==== OPTIONAL: FILTER GENES PRESENT IN BOTH ====
common_genes = tcga_expr.index.intersection(gtex_expr.index)
tcga_expr = tcga_expr.loc[common_genes]
gtex_expr = gtex_expr.loc[common_genes]

print(f"Filtered TCGA shape: {tcga_expr.shape}")
print(f"Filtered GTEx shape: {gtex_expr.shape}")

# ==== SAVE PROCESSED DATA ====
tcga_expr.to_csv(os.path.join(OUTPUT_DIR, "tcga_expression_clean.csv"))
gtex_expr.to_csv(os.path.join(OUTPUT_DIR, "gtex_expression_clean.csv"))

print("Processing complete. Cleaned files saved to:", OUTPUT_DIR)
