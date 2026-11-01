#!/usr/bin/env python3
"""
Download GTEx normal brain tissue expression data
"""

import os
import sys

sys.path.append(os.path.expanduser('~/gpr153_project'))
from scripts.utils import log_message


def download_gtex_brain():
    log_message("Starting GTEx brain data instructions")

    output_dir = os.path.join(
        os.path.expanduser('~/gpr153_project'),
        'data', 'raw', 'gtex'
    )
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print("GTEx DATA DOWNLOAD (REQUIRED)")
    print("=" * 70)
    print("GTEx data must be downloaded manually:")
    print("\n1. Go to: https://gtexportal.org/home/datasets")
    print("2. Register (free)")
    print("3. Download:")
    print("   GTEx_Analysis_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
    print("4. Gunzip the file")
    print("5. Place it in:")
    print(f"   {output_dir}/")
    print("=" * 70)

    log_message("GTEx instructions printed")


if __name__ == "__main__":
    download_gtex_brain()
