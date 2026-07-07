#!/usr/bin/env python3
"""
Download GEO GBM datasets
"""

import os
import sys
import GEOparse

sys.path.append(os.path.expanduser('~/gpr153_project'))
from scripts.utils import log_message


def main():
    log_message("Starting GEO downloads")

    output_dir = os.path.join(
        os.path.expanduser('~/gpr153_project'),
        'data', 'raw', 'geo'
    )
    os.makedirs(output_dir, exist_ok=True)

    geo_series = ["GSE13041", "GSE107850", "GSE4290"]

    for series in geo_series:
        print(f"\nDownloading {series}...")
        try:
            gse = GEOparse.get_GEO(series, destdir=output_dir)
            series_dir = os.path.join(output_dir, series)
            os.makedirs(series_dir, exist_ok=True)

            expr = gse.pivot_samples("VALUE")
            expr.to_csv(os.path.join(series_dir, "expression.csv"))
            gse.phenotype_data.to_csv(
                os.path.join(series_dir, "clinical.csv")
            )

            print(f"✓ {series} complete")
        except Exception as e:
            print(f"✗ {series} failed: {e}")

    log_message("GEO download complete")


if __name__ == "__main__":
    main()
