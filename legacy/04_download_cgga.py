#!/usr/bin/env python3
"""
CGGA download instructions
"""

import os
import sys

sys.path.append(os.path.expanduser('~/gpr153_project'))
from scripts.utils import log_message


def main():
    log_message("Starting CGGA instructions")

    output_dir = os.path.join(
        os.path.expanduser('~/gpr153_project'),
        'data', 'raw', 'cgga'
    )
    os.makedirs(output_dir, exist_ok=True)

    readme = os.path.join(output_dir, "README.txt")
    with open(readme, "w") as f:
        f.write(
            "CGGA DATA (REQUIRED)\n"
            "1. Register at http://www.cgga.org.cn/\n"
            "2. Download mRNAseq_693\n"
            "3. Download clinical data\n"
            "4. Place files here\n"
        )

    print(f"CGGA instructions saved → {readme}")
    log_message("CGGA instructions complete")


if __name__ == "__main__":
    main()
