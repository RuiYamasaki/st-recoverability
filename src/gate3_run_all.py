"""One command to re-derive every Gate 3 result and figure.

    python src/gate3_run_all.py

Downloads the breast and lung Xenium cell-feature matrices if absent (cached under data/,
gitignored), then runs:
  1. gate3_validity.py  -> results/gate3_validity_breast.csv, gate3_pin_breast.csv
  2. gate3_typing.py    -> results/gate3_typing.csv
  3. gate3_replicate.py -> results/gate3_validity_lung.csv, gate3_pin_lung.csv
  4. gate3_figures.py   -> figures/gate3_*.png
Gates 0, 1 and 2 are not touched.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    t0 = time.time()
    print("[1/4] gate3 marker-validity test + admissible-set pin (breast) ...")
    import gate3_validity; gate3_validity.main()
    print("\n[2/4] gate3 sigma robustness to cell-typing ...")
    import gate3_typing; gate3_typing.main()
    print("\n[3/4] gate3 independent Xenium replication (lung) ...")
    import gate3_replicate; gate3_replicate.main()
    print("\n[4/4] gate3 figures ...")
    import gate3_figures; gate3_figures.main()
    print(f"\nall done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
