"""One command to re-derive every Gate 2 result and figure.

    python src/gate2_run_all.py

Downloads the Xenium cell-feature matrix if absent (cached under data/, gitignored),
then runs:
  1. gate2_pin.py         -> results/gate2_xenium_pin.csv
  2. gate2_uncertainty.py -> results/gate2_sigma_uncertainty.csv, results/gate2_design_grid.csv
  3. gate2_nbinom.py      -> results/gate2_nbinom.csv, results/gate2_nbinom_frontier.csv
  4. gate2_figures.py     -> figures/gate2_*.png
Deterministic given the cached data. Gate 0 and Gate 1 are not touched.
"""
from __future__ import annotations

import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

XENIUM_MATRIX_URL = (
    "https://cf.10xgenomics.com/samples/xenium/1.0.1/"
    "Xenium_FFPE_Human_Breast_Cancer_Rep1/"
    "Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5"
)


def _ensure(url, name):
    dest = os.path.join(DATA, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  cached: {name} ({os.path.getsize(dest):,} bytes)")
        return
    os.makedirs(DATA, exist_ok=True)
    print(f"  downloading {name} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "gate2/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as f:
        f.write(r.read())
    print(f"  done: {os.path.getsize(dest):,} bytes")


def main():
    t0 = time.time()
    print("[data] Xenium cell-feature matrix + cells.parquet")
    _ensure(XENIUM_MATRIX_URL, "xenium_breast_rep1_cell_feature_matrix.h5")
    # cells.parquet is fetched by src/realism.py in Gate 0; ensure present too
    from realism import XENIUM_CELLS_URL, _download
    _download(XENIUM_CELLS_URL, os.path.join(DATA, "xenium_breast_rep1_cells.parquet"))

    print("\n[1/4] gate2 pin sigma on Xenium ...")
    import gate2_pin; gate2_pin.main()
    print("\n[2/4] gate2 sigma uncertainty (CI + design) ...")
    import gate2_uncertainty; gate2_uncertainty.main()
    print("\n[3/4] gate2 non-idealized (negative-binomial) emission ...")
    import gate2_nbinom; gate2_nbinom.main()
    print("\n[4/4] gate2 figures ...")
    import gate2_figures; gate2_figures.main()
    print(f"\nall done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
