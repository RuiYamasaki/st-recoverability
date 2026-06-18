"""One command to re-derive every Gate 1 result and figure.

    python src/gate1_run_all.py

Reuses the cached real data under data/ (download with src/realism.py if absent) and
the committed results/realism.csv from Gate 0. Deterministic: same bytes every run.
Steps:
  1. gate1_sweep.py           -> results/gate1_sweep.csv (+ gate1_expression_meta.json)
  2. gate1_leakage_anchor.py  -> results/gate1_leakage.csv
  3. gate1_structural.py      -> results/gate1_structural.csv
  4. gate1_figures.py         -> figures/gate1_*.png
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    t0 = time.time()
    print("[1/4] gate1 sweep (realistic expression) ...")
    import gate1_sweep
    gate1_sweep.main()
    print("\n[2/4] gate1 leakage + data-pinned sigma ...")
    import gate1_leakage_anchor
    gate1_leakage_anchor.main()
    print("\n[3/4] gate1 structural sensitivity ...")
    import gate1_structural
    gate1_structural.main()
    print("\n[4/4] gate1 figures ...")
    import gate1_figures
    gate1_figures.main()
    print(f"\nall done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
