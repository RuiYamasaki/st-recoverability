"""One command to re-derive every Gate-0 result and figure from committed code.

    python src/run_all.py

Steps (all seeded and deterministic; same bytes every run):
  1. realism.py : stream/cache the public FOVs, write results/realism.csv (+meta).
  2. sweep.py   : the (density x packing x sigma) oracle frontier, write results/sweep.csv.
  3. anchor.py  : oracle evaluated at the real-data anchors, write results/realism_oracle.csv.
  4. figures.py : the figures under figures/.

Raw public data is cached under data/ (gitignored). Re-running with data cached is
fast; the first run downloads ~7 MB total. Pass --skip-realism to reuse an existing
results/realism.csv without touching the network.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def main():
    skip_realism = "--skip-realism" in sys.argv
    t0 = time.time()

    if skip_realism and os.path.exists(os.path.join(RESULTS, "realism.csv")):
        print("[1/4] realism: skipped (using existing results/realism.csv)")
    else:
        print("[1/4] realism: streaming public FOVs ...")
        import realism
        realism.main()

    print("\n[2/4] sweep: oracle frontier over the grid ...")
    import sweep
    sweep.main()

    print("\n[3/4] anchor: oracle at the real-data anchors ...")
    import anchor
    anchor.main()

    print("\n[4/4] figures ...")
    import figures
    figures.main()

    print(f"\nall done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
