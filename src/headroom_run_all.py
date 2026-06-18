"""Re-derive the segmenter-headroom analysis:

    python src/headroom_run_all.py

  1. headroom_export.py  -> data/headroom/*, results/headroom_oracle_naive.csv
  2. headroom_methods.py -> results/headroom_methods.csv  (runs pciSeq + ComSeg)
  3. headroom_figures.py -> figures/headroom_*.png
Gates 0 to 3 are not touched. Real methods use documented defaults; failures are recorded.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    t0 = time.time()
    print("[1/3] export known-truth synthetic data + oracle/naive ceiling ...")
    import headroom_export; headroom_export.main()
    print("\n[2/3] run real published methods (pciSeq, ComSeg) ...")
    import headroom_methods; headroom_methods.main()
    print("\n[3/3] figures ...")
    import headroom_figures; headroom_figures.main()
    print(f"\nall done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
