"""Combine the oracle/naive ceiling and the real-method results into the headroom table
and the best-real-method-to-oracle gap. Regen: python src/headroom_report.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def main():
    on = pd.read_csv(os.path.join(RESULTS, "headroom_oracle_naive.csv")).set_index("config")
    mp = os.path.join(RESULTS, "headroom_methods.csv")
    methods = pd.read_csv(mp) if os.path.exists(mp) else pd.DataFrame()

    print("=" * 104)
    print("HEADROOM: real methods vs the oracle ceiling (matched transcript-assignment accuracy vs known truth)")
    print("source: results/headroom_oracle_naive.csv (oracle_acc_matched, naive_acc_matched); "
          "results/headroom_methods.csv (matched_accuracy)")
    print("-" * 104)
    hdr = f"{'config':26s} {'oracle':>7s} {'naive':>7s} {'pciSeq':>7s} {'ComSeg':>7s} {'ComSeg_free':>11s} " \
          f"{'best_real':>9s} {'oracle-best':>11s}"
    print(hdr)
    for cfg in on.index:
        orc = on.loc[cfg, "oracle_acc_matched"]
        nai = on.loc[cfg, "naive_acc_matched"]

        def mval(method, mode):
            if len(methods) == 0:
                return np.nan
            r = methods[(methods.config == cfg) & (methods.method == method) & (methods["mode"] == mode)]
            return float(r.matched_accuracy.iloc[0]) if len(r) and r.status.iloc[0] == "ok" else np.nan

        pci = mval("pciSeq", "nuclei_prior")
        com = mval("ComSeg", "nuclei_prior")
        comf = mval("ComSeg", "free_segmentation")
        # best real method = best of the sophisticated published methods AND the naive baseline
        cands = {"naive": nai, "pciSeq": pci, "ComSeg": com}
        cands = {k: v for k, v in cands.items() if v == v}
        best = max(cands.values()) if cands else np.nan
        gap = orc - best if best == best else np.nan
        def f(v):
            return f"{v:7.3f}" if v == v else f"{'n/a':>7s}"
        print(f"{cfg:26s} {f(orc)} {f(nai)} {f(pci)} {f(com)} {('%11.3f'%comf) if comf==comf else ('%11s'%'n/a')} "
              f"{f(best)} {('%11.3f'%gap) if gap==gap else ('%11s'%'n/a')}")

    if len(methods):
        print("\nmethod run status (source: results/headroom_methods.csv):")
        for _, r in methods.iterrows():
            print(f"  {r.method:8s} [{r['mode']:17s}] {r.config:26s} status={r.status} "
                  f"{'err='+r.error if isinstance(r.error,str) and r.error else ''}")


if __name__ == "__main__":
    main()
