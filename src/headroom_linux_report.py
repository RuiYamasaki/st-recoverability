"""Combine the within-Linux oracle/naive ceiling and the field-standard method results
(Baysor, Proseg; optionally pciSeq, ComSeg) into the headroom table, the best-method-to-
oracle gap, and the beats-naive checks, in both nuclei-prior and free-segmentation modes.

Reports raw numbers only; asserts no conclusion about near-the-wall vs headroom.
Regen: micromamba run -n st python src/headroom_linux_report.py
Source files: results/headroom_linux_oracle_naive.csv, results/headroom_linux_methods.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")

METHOD_MODES = [
    ("Baysor", "nuclei_prior"), ("Baysor", "free_segmentation"),
    ("Proseg", "nuclei_prior"), ("Proseg", "free_segmentation"),
    ("pciSeq", "nuclei_prior"),
    ("ComSeg", "nuclei_prior"), ("ComSeg", "free_segmentation"),
]


def main():
    on = pd.read_csv(os.path.join(RESULTS, "headroom_linux_oracle_naive.csv")).set_index("config")
    mp = os.path.join(RESULTS, "headroom_linux_methods.csv")
    methods = pd.read_csv(mp) if os.path.exists(mp) else pd.DataFrame()

    def cell(cfg, method, mode):
        if len(methods) == 0:
            return (np.nan, np.nan, "missing")
        r = methods[(methods.config == cfg) & (methods.method == method) & (methods["mode"] == mode)]
        if not len(r):
            return (np.nan, np.nan, "missing")
        st = r.status.iloc[0]
        if st != "ok":
            return (np.nan, np.nan, st)
        return (float(r.matched_accuracy.iloc[0]), float(r.frac_assigned.iloc[0]), st)

    print("=" * 120)
    print("HEADROOM (within Linux): field-standard methods vs the oracle ceiling")
    print("matched transcript-assignment accuracy vs known truth; frac = fraction of interior transcripts assigned a cell")
    print("source: results/headroom_linux_oracle_naive.csv (oracle_acc_matched, naive_acc_matched);")
    print("        results/headroom_linux_methods.csv (matched_accuracy, frac_assigned)")
    print("-" * 120)

    def f(v):
        return f"{v:5.3f}" if v == v else f"{'n/a':>5s}"

    for cfg in on.index:
        orc = float(on.loc[cfg, "oracle_acc_matched"])
        nai = float(on.loc[cfg, "naive_acc_matched"])
        print(f"\n[{cfg}]  oracle={f(orc)}  naive={f(nai)}")
        soph = {}
        for method, mode in METHOD_MODES:
            acc, frac, st = cell(cfg, method, mode)
            tag = "" if st == "ok" else f"  ({st})"
            print(f"    {method:8s} {mode:17s}: acc={f(acc)}  frac_assigned={f(frac)}"
                  f"{'  beats_naive=' + ('Y' if (acc==acc and acc>nai) else 'N') if acc==acc else ''}{tag}")
            if st == "ok":
                soph[(method, mode)] = acc
        # gaps
        if soph:
            best_m, best_v = max(soph.items(), key=lambda kv: kv[1])
            print(f"    -> best method = {best_m[0]} ({best_m[1]}) acc={f(best_v)}; "
                  f"oracle-best={f(orc-best_v)}; naive-best={f(nai-best_v)}; "
                  f"oracle-naive={f(orc-nai)}")

    # dense sigma-CI gap summary
    print("\n" + "-" * 120)
    print("Best-method-to-oracle gap across the dense sigma CI (1.43, 1.99, 2.63 um), per mode:")
    dense = ["dense_sigma1.43", "dense_sigma1.99", "dense_sigma2.63"]
    for mode in ["nuclei_prior", "free_segmentation"]:
        line = []
        for cfg in dense:
            orc = float(on.loc[cfg, "oracle_acc_matched"])
            vals = []
            for method in ["Baysor", "Proseg", "pciSeq", "ComSeg"]:
                acc, _, st = cell(cfg, method, mode)
                if st == "ok":
                    vals.append(acc)
            if vals:
                line.append(f"{orc - max(vals):.3f}")
            else:
                line.append("n/a")
        print(f"  {mode:17s}: oracle-best = {' / '.join(line)}")

    if len(methods):
        print("\nrun status (source: results/headroom_linux_methods.csv):")
        for _, r in methods.iterrows():
            err = f"  err={r.error}" if isinstance(r.error, str) and r.error else ""
            print(f"  {r.method:8s} [{r['mode']:17s}] {r.config:26s} "
                  f"acc={f(r.matched_accuracy)} frac={f(r.frac_assigned)} "
                  f"{r.runtime_s:6.0f}s status={r.status}{err}")


if __name__ == "__main__":
    main()
