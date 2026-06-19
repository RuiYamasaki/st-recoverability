"""Reproducibility + over-segmentation diagnostic at the dense operating point.

Baysor and Proseg expose no CLI seed, so their output is not bit-for-bit reproducible. This
script repeats each method/mode N_REP times on the dense, data-pinned-sigma config
(dense_sigma1.99) and records matched accuracy, fraction assigned, and the number of distinct
method-cells over interior transcripts (n_method_cells). The latter shows over- vs
under-segmentation, which the one-to-one matched-accuracy metric penalises (256 true interior
cells exist). Reuses headroom_common.matched_accuracy and the method drivers unchanged.

  micromamba run -n st python src/headroom_linux_repro.py
Writes results/headroom_linux_repro.csv.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import headroom_common as hc  # noqa: E402

DATA_LINUX = os.path.join(hc.DATA, "headroom_linux")
CONFIG = "dense_sigma1.99"
N_REP = 3
RUNS = [("Baysor", "nuclei_prior", True), ("Baysor", "free_segmentation", False),
        ("Proseg", "nuclei_prior", True)]


def main():
    from gate2_pin import build_xenium
    model, _, _ = build_xenium()
    from methods_baysor import run_baysor
    from methods_proseg import run_proseg
    fns = {"Baysor": run_baysor, "Proseg": run_proseg}

    d = os.path.join(DATA_LINUX, CONFIG)
    t = pd.read_csv(os.path.join(d, "transcripts.csv"))
    c = pd.read_csv(os.path.join(d, "cells.csv"))
    inter = t["interior"].to_numpy().astype(bool)
    true = t["true_cell"].to_numpy()
    n_true = int(len(np.unique(true[inter])))
    print(f"[{CONFIG}] interior true cells = {n_true}")

    rows = []
    for method, mode, wp in RUNS:
        for rep in range(N_REP):
            t0 = time.time()
            try:
                assigned, _ = fns[method](t, c, model, with_prior=wp)
                acc = hc.matched_accuracy(assigned, true, inter)
                a = assigned[inter]
                frac = float((a >= 0).mean())
                n_cells = int(len(np.unique(a[a >= 0])))
                st, err = "ok", ""
            except Exception as e:  # noqa: BLE001
                acc, frac, n_cells, st, err = np.nan, np.nan, -1, "failed", repr(e)[:200]
            dt = time.time() - t0
            print(f"  {method:8s} {mode:17s} rep{rep}: acc={acc:.3f} frac={frac:.3f} "
                  f"n_method_cells={n_cells} (true {n_true}) {dt:.0f}s {st}")
            rows.append({"method": method, "mode": mode, "config": CONFIG, "rep": rep,
                         "matched_accuracy": acc, "frac_assigned": frac,
                         "n_method_cells": n_cells, "n_true_interior_cells": n_true,
                         "runtime_s": dt, "status": st, "error": err})
            pd.DataFrame(rows).to_csv(os.path.join(hc.RESULTS, "headroom_linux_repro.csv"), index=False)

    df = pd.DataFrame(rows)
    print("\nsummary (mean +/- range over reps):")
    for (method, mode), g in df[df.status == "ok"].groupby(["method", "mode"]):
        a = g.matched_accuracy
        print(f"  {method:8s} {mode:17s}: acc {a.mean():.3f} [{a.min():.3f}, {a.max():.3f}] "
              f"n_method_cells {g.n_method_cells.min()}-{g.n_method_cells.max()}")
    print(f"\nwrote {os.path.join(hc.RESULTS, 'headroom_linux_repro.csv')}")


if __name__ == "__main__":
    main()
