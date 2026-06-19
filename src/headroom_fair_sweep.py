"""Task 1, method-fairness gate: give each method its best honest configuration.

For Baysor, Proseg, and ComSeg, sweep a small grid over the 1 to 3 parameters each method's
documentation/paper flags as most impactful for high-plex imaging-based spatial
transcriptomics at this data scale, in nuclei-prior mode, at the headline dense operating
point (dense_sigma1.99, the data-pinned sigma). Record matched accuracy and all
over-segmentation-robust metrics per config; the best config per method (by one-to-one matched
accuracy) is then evaluated at every operating point by headroom_fair_eval.py.

Swept parameters (documented as impactful):
  Baysor: scale (expected cell radius, relative to the true mean cell radius), min-molecules-
          per-cell, prior-segmentation-confidence.
  Proseg: voxel-size, prior-seg-reassignment-prob, nuclear-reassignment-prob, diffusion model.
  ComSeg: mean_cell_diameter (the co-expression graph radius, ComSeg's key parameter).

ComSeg is swept at dense_sigma1.99 only because each run is about 210 s (NUMBA_DISABLE_JIT=1
on this aarch64 host); Baysor and Proseg are also swept at dense_sigma1.99 for one consistent
sweep point. Writes results/headroom_fair_sweep.csv (append-only, checkpointed).

  micromamba run -n st python src/headroom_fair_sweep.py
  (ComSeg rows require NUMBA_DISABLE_JIT=1; run that method separately with the env var.)
"""
from __future__ import annotations

import os
import sys
import time
import traceback

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import headroom_common as hc  # noqa: E402
from headroom_fair_common import all_metrics  # noqa: E402

DATA_LINUX = os.path.join(hc.DATA, "headroom_linux")
SWEEP_CONFIG = "dense_sigma1.99"
OUT = os.path.join(hc.RESULTS, "headroom_fair_sweep.csv")


def baysor_grid(cells_df):
    from methods_baysor import _scale_um, run_baysor
    base = _scale_um(cells_df)   # true mean cell radius (um)
    runs = []
    for rel in (1.0, 1.25, 1.5, 2.0):
        for conf in (0.5, 0.8):
            runs.append((f"scale{rel}x_m50_c{conf}",
                         dict(scale=rel * base, min_molecules=50, prior_confidence=conf,
                              scale_rel=rel)))
    for rel in (1.25, 1.5):
        runs.append((f"scale{rel}x_m20_c0.5",
                     dict(scale=rel * base, min_molecules=20, prior_confidence=0.5,
                          scale_rel=rel)))
    return run_baysor, runs


def proseg_grid(cells_df):
    from methods_proseg import run_proseg
    runs = [
        ("defaults", dict(extra_args=None)),
        ("voxel0.5", dict(extra_args=["--voxel-size", "0.5"])),
        ("priorreassign0.3", dict(extra_args=["--prior-seg-reassignment-prob", "0.3"])),
        ("priorreassign0.7", dict(extra_args=["--prior-seg-reassignment-prob", "0.7"])),
        ("nodiffusion", dict(extra_args=["--no-diffusion"])),
        ("nuclearreassign0.4", dict(extra_args=["--nuclear-reassignment-prob", "0.4"])),
    ]
    return run_proseg, runs


def comseg_grid(cells_df):
    from methods_comseg import run_comseg
    runs = [(f"mcd{d}", dict(mean_cell_diameter=float(d))) for d in (8, 10, 12)]
    return run_comseg, runs


GRIDS = {"Baysor": baysor_grid, "Proseg": proseg_grid, "ComSeg": comseg_grid}


def main(methods=("Baysor", "Proseg")):
    from gate2_pin import build_xenium
    model, _, _ = build_xenium()
    d = os.path.join(DATA_LINUX, SWEEP_CONFIG)
    t = pd.read_csv(os.path.join(d, "transcripts.csv"))
    c = pd.read_csv(os.path.join(d, "cells.csv"))
    inter = t["interior"].to_numpy().astype(bool)
    true = t["true_cell"].to_numpy()

    prev = pd.read_csv(OUT).to_dict("records") if os.path.exists(OUT) else []
    rerun = {m for m in methods}
    kept = [r for r in prev if r["method"] not in rerun]
    rows = []
    for method in methods:
        fn, runs = GRIDS[method](c)
        for label, kw in runs:
            run_kw = {k: v for k, v in kw.items() if k != "scale_rel"}
            t0 = time.time()
            try:
                assigned = fn(t, c, model, with_prior=True, **run_kw)
                assigned = assigned[0] if isinstance(assigned, tuple) else assigned
                met = all_metrics(assigned, true, inter)
                status, err = "ok", ""
            except Exception as e:  # noqa: BLE001
                met = {}
                status, err = "failed", repr(e)[:300]
                traceback.print_exc()
            dt = time.time() - t0
            row = {"method": method, "config_label": label, "sweep_at": SWEEP_CONFIG,
                   "scale_rel": kw.get("scale_rel", ""), "params": str(run_kw),
                   "runtime_s": dt, "status": status, "error": err, **met}
            rows.append(row)
            print(f"  {method:8s} {label:22s} acc1to1={met.get('acc_one_to_one', float('nan')):.3f} "
                  f"many2one={met.get('acc_many_to_one', float('nan')):.3f} "
                  f"ari={met.get('ari', float('nan')):.3f} "
                  f"ratio={met.get('pred_true_ratio', float('nan')):.2f} ({dt:.0f}s) {status}")
            pd.DataFrame(kept + rows).to_csv(OUT, index=False)

    df = pd.DataFrame(kept + rows)
    df.to_csv(OUT, index=False)
    # report best per method (by one-to-one matched accuracy)
    print("\nbest config per method (by acc_one_to_one):")
    cur = df[df.status == "ok"]
    for method in cur.method.unique():
        g = cur[cur.method == method]
        b = g.loc[g.acc_one_to_one.idxmax()]
        print(f"  {method:8s} -> {b.config_label}  acc1to1={b.acc_one_to_one:.3f} "
              f"many2one={b.acc_many_to_one:.3f} ari={b.ari:.3f}")
    print(f"\nwrote {OUT}")
    return df


if __name__ == "__main__":
    sel = sys.argv[1:] or ["Baysor", "Proseg"]
    main(tuple(sel))
