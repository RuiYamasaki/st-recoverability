"""Task 2, method-fairness gate: evaluate each method's BEST swept configuration at every
operating point under all three metrics.

Reads the best configuration per method from results/headroom_fair_sweep.csv (the row with the
highest one-to-one matched accuracy), then runs that configuration in nuclei-prior mode at the
dense (1.43, 1.99, 2.63 um), sparse, and representative operating points, computing the full
metric set (one-to-one matched, many-to-one homogeneity, transcript co-assignment ARI,
predicted/true cell-count ratio, frac_assigned) for every run. The oracle and naive baselines
are scored under the same metrics for reference. Writes results/headroom_fair_metrics.csv
(append-only, merge-safe per method, checkpointed).

  micromamba run -n st python src/headroom_fair_eval.py Baysor Proseg
  NUMBA_DISABLE_JIT=1 micromamba run -n st python src/headroom_fair_eval.py ComSeg
"""
from __future__ import annotations

import os
import sys
import time
import traceback

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import headroom_common as hc  # noqa: E402
from headroom_fair_common import all_metrics  # noqa: E402
from headroom_fair_sweep import GRIDS  # noqa: E402  (label -> kwargs grids per method)

DATA_LINUX = os.path.join(hc.DATA, "headroom_linux")
SWEEP_CSV = os.path.join(hc.RESULTS, "headroom_fair_sweep.csv")
OUT = os.path.join(hc.RESULTS, "headroom_fair_metrics.csv")


def best_label(method):
    df = pd.read_csv(SWEEP_CSV)
    g = df[(df.method == method) & (df.status == "ok")]
    if not len(g):
        return None
    return str(g.loc[g.acc_one_to_one.idxmax(), "config_label"])


def _runner(method):
    if method == "Baysor":
        from methods_baysor import run_baysor; return run_baysor
    if method == "Proseg":
        from methods_proseg import run_proseg; return run_proseg
    if method == "ComSeg":
        from methods_comseg import run_comseg; return run_comseg
    raise ValueError(method)


def _kwargs_for(method, cells_df, label):
    """Look up the kwargs for `label` from the method's grid evaluated on THIS config's cells
    (so Baysor's scale, relative to the true radius, is recomputed per operating point)."""
    _, runs = GRIDS[method](cells_df)
    d = dict(runs)
    kw = dict(d[label])
    kw.pop("scale_rel", None)
    return kw


def main(methods=("Baysor", "Proseg")):
    from gate2_pin import build_xenium
    model, _, _ = build_xenium()

    prev = pd.read_csv(OUT).to_dict("records") if os.path.exists(OUT) else []
    rerun = set(methods) | {"oracle", "naive"}
    kept = [r for r in prev if r["method"] not in rerun]
    rows = []

    for name, packing, sigma in hc.CONFIGS:
        d = os.path.join(DATA_LINUX, name)
        t = pd.read_csv(os.path.join(d, "transcripts.csv"))
        c = pd.read_csv(os.path.join(d, "cells.csv"))
        inter = t["interior"].to_numpy().astype(bool)
        true = t["true_cell"].to_numpy()
        print(f"[{name}]")

        # reference baselines (recomputed every invocation, merged)
        for ref in ("oracle", "naive"):
            met = all_metrics(t[f"{ref}_cell"].to_numpy(), true, inter)
            rows.append({"method": ref, "mode": "reference", "config": name,
                         "best_config": "", "runtime_s": 0.0, "status": "ok", "error": "", **met})

        for method in methods:
            lbl = best_label(method)
            if lbl is None:
                rows.append({"method": method, "mode": "nuclei_prior", "config": name,
                             "best_config": "(no sweep row)", "status": "missing", "error": "",
                             "runtime_s": 0.0})
                continue
            kw = _kwargs_for(method, c, lbl)
            t0 = time.time()
            try:
                out = _runner(method)(t, c, model, with_prior=True, **kw)
                assigned = out[0] if isinstance(out, tuple) else out
                met = all_metrics(assigned, true, inter)
                status, err = "ok", ""
            except Exception as e:  # noqa: BLE001
                met, status, err = {}, "failed", repr(e)[:300]
                traceback.print_exc()
            dt = time.time() - t0
            rows.append({"method": method, "mode": "nuclei_prior", "config": name,
                         "best_config": lbl, "runtime_s": dt, "status": status, "error": err, **met})
            print(f"  {method:8s} [{lbl}] one2one={met.get('acc_one_to_one', float('nan')):.3f} "
                  f"many2one={met.get('acc_many_to_one', float('nan')):.3f} "
                  f"ari={met.get('ari', float('nan')):.3f} ratio={met.get('pred_true_ratio', float('nan')):.2f} "
                  f"frac={met.get('frac_assigned', float('nan')):.3f} ({dt:.0f}s) {status}")
            pd.DataFrame(kept + rows).to_csv(OUT, index=False)

    pd.DataFrame(kept + rows).to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")
    return pd.DataFrame(kept + rows)


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("Baysor", "Proseg"))
