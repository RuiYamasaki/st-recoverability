"""Task 3, method-fairness gate: free-segmentation (no per-cell nucleus prior) for each method.

Baysor: native de-novo joint segmentation (no prior). Run at every operating point using the
best scale / min-molecules found by the nuclei-prior sweep (scale is a cell-size geometry
input that transfers; prior-confidence is unused without a prior).

Proseg: requires an initial segmentation (--cell-id-column from a nucleus or cellpose mask)
and aborts without one. It has no documented prior-free mode. Recorded as unsupported.

ComSeg: the valid prior-free graph-partition method is "louvain" (prior_name=None); the
earlier "without_prior" failure was an invalid argument (valid options are with_prior /
louvain). ComSeg's run_all still associates RNA communities to the centroid landmarks
(classify_centroid + associate_rna2landmark), so the true centroids are retained: this is "no
per-RNA nucleus prior", NOT fully landmark-free. Run at the dense operating points and label
accordingly. Run this driver with NUMBA_DISABLE_JIT=1 (ComSeg) on aarch64.

Writes results/headroom_fair_free.csv. All metrics from headroom_fair_common.all_metrics.
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
from headroom_fair_sweep import GRIDS  # noqa: E402

DATA_LINUX = os.path.join(hc.DATA, "headroom_linux")
SWEEP_CSV = os.path.join(hc.RESULTS, "headroom_fair_sweep.csv")
OUT = os.path.join(hc.RESULTS, "headroom_fair_free.csv")
DENSE = ["dense_sigma1.43", "dense_sigma1.99", "dense_sigma2.63"]


def best_baysor_kwargs(cells_df):
    df = pd.read_csv(SWEEP_CSV)
    g = df[(df.method == "Baysor") & (df.status == "ok")]
    lbl = str(g.loc[g.acc_one_to_one.idxmax(), "config_label"])
    _, runs = GRIDS["Baysor"](cells_df)
    kw = dict(dict(runs)[lbl]); kw.pop("scale_rel", None)
    return lbl, kw


def _row(method, mode, config, label, met, dt, status, err):
    return {"method": method, "mode": mode, "config": config, "config_label": label,
            "runtime_s": dt, "status": status, "error": err, **met}


def main(configs=None):
    from gate2_pin import build_xenium
    model, _, _ = build_xenium()
    configs = configs or [name for name, _, _ in hc.CONFIGS]
    rows = []

    for name, packing, sigma in hc.CONFIGS:
        if name not in configs:
            continue
        d = os.path.join(DATA_LINUX, name)
        t = pd.read_csv(os.path.join(d, "transcripts.csv"))
        c = pd.read_csv(os.path.join(d, "cells.csv"))
        inter = t["interior"].to_numpy().astype(bool)
        true = t["true_cell"].to_numpy()
        print(f"[{name}]")

        # --- Baysor free (best scale/min_mol) ---
        from methods_baysor import run_baysor
        lbl, kw = best_baysor_kwargs(c)
        kw.pop("prior_confidence", None)
        t0 = time.time()
        try:
            a, _ = run_baysor(t, c, model, with_prior=False, **kw)
            met, status, err = all_metrics(a, true, inter), "ok", ""
        except Exception as e:  # noqa: BLE001
            met, status, err = {}, "failed", repr(e)[:300]; traceback.print_exc()
        rows.append(_row("Baysor", "free_segmentation", name, lbl + "_free", met,
                         time.time() - t0, status, err))
        print(f"  Baysor   free [{lbl}] one2one={met.get('acc_one_to_one', float('nan')):.3f} "
              f"ari={met.get('ari', float('nan')):.3f} ({rows[-1]['runtime_s']:.0f}s) {status}")

        # --- Proseg free (expected: unsupported) ---
        from methods_proseg import run_proseg
        t0 = time.time()
        try:
            a, _ = run_proseg(t, c, model, with_prior=False)
            met, status, err = all_metrics(a, true, inter), "ok", ""
        except RuntimeError as e:
            met, status, err = {}, "unsupported", repr(e)[:300]
        except Exception as e:  # noqa: BLE001
            met, status, err = {}, "failed", repr(e)[:300]
        rows.append(_row("Proseg", "free_segmentation", name, "native", met,
                         time.time() - t0, status, err))
        print(f"  Proseg   free: {status}")

        # --- ComSeg free (louvain; centroids retained, no per-RNA prior) at dense only ---
        if name in DENSE:
            from methods_comseg import run_comseg
            t0 = time.time()
            try:
                out = run_comseg(t, c, model, with_prior=False)
                a = out[0] if isinstance(out, tuple) else out
                met, status, err = all_metrics(a, true, inter), "ok", ""
            except Exception as e:  # noqa: BLE001
                met, status, err = {}, "failed", repr(e)[:300]; traceback.print_exc()
            rows.append(_row("ComSeg", "free_louvain_centroids", name, "louvain", met,
                             time.time() - t0, status, err))
            print(f"  ComSeg   free(louvain,centroids): "
                  f"one2one={met.get('acc_one_to_one', float('nan')):.3f} "
                  f"ari={met.get('ari', float('nan')):.3f} ({rows[-1]['runtime_s']:.0f}s) {status}")

        pd.DataFrame(rows).to_csv(OUT, index=False)

    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
