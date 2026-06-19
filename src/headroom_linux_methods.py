"""Run the field-standard Xenium methods (Baysor, Proseg) against the within-Linux oracle
ceiling on every exported headroom config, in nuclei-prior and free-segmentation modes,
score each against the known truth with the same matched-accuracy scorer used for the
Windows run, and write results/headroom_linux_methods.csv (append-only, namespaced).

Reuses headroom_common.matched_accuracy (the scorer) and CONFIGS unchanged. Reads the raw
tables exported by headroom_linux_export.py from data/headroom_linux/<config>/. Optionally
re-runs pciSeq/ComSeg for one consistent within-Linux table.

  micromamba run -n st python src/headroom_linux_methods.py            # Baysor + Proseg
  micromamba run -n st python src/headroom_linux_methods.py pciseq comseg baysor proseg

Methods/modes that fail or are unsupported are recorded with their error/status, never
silently dropped. frac_assigned (fraction of interior transcripts given a real cell, not
background) is recorded alongside matched accuracy so background-dumping is visible.
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

DATA_LINUX = os.path.join(hc.DATA, "headroom_linux")

# (method, mode, with_prior). Proseg free is unsupported (no prior-free mode) and recorded.
RUNS = [
    ("Baysor", "nuclei_prior", True),
    ("Baysor", "free_segmentation", False),
    ("Proseg", "nuclei_prior", True),
    ("Proseg", "free_segmentation", False),
]
OPTIONAL_RUNS = [
    ("pciSeq", "nuclei_prior", True),
    ("ComSeg", "nuclei_prior", True),
    ("ComSeg", "free_segmentation", False),
]


def _runner(method):
    if method == "Baysor":
        from methods_baysor import run_baysor
        return run_baysor
    if method == "Proseg":
        from methods_proseg import run_proseg
        return run_proseg
    if method == "pciSeq":
        from methods_pciseq import run_pciseq
        return lambda t, c, m, with_prior=True: run_pciseq(t, c, m)
    if method == "ComSeg":
        from methods_comseg import run_comseg
        return run_comseg
    raise ValueError(method)


def _score(method, mode, config_name, with_prior, t, c, model):
    inter = t["interior"].to_numpy().astype(bool)
    true = t["true_cell"].to_numpy()
    t0 = time.time()
    try:
        fn = _runner(method)
        out = fn(t, c, model, with_prior=with_prior)
        assigned = out[0] if isinstance(out, tuple) else out
        acc = hc.matched_accuracy(assigned, true, inter)
        frac = float((assigned[inter] >= 0).mean())
        status, err = "ok", ""
    except RuntimeError as e:
        # method genuinely does not support this mode (e.g. Proseg has no prior-free mode)
        acc, frac, status, err = float("nan"), float("nan"), "unsupported", repr(e)[:300]
    except Exception as e:  # noqa: BLE001
        acc, frac, status, err = float("nan"), float("nan"), "failed", repr(e)[:300]
        traceback.print_exc()
    dt = time.time() - t0
    aa = acc if acc == acc else float("nan")
    ff = frac if frac == frac else float("nan")
    print(f"  {method:8s} [{mode:17s}] {config_name:26s} acc={aa:.3f} "
          f"frac_assigned={ff:.3f} ({dt:.0f}s) {status}")
    return {"method": method, "mode": mode, "config": config_name,
            "matched_accuracy": acc, "frac_assigned": frac, "runtime_s": dt,
            "status": status, "error": err}


def main(methods=None):
    from gate2_pin import build_xenium
    model, real, abundance = build_xenium()

    runs = list(RUNS)
    if methods:
        want = {m.lower() for m in methods}
        runs = [r for r in (RUNS + OPTIONAL_RUNS) if r[0].lower() in want]

    out = os.path.join(hc.RESULTS, "headroom_linux_methods.csv")
    # merge-safe: keep rows for methods not being (re-)run, replace the ones we run now
    rerun = {m for m, _, _ in runs}
    if os.path.exists(out):
        prev = pd.read_csv(out)
        kept = prev[~prev["method"].isin(rerun)].to_dict("records")
    else:
        kept = []

    rows = []
    for name, packing, sigma in hc.CONFIGS:
        d = os.path.join(DATA_LINUX, name)
        t = pd.read_csv(os.path.join(d, "transcripts.csv"))
        c = pd.read_csv(os.path.join(d, "cells.csv"))
        print(f"[{name}] T={len(t)} cells={len(c)}")
        for method, mode, with_prior in runs:
            rows.append(_score(method, mode, name, with_prior, t, c, model))
            # checkpoint after every run so a long sweep is never lost
            pd.DataFrame(kept + rows).to_csv(out, index=False)

    pd.DataFrame(kept + rows).to_csv(out, index=False)
    print(f"\nwrote {out}")
    return pd.DataFrame(kept + rows)


if __name__ == "__main__":
    main(sys.argv[1:] or None)
