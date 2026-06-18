"""Run the real published methods (pciSeq nuclei-prior; ComSeg nuclei-prior and free) on
every exported headroom config, score each against the known truth with the matched
accuracy, and write results/headroom_methods.csv. The oracle ceiling and the naive
baseline come from results/headroom_oracle_naive.csv (headroom_export.py).

Methods that fail to run are recorded with their error, not silently dropped.
"""
from __future__ import annotations

import os
import sys
import time
import traceback

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate2_pin import build_xenium  # noqa: E402
from headroom_common import CONFIGS, EXPORT_DIR, RESULTS, matched_accuracy  # noqa: E402


def _score(name, config_name, mode, fn, t, c, model):
    inter = t["interior"].to_numpy().astype(bool)
    true = t["true_cell"].to_numpy()
    t0 = time.time()
    try:
        out = fn(t, c, model)
        assigned = out[0] if isinstance(out, tuple) else out
        acc = matched_accuracy(assigned, true, inter)
        frac = float((assigned[inter] >= 0).mean())
        status, err = "ok", ""
    except Exception as e:
        acc, frac, status, err = float("nan"), float("nan"), "failed", repr(e)[:200]
        traceback.print_exc()
    dt = time.time() - t0
    print(f"  {name:16s} [{mode}] {config_name:26s} acc={acc if acc==acc else float('nan'):.3f} "
          f"frac_assigned={frac if frac==frac else float('nan'):.3f} ({dt:.0f}s) {status}")
    return {"method": name, "mode": mode, "config": config_name,
            "matched_accuracy": acc, "frac_assigned": frac, "runtime_s": dt,
            "status": status, "error": err}


def main(methods=("pciseq", "comseg")):
    model, real, abundance = build_xenium()
    from methods_pciseq import run_pciseq
    from methods_comseg import run_comseg

    rows = []
    for name, packing, sigma in CONFIGS:
        d = os.path.join(EXPORT_DIR, name)
        t = pd.read_csv(os.path.join(d, "transcripts.csv"))
        c = pd.read_csv(os.path.join(d, "cells.csv"))
        print(f"[{name}] T={len(t)} cells={len(c)}")
        if "pciseq" in methods:
            rows.append(_score("pciSeq", name, "nuclei_prior",
                               lambda t, c, m: run_pciseq(t, c, m), t, c, model))
        if "comseg" in methods:
            rows.append(_score("ComSeg", name, "nuclei_prior",
                               lambda t, c, m: run_comseg(t, c, m, with_prior=True), t, c, model))
    # one free-segmentation run (context): ComSeg without the nuclei prior, dense point
    if "comseg" in methods:
        name = "dense_sigma1.99"
        d = os.path.join(EXPORT_DIR, name)
        t = pd.read_csv(os.path.join(d, "transcripts.csv"))
        c = pd.read_csv(os.path.join(d, "cells.csv"))
        rows.append(_score("ComSeg", name, "free_segmentation",
                           lambda t, c, m: run_comseg(t, c, m, with_prior=False), t, c, model))

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "headroom_methods.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")
    return df


if __name__ == "__main__":
    main()
