"""Consolidated method-fairness report from the committed result CSVs. Raw numbers only.
Regen: micromamba run -n st python src/headroom_fair_report.py

Sources:
  results/headroom_fair_sweep.csv    (Task 1 parameter sweep, best config per method)
  results/headroom_fair_metrics.csv  (best config at every operating point, all 3 metrics)
  results/headroom_fair_free.csv     (Task 3 free-segmentation)
  results/headroom_linux_oracle_naive.csv (within-Linux oracle/naive, for reference)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DENSE = ["dense_sigma1.43", "dense_sigma1.99", "dense_sigma2.63"]
METRICS = ["acc_one_to_one", "acc_many_to_one", "ari"]


def _load(name):
    p = os.path.join(RESULTS, name)
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()


def f(v, w=6, p=3):
    return f"{v:{w}.{p}f}" if isinstance(v, (int, float)) and v == v else f"{'n/a':>{w}s}"


def sweep_section(sw):
    print("=" * 110)
    print("TASK 1 - fair configuration sweep (nuclei-prior, at dense_sigma1.99)")
    print("metric columns: one2one (Hungarian), many2one (homogeneity), ari (co-assignment), ratio (pred/true cells)")
    print("source: results/headroom_fair_sweep.csv")
    for method in sw.method.unique():
        g = sw[(sw.method == method) & (sw.status == "ok")]
        if not len(g):
            continue
        best = g.loc[g.acc_one_to_one.idxmax(), "config_label"]
        print(f"\n  {method}:")
        for _, r in g.iterrows():
            star = " <- BEST (one2one)" if r.config_label == best else ""
            print(f"    {r.config_label:24s} one2one={f(r.acc_one_to_one)} many2one={f(r.acc_many_to_one)} "
                  f"ari={f(r.ari)} ratio={f(r.pred_true_ratio,5,2)} frac={f(r.frac_assigned,5,3)}{star}")


def metrics_table(mt, on):
    print("\n" + "=" * 110)
    print("TASK 2 - best config per method at every operating point, all three metrics")
    print("source: results/headroom_fair_metrics.csv (best config from the sweep); oracle/naive recomputed within Linux")
    order = ["oracle", "naive", "Baysor", "Proseg", "ComSeg"]
    for cfg in on.config.tolist():
        print(f"\n[{cfg}]")
        sub = mt[mt.config == cfg]
        for method in order:
            r = sub[sub.method == method]
            if not len(r):
                continue
            r = r.iloc[0]
            tag = "" if r.get("status", "ok") == "ok" else f" ({r['status']})"
            lbl = f" [{r['best_config']}]" if isinstance(r.get("best_config"), str) and r["best_config"] else ""
            print(f"  {method:8s}{lbl:26s} one2one={f(r.get('acc_one_to_one'))} "
                  f"many2one={f(r.get('acc_many_to_one'))} ari={f(r.get('ari'))} "
                  f"ratio={f(r.get('pred_true_ratio'),5,2)} frac={f(r.get('frac_assigned'),5,3)}{tag}")


def margins(mt):
    print("\n" + "=" * 110)
    print("Best-fair-method margins in the DENSE regime, under each metric (source: headroom_fair_metrics.csv)")
    sophisticated = ["Baysor", "Proseg", "ComSeg"]
    for metric in METRICS:
        print(f"\n  metric = {metric}")
        for cfg in DENSE:
            sub = mt[mt.config == cfg]
            if not len(sub):
                continue
            orc = sub[sub.method == "oracle"][metric]
            nai = sub[sub.method == "naive"][metric]
            if not len(orc) or not len(nai):
                continue
            orc, nai = float(orc.iloc[0]), float(nai.iloc[0])
            vals = {m: float(sub[sub.method == m][metric].iloc[0])
                    for m in sophisticated
                    if len(sub[sub.method == m]) and sub[sub.method == m].status.iloc[0] == "ok"}
            if not vals:
                continue
            bm = max(vals, key=vals.get); bv = vals[bm]
            print(f"    {cfg:18s} oracle={f(orc)} naive={f(nai)} best={bm}({f(bv)}) "
                  f"best-minus-naive={f(bv-nai)} oracle-minus-best={f(orc-bv)}")


def free_section(fr):
    print("\n" + "=" * 110)
    print("TASK 3 - free-segmentation (source: results/headroom_fair_free.csv)")
    if not len(fr):
        print("  (no free results yet)")
        return
    for _, r in fr.iterrows():
        if r.status == "ok":
            print(f"  {r.method:8s} [{r['mode']}] {r.config:22s} one2one={f(r.get('acc_one_to_one'))} "
                  f"many2one={f(r.get('acc_many_to_one'))} ari={f(r.get('ari'))} "
                  f"ratio={f(r.get('pred_true_ratio'),5,2)} frac={f(r.get('frac_assigned'),5,3)}")
        else:
            print(f"  {r.method:8s} [{r['mode']}] {r.config:22s} {r.status}: {str(r.error)[:90]}")


def main():
    sw = _load("headroom_fair_sweep.csv")
    mt = _load("headroom_fair_metrics.csv")
    fr = _load("headroom_fair_free.csv")
    on = _load("headroom_linux_oracle_naive.csv")
    if len(sw):
        sweep_section(sw)
    if len(mt):
        metrics_table(mt, on)
        margins(mt)
    if len(fr):
        free_section(fr)


if __name__ == "__main__":
    main()
