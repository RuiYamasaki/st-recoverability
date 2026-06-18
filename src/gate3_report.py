"""Re-derive the Gate 3 report numbers from the committed CSVs. Each number cites its
source file. Regen command:  python src/gate3_report.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def _validity_summary(tag):
    v = pd.read_csv(os.path.join(RESULTS, f"gate3_validity_{tag}.csv"))
    excl = v[v.exclusivity >= 0.7]
    loose = v[(v.exclusivity >= 0.6) & (v.exclusivity < 0.7)]
    print(f"  [{tag}] candidates={len(v)} admissible={int(v.admissible.sum())} (p<0.05, 1000 perms)")
    print(f"  [{tag}] mean displacement signal: admissible={v[v.admissible].real_signal.mean():+.4f}  "
          f"rejected={v[~v.admissible].real_signal.mean():+.4f}")
    print(f"  [{tag}] exclusive>=0.7: {int(excl.admissible.sum())}/{len(excl)} admissible; "
          f"loose[0.6,0.7): {int(loose.admissible.sum())}/{len(loose)} admissible "
          f"(rejected loose are the sigma-destabilisers)")
    if len(loose):
        rej = loose[~loose.admissible]
        print(f"  [{tag}] rejected loose markers (p-value): "
              f"{[f'{n}:{p:.3f}' for n, p in zip(rej.gene_name, rej.pval)]}")


def _pin_summary(tag):
    p = pd.read_csv(os.path.join(RESULTS, f"gate3_pin_{tag}.csv")).iloc[0]
    lo = min(p.dense_oracle_at_sigma_ci_lo, p.dense_oracle_at_sigma_ci_hi)
    hi = max(p.dense_oracle_at_sigma_ci_lo, p.dense_oracle_at_sigma_ci_hi)
    print(f"  [{tag}] re-pinned sigma={p.sigma_point_um:.2f} um, bootstrap CI "
          f"[{p.sigma_ci_lo_um:.2f}, {p.sigma_ci_hi_um:.2f}] (admissible markers={int(p.n_admissible_markers)})")
    print(f"  [{tag}] dense oracle accuracy: point={p.dense_oracle_point:.3f}; across sigma CI "
          f"[{lo:.3f}, {hi:.3f}]  -> below 0.95: {bool(hi < 0.95)}; below 0.9: {bool(hi < 0.9)} "
          f"(naive {p.naive_dense_point:.3f}); sparse point={p.sparse_oracle_point:.3f}")


def main():
    print("=" * 92)
    print("TASK 1  marker-validity test on BREAST  (results/gate3_validity_breast.csv, gate3_pin_breast.csv)")
    print("-" * 92)
    _validity_summary("breast")
    _pin_summary("breast")

    print("\n" + "=" * 92)
    print("TASK 2  sigma robustness to cell-typing  (results/gate3_typing.csv)")
    print("-" * 92)
    t = pd.read_csv(os.path.join(RESULTS, "gate3_typing.csv"))
    for _, r in t.iterrows():
        print(f"  {r.variant:20s} K={int(r.n_types):2d}: admissible={int(r.n_admissible):2d}  "
              f"sigma={r.sigma_point_um:.2f} CI[{r.sigma_ci_lo_um:.2f},{r.sigma_ci_hi_um:.2f}]  "
              f"dense point={r.dense_oracle_point:.3f} CI-optimistic={r.dense_oracle_at_sigma_ci_lo:.3f}")
    print(f"  sigma spread across typing: [{t.sigma_point_um.min():.2f}, {t.sigma_point_um.max():.2f}] um")
    print(f"  dense point spread: [{t.dense_oracle_point.min():.3f}, {t.dense_oracle_point.max():.3f}]; "
          f"max dense at any sigma-CI optimistic end: {t.dense_oracle_at_sigma_ci_lo.max():.3f}")

    print("\n" + "=" * 92)
    print("TASK 3  independent Xenium replication LUNG  (results/gate3_validity_lung.csv, gate3_pin_lung.csv)")
    print("-" * 92)
    _validity_summary("lung")
    _pin_summary("lung")


if __name__ == "__main__":
    main()
