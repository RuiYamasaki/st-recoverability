"""Re-derive the Gate-0 report table and headline numbers from the committed CSVs.

    python src/report_table.py

Every printed number cites its source file and column. This is the regen command
for the representative grid-point table in the Gate-0 report.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def _row(df, packing, sigma, density):
    r = df[(df.packing_cells_per_mm2 == packing) & (df.sigma_um == sigma) &
           (df.mean_tx_per_cell == density)]
    return r.iloc[0]


def main():
    sweep = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    anchor = pd.read_csv(os.path.join(RESULTS, "realism_oracle.csv"))
    real = pd.read_csv(os.path.join(RESULTS, "realism.csv"))

    print("=" * 100)
    print("REPRESENTATIVE GRID POINTS  (source: results/sweep.csv)")
    print("cols: oracle_acc | profile_ari_vs_truetype | assign_ari_vs_truecells_oracle | naive_acc")
    print("-" * 100)
    reps = [
        ("low-difficulty extreme", 2000.0, 0.5, 160.0),
        ("MERFISH-like packing, realistic sigma", 2000.0, 2.0, 160.0),
        ("mid", 5700.0, 4.0, 160.0),
        ("mid", 9500.0, 2.0, 160.0),
        ("Xenium-like packing, realistic sigma", 16000.0, 2.0, 160.0),
        ("high-difficulty extreme", 16000.0, 8.0, 160.0),
    ]
    print(f"{'label':38s} {'pack':>6s} {'sig':>4s} {'dens':>5s} | "
          f"{'oracle_acc':>10s} {'prof_ARI':>9s} {'assignARI':>9s} {'naive_acc':>9s}")
    for lab, p, s, d in reps:
        r = _row(sweep, p, s, d)
        print(f"{lab:38s} {p:6.0f} {s:4.1f} {d:5.0f} | "
              f"{r.oracle_acc:10.3f} {r.profile_ari_vs_truetype:9.3f} "
              f"{r.assign_ari_vs_truecells_oracle:9.3f} {r.naive_acc:9.3f}")

    print("\n" + "=" * 100)
    print("HEADLINE SCALARS")
    print("-" * 100)
    print(f"[sweep.csv] full-grid oracle_acc range            : "
          f"[{sweep.oracle_acc.min():.3f}, {sweep.oracle_acc.max():.3f}]")
    print(f"[sweep.csv] full-grid oracle assign-ARI range     : "
          f"[{sweep.assign_ari_vs_truecells_oracle.min():.3f}, {sweep.assign_ari_vs_truecells_oracle.max():.3f}]")
    print(f"[sweep.csv] full-grid oracle profile-ARI range    : "
          f"[{sweep.profile_ari_vs_truetype.min():.3f}, {sweep.profile_ari_vs_truetype.max():.3f}]")
    print(f"[sweep.csv] full-grid naive profile-ARI range     : "
          f"[{sweep.naive_profile_ari_vs_truetype.min():.3f}, {sweep.naive_profile_ari_vs_truetype.max():.3f}]")
    print(f"[sweep.csv] oracle_acc >= naive_acc at all 125 pts : {bool(sweep.oracle_ge_naive.all())}")
    print(f"[sweep.csv] min(oracle_acc - naive_acc)            : {sweep.oracle_minus_naive.min():.4f}")

    # density-independence of assignment accuracy: spread across densities at fixed (packing,sigma)
    g = sweep.groupby(["packing_cells_per_mm2", "sigma_um"])["oracle_acc"]
    spread = (g.max() - g.min())
    print(f"[sweep.csv] max spread of oracle_acc across the 5 densities (fixed pack,sig): "
          f"{spread.max():.4f}  (mean {spread.mean():.4f})  -> assignment accuracy is ~density-independent")

    # monotonicity: oracle_acc non-increasing in sigma, and in packing
    viol_s = 0
    for (p, d), sub in sweep.groupby(["packing_cells_per_mm2", "mean_tx_per_cell"]):
        a = sub.sort_values("sigma_um").oracle_acc.to_numpy()
        viol_s += int((np.diff(a) > 1e-3).sum())
    viol_p = 0
    for (s, d), sub in sweep.groupby(["sigma_um", "mean_tx_per_cell"]):
        a = sub.sort_values("packing_cells_per_mm2").oracle_acc.to_numpy()
        viol_p += int((np.diff(a) > 1e-3).sum())
    print(f"[sweep.csv] monotonic decrease in sigma  : violations={viol_s} (of {5*5} curves x4 steps)")
    print(f"[sweep.csv] monotonic decrease in packing: violations={viol_p} (of {5*5} curves x4 steps)")

    print("\n" + "=" * 100)
    print("REALISTIC BAND  (source: results/realism_oracle.csv, in_realistic_sigma_band == True)")
    print("-" * 100)
    band = anchor[anchor.in_realistic_sigma_band]
    print(f"oracle_acc in realistic band (sigma<=2um)         : "
          f"[{band.oracle_acc.min():.3f}, {band.oracle_acc.max():.3f}]")
    print(f"  Xenium-like FOVs (dense) oracle_acc in band     : "
          f"[{band[band.dataset=='xenium_breast'].oracle_acc.min():.3f}, "
          f"{band[band.dataset=='xenium_breast'].oracle_acc.max():.3f}]")
    print(f"  MERFISH-like FOVs (sparse) oracle_acc in band   : "
          f"[{band[band.dataset=='merfish_hypothal'].oracle_acc.min():.3f}, "
          f"{band[band.dataset=='merfish_hypothal'].oracle_acc.max():.3f}]")
    print(f"oracle profile-ARI in realistic band              : "
          f"[{band.profile_ari_vs_truetype.min():.3f}, {band.profile_ari_vs_truetype.max():.3f}]")
    print(f"oracle assign-ARI in realistic band               : "
          f"[{band.assign_ari_vs_truecells_oracle.min():.3f}, {band.assign_ari_vs_truecells_oracle.max():.3f}]")
    print(f"max tx-match rel err at any anchor                : {anchor.tx_match_rel_err.max():.4f}  "
          f"(within 15%: {bool(anchor.tx_match_within_15pct.all())})")
    print(f"max packing-match rel err at any anchor           : {anchor.packing_match_rel_err.max():.2e}  "
          f"(within 15%: {bool(anchor.packing_match_within_15pct.all())})")

    print("\n" + "=" * 100)
    print("REAL-DATA SUMMARY STATISTICS  (source: results/realism.csv)")
    print("-" * 100)
    cols = ["dataset", "fov_id", "n_cells", "median_tx_per_cell",
            "packing_cells_per_mm2", "median_nn_distance_um"]
    with pd.option_context("display.width", 160):
        print(real[cols].to_string(index=False))


if __name__ == "__main__":
    main()
