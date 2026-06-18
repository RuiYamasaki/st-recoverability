"""Re-derive the Gate 2 report numbers from the committed CSVs (each cites its source).
Regen command:  python src/gate2_report.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def main():
    pin = pd.read_csv(os.path.join(RESULTS, "gate2_xenium_pin.csv")).iloc[0]
    unc = pd.read_csv(os.path.join(RESULTS, "gate2_sigma_uncertainty.csv")).iloc[0]
    dg = pd.read_csv(os.path.join(RESULTS, "gate2_design_grid.csv"))
    nb = pd.read_csv(os.path.join(RESULTS, "gate2_nbinom.csv")).iloc[0]

    print("=" * 92)
    print("TASK 1  pin sigma directly on Xenium  (source: results/gate2_xenium_pin.csv)")
    print("-" * 92)
    print(f"  Xenium types: MiniBatchKMeans, K={int(pin.n_types)}; exclusive markers={int(pin.n_markers)}; "
          f"median NN={pin.median_nn_um:.2f} um")
    print(f"  real spatial leakage={pin.real_spatial_leakage_mean:.4f} (literal diagnostic={pin.real_literal_leakage_mean:.3f})")
    print(f"  Xenium-pinned sigma={pin.xenium_pinned_sigma_um:.2f} um  (MERFISH-pinned was {pin.merfish_pinned_sigma_um} um)")
    print(f"  oracle @ pinned: dense(pack={pin.dense_packing_cells_per_mm2:.0f})={pin.oracle_acc_dense_at_pinned:.3f} "
          f"(naive {pin.naive_acc_dense_at_pinned:.3f}); sparse={pin.oracle_acc_sparse_at_pinned:.3f}")

    print("\n" + "=" * 92)
    print("TASK 2  statistical CI + design-sensitivity on sigma  (results/gate2_sigma_uncertainty.csv, gate2_design_grid.csv)")
    print("-" * 92)
    print(f"  leakage point={unc.leakage_point:.4f}  bootstrap 95% CI=[{unc.leakage_ci_lo:.4f}, {unc.leakage_ci_hi:.4f}] "
          f"(B={int(unc.n_boot)})")
    print(f"  sigma point={unc.sigma_point_um:.2f}  statistical 95% CI=[{unc.sigma_stat_ci_lo_um:.2f}, {unc.sigma_stat_ci_hi_um:.2f}] um")
    print(f"  design-sensitivity sigma range=[{unc.sigma_design_lo_um:.2f}, {unc.sigma_design_hi_um:.2f}] um")
    print(f"  combined band=[{unc.sigma_combined_lo_um:.2f}, {unc.sigma_combined_hi_um:.2f}] um")
    print(f"  DENSE oracle (converged grid res_cell={unc.oracle_res_cell:.0f}): point={unc.oracle_dense_at_point:.3f}")
    print(f"    at statistical-CI ends: low-sigma={unc.oracle_dense_at_stat_ci_lo:.3f}  high-sigma={unc.oracle_dense_at_stat_ci_hi:.3f}")
    print(f"    at combined-band ends:  low-sigma(optimistic)={unc.oracle_dense_at_combined_lo:.3f}  "
          f"high-sigma={unc.oracle_dense_at_combined_hi:.3f}")
    print(f"  SPARSE oracle: point={unc.oracle_sparse_at_point:.3f}  combined ends "
          f"[{unc.oracle_sparse_at_combined_hi:.3f}, {unc.oracle_sparse_at_combined_lo:.3f}]")
    print("  design grid sigma by marker cutoff (the dominant design lever):")
    for thr, sub in dg.groupby("excl_threshold"):
        ok = sub[np.isfinite(sub.sigma_um)]
        if len(ok):
            print(f"    cutoff {thr}: {int(sub.n_markers.iloc[0])} markers, sigma range "
                  f"[{ok.sigma_um.min():.2f}, {ok.sigma_um.max():.2f}] um")
        else:
            print(f"    cutoff {thr}: 0 markers (degenerate, excluded)")

    print("\n" + "=" * 92)
    print("TASK 3  non-idealized negative-binomial emission  (source: results/gate2_nbinom.csv)")
    print("-" * 92)
    print(f"  NB dispersion phi median={nb.nb_dispersion_median:.3f} (Var = mu + phi*mu^2)")
    print(f"  re-pinned sigma: poisson={nb.sigma_poisson_um:.2f} -> nbinom={nb.sigma_nbinom_um:.2f} um "
          f"(shift {nb.sigma_shift_um:+.2f})")
    print(f"  dense oracle: poisson={nb.dense_oracle_poisson:.3f} -> nbinom={nb.dense_oracle_nbinom:.3f} "
          f"(shift {nb.dense_oracle_shift:+.3f}); naive(nb)={nb.dense_naive_nbinom:.3f}")
    print(f"  sparse oracle: poisson={nb.sparse_oracle_poisson:.3f} -> nbinom={nb.sparse_oracle_nbinom:.3f}")
    print(f"  NB frontier: oracle range [{nb.nb_frontier_oracle_min:.3f}, {nb.nb_frontier_oracle_max:.3f}]; "
          f"mono(sig,pack)=({int(nb.nb_frontier_mono_viol_sigma)},{int(nb.nb_frontier_mono_viol_packing)}); "
          f"oracle>=naive all={bool(nb.nb_frontier_oracle_ge_naive_all)} (min gap {nb.nb_frontier_min_gap:+.3f})")


if __name__ == "__main__":
    main()
