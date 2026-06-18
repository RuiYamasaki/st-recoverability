"""Gate 1 Task 3: structural-assumption sensitivity (a sensitivity check, not the
primary decision).

Re-run a reduced (packing x sigma) grid under realistic expression for three conditions
and compare the frontier:
  - baseline : Voronoi cells, Gaussian displacement (the Task 1 model).
  - aniso    : non-Voronoi geometry (irregular, elongated cells via anisotropic nucleus
               expansion), Gaussian displacement.
  - mixture  : Voronoi cells, non-Gaussian heavy-tailed displacement (a fraction
               MIX_EPSILON of transcripts land uniformly in the domain). The oracle is
               made Bayes-optimal for this mixture.

For each condition we report the oracle-accuracy frontier range, monotonicity-violation
counts, whether oracle stays at or above naive, and the dense-regime value, each versus
the baseline. Output: results/gate1_structural.csv. Seeds recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import build_realistic_model_from_merfish  # noqa: E402
from generator import build_field  # noqa: E402
from oracle import build_oracle_maps  # noqa: E402
from evaluate import eval_point  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")

PACKINGS = [2575.0, 6000.0, 13625.0]   # MERFISH-like, mid, Xenium-like
SIGMAS = [1.0, 2.0, 4.0, 8.0]
DENSITY = 160.0
MIX_EPSILON = 0.10                      # heavy-tail background fraction (fixed, plausible)
GATE1_SEED = config.MASTER_SEED + 300000
CONDITIONS = ["baseline", "aniso", "mixture"]


def run_condition(model, cond):
    rows = []
    for pi, packing in enumerate(PACKINGS):
        for si, sigma in enumerate(SIGMAS):
            seed = GATE1_SEED + 1000 * (pi + 1) + 17 * (si + 1) + (CONDITIONS.index(cond) + 1) * 100000
            geometry = "aniso" if cond == "aniso" else "voronoi"
            disp = "gauss_uniform" if cond == "mixture" else "gaussian"
            eps = MIX_EPSILON if cond == "mixture" else 0.0
            field = build_field(packing, sigma, seed, model=model, geometry=geometry)
            Dmax, argcell = build_oracle_maps(field)
            r = eval_point(field, Dmax, argcell, DENSITY, seed + 1,
                           displacement=disp, disp_epsilon=eps)
            rows.append({
                "condition": cond, "packing_cells_per_mm2": packing, "sigma_um": sigma,
                "mean_tx_per_cell": DENSITY, "mix_epsilon": eps, "geometry": geometry,
                "oracle_acc": r["oracle_acc"], "naive_acc": r["naive_acc"],
                "oracle_minus_naive": r["oracle_minus_naive"], "oracle_ge_naive": r["oracle_ge_naive"],
                "n_interior_transcripts": r["n_interior_transcripts"],
                "n_cells": r["n_cells"], "n_interior_cells": r["n_interior_cells"],
                "seed": seed,
            })
    return rows


def _mono_violations(sub):
    vs = vp = 0
    for p, s in sub.groupby("packing_cells_per_mm2"):
        vs += int((np.diff(s.sort_values("sigma_um").oracle_acc.values) > 1e-3).sum())
    for sg, s in sub.groupby("sigma_um"):
        vp += int((np.diff(s.sort_values("packing_cells_per_mm2").oracle_acc.values) > 1e-3).sum())
    return vs, vp


def main():
    model = build_realistic_model_from_merfish()
    rows = []
    for cond in CONDITIONS:
        rows += run_condition(model, cond)
    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "gate1_structural.csv")
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} rows)\n")

    print(f"{'condition':9s} | oracle_acc range | mono(sig,pack) | oracle>=naive all | "
          f"dense(pack={PACKINGS[-1]:.0f}) sig2/sig8")
    for cond in CONDITIONS:
        sub = df[df.condition == cond]
        vs, vp = _mono_violations(sub)
        dense = sub[sub.packing_cells_per_mm2 == PACKINGS[-1]]
        d2 = float(dense[dense.sigma_um == 2.0].oracle_acc.iloc[0])
        d8 = float(dense[dense.sigma_um == 8.0].oracle_acc.iloc[0])
        print(f"{cond:9s} | [{sub.oracle_acc.min():.3f}, {sub.oracle_acc.max():.3f}] | "
              f"({vs},{vp}) | {bool(sub.oracle_ge_naive.all())} (min gap {sub.oracle_minus_naive.min():+.3f}) | "
              f"{d2:.3f}/{d8:.3f}")
    return df


if __name__ == "__main__":
    main()
