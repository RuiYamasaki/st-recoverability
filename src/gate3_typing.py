"""Gate 3 Task 2: sigma robustness to cell-typing.

Re-derive Xenium-breast cell types under several reasonable clustering choices (different
k and different clustering seeds), re-select admissible markers with the Task 1 validity
test for each, re-pin sigma, and report the spread of sigma and dense-regime oracle
accuracy across the typing choices. Stability here is required.

Output: results/gate3_typing.csv. Seeds recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import build_realistic_model_from_xenium  # noqa: E402
from gate3_common import marker_validity, admissible_owner, pin_admissible, candidate_markers  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
GATE3_SEED = config.MASTER_SEED + 510000

# (label, n_types, clustering_seed)
VARIANTS = [
    ("kmeans_K10", 10, config.MASTER_SEED),
    ("kmeans_K14", 14, config.MASTER_SEED),          # the Gate 2 baseline
    ("kmeans_K14_altseed", 14, config.MASTER_SEED + 1),
    ("kmeans_K20", 20, config.MASTER_SEED),
]


def main():
    rows = []
    for vi, (label, K, cseed) in enumerate(VARIANTS):
        model, real = build_realistic_model_from_xenium(n_types=K, seed=cseed,
                                                        name=f"xenium_breast_{label}")
        abundance = model.mean_expr.sum(axis=1); abundance = abundance / abundance.mean()
        owner_cand, candidates = candidate_markers(model)
        seed = GATE3_SEED + 1000 * (vi + 1)
        vdf = marker_validity(real, owner_cand, candidates, seed)
        n_adm = int(vdf.admissible.sum())
        owner_adm = admissible_owner(model, vdf)
        if n_adm == 0:
            print(f"[{label}] K={K} no admissible markers; skipping pin")
            continue
        pin = pin_admissible(model, real, abundance, owner_adm, seed + 100)
        rows.append({
            "variant": label, "n_types": K, "clustering_seed": cseed,
            "n_candidate": len(candidates), "n_admissible": n_adm,
            "real_leakage": pin["real_leakage"],
            "sigma_point_um": pin["sigma_point_um"],
            "sigma_ci_lo_um": pin["sigma_ci_lo_um"], "sigma_ci_hi_um": pin["sigma_ci_hi_um"],
            "dense_oracle_point": pin["dense_oracle_point"],
            "dense_oracle_at_sigma_ci_lo": pin["dense_oracle_at_sigma_ci_lo"],
            "dense_oracle_at_sigma_ci_hi": pin["dense_oracle_at_sigma_ci_hi"],
            "naive_dense_point": pin["naive_dense_point"],
            "sparse_oracle_point": pin["sparse_oracle_point"], "seed": seed,
        })
        print(f"[{label}] K={K}: {n_adm} admissible; sigma={pin['sigma_point_um']:.2f} "
              f"CI[{pin['sigma_ci_lo_um']:.2f},{pin['sigma_ci_hi_um']:.2f}]; "
              f"dense point={pin['dense_oracle_point']:.3f} "
              f"CI[{pin['dense_oracle_at_sigma_ci_hi']:.3f},{pin['dense_oracle_at_sigma_ci_lo']:.3f}]")

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "gate3_typing.csv")
    df.to_csv(out, index=False)
    if len(df):
        print(f"\nsigma spread across typing: [{df.sigma_point_um.min():.2f}, {df.sigma_point_um.max():.2f}] um")
        print(f"dense-accuracy-at-point spread: [{df.dense_oracle_point.min():.3f}, {df.dense_oracle_point.max():.3f}]")
        print(f"max dense accuracy at any sigma-CI optimistic end: {df.dense_oracle_at_sigma_ci_lo.max():.3f}")
    print(f"wrote {out}")
    return df


if __name__ == "__main__":
    main()
