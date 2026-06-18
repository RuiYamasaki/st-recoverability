"""Gate 2 Task 1: pin the displacement sigma DIRECTLY on Xenium (no cross-platform
borrowing), using the same distance-stratified spatial marker-leakage statistic as
Gate 1.

Xenium cell types: MiniBatchKMeans (seed = MASTER_SEED) on median-normalised log1p
counts of the cell-feature matrix (clusters < 50 cells dropped), giving K=14 types over
the 313-gene panel; exclusive markers as in Gate 1 (exclusivity threshold 0.7). Real
spatial leakage is measured on the whole imaged section. Sigma is pinned by inverting the
synthetic leakage curve built at Xenium's own representative (median local) packing and
density, then the same physical sigma is evaluated in the dense Xenium regime.

Dataset: 10x Genomics "Xenium FFPE Human Breast Cancer, Replicate 1" (release 1.0.1),
cell_feature_matrix.h5 and cells.parquet, direct HTTP, no auth.

Output: results/gate2_xenium_pin.csv. Seeds: GATE2_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import build_realistic_model_from_xenium  # noqa: E402
from generator import build_field, generate_transcripts  # noqa: E402
from oracle import naive_assign  # noqa: E402
from gate1_leakage_anchor import (spatial_leakage, literal_leakage,  # noqa: E402
                                  fit_sigma, oracle_acc_at, SIGMA_GRID, _counts_matrix, _median_nn)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
GATE2_SEED = config.MASTER_SEED + 400000

DENSE_PACKING = 13625.0    # Xenium dense-tissue packing (Gate 0/1 realism FOVs; section p90 ~ 12000)
SPARSE_PACKING = 2575.0    # MERFISH-like sparse packing, for context
MERFISH_PINNED_UM = 2.67   # Gate 1 value, for comparison
LEAK_NTARGET = 4000        # cells in the synthetic field for leakage estimation (large -> stable curve)
LEAK_NREP = 3              # replicate fields per sigma, averaged, to denoise the leakage curve
R_NEAR_BASE, R_FAR_BASE = 1.5, 4.0
# higher grid resolution for the dense-regime oracle accuracy: at small sigma the default
# grid under-resolves and UNDER-estimates accuracy; (res_cell, grid_max) below is converged
# enough to decide the 0.95 question (see Gate 2 report, resolution check).
ORACLE_RES_CELL, ORACLE_GRID_MAX = 24.0, 1500


def gen_synth_data(model, packing, density, abundance, seed, n_target=LEAK_NTARGET,
                   n_rep=LEAK_NREP, emission="poisson"):
    """Generate, per sigma on SIGMA_GRID, n_rep large synthetic fields with naive-assigned
    per-cell counts (averaged later to denoise the leakage curve). Returns data[i] = list
    of n_rep (centers, counts_matrix, types, median_nn); reusable for any markers/R."""
    data = []
    for i, s in enumerate(SIGMA_GRID):
        reps = []
        for r in range(n_rep):
            sd = seed + 31 * i + 7919 * r
            f = build_field(packing, s, sd, model=model, n_target=n_target)
            tx = generate_transcripts(f, density, sd + 1, abundance=abundance, emission=emission)
            na = naive_assign(f, tx.obs_xy)
            C = _counts_matrix(f.n_cells, na, tx.gene, model.n_genes)
            reps.append((f.centers, C, f.types, _median_nn(f.centers)))
        data.append(reps)
    return data


def _leak_one(centers, C, types, mnn, owner, r_near, r_far):
    leaks = []
    for g in np.where(owner >= 0)[0]:
        own = types == owner[g]
        if own.sum() < 3:
            continue
        d, _ = cKDTree(centers[own]).query(centers, k=1)
        nonown = ~own
        near = nonown & (d <= r_near * mnn)
        far = nonown & (d >= r_far * mnn)
        mo = C[own, g].mean()
        if near.sum() < 5 or far.sum() < 5 or mo <= 0:
            continue
        leaks.append((C[near, g].mean() - C[far, g].mean()) / mo)
    return float(np.mean(leaks)) if leaks else np.nan


def leak_curve_from_data(data, owner, r_near, r_far):
    """Mean (over replicates) spatial marker leakage per sigma, for given owner set and
    adjacency thresholds."""
    cur = []
    for reps in data:
        vals = [_leak_one(c, C, t, m, owner, r_near, r_far) for (c, C, t, m) in reps]
        vals = [v for v in vals if np.isfinite(v)]
        cur.append(float(np.mean(vals)) if vals else np.nan)
    return np.array(cur)


def fit_sigma_safe(curve, target):
    """fit_sigma that drops non-finite curve points (robust to empty/all-nan curves)."""
    ok = np.isfinite(curve)
    if ok.sum() < 2 or not np.isfinite(target):
        return float("nan"), "degenerate"
    c, g = curve[ok], SIGMA_GRID[ok]
    if target <= c[0]:
        return 0.0, "below_floor"
    if target > c[-1]:
        return float("nan"), "above_15um"
    order = np.argsort(c)
    return float(np.interp(target, c[order], g[order])), "ok"


def build_xenium():
    model, real = build_realistic_model_from_xenium(seed=config.MASTER_SEED)
    abundance = model.mean_expr.sum(axis=1)
    abundance = abundance / abundance.mean()
    return model, real, abundance


def main():
    model, real, abundance = build_xenium()
    # real spatial leakage on the whole Xenium section
    spat = spatial_leakage(model, real["coords"], real["X"], real["types"], real["median_nn_um"])
    lit = literal_leakage(model, real["X"], real["types"])
    sv = np.array(list(spat.values()))
    s_mean = float(sv.mean())
    pack_pin = real["packing_median_cells_per_mm2"]
    dens_pin = real["density_median_tx_per_cell"]
    print(f"Xenium: K={model.n_types} markers={model.n_exclusive} "
          f"median_nn={real['median_nn_um']:.2f}um pack_pin={pack_pin:.0f} dens_pin={dens_pin:.0f}")
    print(f"  real spatial leakage mean={s_mean:.4f} (n_markers={len(sv)}); literal={np.mean(list(lit.values())):.3f}")

    # pin sigma at Xenium's own representative packing/density (large synthetic field)
    data = gen_synth_data(model, pack_pin, dens_pin, abundance, GATE2_SEED + 11)
    curve = leak_curve_from_data(data, model.excl_owner, R_NEAR_BASE, R_FAR_BASE)
    bio = float(curve[0])
    sig_pin, status = fit_sigma_safe(curve, s_mean)
    print(f"  synthetic baseline (sigma=0)={bio:.4f}; Xenium-pinned sigma={sig_pin:.2f}um ({status}) "
          f"[MERFISH-pinned was {MERFISH_PINNED_UM} um]")

    # evaluate oracle assignment accuracy at the Xenium-pinned physical sigma (converged grid)
    oa_dense, na_dense = oracle_acc_at(model, DENSE_PACKING, dens_pin, sig_pin, GATE2_SEED + 21,
                                       res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    oa_sparse, na_sparse = oracle_acc_at(model, SPARSE_PACKING, dens_pin, sig_pin, GATE2_SEED + 23,
                                         res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    print(f"  oracle @ pinned sigma: dense(pack={DENSE_PACKING:.0f})={oa_dense:.3f} (naive {na_dense:.3f}); "
          f"sparse(pack={SPARSE_PACKING:.0f})={oa_sparse:.3f} (naive {na_sparse:.3f})")

    row = {
        "dataset": "xenium_breast_rep1", "n_types": model.n_types, "n_markers": int(len(sv)),
        "median_nn_um": real["median_nn_um"],
        "pin_packing_cells_per_mm2": pack_pin, "pin_density_tx_per_cell": dens_pin,
        "real_spatial_leakage_mean": s_mean,
        "real_literal_leakage_mean": float(np.mean(list(lit.values()))),
        "synthetic_baseline_leakage_sigma0": bio,
        "xenium_pinned_sigma_um": sig_pin, "sigma_fit_status": status,
        "merfish_pinned_sigma_um": MERFISH_PINNED_UM,
        "dense_packing_cells_per_mm2": DENSE_PACKING,
        "oracle_acc_dense_at_pinned": oa_dense, "naive_acc_dense_at_pinned": na_dense,
        "sparse_packing_cells_per_mm2": SPARSE_PACKING,
        "oracle_acc_sparse_at_pinned": oa_sparse, "naive_acc_sparse_at_pinned": na_sparse,
        "seed": GATE2_SEED,
    }
    df = pd.DataFrame([row])
    out = os.path.join(RESULTS, "gate2_xenium_pin.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")
    return df


if __name__ == "__main__":
    main()
