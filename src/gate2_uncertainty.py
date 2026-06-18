"""Gate 2 Task 2: statistical confidence interval and design-sensitivity range on the
Xenium-pinned displacement sigma, and the dense-regime oracle accuracy at the combined
low (optimistic) and high ends.

Statistical CI: bootstrap the real spatial-leakage estimate (B resamples) over markers
(resample the marker set with replacement) AND over cells (resample the near and far cell
sets with replacement), propagated to a CI on sigma through the synthetic leakage curve.

Design-sensitivity: vary the leakage statistic's design choices over a grid:
  - adjacency (near) threshold R_NEAR in {1.0, 1.5, 2.0} (units of median NN distance),
  - distant (far) threshold R_FAR in {3.0, 4.0, 5.0},
  - marker-selection exclusivity cutoff in {0.6, 0.7, 0.8},
recompute the real leakage, the synthetic curve, and the pinned sigma for each, and
report the range of pinned sigma.

Combined band = [min(stat CI low, design min), max(stat CI high, design max)]. Dense and
sparse oracle accuracy are reported at the combined low and high sigma ends. Output:
results/gate2_sigma_uncertainty.csv. Seeds: GATE2_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import _exclusivity  # noqa: E402
from gate1_leakage_anchor import oracle_acc_at, SIGMA_GRID  # noqa: E402
from gate2_pin import (build_xenium, gen_synth_data, leak_curve_from_data, fit_sigma_safe,  # noqa: E402
                       DENSE_PACKING, SPARSE_PACKING, MERFISH_PINNED_UM, R_NEAR_BASE, R_FAR_BASE,
                       ORACLE_RES_CELL, ORACLE_GRID_MAX, GATE2_SEED as PIN_SEED)


def _oa(model, packing, density, sigma, seed):
    return oracle_acc_at(model, packing, density, sigma, seed,
                         res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
GATE2_SEED = config.MASTER_SEED + 410000

N_BOOT = 2000
THR_BASE = 0.7
R_NEAR_GRID = [1.0, 1.5, 2.0]
R_FAR_GRID = [3.0, 4.0, 5.0]
THR_GRID = [0.6, 0.7, 0.8]


def _marker_distance_cache(real, owner):
    """Per marker g: distance of every cell to its nearest owner-type cell, plus the
    owner mask and the gene-g counts. Computed once; reused for any R thresholds."""
    coords, X, types = real["coords"], real["X"], real["types"]
    cache = {}
    for g in np.where(owner >= 0)[0]:
        tau = owner[g]
        own = types == tau
        if own.sum() < 3:
            continue
        d, _ = cKDTree(coords[own]).query(coords, k=1)
        cache[int(g)] = (d, own, X[:, g])
    return cache


def _leak_from_cache(cache, mnn, r_near, r_far):
    """Per-marker spatial leakage from the distance cache, for given R thresholds."""
    out = {}
    for g, (d, own, cg) in cache.items():
        nonown = ~own
        near = nonown & (d <= r_near * mnn)
        far = nonown & (d >= r_far * mnn)
        mo = cg[own].mean()
        if near.sum() < 5 or far.sum() < 5 or mo <= 0:
            continue
        out[g] = (cg[near].mean() - cg[far].mean()) / mo
    return out


def main():
    model, real, abundance = build_xenium()
    mnn = real["median_nn_um"]
    pack_pin = real["packing_median_cells_per_mm2"]
    dens_pin = real["density_median_tx_per_cell"]
    rng = np.random.default_rng(GATE2_SEED)

    # synthetic data once per sigma (same seed as gate2_pin so the point sigma matches;
    # reused for the CI curve and every design combo)
    data = gen_synth_data(model, pack_pin, dens_pin, abundance, PIN_SEED + 11)

    # ---- statistical bootstrap (baseline design) ----
    owner_base = _exclusivity(model.composition, THR_BASE)
    cache = _marker_distance_cache(real, owner_base)
    comp = {}  # per marker: near_counts, far_counts, owner_mean
    for g, (d, own, cg) in cache.items():
        nonown = ~own
        near = cg[nonown & (d <= R_NEAR_BASE * mnn)]
        far = cg[nonown & (d >= R_FAR_BASE * mnn)]
        mo = cg[own].mean()
        if near.size >= 5 and far.size >= 5 and mo > 0:
            comp[g] = (near, far, mo)
    markers = list(comp.keys())
    point_leak = float(np.mean([(comp[g][0].mean() - comp[g][1].mean()) / comp[g][2] for g in markers]))

    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        gs = rng.choice(markers, size=len(markers), replace=True)
        vals = []
        for g in gs:
            near, far, mo = comp[g]
            nb = near[rng.integers(0, near.size, near.size)].mean()
            fb = far[rng.integers(0, far.size, far.size)].mean()
            vals.append((nb - fb) / mo)
        boot[b] = np.mean(vals)
    leak_lo, leak_hi = np.percentile(boot, [2.5, 97.5])

    curve_base = leak_curve_from_data(data, owner_base, R_NEAR_BASE, R_FAR_BASE)
    sig_point, _ = fit_sigma_safe(curve_base, point_leak)
    sig_stat_lo, _ = fit_sigma_safe(curve_base, leak_lo)   # lower leakage -> lower sigma (optimistic)
    sig_stat_hi, _ = fit_sigma_safe(curve_base, leak_hi)
    print(f"point leakage={point_leak:.4f}  95% CI=[{leak_lo:.4f}, {leak_hi:.4f}] (B={N_BOOT})")
    print(f"sigma point={sig_point:.2f}um  statistical 95% CI=[{sig_stat_lo:.2f}, {sig_stat_hi:.2f}]um")

    # ---- design-sensitivity grid ----
    design_sigmas = []
    caches_by_thr = {}
    for thr in THR_GRID:
        caches_by_thr[thr] = (_exclusivity(model.composition, thr),)
        caches_by_thr[thr] += (_marker_distance_cache(real, caches_by_thr[thr][0]),)
    for thr in THR_GRID:
        owner_t, cache_t = caches_by_thr[thr]
        for rn in R_NEAR_GRID:
            for rf in R_FAR_GRID:
                if rf <= rn + 0.5:
                    continue
                real_leak = np.mean(list(_leak_from_cache(cache_t, mnn, rn, rf).values()))
                curve = leak_curve_from_data(data, owner_t, rn, rf)
                sig, status = fit_sigma_safe(curve, real_leak)
                design_sigmas.append({"excl_threshold": thr, "r_near": rn, "r_far": rf,
                                      "real_leakage": float(real_leak), "sigma_um": sig,
                                      "status": status, "n_markers": int((owner_t >= 0).sum())})
    dsig = np.array([d["sigma_um"] for d in design_sigmas if np.isfinite(d["sigma_um"])])
    design_lo, design_hi = float(dsig.min()), float(dsig.max())
    print(f"design-sensitivity sigma range over {len(design_sigmas)} design choices: "
          f"[{design_lo:.2f}, {design_hi:.2f}]um")

    # ---- combined band and dense/sparse oracle accuracy at the ends (converged grid) ----
    combined_lo = float(min(sig_stat_lo, design_lo))
    combined_hi = float(max(sig_stat_hi, design_hi))
    oa_dense_pt, _ = _oa(model, DENSE_PACKING, dens_pin, sig_point, GATE2_SEED + 55)
    # statistical-CI ends (pre-registered threshold-0.7 markers)
    oa_dense_slo, _ = _oa(model, DENSE_PACKING, dens_pin, sig_stat_lo, GATE2_SEED + 41)
    oa_dense_shi, _ = _oa(model, DENSE_PACKING, dens_pin, sig_stat_hi, GATE2_SEED + 43)
    # combined ends (statistical CI + design sensitivity)
    oa_dense_lo, na_dense_lo = _oa(model, DENSE_PACKING, dens_pin, combined_lo, GATE2_SEED + 51)
    oa_dense_hi, na_dense_hi = _oa(model, DENSE_PACKING, dens_pin, combined_hi, GATE2_SEED + 53)
    oa_sparse_pt, _ = _oa(model, SPARSE_PACKING, dens_pin, sig_point, GATE2_SEED + 56)
    oa_sparse_lo, _ = _oa(model, SPARSE_PACKING, dens_pin, combined_lo, GATE2_SEED + 57)
    oa_sparse_hi, _ = _oa(model, SPARSE_PACKING, dens_pin, combined_hi, GATE2_SEED + 59)
    print(f"combined sigma band=[{combined_lo:.2f}, {combined_hi:.2f}]um  "
          f"(statistical CI [{sig_stat_lo:.2f}, {sig_stat_hi:.2f}], design [{design_lo:.2f}, {design_hi:.2f}])")
    print(f"DENSE oracle (converged grid): point({sig_point:.2f})={oa_dense_pt:.3f}; "
          f"stat-CI ends [{oa_dense_shi:.3f}@hi, {oa_dense_slo:.3f}@lo]; "
          f"combined ends [{oa_dense_hi:.3f}@hi, {oa_dense_lo:.3f}@lo]")
    print(f"SPARSE oracle: point={oa_sparse_pt:.3f}; combined ends [{oa_sparse_hi:.3f}, {oa_sparse_lo:.3f}]")

    rows = [{
        "leakage_point": point_leak, "leakage_ci_lo": leak_lo, "leakage_ci_hi": leak_hi, "n_boot": N_BOOT,
        "sigma_point_um": sig_point, "sigma_stat_ci_lo_um": sig_stat_lo, "sigma_stat_ci_hi_um": sig_stat_hi,
        "sigma_design_lo_um": design_lo, "sigma_design_hi_um": design_hi,
        "sigma_combined_lo_um": combined_lo, "sigma_combined_hi_um": combined_hi,
        "dense_packing": DENSE_PACKING, "sparse_packing": SPARSE_PACKING,
        "oracle_res_cell": ORACLE_RES_CELL, "oracle_grid_max": ORACLE_GRID_MAX,
        "oracle_dense_at_point": oa_dense_pt,
        "oracle_dense_at_stat_ci_lo": oa_dense_slo, "oracle_dense_at_stat_ci_hi": oa_dense_shi,
        "oracle_dense_at_combined_lo": oa_dense_lo, "oracle_dense_at_combined_hi": oa_dense_hi,
        "naive_dense_at_combined_lo": na_dense_lo, "naive_dense_at_combined_hi": na_dense_hi,
        "oracle_sparse_at_point": oa_sparse_pt,
        "oracle_sparse_at_combined_lo": oa_sparse_lo, "oracle_sparse_at_combined_hi": oa_sparse_hi,
        "merfish_pinned_sigma_um": MERFISH_PINNED_UM, "seed": GATE2_SEED,
    }]
    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "gate2_sigma_uncertainty.csv")
    df.to_csv(out, index=False)
    pd.DataFrame(design_sigmas).to_csv(os.path.join(RESULTS, "gate2_design_grid.csv"), index=False)
    print(f"\nwrote {out} and results/gate2_design_grid.csv")
    return df


if __name__ == "__main__":
    main()
