"""Gate 3 shared machinery: a principled, non-circular marker-validity test and
admissible-set sigma pinning. Reused by the breast (Task 1), cell-typing (Task 2) and
lung-replication (Task 3) scripts.

Validity test (defined WITHOUT reference to the resulting sigma): a candidate marker is
admissible only if its real adjacent-minus-distant displacement signal exceeds a
position-permutation null at p < 0.05. The null shuffles the marker's per-cell counts
across the non-owner cells (equivalent to scrambling positions relative to counts) and
recomputes the near-minus-far contrast on the FIXED real adjacency geometry, n_perm times.
A marker whose signal is indistinguishable from this null is biology-diluted (no spatial
displacement gradient) and is rejected. Marker selection never uses the sigma it yields.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from expression import _exclusivity  # noqa: E402
from gate2_pin import (gen_synth_data, leak_curve_from_data, fit_sigma_safe, oracle_acc_at,  # noqa: E402
                       R_NEAR_BASE, R_FAR_BASE, DENSE_PACKING, SPARSE_PACKING,
                       ORACLE_RES_CELL, ORACLE_GRID_MAX)

N_PERM = 1000
P_ADMIT = 0.05
CANDIDATE_CUTOFF = 0.5     # broad candidate pool; the validity test (not this cutoff) decides admissibility
LOOSE_LO, LOOSE_HI = 0.6, 0.7   # the markers that destabilised sigma at cutoff 0.6 (exclusivity in [0.6, 0.7))


def _near_far(coords, types, owner_g, mnn, r_near, r_far):
    own = types == owner_g
    d, _ = cKDTree(coords[own]).query(coords, k=1)
    nonown = ~own
    return own, nonown, (nonown & (d <= r_near * mnn)), (nonown & (d >= r_far * mnn))


def marker_validity(real, owner, candidates, seed, r_near=R_NEAR_BASE, r_far=R_FAR_BASE,
                    n_perm=N_PERM):
    """Per-candidate-marker permutation validity test. Returns a DataFrame with the real
    signal, the permutation p-value, and admissibility."""
    coords, X, types, mnn = real["coords"], real["X"], real["types"], real["median_nn_um"]
    rng = np.random.default_rng(seed)
    rows = []
    for g in candidates:
        tau = int(owner[g])
        own, nonown, near, far = _near_far(coords, types, tau, mnn, r_near, r_far)
        cg = X[:, g]
        mo = cg[own].mean() if own.sum() >= 3 else 0.0
        if own.sum() < 3 or near.sum() < 5 or far.sum() < 5 or mo <= 0:
            rows.append({"gene": int(g), "owner": tau, "real_signal": np.nan, "pval": np.nan,
                         "n_near": int(near.sum()), "n_far": int(far.sum()), "admissible": False,
                         "reason": "too_few_cells"})
            continue
        real_sig = (cg[near].mean() - cg[far].mean()) / mo
        pool = cg[nonown]
        near_in = near[nonown]
        far_in = far[nonown]
        null = np.empty(n_perm)
        for p in range(n_perm):
            sh = rng.permutation(pool)
            null[p] = (sh[near_in].mean() - sh[far_in].mean()) / mo
        pval = (np.sum(null >= real_sig) + 1) / (n_perm + 1)
        rows.append({"gene": int(g), "owner": tau, "real_signal": float(real_sig),
                     "pval": float(pval), "n_near": int(near.sum()), "n_far": int(far.sum()),
                     "admissible": bool(pval < P_ADMIT), "reason": "ok"})
    return pd.DataFrame(rows)


def admissible_owner(model, validity_df):
    """Owner array restricted to admissible markers (others -> -1)."""
    owner = np.full(model.n_genes, -1, dtype=int)
    for _, r in validity_df[validity_df.admissible].iterrows():
        owner[int(r.gene)] = int(r.owner)
    return owner


def _marker_components(real, owner_adm, r_near, r_far):
    """Per admissible marker: near and far non-owner counts and owner mean (for bootstrap)."""
    coords, X, types, mnn = real["coords"], real["X"], real["types"], real["median_nn_um"]
    comp = {}
    for g in np.where(owner_adm >= 0)[0]:
        own, nonown, near, far = _near_far(coords, types, int(owner_adm[g]), mnn, r_near, r_far)
        cg = X[:, g]
        mo = cg[own].mean()
        if near.sum() >= 5 and far.sum() >= 5 and mo > 0:
            comp[int(g)] = (cg[near], cg[far], mo)
    return comp


def pin_admissible(model, real, abundance, owner_adm, seed, density=None,
                   r_near=R_NEAR_BASE, r_far=R_FAR_BASE, n_boot=2000):
    """Re-pin sigma on the admissible marker set with a bootstrap CI, and report the
    dense and sparse oracle accuracy at the sigma point and CI ends (converged grid)."""
    if density is None:
        density = real["density_median_tx_per_cell"]
    pack_pin = real["packing_median_cells_per_mm2"]
    comp = _marker_components(real, owner_adm, r_near, r_far)
    markers = list(comp.keys())
    point_leak = float(np.mean([(comp[g][0].mean() - comp[g][1].mean()) / comp[g][2] for g in markers]))

    rng = np.random.default_rng(seed + 1)
    boot = np.empty(n_boot)
    for b in range(n_boot):
        gs = rng.choice(markers, size=len(markers), replace=True)
        vals = []
        for g in gs:
            near, far, mo = comp[g]
            nb = near[rng.integers(0, near.size, near.size)].mean()
            fb = far[rng.integers(0, far.size, far.size)].mean()
            vals.append((nb - fb) / mo)
        boot[b] = np.mean(vals)
    leak_lo, leak_hi = np.percentile(boot, [2.5, 97.5])

    data = gen_synth_data(model, pack_pin, density, abundance, seed + 11)
    curve = leak_curve_from_data(data, owner_adm, r_near, r_far)
    sig_pt, st = fit_sigma_safe(curve, point_leak)
    sig_lo, _ = fit_sigma_safe(curve, leak_lo)
    sig_hi, _ = fit_sigma_safe(curve, leak_hi)

    def oa(packing, sigma, sd):
        return oracle_acc_at(model, packing, density, sigma, sd,
                             res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    dense_pt, naive_pt = oa(DENSE_PACKING, sig_pt, seed + 21)
    dense_lo, _ = oa(DENSE_PACKING, sig_lo, seed + 23)
    dense_hi, _ = oa(DENSE_PACKING, sig_hi, seed + 25)
    sparse_pt, _ = oa(SPARSE_PACKING, sig_pt, seed + 27)
    return {
        "n_admissible_markers": len(markers),
        "real_leakage": point_leak, "leakage_ci_lo": float(leak_lo), "leakage_ci_hi": float(leak_hi),
        "sigma_point_um": sig_pt, "sigma_fit_status": st,
        "sigma_ci_lo_um": sig_lo, "sigma_ci_hi_um": sig_hi,
        "pin_packing": pack_pin, "density": density, "dense_packing": DENSE_PACKING,
        "dense_oracle_point": dense_pt, "naive_dense_point": naive_pt,
        "dense_oracle_at_sigma_ci_lo": dense_lo, "dense_oracle_at_sigma_ci_hi": dense_hi,
        "sparse_oracle_point": sparse_pt, "n_boot": n_boot,
    }


def candidate_markers(model, cutoff=CANDIDATE_CUTOFF):
    owner = _exclusivity(model.composition, cutoff)
    return owner, list(np.where(owner >= 0)[0])
