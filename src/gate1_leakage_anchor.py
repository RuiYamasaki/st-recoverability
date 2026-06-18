"""Gate 1 Task 2: pin the displacement sigma to a measurable real contamination
statistic (marker leakage) instead of assuming it.

Two leakage statistics on the exclusive-marker genes (owner type tau), both
computable on real and synthetic data with no ground truth:

  - literal leakage: fraction of a marker's signal residing in non-owner-type cells.
    Reported as a DIAGNOSTIC. It is biology-dominated (real cell types co-express the
    marker) and the synthetic model reproduces it at sigma=0 by construction, so it
    cannot identify the displacement. We show this explicitly.

  - spatial leakage (the one used to pin sigma): the EXCESS marker signal in non-owner
    cells that are spatially ADJACENT to owner cells, over the level in DISTANT
    non-owner cells, normalised by the owner level. Biological co-expression is
    spatially uniform and cancels in the near-minus-far contrast; only displacement
    spreads a marker into a cell's immediate neighbours, so this isolates the spatial
    contamination and rises monotonically with sigma.

Real statistics use the dataset's cell types and the segmented cell-by-gene matrix and
centroids (MERFISH Moffitt 2018). Synthetic statistics use the naive nearest-nucleus
assignment of displaced transcripts. The fitted sigma reproduces the real spatial
leakage; the 25-75 bracket comes from the per-marker spatial-leakage distribution. The
oracle assignment accuracy is then reported at the pinned sigma (point + bracket) at
each real anchor's packing and density.

Output: results/gate1_leakage.csv. Seeds: GATE1_SEED + offsets, recorded.
Note: the MERFISH X matrix is volume-normalised expression, so "counts" here are
expression mass (stated honestly).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import build_realistic_model_from_merfish  # noqa: E402
from generator import build_field, generate_transcripts  # noqa: E402
from oracle import build_oracle_maps, naive_assign  # noqa: E402
from evaluate import eval_point  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")

GATE1_SEED = config.MASTER_SEED + 200000
SIGMA_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0])
R_NEAR, R_FAR = 1.5, 4.0          # adjacency thresholds in units of median nearest-neighbour distance
SIGMA_MAX_PLAUSIBLE = 15.0


def _load_merfish_typed(model):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(DATA, "merfish_moffitt.h5ad"))
    X = np.asarray(a.X, dtype=np.float64)
    coords = np.asarray(a.obsm["spatial"], dtype=np.float64)
    name_to_idx = {str(n): i for i, n in enumerate(model.type_names)}
    cell_type = np.array([name_to_idx.get(str(c), -1) for c in a.obs["cell_class"].astype(str)])
    return X, coords, cell_type


def _median_nn(coords):
    d, _ = cKDTree(coords).query(coords, k=2)
    return float(np.median(d[:, 1]))


def literal_leakage(model, counts, types):
    """Diagnostic: fraction of each marker's signal in non-owner-type cells."""
    out = {}
    for g in np.where(model.excl_owner >= 0)[0]:
        tau = model.excl_owner[g]
        total = counts[:, g].sum()
        if total > 0:
            out[int(g)] = float(1.0 - counts[types == tau, g].sum() / total)
    return out


def spatial_leakage(model, coords, counts, types, median_nn):
    """Per-marker spatial leakage: (mean signal in adjacent non-owner cells minus mean
    in distant non-owner cells) / mean signal in owner cells."""
    out = {}
    for g in np.where(model.excl_owner >= 0)[0]:
        tau = model.excl_owner[g]
        own = types == tau
        if own.sum() < 3:
            continue
        mean_owner = counts[own, g].mean()
        if mean_owner <= 0:
            continue
        d, _ = cKDTree(coords[own]).query(coords, k=1)
        nonown = ~own
        near = nonown & (d <= R_NEAR * median_nn)
        far = nonown & (d >= R_FAR * median_nn)
        if near.sum() < 5 or far.sum() < 5:
            continue
        out[int(g)] = float((counts[near, g].mean() - counts[far, g].mean()) / mean_owner)
    return out


def _counts_matrix(n_cells, cell_ids, gene, n_genes):
    C = np.zeros((n_cells, n_genes), dtype=np.float64)
    ok = (cell_ids >= 0) & (cell_ids < n_cells)
    np.add.at(C, (cell_ids[ok], gene[ok]), 1.0)
    return C


def synthetic_spatial_leakage(model, packing, density, sigma, seed, abundance):
    f = build_field(packing, sigma, seed, model=model)
    tx = generate_transcripts(f, density, seed + 1, abundance=abundance)
    na = naive_assign(f, tx.obs_xy)
    C = _counts_matrix(f.n_cells, na, tx.gene, model.n_genes)
    return float(np.mean(list(spatial_leakage(
        model, f.centers, C, f.types, _median_nn(f.centers)).values())))


def synthetic_curve(model, packing, density, seed, abundance):
    return np.array([synthetic_spatial_leakage(model, packing, density, s, seed + 31 * i, abundance)
                     for i, s in enumerate(SIGMA_GRID)])


def fit_sigma(curve, target):
    """Invert the monotone curve to sigma. status: ok / below_floor / above_15um."""
    if target <= curve[0]:
        return 0.0, "below_floor"
    if target > curve[-1]:
        return float("nan"), "above_15um"
    order = np.argsort(curve)
    return float(np.interp(target, curve[order], SIGMA_GRID[order])), "ok"


def oracle_acc_at(model, packing, density, sigma, seed, emission="poisson",
                  res_cell=6.0, grid_max=768):
    if not np.isfinite(sigma):
        return float("nan"), float("nan")
    f = build_field(packing, sigma, seed, model=model, res_cell=res_cell, grid_max=grid_max)
    D, A = build_oracle_maps(f)
    r = eval_point(f, D, A, density, seed + 1, emission=emission)  # uniform abundance: matches the sweep
    return float(r["oracle_acc"]), float(r["naive_acc"])


def main():
    model = build_realistic_model_from_merfish()
    abundance = model.mean_expr.sum(axis=1)
    abundance = abundance / abundance.mean()
    X, coords, cell_type = _load_merfish_typed(model)
    mnn = _median_nn(coords)

    # --- real statistics on the MERFISH section (sub-200um FOVs are too sparse for the
    #     per-marker owner/adjacency split, so the whole imaged section is the FOV) ---
    lit = literal_leakage(model, X, cell_type)
    spat = spatial_leakage(model, coords, X, cell_type, mnn)
    spat_vals = np.array(list(spat.values()))
    lit_mean = float(np.mean(list(lit.values())))
    s_mean = float(spat_vals.mean())
    s_p25, s_p50, s_p75 = [float(np.percentile(spat_vals, q)) for q in (25, 50, 75)]
    print(f"real (MERFISH section, n_markers={len(spat)}):")
    print(f"  literal leakage mean = {lit_mean:.3f}  (biology-dominated diagnostic)")
    print(f"  spatial leakage mean = {s_mean:.4f}  p25={s_p25:.4f} p50={s_p50:.4f} p75={s_p75:.4f}")

    # --- pin the PHYSICAL displacement sigma once, from MERFISH (where leakage is
    #     measured), at MERFISH geometry. sigma is a physical micron quantity. ---
    real = pd.read_csv(os.path.join(RESULTS, "realism.csv"))
    mer = real[real.dataset == "merfish_hypothal"]
    mer_pack, mer_dens = float(mer.packing_cells_per_mm2.mean()), float(mer.median_tx_per_cell.mean())
    pin_curve = synthetic_curve(model, mer_pack, mer_dens, GATE1_SEED + 9000, abundance)
    bio = float(pin_curve[0])
    sig_pt, st_pt = fit_sigma(pin_curve, s_mean)
    sig_lo, st_lo = fit_sigma(pin_curve, s_p25)
    sig_hi, st_hi = fit_sigma(pin_curve, s_p75)
    print(f"  synthetic baseline (sigma=0) = {bio:.4f}")
    print(f"  data-pinned sigma = {sig_pt:.2f} um ({st_pt}); 25-75 bracket = "
          f"[{sig_lo:.2f} ({st_lo}), {sig_hi:.2f} ({st_hi})] um")

    # --- evaluate the oracle at that same physical sigma (point + bracket) at each
    #     real anchor's packing and density ---
    anchors = {
        "merfish_hypothal": real[real.dataset == "merfish_hypothal"],
        "xenium_breast": real[real.dataset == "xenium_breast"],
    }
    rows = []
    for ds, sub in anchors.items():
        packing = float(sub.packing_cells_per_mm2.mean())
        density = float(sub.median_tx_per_cell.mean())
        seed = GATE1_SEED + (1 if ds == "merfish_hypothal" else 2) * 50000
        oa_pt, na_pt = oracle_acc_at(model, packing, density, sig_pt, seed + 5)
        oa_lo, _ = oracle_acc_at(model, packing, density, sig_lo, seed + 7)
        oa_hi, _ = oracle_acc_at(model, packing, density, sig_hi, seed + 9)
        leak_source = ("spatial marker leakage measured on the MERFISH section" if ds == "merfish_hypothal"
                       else "MERFISH-measured physical sigma applied at Xenium packing/density; "
                            "Xenium cell-by-gene matrix not downloaded so leakage not independently measured")
        rows.append({
            "anchor_dataset": ds, "packing_cells_per_mm2": packing, "density_tx_per_cell": density,
            "real_spatial_leakage_mean": s_mean, "real_spatial_leakage_p25": s_p25,
            "real_spatial_leakage_p75": s_p75, "real_literal_leakage_mean": lit_mean,
            "synthetic_baseline_leakage_sigma0": bio,
            "sigma_pinned_um": sig_pt, "sigma_fit_status": st_pt,
            "sigma_bracket_lo_um": sig_lo, "sigma_bracket_hi_um": sig_hi,
            "oracle_acc_at_pinned": oa_pt, "naive_acc_at_pinned": na_pt,
            "oracle_acc_at_bracket_lo_sigma": oa_lo, "oracle_acc_at_bracket_hi_sigma": oa_hi,
            "leakage_source": leak_source, "median_nn_um_real": mnn, "seed": seed,
        })
        print(f"  [{ds}] pack={packing:.0f} dens={density:.0f} | oracle@pinned({sig_pt:.2f}um)={oa_pt:.3f} "
              f"naive={na_pt:.3f} | oracle across sigma bracket=[{oa_hi:.3f},{oa_lo:.3f}]")

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "gate1_leakage.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out} ({len(df)} anchors)")
    return df


if __name__ == "__main__":
    main()
