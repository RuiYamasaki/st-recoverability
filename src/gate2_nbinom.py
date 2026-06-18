"""Gate 2 Task 3: non-idealized emission (within-type variance).

Replace the per-type-mean Poisson emission with a negative-binomial (gamma-Poisson) model:
each cell's per-gene counts are Poisson(mu * f_g) with a per-cell, per-gene gamma factor
f_g of mean 1 and variance phi_g, where phi_g is the per-gene NB dispersion estimated by
within-type method of moments from the real Xenium counts (Var = mu + phi*mu^2). The real
spatial marker leakage is unchanged (it is a property of the real data); only the
synthetic emission changes, so sigma is re-pinned under the non-idealized synthetic curve
and the oracle is re-evaluated.

Reports: the dispersion used, the sigma shift and dense-regime accuracy shift versus the
per-type-mean (Poisson) model, and a reduced-grid frontier under NB emission with
monotonicity-violation counts and the oracle-at-or-above-naive check.

Output: results/gate2_nbinom.csv. Seeds: GATE2_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import build_field  # noqa: E402
from oracle import build_oracle_maps  # noqa: E402
from evaluate import eval_point  # noqa: E402
from gate1_leakage_anchor import spatial_leakage  # noqa: E402
from gate2_pin import (build_xenium, gen_synth_data, leak_curve_from_data, fit_sigma_safe,  # noqa: E402
                       oracle_acc_at, DENSE_PACKING, SPARSE_PACKING, R_NEAR_BASE, R_FAR_BASE,
                       ORACLE_RES_CELL, ORACLE_GRID_MAX)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
GATE2_SEED = config.MASTER_SEED + 420000

FRONTIER_PACKINGS = [2575.0, 6000.0, 13625.0]
FRONTIER_SIGMAS = [1.0, 2.0, 4.0, 8.0]


def _frontier(model, abundance, density, emission, seed):
    rows = []
    for pi, packing in enumerate(FRONTIER_PACKINGS):
        for si, sigma in enumerate(FRONTIER_SIGMAS):
            sd = seed + 1000 * (pi + 1) + 17 * (si + 1)
            f = build_field(packing, sigma, sd, model=model)
            D, A = build_oracle_maps(f)
            r = eval_point(f, D, A, density, sd + 1, emission=emission)
            rows.append({"emission": emission, "packing_cells_per_mm2": packing, "sigma_um": sigma,
                         "oracle_acc": r["oracle_acc"], "naive_acc": r["naive_acc"],
                         "oracle_minus_naive": r["oracle_minus_naive"], "oracle_ge_naive": r["oracle_ge_naive"]})
    return rows


def _mono(rows):
    df = pd.DataFrame(rows)
    vs = sum(int((np.diff(s.sort_values("sigma_um").oracle_acc.values) > 1e-3).sum())
             for _, s in df.groupby("packing_cells_per_mm2"))
    vp = sum(int((np.diff(s.sort_values("packing_cells_per_mm2").oracle_acc.values) > 1e-3).sum())
             for _, s in df.groupby("sigma_um"))
    return vs, vp


def main():
    model, real, abundance = build_xenium()
    dens = real["density_median_tx_per_cell"]
    pack_pin = real["packing_median_cells_per_mm2"]
    phi = model.dispersion
    print(f"NB dispersion phi: median={np.median(phi):.3f} (Var = mu + phi*mu^2); "
          f"IQR=[{np.percentile(phi,25):.3f}, {np.percentile(phi,75):.3f}]")

    # real spatial leakage (unchanged) and re-pin under both emissions at the same seed
    s_mean = float(np.mean(list(spatial_leakage(model, real["coords"], real["X"],
                                                real["types"], real["median_nn_um"]).values())))
    data_pois = gen_synth_data(model, pack_pin, dens, abundance, GATE2_SEED + 11, emission="poisson")
    data_nb = gen_synth_data(model, pack_pin, dens, abundance, GATE2_SEED + 11, emission="nbinom")
    sig_pois, _ = fit_sigma_safe(leak_curve_from_data(data_pois, model.excl_owner, R_NEAR_BASE, R_FAR_BASE), s_mean)
    sig_nb, _ = fit_sigma_safe(leak_curve_from_data(data_nb, model.excl_owner, R_NEAR_BASE, R_FAR_BASE), s_mean)
    print(f"re-pinned sigma: poisson={sig_pois:.2f}um  nbinom={sig_nb:.2f}um  (shift {sig_nb-sig_pois:+.2f}um)")

    # dense + sparse oracle accuracy at each emission's own pinned sigma (converged grid)
    oa_d_pois, na_d_pois = oracle_acc_at(model, DENSE_PACKING, dens, sig_pois, GATE2_SEED + 21,
                                         emission="poisson", res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    oa_d_nb, na_d_nb = oracle_acc_at(model, DENSE_PACKING, dens, sig_nb, GATE2_SEED + 23,
                                     emission="nbinom", res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    oa_s_pois, _ = oracle_acc_at(model, SPARSE_PACKING, dens, sig_pois, GATE2_SEED + 25,
                                 emission="poisson", res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    oa_s_nb, _ = oracle_acc_at(model, SPARSE_PACKING, dens, sig_nb, GATE2_SEED + 27,
                               emission="nbinom", res_cell=ORACLE_RES_CELL, grid_max=ORACLE_GRID_MAX)
    print(f"dense oracle: poisson={oa_d_pois:.3f}  nbinom={oa_d_nb:.3f}  (shift {oa_d_nb-oa_d_pois:+.3f})")
    print(f"sparse oracle: poisson={oa_s_pois:.3f}  nbinom={oa_s_nb:.3f}")

    # frontier under NB emission
    fr_nb = _frontier(model, abundance, dens, "nbinom", GATE2_SEED + 100)
    fr_pois = _frontier(model, abundance, dens, "poisson", GATE2_SEED + 200)
    vs_nb, vp_nb = _mono(fr_nb)
    nb_df = pd.DataFrame(fr_nb)
    print(f"NB frontier: oracle range [{nb_df.oracle_acc.min():.3f}, {nb_df.oracle_acc.max():.3f}]; "
          f"mono(sig,pack)=({vs_nb},{vp_nb}); oracle>=naive all={bool(nb_df.oracle_ge_naive.all())} "
          f"(min gap {nb_df.oracle_minus_naive.min():+.3f})")

    summary = {
        "nb_dispersion_median": float(np.median(phi)),
        "real_spatial_leakage_mean": s_mean,
        "sigma_poisson_um": sig_pois, "sigma_nbinom_um": sig_nb, "sigma_shift_um": sig_nb - sig_pois,
        "dense_oracle_poisson": oa_d_pois, "dense_oracle_nbinom": oa_d_nb,
        "dense_oracle_shift": oa_d_nb - oa_d_pois,
        "dense_naive_nbinom": na_d_nb,
        "sparse_oracle_poisson": oa_s_pois, "sparse_oracle_nbinom": oa_s_nb,
        "nb_frontier_oracle_min": float(nb_df.oracle_acc.min()), "nb_frontier_oracle_max": float(nb_df.oracle_acc.max()),
        "nb_frontier_mono_viol_sigma": vs_nb, "nb_frontier_mono_viol_packing": vp_nb,
        "nb_frontier_oracle_ge_naive_all": bool(nb_df.oracle_ge_naive.all()),
        "nb_frontier_min_gap": float(nb_df.oracle_minus_naive.min()),
        "dense_packing": DENSE_PACKING, "seed": GATE2_SEED,
    }
    pd.DataFrame([summary]).to_csv(os.path.join(RESULTS, "gate2_nbinom.csv"), index=False)
    pd.DataFrame(fr_nb + fr_pois).to_csv(os.path.join(RESULTS, "gate2_nbinom_frontier.csv"), index=False)
    print(f"\nwrote results/gate2_nbinom.csv and results/gate2_nbinom_frontier.csv")
    return summary


if __name__ == "__main__":
    main()
