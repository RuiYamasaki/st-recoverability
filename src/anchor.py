"""Realism-anchored oracle evaluation -> results/realism_oracle.csv.

Takes each real field of view from results/realism.csv and evaluates the oracle at
the generator point whose nominal packing and density equal that FOV's measured
packing and median transcripts/cell, across the full displacement sweep. Marks the
realistic displacement sub-band (sigma <= REALISTIC_SIGMA_MAX um) and records the
generator<->reality match (KILL threshold is a >15% mismatch on median tx/cell or
packing anywhere in the generator's range; here we verify the match at the anchors).

displacement is the one axis not observable from the public per-cell tables, so the
realistic sigma band is an explicit assumption, stated as such; the full sweep is
reported so the reader sees sensitivity.
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")

SIGMAS = [0.5, 1.0, 2.0, 4.0, 8.0]       # same displacement sweep as the grid
REALISTIC_SIGMA_MAX = 2.0                 # assumed upper bound of the realistic band (um)
MATCH_TOL = 0.15                          # 15% match tolerance (locked threshold)
ANCHOR_SEED = config.MASTER_SEED + 50000


def main():
    real = pd.read_csv(os.path.join(RESULTS, "realism.csv"))
    rows = []
    for ai, fov in real.iterrows():
        packing = float(fov["packing_cells_per_mm2"])
        density = float(fov["median_tx_per_cell"])   # match generator mean to real median
        for si, sigma in enumerate(SIGMAS):
            fseed = ANCHOR_SEED + 131 * int(ai) + 7 * si
            field = build_field(packing, sigma, fseed)
            Dmax, argcell = build_oracle_maps(field)
            rec = eval_point(field, Dmax, argcell, density, fseed + 1)
            tx_err = abs(rec["realized_median_tx_per_cell"] - density) / density
            pk_err = abs(rec["realized_packing_cells_per_mm2"] - packing) / packing
            rec.update({
                "fov_id": fov["fov_id"], "dataset": fov["dataset"],
                "real_median_tx_per_cell": density,
                "real_packing_cells_per_mm2": packing,
                "real_median_nn_distance_um": float(fov["median_nn_distance_um"]),
                "sigma_um": sigma,
                "in_realistic_sigma_band": bool(sigma <= REALISTIC_SIGMA_MAX),
                "tx_match_rel_err": tx_err,
                "packing_match_rel_err": pk_err,
                "tx_match_within_15pct": bool(tx_err <= MATCH_TOL),
                "packing_match_within_15pct": bool(pk_err <= MATCH_TOL),
                "field_seed": fseed,
            })
            rows.append(rec)
        print(f"anchored {fov['fov_id']:>26s} (pack={packing:.0f}, dens={density:.0f}) done")

    cols = [
        "fov_id", "dataset", "real_packing_cells_per_mm2", "real_median_tx_per_cell",
        "real_median_nn_distance_um", "sigma_um", "in_realistic_sigma_band",
        "oracle_acc", "naive_acc", "oracle_minus_naive", "oracle_ge_naive",
        "assign_ari_vs_truecells_oracle", "assign_ari_vs_truecells_naive",
        "profile_ari_vs_truetype", "naive_profile_ari_vs_truetype", "perfect_ari_vs_truetype",
        "same_type_error_frac_oracle",
        "realized_median_tx_per_cell", "realized_packing_cells_per_mm2",
        "tx_match_rel_err", "packing_match_rel_err",
        "tx_match_within_15pct", "packing_match_within_15pct",
        "n_cells", "n_interior_cells", "grid", "dx_um", "r_mean_um", "L_um",
        "n_interior_transcripts", "field_seed",
    ]
    df = pd.DataFrame(rows)[cols]
    out = os.path.join(RESULTS, "realism_oracle.csv")
    df.to_csv(out, index=False)
    band = df[df.in_realistic_sigma_band]
    print(f"\nwrote {out}  ({len(df)} rows; {len(band)} in realistic sigma band)")
    print(f"15% match holds at every anchor: "
          f"tx={bool(df.tx_match_within_15pct.all())}, "
          f"packing={bool(df.packing_match_within_15pct.all())}")
    print(f"oracle_acc within realistic band (sigma<= {REALISTIC_SIGMA_MAX}um): "
          f"[{band.oracle_acc.min():.3f}, {band.oracle_acc.max():.3f}]")
    print(df[["fov_id", "sigma_um", "oracle_acc", "naive_acc",
              "tx_match_rel_err", "packing_match_rel_err"]].to_string(index=False))
    return df


if __name__ == "__main__":
    main()
