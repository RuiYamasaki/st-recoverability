"""Gate 1 Task 1: re-run the (density x packing x sigma) sweep under realistic,
overlapping per-type expression (MERFISH cell-class profiles) instead of the Gate 0
disjoint-marker model. Writes results/gate1_sweep.csv and the expression-overlap
metrics to results/gate1_expression_meta.json.

Axes are identical to Gate 0 so the frontiers are directly comparable. Seeds are
distinct (GATE1_SEED base) and recorded per row. The oracle map depends only on
(packing, sigma), so it is built once per (packing, sigma) and reused across densities.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import build_realistic_model_from_merfish, overlap_metrics  # noqa: E402
from generator import build_field  # noqa: E402
from oracle import build_oracle_maps  # noqa: E402
from evaluate import eval_point  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")

PACKINGS = [2000.0, 3400.0, 5700.0, 9500.0, 16000.0]
DENSITIES = [50.0, 90.0, 160.0, 280.0, 500.0]
SIGMAS = [0.5, 1.0, 2.0, 4.0, 8.0]
GATE1_SEED = config.MASTER_SEED + 100000

COLS = [
    "packing_idx", "sigma_idx", "density_idx",
    "packing_cells_per_mm2", "sigma_um", "mean_tx_per_cell",
    "field_seed", "transcript_seed", "model",
    "n_cells", "n_interior_cells", "grid", "dx_um", "r_mean_um", "L_um",
    "n_interior_transcripts",
    "oracle_acc", "naive_acc", "oracle_minus_naive", "oracle_ge_naive",
    "assign_ari_vs_truecells_oracle", "assign_ari_vs_truecells_naive",
    "same_type_error_frac_oracle", "same_type_error_frac_naive",
    "profile_ari_vs_truetype", "profile_ari_vs_trueclust",
    "perfect_ari_vs_truetype", "naive_profile_ari_vs_truetype", "n_cells_clustered",
    "realized_median_tx_per_cell", "realized_packing_cells_per_mm2",
]


def main():
    model = build_realistic_model_from_merfish()
    om = overlap_metrics(model)
    meta = {
        "model": model.name, "n_types": model.n_types, "n_genes": model.n_genes,
        "type_names": list(map(str, model.type_names)),
        "proportions": [float(x) for x in model.proportions],
        "excl_threshold": model.excl_threshold, "n_exclusive_markers": model.n_exclusive,
        "overlap_metrics": om,
        "reference": "Moffitt 2018 MERFISH mouse hypothalamic preoptic (squidpy MERFISH_0.24.h5ad); "
                     "per-type mean expression over the 155-gene panel from native cell_class labels",
    }
    with open(os.path.join(RESULTS, "gate1_expression_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"realistic model: {model.name} K={model.n_types} G={model.n_genes}")
    print(f"overlap: {om}")

    rows, t0 = [], time.time()
    for pi, packing in enumerate(PACKINGS):
        for si, sigma in enumerate(SIGMAS):
            fseed = GATE1_SEED + 1000 * (pi + 1) + 17 * (si + 1)
            field = build_field(packing, sigma, fseed, model=model)
            Dmax, argcell = build_oracle_maps(field)
            for di, density in enumerate(DENSITIES):
                tseed = fseed + 1009 * (di + 1)
                rec = eval_point(field, Dmax, argcell, density, tseed)
                rec.update({
                    "packing_idx": pi, "sigma_idx": si, "density_idx": di,
                    "packing_cells_per_mm2": packing, "sigma_um": sigma,
                    "mean_tx_per_cell": density, "field_seed": fseed,
                    "transcript_seed": tseed, "model": model.name,
                })
                rows.append(rec)
            print(f"pack={packing:>7.0f} sig={sigma:>4.1f} oracle@d160={rows[-3]['oracle_acc']:.3f} "
                  f"profARI@d160={rows[-3]['profile_ari_vs_truetype']:.3f} [{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)[COLS]
    out = os.path.join(RESULTS, "gate1_sweep.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out} ({len(df)} rows)")
    print(f"oracle_acc range: [{df.oracle_acc.min():.3f}, {df.oracle_acc.max():.3f}]")
    print(f"oracle>=naive everywhere: {bool(df.oracle_ge_naive.all())}")
    print(f"oracle profile-ARI range: [{df.profile_ari_vs_truetype.min():.3f}, "
          f"{df.profile_ari_vs_truetype.max():.3f}]  (Gate 0 was [1.000, 1.000])")
    return df


if __name__ == "__main__":
    main()
