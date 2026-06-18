"""Sweep the (density x packing x displacement) grid and write results/sweep.csv.

Grid ranges (locked before any oracle output was seen; chosen to bracket and
extend the real-data summary statistics measured in results/realism.csv):
  - packing: 2000-16000 cells/mm^2. Real anchors: MERFISH ~2600, Xenium ~13500.
  - density (mean tx/cell): 50-500. Real anchors: Xenium ~66-86, MERFISH ~267-312.
  - displacement sigma: 0.5-8 um. The realistic sub-band [0.5, 2] um is the
    least-constrained axis (displacement is not directly observable from the public
    per-cell tables); 4 and 8 um extend past plausibility to map the full frontier.

The oracle assignment map depends only on (packing, sigma), so the field and oracle
maps are built once per (packing, sigma) and reused across densities (this also lets
us confirm empirically that oracle assignment accuracy is density-independent).

Seeds: field/oracle seed = config.config_seed(packing_idx, sigma_idx); the transcript
draw adds a per-density offset. Every seed is recorded in the CSV.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from generator import build_field  # noqa: E402
from oracle import build_oracle_maps  # noqa: E402
from evaluate import eval_point  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)

PACKINGS = [2000.0, 3400.0, 5700.0, 9500.0, 16000.0]   # cells/mm^2
DENSITIES = [50.0, 90.0, 160.0, 280.0, 500.0]          # mean transcripts/cell
SIGMAS = [0.5, 1.0, 2.0, 4.0, 8.0]                     # displacement, micrometres


def main():
    rows = []
    t_start = time.time()
    for pi, packing in enumerate(PACKINGS):
        for si, sigma in enumerate(SIGMAS):
            fseed = config.config_seed(pi, si)
            field = build_field(packing, sigma, fseed)
            Dmax, argcell = build_oracle_maps(field)
            for di, density in enumerate(DENSITIES):
                tseed = fseed + 1009 * (di + 1)
                rec = eval_point(field, Dmax, argcell, density, tseed)
                rec.update({
                    "packing_idx": pi, "sigma_idx": si, "density_idx": di,
                    "packing_cells_per_mm2": packing, "sigma_um": sigma,
                    "mean_tx_per_cell": density,
                    "field_seed": fseed, "transcript_seed": tseed,
                })
                rows.append(rec)
            print(f"pack={packing:>7.0f} sig={sigma:>4.1f} done "
                  f"(grid={field.grid}, r={field.r_mean_um:.2f}um, "
                  f"oracle_acc@d160={rows[-3]['oracle_acc']:.3f})  "
                  f"[{time.time()-t_start:.1f}s]")

    cols = [
        "packing_idx", "sigma_idx", "density_idx",
        "packing_cells_per_mm2", "sigma_um", "mean_tx_per_cell",
        "field_seed", "transcript_seed",
        "n_cells", "n_interior_cells", "grid", "dx_um", "r_mean_um", "L_um",
        "n_interior_transcripts",
        "oracle_acc", "naive_acc", "oracle_minus_naive", "oracle_ge_naive",
        "assign_ari_vs_truecells_oracle", "assign_ari_vs_truecells_naive",
        "same_type_error_frac_oracle", "same_type_error_frac_naive",
        "profile_ari_vs_truetype", "profile_ari_vs_trueclust",
        "perfect_ari_vs_truetype", "naive_profile_ari_vs_truetype", "n_cells_clustered",
        "realized_median_tx_per_cell", "realized_packing_cells_per_mm2",
    ]
    df = pd.DataFrame(rows)[cols]
    out = os.path.join(RESULTS, "sweep.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}  ({len(df)} grid points)")
    print(f"oracle_acc range: [{df.oracle_acc.min():.3f}, {df.oracle_acc.max():.3f}]")
    print(f"oracle>=naive at every grid point: {bool(df.oracle_ge_naive.all())}")
    print(f"oracle ARI(vs truetype) range: "
          f"[{df.profile_ari_vs_truetype.min():.3f}, {df.profile_ari_vs_truetype.max():.3f}]")
    return df


if __name__ == "__main__":
    main()
