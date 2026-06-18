"""Export the known-truth synthetic datasets and compute the oracle ceiling and the naive
baseline (matched-accuracy metric) for each headroom configuration.

Writes:
  data/headroom/<config>/transcripts.csv, cells.csv  (gitignored raw tables for methods)
  results/headroom_oracle_naive.csv                   (committed: oracle + naive accuracy)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate2_pin import build_xenium  # noqa: E402
from headroom_common import (CONFIGS, HEADROOM_SEED, DENSITY, N_TARGET, EMISSION,  # noqa: E402
                             build_config, export_config, matched_accuracy, RESULTS)


def main():
    model, real, abundance = build_xenium()
    rows = []
    for ci, (name, packing, sigma) in enumerate(CONFIGS):
        seed = HEADROOM_SEED + 1000 * (ci + 1)
        field, tx, oa, na = build_config(model, abundance, packing, sigma, seed)
        export_config(model, field, tx, oa, na, name)
        acc_o = matched_accuracy(oa, tx.true_cell, tx.interior)
        acc_n = matched_accuracy(na, tx.true_cell, tx.interior)
        # direct id-match for continuity with the gates (oracle uses true ids)
        m = tx.interior
        direct_o = float((oa[m] == tx.true_cell[m]).mean())
        direct_n = float((na[m] == tx.true_cell[m]).mean())
        rows.append({
            "config": name, "packing_cells_per_mm2": packing, "sigma_um": sigma,
            "density_tx_per_cell": DENSITY, "n_target_cells": N_TARGET, "emission": EMISSION,
            "n_cells": field.n_cells, "n_interior_cells": int(field.interior_cell.sum()),
            "n_interior_transcripts": int(m.sum()),
            "oracle_acc_matched": acc_o, "naive_acc_matched": acc_n,
            "oracle_acc_direct": direct_o, "naive_acc_direct": direct_n,
            "seed": seed,
        })
        print(f"{name:26s} pack={packing:6.0f} sig={sigma:.2f} | Tint={int(m.sum()):6d} "
              f"oracle={acc_o:.3f} naive={acc_n:.3f} (direct {direct_o:.3f}/{direct_n:.3f})")
    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "headroom_oracle_naive.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")
    return df


if __name__ == "__main__":
    main()
