"""Within-Linux recompute of the oracle ceiling and the naive nearest-nucleus baseline on
the SAME headroom configurations as the Windows run, so the method comparison is
apples-to-apples within this environment (cross-platform byte-identity is explicitly not
required; the committed Windows results stay canonical).

Append-only and namespaced: raw tables go to data/headroom_linux/<config>/ (gitignored),
the oracle/naive references to results/headroom_linux_oracle_naive.csv. This script reuses
the existing headroom infrastructure unchanged: the scorer (headroom_common.matched_accuracy),
the data build (headroom_common.build_config -> generator + oracle), and the table export
(headroom_common.export_config). No prior-gate or prior-headroom source file is modified.

  micromamba run -n st python src/headroom_linux_export.py
"""
from __future__ import annotations

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import headroom_common as hc  # noqa: E402
from gate2_pin import build_xenium  # noqa: E402

# Namespace the raw export under data/headroom_linux/ (gitignored); leave the Windows
# export at data/headroom/ untouched. export_config reads hc.EXPORT_DIR at call time.
EXPORT_DIR_LINUX = os.path.join(hc.DATA, "headroom_linux")
hc.EXPORT_DIR = EXPORT_DIR_LINUX


def main():
    model, real, abundance = build_xenium()
    rows = []
    for ci, (name, packing, sigma) in enumerate(hc.CONFIGS):
        seed = hc.HEADROOM_SEED + 1000 * (ci + 1)
        field, tx, oa, na = hc.build_config(model, abundance, packing, sigma, seed)
        hc.export_config(model, field, tx, oa, na, name)
        acc_o = hc.matched_accuracy(oa, tx.true_cell, tx.interior)
        acc_n = hc.matched_accuracy(na, tx.true_cell, tx.interior)
        # direct id-match (oracle/naive use the true ids, so this equals matched accuracy)
        m = tx.interior
        direct_o = float((oa[m] == tx.true_cell[m]).mean())
        direct_n = float((na[m] == tx.true_cell[m]).mean())
        rows.append({
            "config": name, "packing_cells_per_mm2": packing, "sigma_um": sigma,
            "density_tx_per_cell": hc.DENSITY, "n_target_cells": hc.N_TARGET,
            "emission": hc.EMISSION,
            "n_cells": field.n_cells, "n_interior_cells": int(field.interior_cell.sum()),
            "n_interior_transcripts": int(m.sum()),
            "oracle_acc_matched": acc_o, "naive_acc_matched": acc_n,
            "oracle_acc_direct": direct_o, "naive_acc_direct": direct_n,
            "seed": seed,
        })
        print(f"{name:26s} pack={packing:6.0f} sig={sigma:.2f} | Tint={int(m.sum()):6d} "
              f"oracle={acc_o:.3f} naive={acc_n:.3f} (direct {direct_o:.3f}/{direct_n:.3f})")
    df = pd.DataFrame(rows)
    out = os.path.join(hc.RESULTS, "headroom_linux_oracle_naive.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}\nraw tables under {EXPORT_DIR_LINUX}/")
    return df


if __name__ == "__main__":
    main()
