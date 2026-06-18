"""Gate 3 Task 1: principled non-circular marker-validity test on the breast dataset.

Builds the broad candidate marker pool (exclusivity >= 0.5), runs the permutation
validity test (gate3_common.marker_validity), reports the admissible set, and shows
explicitly whether the markers that destabilised sigma at cutoff 0.6 (exclusivity in
[0.6, 0.7)) pass or fail, with their permutation p-values. Re-pins sigma on the admissible
set with a bootstrap CI and reports dense-regime oracle accuracy across that CI.

Outputs: results/gate3_validity_breast.csv (per-marker), results/gate3_pin_breast.csv.
Seeds: GATE3_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from gate2_pin import build_xenium  # noqa: E402
from gate3_common import (marker_validity, admissible_owner, pin_admissible,  # noqa: E402
                          candidate_markers, LOOSE_LO, LOOSE_HI, N_PERM, P_ADMIT)
from expression import _exclusivity  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
GATE3_SEED = config.MASTER_SEED + 500000


def run_dataset(model, real, abundance, seed, tag, out_validity, out_pin):
    owner_cand, candidates = candidate_markers(model)
    print(f"[{tag}] candidate markers (exclusivity>=0.5): {len(candidates)}")
    vdf = marker_validity(real, owner_cand, candidates, seed)
    # label each candidate's exclusivity band
    col = model.composition.sum(axis=0)
    excl = np.array([model.composition[int(owner_cand[g]), g] / col[g] if col[g] > 0 else 0
                     for g in vdf.gene])
    vdf["exclusivity"] = excl
    vdf["gene_name"] = [model.gene_names[int(g)] for g in vdf.gene]
    vdf["is_loose_0p6"] = (vdf.exclusivity >= LOOSE_LO) & (vdf.exclusivity < LOOSE_HI)
    vdf["is_exclusive_0p7"] = vdf.exclusivity >= LOOSE_HI
    vdf.to_csv(out_validity, index=False)

    adm = vdf[vdf.admissible]
    print(f"[{tag}] admissible markers (p<{P_ADMIT}, {N_PERM} perms): {len(adm)} of {len(vdf)}")
    loose = vdf[vdf.is_loose_0p6]
    excl07 = vdf[vdf.is_exclusive_0p7]
    print(f"[{tag}] exclusive (>=0.7) markers: {int(excl07.admissible.sum())}/{len(excl07)} admissible; "
          f"loose [0.6,0.7) markers: {int(loose.admissible.sum())}/{len(loose)} admissible")
    if len(loose):
        print(f"[{tag}] loose-marker p-values: "
              f"{[f'{n}:{p:.3f}' for n, p in zip(loose.gene_name, loose.pval)]}")

    owner_adm = admissible_owner(model, vdf)
    pin = pin_admissible(model, real, abundance, owner_adm, seed + 100)
    pin.update({"dataset": tag, "n_candidate": len(candidates), "n_admissible": len(adm),
                "n_loose_total": int(len(loose)), "n_loose_admissible": int(loose.admissible.sum()),
                "n_exclusive_total": int(len(excl07)), "n_exclusive_admissible": int(excl07.admissible.sum()),
                "seed": seed})
    pd.DataFrame([pin]).to_csv(out_pin, index=False)
    print(f"[{tag}] re-pinned sigma (admissible set) = {pin['sigma_point_um']:.2f} um "
          f"CI [{pin['sigma_ci_lo_um']:.2f}, {pin['sigma_ci_hi_um']:.2f}]")
    print(f"[{tag}] dense oracle: point={pin['dense_oracle_point']:.3f} "
          f"CI ends [{pin['dense_oracle_at_sigma_ci_hi']:.3f}, {pin['dense_oracle_at_sigma_ci_lo']:.3f}]; "
          f"sparse={pin['sparse_oracle_point']:.3f}")
    return vdf, pin


def main():
    model, real, abundance = build_xenium()
    run_dataset(model, real, abundance, GATE3_SEED, "breast",
                os.path.join(RESULTS, "gate3_validity_breast.csv"),
                os.path.join(RESULTS, "gate3_pin_breast.csv"))


if __name__ == "__main__":
    main()
