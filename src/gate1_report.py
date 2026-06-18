"""Re-derive the Gate 1 report numbers from the committed CSVs (every number cites its
source file and column).  Regen command:  python src/gate1_report.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")


def _mono(df):
    vs = vp = 0
    for (p, d), s in df.groupby(["packing_cells_per_mm2", "mean_tx_per_cell"]):
        vs += int((np.diff(s.sort_values("sigma_um").oracle_acc.values) > 1e-3).sum())
    for (sg, d), s in df.groupby(["sigma_um", "mean_tx_per_cell"]):
        vp += int((np.diff(s.sort_values("packing_cells_per_mm2").oracle_acc.values) > 1e-3).sum())
    return vs, vp


def main():
    g0 = pd.read_csv(os.path.join(RESULTS, "sweep.csv"))
    g1 = pd.read_csv(os.path.join(RESULTS, "gate1_sweep.csv"))
    leak = pd.read_csv(os.path.join(RESULTS, "gate1_leakage.csv"))
    st = pd.read_csv(os.path.join(RESULTS, "gate1_structural.csv"))
    meta = json.load(open(os.path.join(RESULTS, "gate1_expression_meta.json")))

    print("=" * 96)
    print("TASK 1  realistic overlapping expression  (source: results/gate1_sweep.csv, gate1_expression_meta.json)")
    print("-" * 96)
    om = meta["overlap_metrics"]
    print(f"  model: {meta['model']}  K={meta['n_types']} G={meta['n_genes']}  ref: MERFISH cell_class profiles")
    print(f"  overlap: mean pairwise cosine={om['mean_pairwise_cosine']:.3f} (max={om['max_pairwise_cosine']:.3f}); "
          f"n_exclusive_markers={om['n_exclusive_markers']}; "
          f"mass on exclusive genes={om['frac_mass_exclusive_genes']:.3f}, shared={om['frac_mass_shared_genes']:.3f}")
    vs, vp = _mono(g1)
    print(f"  oracle_acc frontier range: [{g1.oracle_acc.min():.3f}, {g1.oracle_acc.max():.3f}]  "
          f"(Gate 0: [{g0.oracle_acc.min():.3f}, {g0.oracle_acc.max():.3f}])")
    print(f"  monotonicity violations: sigma={vs} packing={vp}")
    print(f"  oracle>=naive at all {len(g1)} points: {bool(g1.oracle_ge_naive.all())} (min gap {g1.oracle_minus_naive.min():+.4f})")
    print(f"  oracle profile-ARI range: [{g1.profile_ari_vs_truetype.min():.3f}, {g1.profile_ari_vs_truetype.max():.3f}] "
          f"median={g1.profile_ari_vs_truetype.median():.3f}  (Gate 0 was [1.000, 1.000]) -> DE-SATURATED")
    print(f"  naive profile-ARI range: [{g1.naive_profile_ari_vs_truetype.min():.3f}, {g1.naive_profile_ari_vs_truetype.max():.3f}] "
          f"median={g1.naive_profile_ari_vs_truetype.median():.3f}")
    key = ["packing_cells_per_mm2", "sigma_um", "mean_tx_per_cell"]
    m = g0.merge(g1, on=key, suffixes=("_g0", "_g1"))
    print(f"  oracle_acc vs Gate 0 at matched points: mean abs diff={ (m.oracle_acc_g1-m.oracle_acc_g0).abs().mean():.4f} "
          f"max abs diff={ (m.oracle_acc_g1-m.oracle_acc_g0).abs().max():.4f}")

    print("\n" + "=" * 96)
    print("TASK 2  data-pinned displacement sigma  (source: results/gate1_leakage.csv)")
    print("-" * 96)
    r0 = leak.iloc[0]
    print(f"  real spatial marker leakage (MERFISH section): mean={r0.real_spatial_leakage_mean:.4f} "
          f"p25={r0.real_spatial_leakage_p25:.4f} p75={r0.real_spatial_leakage_p75:.4f}")
    print(f"  literal leakage (diagnostic, biology-dominated): {r0.real_literal_leakage_mean:.3f}; "
          f"synthetic baseline (sigma=0)={r0.synthetic_baseline_leakage_sigma0:.4f}")
    print(f"  data-pinned sigma = {r0.sigma_pinned_um:.2f} um ({r0.sigma_fit_status}); "
          f"25-75 bracket = [{r0.sigma_bracket_lo_um:.2f}, {r0.sigma_bracket_hi_um:.2f}] um "
          f"(reproducible within [0,15] um)")
    for _, r in leak.iterrows():
        print(f"  oracle @ pinned sigma at {r.anchor_dataset} (pack={r.packing_cells_per_mm2:.0f}, "
              f"dens={r.density_tx_per_cell:.0f}): {r.oracle_acc_at_pinned:.3f} "
              f"(naive {r.naive_acc_at_pinned:.3f}); across bracket "
              f"[{r.oracle_acc_at_bracket_hi_sigma:.3f}, {r.oracle_acc_at_bracket_lo_sigma:.3f}]")

    print("\n" + "=" * 96)
    print("TASK 3  structural sensitivity  (source: results/gate1_structural.csv)")
    print("-" * 96)
    pk = st.packing_cells_per_mm2.max()
    for cond in ["baseline", "aniso", "mixture"]:
        sub = st[st.condition == cond]
        vs2 = sum(int((np.diff(s.sort_values('sigma_um').oracle_acc.values) > 1e-3).sum())
                  for _, s in sub.groupby('packing_cells_per_mm2'))
        vp2 = sum(int((np.diff(s.sort_values('packing_cells_per_mm2').oracle_acc.values) > 1e-3).sum())
                  for _, s in sub.groupby('sigma_um'))
        dense = sub[sub.packing_cells_per_mm2 == pk]
        d2 = float(dense[dense.sigma_um == 2.0].oracle_acc.iloc[0])
        d8 = float(dense[dense.sigma_um == 8.0].oracle_acc.iloc[0])
        print(f"  {cond:9s}: range [{sub.oracle_acc.min():.3f}, {sub.oracle_acc.max():.3f}] | "
              f"mono(sig,pack)=({vs2},{vp2}) | oracle>=naive all={bool(sub.oracle_ge_naive.all())} "
              f"(min gap {sub.oracle_minus_naive.min():+.3f}) | dense sig2/sig8={d2:.3f}/{d8:.3f}")


if __name__ == "__main__":
    main()
