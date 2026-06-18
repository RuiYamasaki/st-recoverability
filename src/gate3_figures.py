"""Gate 3 figures from the committed result CSVs.

  figures/gate3_validity_separation.png : per-marker displacement signal vs permutation
      p-value on the breast dataset, admissible vs rejected (shows the validity test
      cleanly separates real-signal markers from biology-diluted ones).
  figures/gate3_dense_limit_forest.png  : dense-regime oracle accuracy (point + sigma-CI
      range) under the principled admissible-marker set, across cell-typing choices and
      the independent lung dataset, with the 0.95 and 0.9 lines.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")


def fig_validity():
    v = pd.read_csv(os.path.join(RESULTS, "gate3_validity_breast.csv"))
    v = v[np.isfinite(v.pval)]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    nl = -np.log10(v.pval.clip(lower=1e-4))
    for adm, c, lab in [(True, "C2", "admissible (p<0.05)"), (False, "C3", "rejected")]:
        s = v[v.admissible == adm]
        ax.scatter(s.real_signal, -np.log10(s.pval.clip(lower=1e-4)), c=c, s=45,
                   edgecolor="k", linewidth=0.3, label=lab, alpha=0.85)
    ax.axhline(-np.log10(0.05), color="0.4", ls=":", lw=1, label="p = 0.05")
    ax.axvline(0.0, color="0.7", lw=0.8)
    ax.set_xlabel("real adjacent-minus-distant displacement signal")
    ax.set_ylabel("-log10(permutation p-value)")
    ax.set_title("Gate 3 marker validity (breast): displacement signal vs position-permutation null")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGURES, "gate3_validity_separation.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def fig_forest():
    rows = []
    for tag in ["breast", "lung"]:
        f = os.path.join(RESULTS, f"gate3_pin_{tag}.csv")
        if os.path.exists(f):
            p = pd.read_csv(f).iloc[0]
            rows.append((f"{tag} (admissible set)", p.dense_oracle_point,
                         p.dense_oracle_at_sigma_ci_lo, p.dense_oracle_at_sigma_ci_hi, p.sigma_point_um))
    typ = os.path.join(RESULTS, "gate3_typing.csv")
    if os.path.exists(typ):
        for _, r in pd.read_csv(typ).iterrows():
            rows.append((f"breast {r.variant}", r.dense_oracle_point,
                         r.dense_oracle_at_sigma_ci_lo, r.dense_oracle_at_sigma_ci_hi, r.sigma_point_um))

    labels = [r[0] for r in rows]
    pts = np.array([r[1] for r in rows])
    los = np.array([min(r[2], r[3]) for r in rows])
    his = np.array([max(r[2], r[3]) for r in rows])
    y = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 0.55 * len(rows) + 1.8))
    ax.hlines(y, los, his, color="C0", lw=3, alpha=0.6)
    ax.plot(pts, y, "o", color="C0", ms=7)
    for yi, r in zip(y, rows):
        ax.annotate(f"sigma {r[4]:.2f}um", (max(r[2], r[3]) + 0.005, yi), va="center", fontsize=7.5, color="0.3")
    ax.axvline(0.95, color="C3", ls="--", lw=1.3, label="0.95 (Gate 2/3 kill line)")
    ax.axvline(0.90, color="0.5", ls=":", lw=1, label="0.90 (informative bar)")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("dense-regime oracle assignment accuracy (point and sigma-CI range)")
    ax.set_title("Gate 3: dense limit under the principled admissible-marker set")
    ax.set_xlim(0.55, 1.0)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    out = os.path.join(FIGURES, "gate3_dense_limit_forest.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_validity()
    fig_forest()


if __name__ == "__main__":
    main()
