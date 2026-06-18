"""Headroom figures from the committed result CSVs.

  figures/headroom_methods_vs_oracle.png : per-config matched transcript-assignment
      accuracy for the oracle ceiling, the naive baseline, and each real method.
  figures/headroom_dense_vs_sigma.png    : the same across the dense-regime sigma CI.
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


def _load():
    on = pd.read_csv(os.path.join(RESULTS, "headroom_oracle_naive.csv")).set_index("config")
    mp = os.path.join(RESULTS, "headroom_methods.csv")
    methods = pd.read_csv(mp) if os.path.exists(mp) else pd.DataFrame()
    return on, methods


def _mval(methods, cfg, method, mode):
    if len(methods) == 0:
        return np.nan
    r = methods[(methods.config == cfg) & (methods.method == method) & (methods["mode"] == mode)]
    return float(r.matched_accuracy.iloc[0]) if len(r) and r.status.iloc[0] == "ok" else np.nan


def fig_bars():
    on, methods = _load()
    cfgs = list(on.index)
    series = {
        "oracle (ceiling)": [on.loc[c, "oracle_acc_matched"] for c in cfgs],
        "naive nearest-nucleus": [on.loc[c, "naive_acc_matched"] for c in cfgs],
        "pciSeq (nuclei prior)": [_mval(methods, c, "pciSeq", "nuclei_prior") for c in cfgs],
        "ComSeg (nuclei prior)": [_mval(methods, c, "ComSeg", "nuclei_prior") for c in cfgs],
    }
    colors = {"oracle (ceiling)": "C2", "naive nearest-nucleus": "0.5",
              "pciSeq (nuclei prior)": "C0", "ComSeg (nuclei prior)": "C1"}
    fig, ax = plt.subplots(figsize=(11, 5.2))
    x = np.arange(len(cfgs)); w = 0.2
    for i, (lab, vals) in enumerate(series.items()):
        ax.bar(x + (i - 1.5) * w, np.nan_to_num(vals, nan=0), w, label=lab, color=colors[lab])
    ax.axhline(0.95, color="C3", ls="--", lw=1, label="0.95")
    ax.axhline(0.90, color="0.4", ls=":", lw=1, label="0.90")
    ax.set_xticks(x); ax.set_xticklabels(cfgs, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("matched transcript-assignment accuracy")
    ax.set_title("Segmenter headroom: real methods vs the oracle ceiling (known-truth synthetic)")
    ax.set_ylim(0, 1.0); ax.legend(fontsize=8, ncol=3, loc="upper center")
    fig.tight_layout()
    out = os.path.join(FIGURES, "headroom_methods_vs_oracle.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def fig_dense():
    on, methods = _load()
    dense = [(c, on.loc[c, "sigma_um"]) for c in on.index if c.startswith("dense")]
    dense.sort(key=lambda t: t[1])
    sig = [s for _, s in dense]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(sig, [on.loc[c, "oracle_acc_matched"] for c, _ in dense], "-o", color="C2", label="oracle (ceiling)")
    ax.plot(sig, [on.loc[c, "naive_acc_matched"] for c, _ in dense], "-o", color="0.5", label="naive nearest-nucleus")
    ax.plot(sig, [_mval(methods, c, "pciSeq", "nuclei_prior") for c, _ in dense], "-s", color="C0", label="pciSeq")
    ax.plot(sig, [_mval(methods, c, "ComSeg", "nuclei_prior") for c, _ in dense], "-^", color="C1", label="ComSeg")
    ax.axhline(0.95, color="C3", ls="--", lw=1); ax.axhline(0.90, color="0.4", ls=":", lw=1)
    ax.set_xlabel("displacement sigma (um), dense regime")
    ax.set_ylabel("matched transcript-assignment accuracy")
    ax.set_title("Dense-regime headroom across the sigma CI")
    ax.set_ylim(0.4, 1.0); ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGURES, "headroom_dense_vs_sigma.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_bars()
    fig_dense()


if __name__ == "__main__":
    main()
