"""Within-Linux headroom figures from the committed result CSVs.

  figures/headroom_linux_methods_vs_oracle.png : per-config matched transcript-assignment
      accuracy for the oracle ceiling, naive baseline, and each field-standard method/mode.
  figures/headroom_linux_dense_vs_sigma.png    : the same across the dense-regime sigma CI.

Regen: micromamba run -n st python src/headroom_linux_figures.py
Source: results/headroom_linux_oracle_naive.csv, results/headroom_linux_methods.csv.
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
    on = pd.read_csv(os.path.join(RESULTS, "headroom_linux_oracle_naive.csv")).set_index("config")
    mp = os.path.join(RESULTS, "headroom_linux_methods.csv")
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
    series = [
        ("oracle (ceiling)", [on.loc[c, "oracle_acc_matched"] for c in cfgs], "C2"),
        ("naive nearest-nucleus", [on.loc[c, "naive_acc_matched"] for c in cfgs], "0.5"),
        ("Baysor (nuclei prior)", [_mval(methods, c, "Baysor", "nuclei_prior") for c in cfgs], "C0"),
        ("Baysor (free)", [_mval(methods, c, "Baysor", "free_segmentation") for c in cfgs], "C9"),
        ("Proseg (nuclei prior)", [_mval(methods, c, "Proseg", "nuclei_prior") for c in cfgs], "C1"),
        ("ComSeg (nuclei prior)", [_mval(methods, c, "ComSeg", "nuclei_prior") for c in cfgs], "C4"),
    ]
    fig, ax = plt.subplots(figsize=(13, 5.4))
    x = np.arange(len(cfgs)); w = 0.14
    for i, (lab, vals, col) in enumerate(series):
        ax.bar(x + (i - 2.5) * w, np.nan_to_num(vals, nan=0), w, label=lab, color=col)
    ax.set_xticks(x); ax.set_xticklabels(cfgs, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("matched transcript-assignment accuracy")
    ax.set_title("Headroom redo on Linux (aarch64): field-standard methods vs the oracle ceiling")
    ax.set_ylim(0, 1.0); ax.legend(fontsize=8, ncol=6, loc="upper center")
    fig.tight_layout()
    out = os.path.join(FIGURES, "headroom_linux_methods_vs_oracle.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def fig_dense():
    on, methods = _load()
    dense = [(c, on.loc[c, "sigma_um"]) for c in on.index if c.startswith("dense")]
    dense.sort(key=lambda t: t[1])
    sig = [s for _, s in dense]
    fig, ax = plt.subplots(figsize=(7.8, 5))
    ax.plot(sig, [on.loc[c, "oracle_acc_matched"] for c, _ in dense], "-o", color="C2", label="oracle (ceiling)")
    ax.plot(sig, [on.loc[c, "naive_acc_matched"] for c, _ in dense], "-o", color="0.5", label="naive nearest-nucleus")
    ax.plot(sig, [_mval(methods, c, "Baysor", "nuclei_prior") for c, _ in dense], "-s", color="C0", label="Baysor (nuclei prior)")
    ax.plot(sig, [_mval(methods, c, "Baysor", "free_segmentation") for c, _ in dense], "--s", color="C9", label="Baysor (free)")
    ax.plot(sig, [_mval(methods, c, "Proseg", "nuclei_prior") for c, _ in dense], "-^", color="C1", label="Proseg (nuclei prior)")
    ax.plot(sig, [_mval(methods, c, "ComSeg", "nuclei_prior") for c, _ in dense], "-D", color="C4", label="ComSeg (nuclei prior)")
    ax.set_xlabel("displacement sigma (um), dense regime")
    ax.set_ylabel("matched transcript-assignment accuracy")
    ax.set_title("Dense-regime headroom across the sigma CI (within Linux)")
    ax.set_ylim(0.4, 0.95); ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGURES, "headroom_linux_dense_vs_sigma.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"wrote {out}")


def main():
    fig_bars()
    fig_dense()


if __name__ == "__main__":
    main()
