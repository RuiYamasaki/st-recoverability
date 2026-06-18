"""Gate 3 Task 3: independent Xenium replication on a different tissue (human lung cancer).

Dataset: 10x Genomics "Human Lung Cancer Preview Data, Xenium Human Multi-Tissue and
Cancer Panel" (sample Xenium_V1_hLung_cancer_section, release 1.5.0; 150,365 cells).
Files (direct HTTP, no auth): cell_feature_matrix.h5 and cells.parquet under
cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hLung_cancer_section/.

Applies the identical pipeline as the breast dataset (gate3_validity.run_dataset):
cluster, validity-test markers, pin sigma with a bootstrap CI, evaluate dense-regime
oracle accuracy. Output: results/gate3_validity_lung.csv, results/gate3_pin_lung.csv.
"""
from __future__ import annotations

import os
import sys
import urllib.request

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config  # noqa: E402
from expression import build_realistic_model_from_xenium  # noqa: E402
from gate3_validity import run_dataset  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")
GATE3_SEED = config.MASTER_SEED + 520000

BASE = "https://cf.10xgenomics.com/samples/xenium/1.5.0/Xenium_V1_hLung_cancer_section/Xenium_V1_hLung_cancer_section"
LUNG_H5_URL = BASE + "_cell_feature_matrix.h5"
LUNG_CELLS_URL = BASE + "_cells.parquet"
LUNG_H5 = os.path.join(DATA, "xenium_lung_cancer_cell_feature_matrix.h5")
LUNG_CELLS = os.path.join(DATA, "xenium_lung_cancer_cells.parquet")


def _download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  cached: {os.path.basename(dest)} ({os.path.getsize(dest):,} bytes)")
        return
    os.makedirs(DATA, exist_ok=True)
    print(f"  downloading {os.path.basename(dest)} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "gate3/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as f:
        f.write(r.read())
    print(f"  done: {os.path.getsize(dest):,} bytes")


def main():
    _download(LUNG_H5_URL, LUNG_H5)
    _download(LUNG_CELLS_URL, LUNG_CELLS)
    model, real = build_realistic_model_from_xenium(
        n_types=15, seed=config.MASTER_SEED, h5_path=LUNG_H5, cells_parquet=LUNG_CELLS,
        name="realistic_xenium_lung")
    abundance = model.mean_expr.sum(axis=1); abundance = abundance / abundance.mean()
    print(f"[lung] K={model.n_types} median_nn={real['median_nn_um']:.2f}um "
          f"pack_median={real['packing_median_cells_per_mm2']:.0f} "
          f"pack_p90={real['packing_p90_cells_per_mm2']:.0f} dens={real['density_median_tx_per_cell']:.0f}")
    run_dataset(model, real, abundance, GATE3_SEED, "lung",
                os.path.join(RESULTS, "gate3_validity_lung.csv"),
                os.path.join(RESULTS, "gate3_pin_lung.csv"))


if __name__ == "__main__":
    main()
