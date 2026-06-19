"""Lever analysis Task 2: synthetic mechanism demonstration with known ground truth.

Question (sufficiency, NOT adjudication of the real finding): at the data-pinned
displacement, is transcript misassignment ALONE enough to manufacture a spurious
triple-positive (ESR1+/ERBB2+/PGR+) signal of the magnitude observed, in a region that
by construction contains NO truly triple-positive cell?

This runs entirely on synthetic data with known truth. It does NOT run on the real
transcripts and makes no claim that the real triple-positive is false; it only shows the
artifact is sufficient. The synthetic region is built to match the ROI's measured density
and three-marker cell-type composition (Task 1 / the released matrix):
  - DCIS type:  ER+/HER2+ cells (realistic ESR1 and ERBB2 means), PGR mean = 0.
  - PR-source:  genuinely PGR+ cells (realistic PGR mean), ESR1 = ERBB2 = 0 (matching the
                authors' scFFPE-seq: PGR+ cells do not co-express ESR1/ERBB2).
  - other:      background cells, low/zero on all three markers.
No cell type expresses all three markers, so the TRUE triple-positive count is zero.

We emit transcripts, apply isotropic Gaussian displacement at sigma, assign each
transcript to its nearest nucleus (the standard pipeline that produced the released
matrix), and count cells that then appear positive for all three markers (>= threshold).
sigma = 0 is the control (perfect assignment -> essentially zero spurious). The data-
pinned sigma and its CI ends are the test. Nothing is tuned to the observed number:
composition, marker means, density and packing are read from the data; sigma is the
data-pinned value.

Output: results/lever/lever_mechanism.csv. Seeds: LEVER_SEED + offsets, recorded.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config  # noqa: E402
from expression import ExpressionModel, _exclusivity, load_xenium_cells_genes  # noqa: E402
from generator import build_field, generate_transcripts  # noqa: E402
from oracle import naive_assign  # noqa: E402
from gate1_leakage_anchor import _counts_matrix  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = os.path.join(ROOT, "results", "lever")
DATA = os.path.join(ROOT, "data")
LEVER_SEED = config.MASTER_SEED + 920000

H5 = os.path.join(DATA, "xenium_breast_rep1_cell_feature_matrix.h5")
PARQ = os.path.join(DATA, "xenium_breast_rep1_cells.parquet")
MARKERS = ["ESR1", "PGR", "ERBB2"]   # gene order in the synthetic 4-gene model: ESR1, PGR, ERBB2, BG
POS_T = 2            # primary positivity threshold (counts)
THRESH_SWEEP = [1, 2, 3]
SIGMAS = {"sigma0": 0.0, "ci_lo": 1.43, "pinned": 2.10, "ci_hi": 2.50}
N_TARGET = 6000      # cells per synthetic field (large -> stable rare-event rate)
N_REP = 6            # replicate fields per sigma
KNN = 10


def local_packing(coords, k=KNN):
    d, _ = cKDTree(coords).query(coords, k=k + 1)
    return k / (np.pi * d[:, k] ** 2) * 1e6


def measure_real_params():
    """Read the released matrix to set the synthetic null's data-driven parameters:
    type fractions and per-type marker means in the triple-positive ROI."""
    X, genes, coords, total = load_xenium_cells_genes(H5, PARQ)
    genes = list(genes)
    iE, iP, iH = (genes.index(m) for m in MARKERS)
    esr1, pgr, erbb2 = X[:, iE], X[:, iP], X[:, iH]
    lpack = local_packing(coords)
    t = POS_T

    tp = (esr1 >= t) & (erbb2 >= t) & (pgr >= t)
    # ROI = tissue within 50um of an apparent triple-positive cell
    nbr = cKDTree(coords).query_ball_point(coords[tp], 50.0)
    roi = np.zeros(X.shape[0], bool)
    roi[np.unique(np.concatenate([np.array(ix, int) for ix in nbr]))] = True

    # classify ROI cells into DCIS / PR-source / other on the three markers
    dcis = roi & (esr1 >= t) & (erbb2 >= t) & (pgr < t)
    prsrc = roi & (pgr >= t) & ~((esr1 >= t) & (erbb2 >= t))
    other = roi & ~dcis & ~prsrc
    nroi = int(roi.sum())
    frac = np.array([dcis.sum(), prsrc.sum(), other.sum()], float) / nroi

    params = {
        "roi_n": nroi, "roi_pack_median": float(np.median(lpack[roi])),
        "roi_density_median": float(np.median(total[roi])),
        "frac_dcis": frac[0], "frac_pr": frac[1], "frac_other": frac[2],
        # per-type marker means (the displacement spillover magnitude comes from these)
        "dcis_mean_ESR1": float(esr1[dcis].mean()), "dcis_mean_ERBB2": float(erbb2[dcis].mean()),
        "dcis_total": float(np.median(total[dcis])),
        "pr_mean_PGR": float(pgr[prsrc].mean()) if prsrc.sum() else float(pgr[pgr >= t].mean()),
        "pr_total": float(np.median(total[prsrc])) if prsrc.sum() else float(np.median(total)),
        "other_total": float(np.median(total[other])) if other.sum() else float(np.median(total)),
        # observed apparent-triple rate among real ER+/HER2+ DCIS cells (section-wide), for comparison
        "obs_triple_rate_among_dcis": float(tp.sum()) / float(((esr1 >= t) & (erbb2 >= t)).sum()),
        "obs_n_triple": int(tp.sum()),
        "obs_n_dcis_section": int(((esr1 >= t) & (erbb2 >= t)).sum()),
    }
    # observed apparent-triple rate among real DCIS cells at each threshold (section-wide)
    for tt in (1, 2, 3):
        dcis_tt = (esr1 >= tt) & (erbb2 >= tt)
        trip_tt = dcis_tt & (pgr >= tt)
        params[f"obs_triple_rate_among_dcis_t{tt}"] = float(trip_tt.sum()) / float(dcis_tt.sum())
        params[f"obs_n_triple_t{tt}"] = int(trip_tt.sum())
        params[f"obs_n_dcis_t{tt}"] = int(dcis_tt.sum())
    return params


def build_null_model(p):
    """A 4-gene (ESR1, PGR, ERBB2, BG) 3-type model with NO triple-positive type.
    mu[type, gene] is set exactly to the data-derived per-cell mean counts."""
    gene_names = ["ESR1", "PGR", "ERBB2", "BG"]
    # rows: DCIS, PR-source, other ; cols: ESR1, PGR, ERBB2, BG
    me = np.zeros((3, 4))
    me[0, 0] = p["dcis_mean_ESR1"]; me[0, 2] = p["dcis_mean_ERBB2"]
    me[0, 3] = max(p["dcis_total"] - me[0, 0] - me[0, 2], 0.0)            # BG fills to real total
    me[1, 1] = p["pr_mean_PGR"]
    me[1, 3] = max(p["pr_total"] - me[1, 1], 0.0)
    me[2, 3] = max(p["other_total"], 0.0)
    proportions = np.array([p["frac_dcis"], p["frac_pr"], p["frac_other"]], float)
    proportions = proportions / proportions.sum()
    row = me.sum(axis=1, keepdims=True)
    comp = np.divide(me, row, out=np.zeros_like(me), where=row > 0)
    owner = _exclusivity(comp, 0.7)
    model = ExpressionModel(
        name="lever_null_no_triple", n_types=3, n_genes=4,
        type_names=["DCIS_ERpHER2p", "PR_source", "other"],
        proportions=proportions, composition=comp, gene_names=gene_names,
        excl_threshold=0.7, excl_owner=owner, mean_expr=me)
    abundance = row.ravel() / row.ravel().mean()
    mean_tx = float(row.ravel().mean())   # mu[t,g] = mean_tx*abundance[t]*comp[t,g] = me[t,g]
    return model, abundance, mean_tx


def apparent_triples(model, abundance, mean_tx, packing, sigma, seed, thresholds):
    """Generate one field, displace by sigma, naive-assign, count apparent-triple cells.
    Reports the apparent-triple rate among all DCIS-type cells (field average, primary),
    among all interior cells, and CONDITIONAL on a DCIS cell being spatially adjacent to a
    PR-source cell (the micro-environment the paper describes: DCIS abutting PGR sources).
    Adjacency = nearest PR-source within 1.5x the field median nearest-neighbour distance,
    the same near-band the displacement pin uses. True triple count is the ground truth."""
    f = build_field(packing, sigma, seed, model=model, n_target=N_TARGET)
    tx = generate_transcripts(f, mean_tx, seed + 1, abundance=abundance)
    na = naive_assign(f, tx.obs_xy)
    C = _counts_matrix(f.n_cells, na, tx.gene, model.n_genes)          # apparent (post-displacement)
    Ctrue = _counts_matrix(f.n_cells, tx.true_cell, tx.gene, model.n_genes)  # ground truth
    iE, iP, iH = 0, 1, 2
    interior = f.interior_cell
    is_dcis = (f.types == 0) & interior
    # adjacency of each DCIS cell to a PR-source cell (type 1)
    d_nn, _ = cKDTree(f.centers).query(f.centers, k=2)
    mnn = float(np.median(d_nn[:, 1]))
    pr_centers = f.centers[f.types == 1]
    if pr_centers.shape[0] > 0:
        d_pr, _ = cKDTree(pr_centers).query(f.centers, k=1)
        adj_pr = d_pr <= 1.5 * mnn
    else:
        adj_pr = np.zeros(f.n_cells, bool)
    dcis_adj = is_dcis & adj_pr
    out = {"n_dcis": int(is_dcis.sum()), "n_interior": int(interior.sum()),
           "n_dcis_adj_pr": int(dcis_adj.sum())}
    for thr in thresholds:
        appt = (C[:, iE] >= thr) & (C[:, iH] >= thr) & (C[:, iP] >= thr)
        truet = (Ctrue[:, iE] >= thr) & (Ctrue[:, iH] >= thr) & (Ctrue[:, iP] >= thr)
        out[f"app_triple_dcis_t{thr}"] = int((appt & is_dcis).sum())
        out[f"app_triple_all_t{thr}"] = int((appt & interior).sum())
        out[f"app_triple_dcis_adj_t{thr}"] = int((appt & dcis_adj).sum())
        out[f"true_triple_t{thr}"] = int((truet & interior).sum())
    return out


def main():
    os.makedirs(RESULTS, exist_ok=True)
    p = measure_real_params()
    print("real ROI params:")
    for k, v in p.items():
        print(f"  {k} = {v}")
    model, abundance, mean_tx = build_null_model(p)
    packing = p["roi_pack_median"]
    print(f"\nnull model proportions={np.round(model.proportions,4)} mean_tx={mean_tx:.1f} "
          f"packing={packing:.0f}; per-type marker means mu=\n{np.round(model.mean_expr,3)}")
    print(f"(true triple-positive cells by construction: 0)\n")

    rows = []
    for si, (sname, sig) in enumerate(SIGMAS.items()):
        agg = {f"app_triple_dcis_t{thr}": [] for thr in THRESH_SWEEP}
        agg.update({f"app_triple_all_t{thr}": [] for thr in THRESH_SWEEP})
        agg.update({f"app_triple_dcis_adj_t{thr}": [] for thr in THRESH_SWEEP})
        agg.update({f"true_triple_t{thr}": [] for thr in THRESH_SWEEP})
        ndcis, nint, nadj = [], [], []
        for r in range(N_REP):
            seed = LEVER_SEED + 1000 * si + 17 * r
            o = apparent_triples(model, abundance, mean_tx, packing, sig, seed, THRESH_SWEEP)
            ndcis.append(o["n_dcis"]); nint.append(o["n_interior"]); nadj.append(o["n_dcis_adj_pr"])
            for thr in THRESH_SWEEP:
                agg[f"app_triple_dcis_t{thr}"].append(o[f"app_triple_dcis_t{thr}"])
                agg[f"app_triple_all_t{thr}"].append(o[f"app_triple_all_t{thr}"])
                agg[f"app_triple_dcis_adj_t{thr}"].append(o[f"app_triple_dcis_adj_t{thr}"])
                agg[f"true_triple_t{thr}"].append(o[f"true_triple_t{thr}"])
        nd, ni, naj = float(np.sum(ndcis)), float(np.sum(nint)), float(np.sum(nadj))
        row = {"sigma_name": sname, "sigma_um": sig, "n_rep": N_REP,
               "n_dcis_total": nd, "n_interior_total": ni, "n_dcis_adj_pr_total": naj,
               "packing": packing, "obs_triple_rate_among_dcis": p["obs_triple_rate_among_dcis"]}
        for thr in THRESH_SWEEP:
            adc = float(np.sum(agg[f"app_triple_dcis_t{thr}"]))
            aal = float(np.sum(agg[f"app_triple_all_t{thr}"]))
            aaj = float(np.sum(agg[f"app_triple_dcis_adj_t{thr}"]))
            tru = float(np.sum(agg[f"true_triple_t{thr}"]))
            row[f"spurious_triple_rate_among_dcis_t{thr}"] = adc / nd if nd else np.nan
            row[f"spurious_triple_rate_among_all_t{thr}"] = aal / ni if ni else np.nan
            row[f"spurious_triple_rate_dcis_adj_pr_t{thr}"] = aaj / naj if naj else np.nan
            row[f"true_triple_count_t{thr}"] = tru
        rows.append(row)
        print(f"[{sname} sigma={sig}] DCIS apparent-triple rate (field avg): "
              f"t1={row['spurious_triple_rate_among_dcis_t1']:.4f} "
              f"t2={row['spurious_triple_rate_among_dcis_t2']:.4f} "
              f"t3={row['spurious_triple_rate_among_dcis_t3']:.4f} | "
              f"adjacent-to-PR-source: t1={row['spurious_triple_rate_dcis_adj_pr_t1']:.4f} "
              f"t2={row['spurious_triple_rate_dcis_adj_pr_t2']:.4f} "
              f"t3={row['spurious_triple_rate_dcis_adj_pr_t3']:.4f}  "
              f"(true t2={row['true_triple_count_t2']:.0f}) [obs t2={p['obs_triple_rate_among_dcis']:.4f}]")

    df = pd.DataFrame(rows)
    out = os.path.join(RESULTS, "lever_mechanism.csv")
    df.to_csv(out, index=False)
    # also persist the data-derived params for provenance
    pd.DataFrame([p]).to_csv(os.path.join(RESULTS, "lever_mechanism_params.csv"), index=False)
    print(f"\nwrote {out} and lever_mechanism_params.csv")
    return df


if __name__ == "__main__":
    main()
