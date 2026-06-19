# Extended Data items: value provenance and flags

**Branch:** `figures` (off `lever-analysis`), decision D14. Append-only: no prior result, source, or figure file modified or deleted. Regenerate all: `python src/figures/ed_render.py` (single item: `python src/figures/ed_render.py 7`). Outputs: `figures/ED01_*..ED10_*.{png,pdf}` (vector PDF + 400 dpi PNG, Okabe-Ito palette). No new analysis: every value is read from a committed results CSV (or a committed report, for two cell counts). The worked-example item uses the synthetic-field ceiling, the known-truth synthetic null, and the released matrix read only; no re-assignment of real transcripts.

## Per-item source and key values

| ED | item | source file(s) | key values (all match committed / manuscript) |
|---|---|---|---|
| 1 | Full answerability frontier | `results/sweep.csv` | oracle accuracy over packing x sigma per density slice; full-grid range 0.234 to 0.957 (`oracle_acc` min/max) |
| 2 | Generator calibration | `results/realism_oracle.csv` | transcripts/cell max rel. err 1.52% (`tx_match_rel_err`), packing max rel. err 1.4e-16 (`packing_match_rel_err`) |
| 3 | Overlapping expression | `results/gate1_sweep.csv` vs `results/sweep.csv` | frontier mean abs diff 0.006; oracle-minus-naive gap 0.071 (Gate 1) vs 0.066 (Gate 0), means of `oracle_minus_naive` |
| 4 | Model shape | `results/gate1_structural.csv`, `gate2_nbinom_frontier.csv`, `gate2_nbinom.csv` | baseline/aniso/mixture frontiers, oracle>=naive throughout; NB dispersion 4.19, dense shift -0.016 |
| 5 | Marker validity | `results/gate3_validity_breast.csv`, `gate3_validity_lung.csv` | per-marker signal + permutation p; breast 28/38 admissible, lung 31/32; failing markers shaded |
| 6 | Cell-typing stability | `results/gate3_typing.csv` | sigma 1.99 to 2.65 um; dense optimistic oracle at most 0.803 across K=10/14/14alt/20 |
| 7 | Selection-on-significance margin | `results/gate2_sigma_uncertainty.csv`, `gate3_pin_breast.csv`, `gate3_pin_lung.csv` | dense oracle vs sigma; 0.95 crossing ~0.36 um, ~5.5x below the 1.99 um pin; lung pin 1.83 um (31/32) control |
| 8 | Method-headroom metrics | `results/headroom_fair_metrics.csv`, `headroom_fair_free.csv` | one-to-one, many-to-one, ARI for Baysor/Proseg/ComSeg vs oracle/naive; free-segmentation rows |
| 9 | Per-dataset summary | `results/lever/lever_ceiling.csv`, `gate3_pin_lung.csv`, `realism.csv` + committed reports | packing, NN, tx/cell per dataset; cell counts (see flags) |
| 10 | Worked-example detail | `results/lever/lever_roi_pin199.csv`, `lever_mechanism.csv`, `leg2_local.csv`, `lever_mechanism_params.csv` | ROI ceiling at pin 1.99; full Leg 2 table; observed adjacency 24.2% vs synthetic 1.26-1.96% (19x) |

## Flags for the orchestrator

1. **ED9 cell counts and lung NN are partly outside committed CSVs (flagged in the figure note).** Breast (167,780) and lung (150,365) cell counts come from committed reports (GATE2_REPORT.md, GATE3_REPORT.md). The colorectal count (388,175) is from the 10x portal (dataset provenance), not a committed CSV. MERFISH has only committed per-FOV counts (realism.csv), not a section total. Lung median nearest-neighbour distance was computed during Gate 3 but not written to a committed CSV, so it is shown as n/a. All packing, NN (breast/CRC/MERFISH), and tx/cell values are committed-CSV values.

2. **ED7 0.95-crossing is interpolated from committed points.** The dense-oracle-vs-sigma curve uses the four committed points in `gate2_sigma_uncertainty.csv` (sigma 0.20, 0.79, 1.61, 2.53 um -> dense oracle 0.972, 0.891, 0.787, 0.672). The 0.95 crossing (~0.36 um) is read off this committed curve by linear interpolation, consistent with the Gate 2 report's stated "~0.35 um" crossing. The fall factor (~5.5x = 1.99/0.36) is "roughly fivefold." The dedicated selection-on-significance analysis was a residual item in Gate 3; ED7 renders it from the committed sigma-uncertainty curve plus the lung near-zero-selection (31/32) control, which is the committed basis available.

3. **ED6 per-variant admissible counts are not the primary pin.** ED6 shows each cell-typing variant's own candidate pool (e.g. K=14 here is 27/35), which differs from the primary breast pin (28/38, `gate3_pin_breast.csv`) because each variant re-clusters and re-derives candidate markers. This is the committed Gate 3 typing-robustness analysis, not a discrepancy with the locked 28/38.

**No printed value differs from the manuscript's locked numbers.** The three notes above are provenance/scope clarifications, not value mismatches. Spot-checked locked values reproduced exactly: 0.234-0.957 (ED1), 1.52%/1.4e-16 (ED2), 0.006 and 0.071/0.066 (ED3), 4.19/-0.016 (ED4), 28/38 and 31/32 (ED5), 1.99-2.65 and 0.803 (ED6), method gaps and free<naive (ED8), 24.2% vs 1.26-1.96% (ED10).

## Provenance
- New code: `src/figures/ed_render.py` (own table helper; imports `figstyle` read-only).
- New figures: `figures/ED01_*..ED10_*.{png,pdf}` (20 files). No new result CSVs (ED render reads existing committed results).
- Seeds: none needed (all items read committed values; no field is regenerated for ED).
- No prior result, source, or figure file modified (`git status --porcelain`: additions only, plus the Task-0 D14 STATUS.md append).
