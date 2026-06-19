# Worked lever analysis: the Janesick triple-positive as a Nature-Methods example

**Branch:** `lever-analysis` (off `lever-scout`). Append-only: no prior result or source file modified or deleted (verified by `git status --porcelain`; raw data gitignored). Regenerate: `python src/lever/lever_roi.py`, `python src/lever/lever_mechanism.py`, `python src/lever/lever_figures.py`.

**Integrity linchpin, honored.** The real Xenium triple-positive is never declared false using our own re-assignment of the real transcripts. The bound proves assignment is unrecoverable in this regime, so any re-assignment we perform carries no authority. The only evidence that the real co-expression is spurious is the orthogonal data (the authors' dissociated scFFPE-seq and the clinical PR-negative pathology). The synthetic mechanism test uses known-truth synthetic data only; it shows the artifact is sufficient, it does not adjudicate the real finding. No oracle or segmenter was run on the real transcripts to "correct" the triple-positive.

**What this reports.** Three legs and the artifact-consistent vs corroborated comparison. The executor reports the numbers; the orchestrator frames the claim and applies the go/no-go.

---

## Summary of the three legs (raw outcome)

1. **The ROI is above the ceiling at its own local density.** The apparent triple-positive cells sit in denser tissue than the section median (local packing 7,985 vs 6,534 cells/mm^2). The oracle best-possible assignment accuracy at the ROI is 0.77 to 0.80 (point), with the optimistic end of the sigma CI at 0.84 to 0.86: below 0.90 AND below 0.95 across the entire sigma CI, at all three ROI definitions. Leg holds.
2. **Misassignment at the data-pinned sigma is sufficient to reproduce a spurious triple-positive of the observed magnitude, in the micro-environment the paper describes.** With zero truly triple-positive cells by construction, the sigma=0 control gives 0% spurious (threshold >=2); at the data-pinned sigma, the spurious rate among DCIS cells adjacent to a PR source (the "small region in proximity to adipocytes" configuration) is 1.26%, rising to 1.96% at the upper sigma CI, comparable to the observed 1.85%. Leg holds, with the honest caveat that the effect is local (the field-wide average is smaller, 0.12%, because the artifact concentrates where PGR sources abut DCIS).
3. **The orthogonal contradiction holds and survives the dropout alternative.** The authors' own dissociated scFFPE-seq finds about 30 PGR-positive cells, none co-expressing ESR1 or ERBB2; the block is clinically PR-negative. Both under-report the co-expression that the spatial assay over-reports, the direction misassignment predicts. The dissociation-dropout alternative does not overturn it (below). Leg holds.

All three legs hold on the raw numbers. The go/no-go decision is the orchestrator's.

Figure: `figures/lever_worked_example.png` (Panel A = Leg 1 ceiling; Panel B = Leg 2 mechanism).

---

## Leg 1: the ceiling at the triple-positive ROI's own local density (Task 1)

**How the ROI was located.** Reading the released segmented cell-by-gene matrix and centroids (not a re-assignment of transcripts): the markers ESR1, PGR, ERBB2 are all on the 313-gene panel. ERBB2 is near-universal (mean 11.7 counts/cell, 87% of cells positive), ESR1 is moderate (15% positive), and PGR is rare (2.9% positive), consistent with the PR-negative clinical status. PGR is the limiting marker, so the relevant artifact is a rare PGR transcript spilling into the dense ER+/HER2+ field. We take the apparent triple-positive cells (ESR1>=2, ERBB2>=2, PGR>=2: 165 cells, 0.098% of the section) and measure local packing at their own locations and in their 50 um neighbourhood. Because the exact ROI polygon is not recoverable from the paper, we also report the conservative densest-DCIS-nest proxy.

**Result** (`results/lever/lever_roi.csv`; sigma reused from the lever diagnostic, `results/lever/lever_ceiling.csv` breast row: 2.10 um, CI [1.43, 2.50]; oracle on the converged grid):

| ROI definition | packing (cells/mm^2) | density (tx/cell) | oracle point | oracle optimistic / pessimistic CI end | naive | below 0.90 across CI? | below 0.95 across CI? |
|---|---|---|---|---|---|---|---|
| triple-positive cells' 50 um neighbourhood | 7,103 | 194 | 0.798 | 0.860 / 0.768 | 0.720 | yes | yes |
| triple-positive cells' own local packing | 7,985 | 295 | 0.787 | 0.850 / 0.748 | 0.705 | yes | yes |
| densest DCIS nest (p90 proxy) | 9,329 | 194 | 0.768 | 0.841 / 0.732 | 0.688 | yes | yes |

The triple-positive cells sit at higher local packing than the section median (7,985 vs 6,534), so the ROI is deeper above the ceiling than the section as a whole. Best-possible transcript-to-cell assignment in the ROI is below 0.95 (and below 0.90) across the entire sigma CI: the co-expression call cannot be established by any method in this regime, independent of segmenter.

## Leg 2: synthetic mechanism demonstration with known truth (Task 2)

**Construction** (`results/lever/lever_mechanism_params.csv`, all parameters read from the data; nothing tuned). A synthetic region at the ROI's measured packing (7,103 cells/mm^2) and three-marker cell-type composition: DCIS type (ER+/HER2+: mean ESR1 3.42, mean ERBB2 24.4, PGR 0), PR-source type (mean PGR 2.52, ESR1 0, ERBB2 0, matching the scFFPE-seq finding that PGR+ cells do not co-express the tumor markers), and other/background. No type expresses all three markers, so the true triple-positive count is zero by construction. We emit transcripts, apply Gaussian displacement at sigma, assign each transcript to its nearest nucleus (the standard pipeline), and count apparent triple-positives. Six replicate fields of 6,000 cells per sigma.

**Result** (`results/lever/lever_mechanism.csv`; spurious triple-positive rate among DCIS cells, percent; true triples = 0 at every sigma and threshold):

| sigma (um) | threshold | field average | adjacent to a PR source | observed real |
|---|---|---|---|---|
| 0.0 (control) | >=2 | 0.00 | 0.00 | 1.85 |
| 1.43 (CI low) | >=2 | 0.14 | 1.20 | 1.85 |
| 2.10 (pinned) | >=2 | 0.12 | 1.26 | 1.85 |
| 2.50 (CI high) | >=2 | 0.20 | 1.96 | 1.85 |
| 0.0 (control) | >=1 | 0.38 | 2.62 | 4.24 |
| 2.10 (pinned) | >=1 | 1.65 | 17.23 | 4.24 |

Reading: at the standard positivity threshold (>=2), the sigma=0 control produces zero spurious triple-positives, confirming the effect is entirely displacement-driven. At the data-pinned sigma, in the micro-environment the paper describes (a DCIS cell abutting a PR source, "in proximity to adipocytes"), the pure-displacement artifact produces 1.26% spurious triple-positives, rising to 1.96% at the upper sigma CI: comparable to the observed 1.85%. The adjacency band (within 1.5x the median nearest-neighbour distance) is the same near-band the displacement pin uses, not a tuned choice.

**Honest caveats (reported, not hidden).**
- The effect is local. Averaged over all DCIS cells in the field (including those far from any PR source) the spurious rate at the pinned sigma is only 0.12% at threshold >=2, well below observed. The artifact concentrates where PGR sources abut DCIS, which is exactly why the observed signal is a localized region rather than a section-wide population. The field-average is the conservative number; the adjacency-conditional is the configuration-matched number; both are reported.
- The magnitude is threshold-dependent. At the permissive >=1 threshold the artifact over-produces in the adjacency micro-environment (17.2% vs observed 4.24%) and reproduces the same order field-wide (1.65% vs 4.24%); at >=2 it matches observed in the micro-environment; at >=3 it falls to near zero (manufacturing three spilled PGR transcripts in one cell is very rare at this sigma). The artifact brackets the observed magnitude rather than hitting a single number, which is the expected behaviour and was not tuned.
- Nothing was tuned to the observed count: composition fractions, per-type marker means, density and packing are read from the data; sigma is the data-pinned value; the threshold is reported as a sweep.

Conclusion of Leg 2: transcript misassignment at the data-pinned displacement is a sufficient mechanism to manufacture a spurious triple-positive population of the observed magnitude in the described configuration, with a clean zero control. It does not, and is not claimed to, prove the real triple-positive is entirely artifactual; it shows the artifact is sufficient, so the spatial co-expression alone cannot establish a genuine triple-positive population in this regime.

## Leg 3: the orthogonal contradiction and the dropout alternative (Task 3)

**Quantified orthogonal evidence (primary source: Janesick et al., Nat Commun 14:8353, 2023; PMC10730913).**
- Dissociated single-cell FFPE RNA-seq (Chromium, same tissue, no spatial misassignment): verbatim, "Chromium scFFPE-seq yields only about 30 cells that are positive for PGR, but these cells do not express ERBB2 or ESR1." So the dissociated assay finds PGR+ cells but none co-expressing the tumor markers: it does not reproduce the triple-positive.
- Clinical pathology: the primary block is "TNM stage T2N1M0, ER+/HER2+/PR-." A PR-negative tumor should not contain a genuine PR+ triple-positive tumor population.

**Direction.** The spatial assay (Xenium) over-reports the co-expression (165 apparent triple-positive cells at threshold >=2, 1,053 at >=1); both orthogonal sources under-report it (scFFPE-seq: ~30 PGR+ cells, none triple; IHC: PR-negative). Over-reporting of co-expression in dense tissue is exactly the direction transcript misassignment predicts.

**The dropout alternative, argued both sides.**
- For dropout: the paper reports scFFPE-seq is less sensitive per cell than Xenium, "a median of 34 genes per cell for scFFPE-seq compared to a median of 62 genes per cell in the Xenium data" (downsampled to the 313-gene panel). Lower per-cell sensitivity could in principle cause a genuinely triple-positive cell to register as PGR-only if ESR1/ERBB2 dropped out.
- Against dropout (why the contradiction holds): (1) ERBB2 is the highest-expressed marker in this HER2+ tumor (mean ~12 counts/cell in Xenium, detected in 87% of cells). ERBB2 dropping out of a true tumor cell is very unlikely even at reduced sensitivity, so the scFFPE-seq PGR+ cells lacking ERBB2 indicates they are not tumor cells at all, not that ERBB2 was missed. (2) The PR-negative clinical pathology is an independent second leg on a different modality (protein/IHC) that dropout cannot explain: scRNA dropout does not make a PR-positive tumor read as PR-negative by antibody. (3) Dropout can only cause the orthogonal assay to fail to confirm the co-expression; it cannot manufacture the spatial co-expression, which still requires an explanation, and misassignment (proven unrecoverable here, shown sufficient in Leg 2) provides one. The dropout caveat weakens a strong "the dissociated data proves absence" reading, but the combination of the abundant-ERBB2 argument and the dropout-immune PR-negative pathology keeps the contradiction standing.

---

## Task 4: artifact-consistent vs corroborated comparison

The diagnostic plus orthogonal data separates artifact-consistent from corroborated on principled grounds, so the method is not a flag-everything machine. Note that CRC is more dense, and further above the ceiling, than the breast ROI, yet it is corroborated, not flagged: being above the ceiling is necessary but not sufficient for the artifact-consistent label; the orthogonal data decides.

| case | ceiling position (computed, oracle accuracy) | mechanism status | orthogonal data direction | classification |
|---|---|---|---|---|
| Janesick breast triple-positive (HER2+/ER+/PR+ DCIS) | above ceiling at the ROI: 0.77 to 0.80 (below 0.95 across the sigma CI) | misassignment at the data-pinned sigma reproduces the observed magnitude in the described micro-environment; sigma=0 control = 0 | DISAGREES: scFFPE-seq ~30 PGR+ cells, none triple; PR-negative IHC; spatial over-reports | artifact-consistent |
| CRC tumor (Oliveira, SPP1+ macrophage / tumor niche) | above ceiling, deeper: 0.77 (own median) / 0.68 (dense tail) | not run (corroborated; no need) | AGREES: matched dissociated Chromium scRNA-seq + Visium HD | corroborated |
| MERFISH hypothalamus (Moffitt, hybrid glut/GABA neurons) | sparser regime: 0.836 at the pinned sigma | not run (corroborated) | AGREES: matched dissociated scRNA-seq; co-transmission established orthogonally | corroborated |

The point of the table: the same class of claim (co-expression / double-positive) in above-ceiling or hard regimes is separated into artifact-consistent (Janesick) versus corroborated (CRC, MERFISH) entirely by the orthogonal data and the synthetic null, never by our re-assignment of anyone's real transcripts.

---

## Provenance

Every number traces to a committed result file and column; every Janesick quote traces to the primary source (Nat Commun 14:8353, PMC10730913); seeds recorded.

- **Leg 1**: `results/lever/lever_roi.csv` (columns `packing_cells_per_mm2`, `density_tx_per_cell`, `oracle_point`, `oracle_opt_ci_lo`, `oracle_pess_ci_hi`, `above_ceiling_0p90_whole_ci`, `above_ceiling_0p95_whole_ci`). Sigma source: `results/lever/lever_ceiling.csv` breast row. Code `src/lever/lever_roi.py`, seed LEVER_SEED = 20260618 + 910000 = 21170618 (per-call offsets in the CSV).
- **Leg 2**: `results/lever/lever_mechanism.csv` (columns `spurious_triple_rate_among_dcis_t{1,2,3}`, `spurious_triple_rate_dcis_adj_pr_t{1,2,3}`, `true_triple_count_t{1,2,3}`) and `results/lever/lever_mechanism_params.csv` (data-derived parameters and observed rates). Code `src/lever/lever_mechanism.py`, seed 20260618 + 920000 = 21180618; 6 replicates per sigma.
- **Leg 3**: Janesick et al., Nat Commun 14:8353 (2023), DOI 10.1038/s41467-023-43458-x; quotes from PMC10730913.
- **Figure**: `figures/lever_worked_example.png` from `src/lever/lever_figures.py`.
- **No prior file modified**: `git status --porcelain` shows only additions under `src/lever/`, `results/lever/`, `figures/`, `outputs/`, and the Task-0 spine appends. No prior-gate result or source file changed.

**Do not assert the real triple-positive is false.** The three legs are reported; the orchestrator frames the claim. The defensible statement is: the triple-positive sits in a regime where assignment is unrecoverable (so the spatial co-expression is not separable from a misassignment artifact by any method), a known-truth synthetic null shows misassignment at the data-pinned displacement is sufficient to manufacture the observed signature, and the authors' own orthogonal scFFPE-seq plus the PR-negative pathology independently indicate the co-expression is not genuine.
