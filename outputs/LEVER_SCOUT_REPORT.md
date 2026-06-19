# Lever-scout report: highest-tier biological lever above the answerability ceiling

**Branch:** `lever-scout` (off HEAD `f963190`, headroom-fair). Append-only: no prior-gate result or source file was modified or deleted (verified by `git status --porcelain`: only `results/lever/` and `src/lever/` are new; raw data under `data/lever/` is gitignored).

**What this is.** A scouting pass, not a verdict. It (1) finds notable published dense-tissue spatial-transcriptomics conclusions that depend on transcript-to-cell assignment with open data, (2) checks each for independent orthogonal ground truth, and (3) runs the existing displacement-pinning + oracle diagnostic on the real data to read off where each sits relative to the answerability ceiling. The executor classifies; it does not assert which lever to pursue, and it asserts no published finding is false.

**The hard logical constraint, honored throughout.** The bound proves transcript assignment is unrecoverable above the ceiling, so our own re-assignment of real data carries no authority. No "this conclusion is wrong" claim is ever made from our re-assignment. Every disagreement reported below rests entirely on independent orthogonal data (matched dissociated sequencing, clinical pathology, or an orthogonal assay), never on us re-assigning anyone's transcripts. A Tier A entry means the orthogonal data disagrees; the diagnostic only certifies that the regime is unrecoverable.

---

## 1. Diagnostic results (computed this session on real data)

Pipeline per dataset: load the standard 10x Xenium bundle (`cell_feature_matrix.h5` + `cells.parquet`); measure local packing (median and 90th-percentile dense tail) from cell centroids and median transcripts per cell; pin the displacement sigma with the Gate 3 principled non-circular marker-validity test plus admissible-set bootstrap; evaluate the oracle (Bayes-optimal, method-independent) assignment accuracy at the dataset's own packing and density on the converged grid. Source: `src/lever/lever_diagnostic.py`; output `results/lever/lever_ceiling.csv`. Master seed 20260618; LEVER_SEED 21160618 (breast row seed 21160618, CRC row seed 21161618).

| dataset (real data run) | median packing (cells/mm^2) | dense-tail p90 packing | median tx/cell | admissible markers | pinned sigma (um), 95% CI | oracle @ own median packing (point; optimistic CI end) | oracle @ own dense tail (point; optimistic CI end) | above 0.9 ceiling? |
|---|---|---|---|---|---|---|---|---|
| Xenium breast Rep1 (Janesick 2023) | 6,534 | 11,956 | 164 | 28 of 38 | 2.10, [1.43, 2.50] | 0.802; 0.863 | 0.736; 0.819 | YES (both regimes) |
| Xenium colorectal cancer (10x public, Onboard 2.0.0) | 17,424 | 33,607 | 78 | 38 of 40 | 1.46, [1.34, 2.05] | 0.769; 0.787 | 0.679; 0.706 | YES (both regimes) |
| MERFISH hypothalamus (Moffitt 2018) [from committed Gate 1] | 2,575 | n/a (sparse) | ~286 | n/a | 2.67, bracket [0.10, 3.03] | 0.836 at pin (sparse) | n/a | central pin below 0.9; optimistic end recoverable |

"Optimistic CI end" is the oracle accuracy at the low-sigma end of the bootstrap CI (the most generous case); if even that is below 0.9 the regime is unrecoverable across the whole CI. Both Xenium tumor datasets are above the ceiling at every regime and across the entire sigma CI. The independently pinned sigma on three tissues (breast 2.10, CRC 1.46, MERFISH 2.67) stays in the same 1.4 to 2.7 um band the gates established, so the ceiling position is not an artifact of a borrowed displacement value. The CRC tumor is the densest tissue measured in this project (median 17,424 cells/mm^2, nearest-neighbour 5.6 um), and it sits deepest below the ceiling.

Cross-check: the breast run reproduces Gate 3 (sigma 1.99 to 2.10 um, 28 of 38 admissible markers identical, dense-regime oracle 0.66 to 0.82), confirming the lever pipeline is the same machinery.

---

## 2. Candidate table

Tier key: **A** = above the ceiling (diagnostic-confirmed) AND independent orthogonal ground truth disagrees in the ceiling-predicted direction (spatial over-reports co-expression/interaction). **B** = above the ceiling AND provably assignment-dependent AND no orthogonal contradiction available (supports "not separable from a misassignment artifact," never "false"). **None** = at or below the ceiling, or orthogonal evidence agrees (claim corroborated), or the diagnostic could not be run on the exact data.

| # | claim (assignment-dependent) | paper + DOI | dataset + access | tissue density (diagnostic) | ceiling position (computed) | assignment-dependent? | orthogonal ground truth | tier |
|---|---|---|---|---|---|---|---|---|
| 1 | A triple-positive ERBB2+/ESR1+/PGR+ (HER2+/ER+/PR+) co-expressing DCIS region, seen only by Xenium in situ | Janesick et al., *Nat Commun* 14:8353 (2023), 10.1038/s41467-023-43458-x | Xenium FFPE Human Breast Cancer Rep1, 10x release 1.0.1; direct HTTP h5+parquet (also GEO GSE243280) | dense DCIS nests; section median 6,534, dense tail 11,956 cells/mm^2 | oracle 0.736 (CI 0.70-0.82) dense tail; 0.802 (opt 0.863) at median. Below 0.9 across the whole sigma CI -> ABOVE | YES: a triple-positive call requires ERBB2, ESR1 and PGR transcripts to land in the same segmented cell; boundary spillover in dense DCIS manufactures the co-expression | **DISAGREES.** (a) Matched dissociated scFFPE-seq, same tissue: "While few PGR+ cells were found in the scFFPE-seq data, these cells did not coexpress ESR1 or ERBB2." (b) Sample clinical status ER+/HER2+/**PR-**. Both contradict a real PR+ triple-positive, in the over-reporting direction the ceiling predicts | **A** |
| 2 | SPP1+ macrophages localised near TGFBI+PERP+ tumor cells; SELENOP+STAB1+ macrophages near REG1A+LCN2+ tumor/goblet cells (spatial niche co-localisation) | de Oliveira et al., *Nat Genet* 57(6):1512-1523 (2025), 10.1038/s41588-025-02193-3 | GEO GSE280318 (P1/P2/P5 CRC). Diagnostic run on the 10x public CRC Xenium (Xenium V1 Human Colorectal Cancer Add-on FFPE, Onboard 2.0.0), same tissue and panel family, NOT byte-identical to the P-samples | extreme dense: median 17,424, dense tail 33,607 cells/mm^2 | oracle 0.769 (median); 0.679 (dense tail). Below 0.9 across the sigma CI -> ABOVE | YES: macrophage-versus-tumor identity and their proximity both require correct per-cell assignment of densely packed cells | **AGREES / no contradiction.** Matched dissociated Chromium Flex scRNA-seq and Visium HD on the same patients corroborate the cell identities and heterogeneity (the proximity component is supported by Visium HD, a lower-resolution spatial platform, not by dissociated data) | **None** (corroborated; dense but not separable as artifact because orthogonal data supports it) |
| 3 | Hybrid glutamatergic/GABAergic neurons co-expressing Slc17a6/8 and Slc32a1 in individual cells (GABA/glutamate co-release) | Moffitt et al., *Science* 362:eaau5324 (2018), 10.1126/science.aau5324 | Dryad 10.5061/dryad.8t8s248 (also cached `merfish_moffitt.h5ad`) | sparse: 2,575 cells/mm^2 | oracle 0.836 at pinned sigma 2.67 um (sparse). Moderate-hard regime, not the dense tail | YES: single-cell co-expression; misassignment between adjacent excitatory/inhibitory neurons could fabricate double-positives | **AGREES.** The hybrid co-expressing clusters were independently identified in matched dissociated scRNA-seq (same study); glut/GABA co-transmission is established by electrophysiology/optogenetics in other regions | **None** (corroborated by orthogonal dissociated data) |
| 4 | 16 ligand-receptor pairs significantly enriched at the tumor-T-cell interface (NSCLC) | He et al., *Nat Biotechnol* 40:1794-1806 (2022), 10.1038/s41587-022-01483-z | CosMx SMI NSCLC FFPE (Bruker/NanoString portal; exprMat + metadata + tx CSVs); ~6,800 cells/mm^2 reported | moderate (~6,800 cells/mm^2 reported) | **NOT computed.** The CosMx flat-file format was not ingested this session; ceiling position is inferred (a ~6,500 cells/mm^2 Xenium-like regime gives oracle ~0.8) but not measured on its data | YES: interface L-R enrichment needs ligand-to-tumor and receptor-to-T-cell correct assignment; spillover inflates apparent interface signalling | **NONE available.** Only count-level scRNA-seq concordance is reported; the spatial interface claim is not orthogonally validated, and no contradicting source was found | **B-shaped but UNCONFIRMED** (no orthogonal contradiction; ceiling not computed on its exact data) |
| 5 | Ligand-receptor pairs across neighbouring cells (e.g. Tgfb1 in microglia, Eng/Acvrl1 in adjacent endothelium) | Eng et al., *Nature* 568:235-239 (2019), 10.1038/s41586-019-1049-y | github.com/CaiGroup/seqFISH-PLUS; Zenodo 10.5281/zenodo.2669683 | dense cortex, small fields | not computed (custom format) | YES: the neighbour-pair call needs ligand and receptor assigned to the correct adjacent cells | **AGREES.** The co-localised transcripts were validated by smFISH, an orthogonal higher-specificity imaging assay (same study) | **None** (corroborated by orthogonal assay) |

---

## 3. The Tier A candidate, laid out for independent verification

**Candidate 1: the Xenium triple-positive HER2+/ER+/PR+ DCIS co-expression claim (Janesick et al. 2023).** This is the only candidate that is both diagnostic-confirmed above the ceiling and contradicted by clean independent orthogonal data. The orchestrator can verify every piece from primary sources:

1. **The claim is real and assignment-dependent.** Janesick et al., *Nature Communications* 14:8353 (2023), DOI 10.1038/s41467-023-43458-x (verified via Crossref). Verbatim from the full text (PMC10730913, Figure 5 region): "we identified a small triple positive (*ERBB2*+/*ESR1*+/*PGR*+ (HER2+/ER+/PR+)) DCIS region." The authors state this triple-positive population was detected by Xenium in situ and not by the two whole-transcriptome methods. A triple-positive cell exists only if all three transcripts are assigned to the same segmented cell, which is exactly the quantity transcript misassignment corrupts.

2. **The diagnostic places this dataset above the ceiling.** The dataset is the Xenium FFPE Human Breast Cancer Rep1 release (10x 1.0.1), the same data analysed in the paper. Run this session: pinned displacement sigma 2.10 um (95% CI [1.43, 2.50]) on 28 admissible markers; oracle best-possible assignment accuracy 0.736 at the dense tumor-nest tail (11,956 cells/mm^2, CI 0.70 to 0.82) and 0.802 at the section median (optimistic end 0.863). Every value is below 0.9 across the entire sigma CI, so transcript-to-cell assignment is provably unrecoverable in the regime where the DCIS nests sit. `results/lever/lever_ceiling.csv`, row `breast`.

3. **Independent orthogonal ground truth disagrees, in the ceiling-predicted direction.** Two independent sources, neither of which involves us re-assigning any transcript:
   - **Matched dissociated single-cell FFPE RNA-seq (same tissue block).** Verbatim: "While few *PGR*+ cells were found in the scFFPE-seq data, these cells did not coexpress *ESR1* or *ERBB2*." Dissociated sequencing has no spatial misassignment, so it is the orthogonal truth for whether a genuine triple-positive cell population exists. It does not reproduce the co-expression.
   - **Clinical pathology of the sample.** Verbatim: the primary block is "TNM stage T2N1M0, ER+/HER2+/**PR-**." A PR-negative tumor should not contain a real PR+ triple-positive population.
   - **Direction.** Misassignment in dense tissue is predicted to over-report co-expression (boundary or background PGR transcripts bleed into ERBB2+/ESR1+ cells, fabricating triple-positives). The spatial data reports the co-expression; both orthogonal sources do not. The disagreement points exactly the way the ceiling predicts.

**Honest weaknesses the orchestrator must weigh (these are real and not hidden):**
- The discordance is **visible in the original paper**, not a hidden refutation. The authors themselves reported that the triple-positive appears only in Xenium and not in scFFPE-seq. This makes the Tier A claim safe and easy to verify, but it is not an "overturning" of a result the field relied on.
- The triple-positive is a **minor sub-finding** in a platform and methods paper, not a headline biological discovery. Its prominence is modest.
- A **separate** double-positive claim in the same paper (a tumor + myoepithelial boundary subcluster co-expressing ERBB2/ABCC11 and MYLK/DST, Figure 6) is **NOT** a clean Tier A target: the authors pre-empt the artifact interpretation with a normal-duct control ("cell type-specific markers are not commingled ... indicating that our observations are not an artifact of gross segmentation errors"). Only the Figure 5 triple-positive is contradicted by orthogonal data; the Figure 6 subcluster is defended as real. Do not conflate the two.

**Independent corroboration that this failure mode is accepted (peer-reviewed, strengthens the Tier B floor under the Tier A claim):**
- SPLIT: Bilous et al., *Nat Methods* 23:1152-1162 (2026), 10.1038/s41592-026-03089-8. Quantifies Xenium transcript spillover against matched snRNA-seq + IHC in breast and lung tumor; double-positive/co-expression calls in dense tumor are frequently contamination.
- cellAdmix: Mitchel et al., *Nat Genet* 58:434-444 (2026), 10.1038/s41588-025-02497-4. Segmentation/misassignment errors "frequently dominate" differential-expression, neighbor-influence, and ligand-receptor results. (Its correction re-models the same data, so it is NOT itself orthogonal truth; the dissociated references it uses are.)
- ovrlpy: Tiesmeyer et al., *Nat Biotechnol* (2026), 10.1038/s41587-026-03004-8. 3D signal overlap produces mixed/co-expression profiles.
- MisTIC (preprint, not peer-reviewed): Yang et al., bioRxiv 2025, DOI 10.64898/2025.12.11.693759 (PMC12724677). In LUAD Xenium/MERSCOPE, 35.4% of "B cells" express the T-cell gene CD3E and 29.5% of "macrophages" express the epithelial gene EPCAM, contradicted by dissociated LUAD scRNA-seq where those genes are mutually exclusive. Direct documentation of spurious double-positives; cited as mechanism support, not as a primary candidate.

---

## 4. What the contrast cases establish

The "None" candidates are not failures; they are the controls that keep the framework honest.
- **Oliveira CRC (Candidate 2)** is the densest tissue measured and sits deepest below the ceiling, yet its orthogonal data agrees. This is the proof that "above the ceiling" alone does NOT make a conclusion wrong. Being above the ceiling means assignment is unrecoverable, so the spatial data cannot by itself establish the claim; when independent orthogonal data confirms it, the claim stands. The framework correctly does not flag it.
- **Moffitt hybrid neurons (Candidate 3)** is a co-expression / double-positive claim in a moderately hard regime (oracle 0.836), exactly the kind of call misassignment can fabricate, but matched dissociated scRNA-seq independently found the same co-expressing population. Corroborated, so None.
- Together, Candidates 1, 2 and 3 are a clean triplet: same class of claim (co-expression / double-positive), all in imperfect-assignment regimes, separated into artifact-consistent (1) versus corroborated (2, 3) purely by the orthogonal data, never by our re-assignment.

---

## 5. Classification summary (executor does not choose; orchestrator judges)

- **Tier A, confirmed:** Candidate 1 (Janesick triple-positive HER2+/ER+/PR+ DCIS). Above the ceiling by the diagnostic on the exact dataset, with clean independent orthogonal disagreement (matched dissociated scFFPE-seq plus PR-negative clinical status) in the ceiling-predicted direction. Weaknesses: the discordance is already visible in the source paper and the claim is a minor sub-finding.
- **Tier B, shaped but unconfirmed:** Candidate 4 (He CosMx NSCLC tumor-T-cell interface L-R). Assignment-dependent, no orthogonal contradiction available, but the diagnostic was not run on the CosMx flat-file data, so the ceiling position is inferred, not measured. A future run with a CosMx loader would be needed to confirm.
- **None (corroborated or recoverable):** Candidates 2 (Oliveira CRC, above ceiling but orthogonal agrees), 3 (Moffitt, orthogonal agrees), 5 (Eng seqFISH+, orthogonal smFISH agrees).

**Against the pre-committed go/no-go (for the orchestrator's audit, not an executor decision):** a Tier A candidate exists that is diagnostic-confirmed above the ceiling with independent orthogonal data contradicting it in the predicted direction, verifiable from primary sources. Its weaknesses (visible-in-paper discordance, minor sub-finding) are the orchestrator's to weigh against the alternative of shipping the bound plus tool. The diagnostic did not return any dense tumor dataset at or below the ceiling: both Xenium tumors are firmly above it.

---

## 6. Provenance

- New code: `src/lever/lever_diagnostic.py` (reuses Gate 0-3 machinery unchanged).
- New result: `results/lever/lever_ceiling.csv` (breast and CRC rows).
- Raw data (gitignored): `data/xenium_breast_rep1_*` (cached), `data/lever/xenium_crc_addon_*` (downloaded this session from `cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE/`).
- Seeds: MASTER_SEED 20260618; LEVER_SEED 21160618; per-dataset and per-call offsets recorded in the CSV and the script.
- Prior gates untouched: `git status --porcelain` shows only additions; no tracked file modified or deleted.
- Every DOI in this report was confirmed via Crossref; every verbatim quote was pulled from the primary full text (PMC10730913 for Janesick); dataset IDs were confirmed on the 10x, GEO and Dryad pages. Anything that could not be verified was excluded (see below).

**Excluded for non-verifiability or scope:** a public Xenium-bundle accession for the SPLIT breast+lung data (data-availability statement not readable); per-sample h5/parquet for the Oliveira P-samples (only H&E images are individually downloadable; the matrices are inside a 56 GB combined archive, so the diagnostic used the 10x public CRC dataset instead and this is stated, not hidden); a precise published cells/mm^2 figure for any candidate (the diagnostic measures packing directly rather than quoting one); any specific Moffitt excitatory/inhibitory spurious-co-expression artifact claim (not found in the source; not asserted).
