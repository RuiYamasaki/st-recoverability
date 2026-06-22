# A method-independent recoverability ceiling for transcript-to-cell assignment in imaging-based spatial transcriptomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20760704.svg)](https://doi.org/10.5281/zenodo.20760704)
https://doi.org/10.21203/rs.3.rs-10093252/v1 

This repository maps a **method-independent recoverability ceiling** for transcript-to-cell assignment in imaging-based spatial transcriptomics (Xenium, MERFISH, and similar). Using a by-construction synthetic ground-truth generator and a Bayes-optimal ("oracle") assignment, it characterizes when assigning transcripts to their cell of origin is recoverable at all, as a function of transcript density, cell packing, and molecular displacement. The displacement is pinned to real Xenium and MERFISH data, and the resulting ceiling shows that in densely packed tissue (around 13,600 cells/mm^2, realistic for breast, lung, and colorectal tumor) the best-possible assignment accuracy falls below 0.9: assignment there is unrecoverable by any method. Real published segmenters (Baysor, Proseg, ComSeg) given true nuclei sit below this ceiling in dense tissue and do not exceed a naive nearest-nucleus baseline by more than the measured run-to-run noise. A worked diagnostic applies the ceiling to a published Xenium triple-positive co-expression claim and places it, together with the orthogonal data that bears on it, on principled footing.

## Repository structure
```
.
├── README.md
├── LICENSE                            # MIT
├── CITATION.cff
├── STATUS.md                          # decision log and project history (D1 to D14)
├── LOCKED_NUMBERS_AND_GUARDRAILS.md   # verified headline numbers and guardrails
├── requirements.txt                   # direct dependencies
├── requirements.lock.txt              # exact pinned dependency closure
├── src/                               # generator, oracle, gates, headroom, lever, figures
│   ├── generator.py  oracle.py  config.py  expression.py  metrics.py  evaluate.py
│   ├── sweep.py  realism.py  anchor.py  run_all.py        # Gate 0
│   ├── gate1_*.py  gate2_*.py  gate3_*.py                 # robustness gates
│   ├── headroom_*.py  methods_*.py                        # real-segmenter comparison
│   ├── lever/                                             # worked diagnostic
│   ├── figures/                                           # main + Extended Data figures
│   └── demo.py                                            # quick self-contained demo
├── results/                           # committed: small derived CSV / JSON only
├── figures/                           # committed: generated figures (PNG + PDF)
├── outputs/                           # committed: per-stage reports
└── data/                              # NOT committed: raw public data, streamed and cached
```

## Installation
Python 3.11 or newer.
```bash
git clone https://github.com/RuiYamasaki/st-recoverability.git
cd st-recoverability
pip install -r requirements.lock.txt   # exact pinned versions used for the results
# or: pip install -r requirements.txt   # latest compatible versions
```
External segmenters used in the headroom comparison are not pip-installable and are optional (only needed to re-run that one stage); see "Segmenter versions" below.

## Quick demo (under 1 minute, no data download)
```bash
python src/demo.py
```
Builds two synthetic fields (sparse and dense) at a realistic displacement and prints the oracle versus naive assignment accuracy, writing `figures/demo_oracle_vs_naive.png`. Expected: oracle accuracy near 0.87 in the sparse regime and around 0.72 in the dense regime, illustrating that the dense-tissue ceiling sits below 0.9.

## Reproducing the results
All randomness uses fixed, recorded seeds. Raw data is streamed and cached under `data/` (never committed); only small derived files live in `results/`.

- **Gate 0 (the frontier), one command:** `python src/run_all.py` re-derives `results/sweep.csv`, `results/realism.csv`, `results/realism_oracle.csv` and the Gate 0 figures from code and cached data. Use `--skip-realism` to avoid the network if `results/realism.csv` is present.
- **Robustness gates:** `python src/gate1_run_all.py`, `gate2_run_all.py`, `gate3_run_all.py` (overlapping expression, data-pinned displacement with bootstrap CIs, marker-validity permutation test, independent-tissue replication).
- **Worked diagnostic:** `python src/lever/lever_diagnostic.py`, `lever_roi.py`, `lever_mechanism.py`.
- **Figures:** `python src/figures/fig1_frontier.py` (and `fig2_bound.py`, `fig3_methods.py`, `fig4_worked.py`, `ed_render.py`, `leg2_local.py`).

The analysis stages are deterministic and reproduce the committed CSVs byte-for-byte on the pinned environment (verified for the Gate 0 sweep and anchor). The real-segmenter headroom stage requires the external tools below and is the one stage not reproducible without them.

## Data sources (exact identifiers)
Only small derived summary statistics are committed; raw data is downloaded from the original repositories.

| dataset | identifier | reference |
|---|---|---|
| Xenium human breast cancer (Rep 1) | 10x Genomics `Xenium_FFPE_Human_Breast_Cancer_Rep1`, release 1.0.1; GEO **GSE243280** | Janesick et al., Nat. Commun. 14:8353 (2023), doi:10.1038/s41467-023-43458-x |
| Xenium human lung cancer | 10x Genomics sample `Xenium_V1_hLung_cancer_section`, release 1.5.0 | 10x Genomics public dataset |
| Xenium human colorectal cancer | 10x Genomics sample `Xenium_V1_Human_Colorectal_Cancer_Addon_FFPE`, Xenium Onboard Analysis 2.0.0 | 10x Genomics public dataset |
| MERFISH mouse hypothalamic preoptic region | Dryad **doi:10.5061/dryad.8t8s248** (cached via the squidpy `MERFISH_0.24.h5ad` mirror) | Moffitt et al., Science 362:eaau5324 (2018) |

## Segmenter versions (external tools, not pip-installable)
The headroom comparison evaluated real published methods against the oracle on known-truth synthetic data. The fairness-gate (Linux) run used:
- **Baysor 0.7.1** (Petukhov et al., Nat. Biotechnol. 2022)
- **Proseg 3.1.1** (Newell lab, Nat. Methods 2025; via bioconda)
- **ComSeg 1.8.5** (Defard et al. 2024; pip)

An earlier Windows run additionally used pciSeq 0.0.61. Installation specifics and outcomes are recorded in `outputs/HEADROOM_LINUX_REPORT.md` and `outputs/HEADROOM_FAIR_REPORT.md`.

## Expected outputs
- `results/*.csv`, `results/lever/*.csv`: the frontier grid, data-pinned displacement and CIs, dense-regime ceiling, method metrics, and the worked-example numbers.
- `figures/fig1_frontier.*` to `fig4_worked.*`: the four main figures (PNG + vector PDF).
- `figures/ED01_*` to `ED10_*`: the ten Extended Data items.
- `outputs/*REPORT.md`: a report per stage, with every headline number traced to its source file.

## License
MIT (see [LICENSE](LICENSE)). MIT was chosen as a permissive default; it can be swapped for BSD-3-Clause or Apache-2.0 if preferred.

## Citation
If you use this software or its results, please cite it (see [CITATION.cff](CITATION.cff)):

> Yamasaki, R. (2026). A method-independent recoverability ceiling for transcript-to-cell assignment in imaging-based spatial transcriptomics. Zenodo. https://doi.org/10.5281/zenodo.20760704

The DOI [10.5281/zenodo.20760704](https://doi.org/10.5281/zenodo.20760704) is the Zenodo concept DOI: it always resolves to the latest version. A preprint/article DOI: https://doi.org/10.21203/rs.3.rs-10093252/v1 
