# Prolonged gaps, more than the missingness rate, change which imputation method to recommend

Code for a simulation study of ten imputation methods on weekly confirmed malaria
cases from the eleven health districts of the Adamaoua region, Cameroon
(417 weeks, January 2018 to December 2025).

Authors: A. Bagmbaye, A. Batoure Bamana, L. Zhang.
Article: submitted to *BMC Medical Research Methodology* (reference and DOI to be added).

---

## Data availability

The weekly counts come from the DHIS2 system of the Ministry of Public Health of
Cameroon. Two extraction campaigns were used, each made of yearly exports merged
afterwards: the first, in January 2023, covers 2018-2022; the second, in October 2026,
covers 2021-2025. **The counts are not included in this
repository**, under the access conditions set by the Ministry. They are available from the corresponding author on reasonable
request, with the authorisation of the Ministry of Public Health.

For the same reason, the repository leaves out every file that contains the observed
counts or allows them to be rebuilt:

| Not included | Why |
|---|---|
| `data/raw/` | raw DHIS2 extracts |
| `data/derived/reference_adamaoua_weekly_*.csv` | reference series (417 weeks x 11 districts) |
| `data/derived/reference_2018_2025_overlap_differences.csv` | individual observed cells |
| `output/imputed/`, `output/imputed_examples/` | imputed datasets keep the observed cells |
| `output/metrics/metrics_by_district.csv`, `output_mechanism/metrics/district_mechanism.csv`, `output_mechanism/metrics/district_mar_strong.csv` | per-district, per-mask signed errors: combined with the masks, they give linear equations in the true counts |

Everything else is included: the analysis code, the 240 simulated masks, the
aggregated results and the figures produced by R.

---

## Study design

| Factor | Levels |
|---|---|
| Missingness rate | 5 %, 10 %, 20 % |
| Gap length | 1, 4 and 8 consecutive weeks |
| Mechanism | MCAR (main design); MAR, strong MAR and MNAR (reduced design: 10 %, gaps of 1 and 8 weeks, 10 repetitions) |
| Predictor block | `NoTime`, `Seasonal`, `SeasonalTrend` |
| Repetitions | 20 per configuration in the main design (180 masks) |

Ten methods, 24 method-predictor combinations:

- naive: `Mean`;
- cross-sectional: `MICE_inference` (MICE-pmm, m = 5), `MICE_precision` (MICE-pmm, m = 20), `MICE_deterministic` (MICE-norm.predict), `KNN`, `MissForest`;
- temporal: `LinearInterp`, `Kalman`, `SeasonalDecomp`;
- hybrid: `MissForest_ST`.

Criteria: MAE, RMSE, MAPE, bias, variance ratio, R², and the detection of
high-incidence weeks (sensitivity, specificity, Youden index, threshold at each
district's 90th percentile).

---

## Repository layout

```
malaria-imputation-benchmark/
├── run_all.R                       master script
├── R/
│   ├── 00_install_packages.R       run once
│   ├── functions.R                 shared setup, masking, the ten methods, metrics
│   ├── functions_mechanism.R       shared machinery for the mechanism arms
│   ├── 01_build_reference.R        builds and certifies the reference series
│   ├── 02_run_benchmark.R          masks, imputations, metrics (main design)
│   ├── 03_benchmark_analysis.R     tables, Friedman tests, forecasting analysis
│   ├── 04_final_analysis.R         districts, decision guide, peak detection by gap length
│   ├── 05_mechanism_sensitivity.R  MAR and MNAR arms
│   ├── 06_mechanism_diagnostics.R  mechanism diagnostics
│   ├── 07_mar_strong_arm.R         strong MAR arm
│   ├── 10_example_series.R         imputed series of the example figure
│   ├── 11_supplementary_statistics.R  supplementary statistics
│   └── 99_session_info.R           records the software environment
├── data/
│   ├── extend_reference_2018_2025.py  builds the 2018-2025 reference from the raw extracts
│   ├── check_peaks_2023.py            checks the autumn 2023 peaks
│   ├── raw/                           (not included) raw DHIS2 extracts
│   └── derived/                       district context and aggregated checks
├── output/                         main design: masks, aggregated metrics, figures
├── output_mechanism/               mechanism arms: masks, aggregated metrics, figures
├── redaction/
│   ├── figures/make_figures.py         article figures
│   └── supplements/build_supplements.py  additional files 1 to 6
└── docs/                           software environment
```

---

## How to reproduce

1. Obtain the data (see *Data availability*) and place the raw extracts in `data/raw/`
   (file names are given at the top of `R/01_build_reference.R` and
   `data/extend_reference_2018_2025.py`).
2. Build the reference series `data/derived/reference_adamaoua_weekly_2018_2025.csv`
   with `R/01_build_reference.R` and `data/extend_reference_2018_2025.py`.
3. Install the R packages once, then run the pipeline:
   ```r
   source("R/00_install_packages.R")
   source("run_all.R")
   ```
4. Draw the figures and build the additional files:
   ```
   python redaction/figures/make_figures.py
   python redaction/supplements/build_supplements.py
   ```

Every mask of the main design is drawn with a seed derived from its position in
the design (`1000 * rep + rate + 7 * gap`). All masks are saved in `output/masks/` and
`output_mechanism/masks/`, so the imputed datasets can be regenerated exactly once
the reference series is in place. Software versions are listed in
`docs/ENVIRONMENT.md` and `docs/package_versions.csv` (R 4.6.1).

---

## What can be checked without the data

The scripts read the reference series when they start, so they do not run without it.
The aggregated outputs are enough to recompute the rankings, the Friedman and Nemenyi
tests, the excess-change measure and the paired comparisons:
`output/metrics/metrics_long.csv` (one row per repetition, rate, gap length and
method-predictor combination) and `output/metrics/selected_variant.csv`.
`output/metrics/peak_weeks_per_mask.csv` gives only a summary by rate and gap length.

---

## Licence

Code: MIT (see `LICENSE`). The licence does not cover the DHIS2 data, which are
not distributed here.

## Citation

To be added after publication.
