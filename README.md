# pkbm-4c-plssem

Analysis code for the article:

> Ansar A, Hidayatullah F, Puspitasari SW, Latang L, Gaffar F. **Tutor competence as the pathway from learning service quality to 4C skills: Evidence from urban, coastal, and highland community learning centers in Indonesia.** *F1000Research* (2026).

The published estimates were produced with SmartPLS 4 (proprietary). This repository provides an open-source R script, based on the CRAN package [seminr](https://CRAN.R-project.org/package=seminr), that implements the same model and analysis steps so that every step can be inspected and rerun.

## What the script does

`R/analysis.R` estimates the model learning service quality (X, `LSQ`) → tutor competence (Z, `TC`) → 4C skills (Y, `FourC`), with a direct path X → Y, and reports:

| Step | Output file (in `output/`) | Article |
|---|---|---|
| Initial model with all 59 items (basis for item purification) | `01_loadings_initial_59items.csv` | Methods, Table 1 |
| Outer loadings of the 35 retained indicators | `02_outer_loadings.csv` | Table 3 |
| Cronbach's α, ρA, ρC, AVE | `03_reliability_validity.csv` | Table 3 |
| Fornell–Larcker criterion and HTMT | `04a_fornell_larcker.csv`, `04b_htmt.csv` | Table 4 |
| Full collinearity VIF (common method bias; Kock, 2015) | `05_full_collinearity_vif.csv` | Results |
| Inner VIF, R², adjusted R², f² | `06a`–`06c` | Tables 5–6 |
| Bootstrapping, 5,000 resamples; 95% CI = β ± 1.96 SD (percentile CI also given) | `07_hypothesis_tests_direct.csv` | Table 6 (H1–H3) |
| Indirect effect X → Z → Y, mediation type, upsilon (υ) | `08_indirect_effect.csv` | Table 6 (H4) |
| Q²predict, RMSE, MAE (10-fold, 10 repetitions) | `09_explanatory_predictive_power.csv` | Table 5 |

## Data

The respondents included learners younger than 18 years, and the ethics approval and consent did not cover public release of the data. The real item-level dataset is therefore **not** included. It is available to qualified researchers from the corresponding author (Fauzan Hidayatullah, fauzan.hidayatullah@unm.ac.id) under a data-use agreement; see the Data availability statement of the article.

To let anyone run the code, `data/example_data_simulated.csv` contains **simulated** data (354 rows, 59 items, 1–10 scale) with the same column names as the real dataset (`X1`–`X24`, `Z1`–`Z20`, `Y1`–`Y15`). Results from the simulated data do **not** reproduce the published values.

## How to run

Requirements: R ≥ 4.1 and the package `seminr` (installed automatically if missing).

```sh
# simulated example data
Rscript R/analysis.R

# real data (after obtaining it from the corresponding author)
Rscript R/analysis.R path/to/pkbm_data.csv
```

Because seminr and SmartPLS use different random resampling and cross-validation routines, bootstrap standard errors and Q²predict may differ from the published values in the third decimal place. Point estimates (loadings, path coefficients, R², f²) should match.

## License

MIT License (OSI-approved), see [`LICENSE`](LICENSE).

## How to cite

Please cite the article above and the archived version of this code (see `CITATION.cff`; DOI from Zenodo).
