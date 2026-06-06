# AGB-MetaAnalysis-2026

**Advances in forest Above-Ground Biomass estimation from Remote Sensing:  
A Meta-Analysis of GEDI LiDAR and Multi-Sensor Approaches**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20575421.svg)](https://doi.org/10.5281/zenodo.20575421)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.5.3-blue)](https://www.r-project.org/)

---

## Overview

This repository contains all data, analysis scripts, and outputs supporting the systematic review and meta-analysis:

> **Musumba Teso P., Okole N., Birindwa J., Muhindo Sahani W., Mwanjalolo M.J.G., Karume K. (2026).**  
> *Advances in forest Above-Ground Biomass estimation from Remote Sensing: A Meta-Analysis of GEDI LiDAR and Multi-Sensor Approaches.*  
> ISPRS Journal of Photogrammetry and Remote Sensing. [Under review]

### Key findings

| Metric | Value |
|--------|-------|
| Studies included | 77 (70 Fisher's Z · 77 RMSE) |
| Observations analysed | 314 Fisher's Z · 308 RMSE |
| Countries represented | 25 |
| Publication period | 2008–2026 |
| Pooled Fisher's Z | 1.234 (95% CI: 1.161–1.308) |
| Equivalent R² | ≈ 0.71 |
| Between-study heterogeneity | I² > 99% |
| Strongest moderator | Sensor combination (ε² = 0.181) |
| Best sensor configuration | LiDAR + Optical (R² = 0.77) |
| Publication bias | None detected (Egger p = 0.559) |

---

## Repository structure

AGB-MetaAnalysis-2026/
│
├── R/                                  # Analysis scripts (run in order)
│   ├── 00_setup.R                      # Package loading and global seed (set.seed(2024))
│   ├── 01_data_import.R                # Data loading and validation
│   ├── 02_descriptive_stats.R          # Summary statistics and dataset overview
│   ├── 03_fisher_z_meta.R              # Primary meta-analysis (rma.mv, REML) + LOSO
│   ├── 04_nonparametric_tests.R        # Kruskal-Wallis and Mann-Whitney tests
│   ├── 04b_VIF_diagnostics.R           # VIF, Cramér's V, η² moderator screening
│   ├── 04c_distributional_diagnostics.R # Shapiro-Wilk, Levene, skewness/kurtosis
│   ├── 05_meta_regression.R            # Multilevel meta-regression (5 moderators)
│   ├── 05b_dharma_diagnostics.R        # DHARMa residual diagnostics (1,000 simulations)
│   ├── 05c_random_structure_comparison.R # Random-effects structure comparison
│   ├── 06_publication_bias.R           # Egger's test + trim-and-fill
│   ├── 07_figures.R                    # All manuscript and supplementary figures
│   ├── 08_tables.R                     # All manuscript and supplementary tables
│   └── 09_sensitivity.R               # Scenario-based sensitivity analyses (S1–S8)
│
├── data/                               # Input datasets
│   ├── DATA_FisherZ_latest.csv         # Fisher's Z dataset (314 obs, 70 studies)
│   ├── DATA_RMS_latest.csv             # RMSE dataset (308 obs, 77 studies)
│   └── AGB_master_database.csv         # Full master database (342 obs, 38 variables)
│
├── output/
│   ├── tables/                         # CSV outputs from all analysis scripts
│   └── figures/                        # PNG outputs (manuscript + supplementary)
│
├── CHANGES.md                          # Data correction log
├── LICENSE                             # CC BY 4.0
└── README.md                           # This file
---

## Data description

### `DATA_FisherZ_latest.csv` — primary effect-size dataset

| Column | Type | Description |
|--------|------|-------------|
| `Obs_ID` | character | Unique observation identifier (M1, M2, …) |
| `Study_ID` | character | Unique study identifier (S1, S2, …) |
| `Reference` | character | First author and year |
| `Year` | integer | Publication year |
| `Country` | character | Country of study |
| `Climate_zone` | character | Köppen–Geiger climate zone |
| `Forest_type` | character | Dominant forest type |
| `Sensor_combination` | character | Sensor modalities used as model inputs |
| `Model_category` | character | Statistical / Machine-learning / Deep-learning |
| `Model_algorithm` | character | Specific algorithm (e.g., Random Forest, CNN) |
| `GEDI_used` | character | Yes / No — whether GEDI data were used as predictor |
| `GEDI_product` | character | Specific GEDI product (L2A, L2B, L4A, L4B, etc.) |
| `Study_scale` | character | Local / Regional / National / Multi-country |
| `Training_samples` | integer | Number of ground reference plots |
| `Spatial_resolution_m` | numeric | Spatial resolution of primary input sensor (metres) |
| `R2` | numeric | Reported coefficient of determination |
| `r` | numeric | Pearson correlation coefficient |
| `Fisher_Z` | numeric | Fisher's Z transformed correlation |
| `vi` | numeric | Sampling variance = 1 / (Training_samples − 3) |

### Sensor combination distribution (Fisher's Z dataset)

| Sensor combination | n obs | n studies |
|--------------------|-------|-----------|
| Optical-only | 75 | 14 |
| Optical + SAR + LiDAR | 68 | 13 |
| LiDAR + Optical | 58 | 15 |
| LiDAR-only | 58 | 13 |
| Optical + SAR | 28 | 9 |
| SAR-only | 21 | 4 |
| SAR + LiDAR | 6 | 2 |

---

## How to reproduce the analysis

### Requirements

```r
install.packages(c(
  "metafor", "robumeta", "DHARMa", "dplyr", "tidyr",
  "ggplot2", "ggdist", "patchwork", "dunn.test", "rstatix",
  "car", "lmtest", "readr", "openxlsx", "viridis"
))
```

### Run order

```r
source("R/00_setup.R")
source("R/01_data_import.R")
source("R/02_descriptive_stats.R")
source("R/03_fisher_z_meta.R")
source("R/04_nonparametric_tests.R")
source("R/04b_VIF_diagnostics.R")
source("R/04c_distributional_diagnostics.R")
source("R/05_meta_regression.R")
source("R/05b_dharma_diagnostics.R")
source("R/05c_random_structure_comparison.R")
source("R/06_publication_bias.R")
source("R/07_figures.R")
source("R/08_tables.R")
source("R/09_sensitivity.R")
```

> A global random seed `set.seed(2024)` is applied in `00_setup.R` to ensure full reproducibility.

---

## Citation

```bibtex
@article{musumbateso2026agb,
  author  = {Musumba Teso, Philippe and Okole, Nathan and Birindwa, Jovianne
             and Muhindo Sahani, Walere and Mwanjalolo, Majaliwa Jackson Gilbert
             and Karume, Katcho},
  title   = {Advances in forest Above-Ground Biomass estimation from Remote Sensing:
             A Meta-Analysis of {GEDI} {LiDAR} and Multi-Sensor Approaches},
  journal = {ISPRS Journal of Photogrammetry and Remote Sensing},
  year    = {2026},
  note    = {Under review},
  doi     = {10.5281/zenodo.20575421}
}
```

---

## Data corrections log

See [`CHANGES.md`](CHANGES.md) for corrections applied:
- 9 duplicate observations removed (S6, S89, S90, S101, S111, S308–S311)
- Yang et al. 2025 disambiguated into 2025a, 2025b, 2025c
- Final corrected dataset: v_corrected (2026-06-06)

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — free to share and adapt with attribution.

---

## Contact

**Musumba Teso Philippe**  
Doctoral Researcher — École Doctorale d'Agroécologie et Sciences du Climat  
Executive Director — Baraka Climat ASBL, Goma/Bukavu, DR Congo  
✉ musumbateso@yahoo.fr  
🔗 https://github.com/musumbateso1975

---

*Archived on Zenodo: https://doi.org/10.5281/zenodo.20575421*
