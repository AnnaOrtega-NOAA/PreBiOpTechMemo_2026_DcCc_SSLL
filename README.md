# PreBiOpTechMemo_2026_DcCc_SSLL

Analytical pipeline and codebase supporting the NOAA Fisheries 2026 Pre-Biological Opinion (Pre-BiOp) Technical Memorandum evaluating loggerhead (*Caretta caretta*) and leatherback (*Dermochelys coriacea*) sea turtle interactions in the Hawaii-based **shallow-set longline (SSLL)** fishery.

> **Methodological Note:** The analytical workflows, state-space trend models, and population viability analyses in this repository recreate and extend the baseline methods established in **Martin et al. (2020)** (*NOAA Technical Memorandum NMFS-PIFSC-95*) using updated shallow-set longline interaction records and nesting observations.

---

## Overview

This repository contains the complete R analytical workflow for data preparation, nesting trend state-space modeling, historical shallow-set longline interaction and Adult Nester Equivalent (ANE) estimation, Conway-Maxwell-Poisson (CMP) Anticipated Take Level (ATL) fitting, Stan-based take-demographic modeling, 100-year Population Viability Analysis (PVA), and production figure/table generation.

---

## Repository Structure

```text
PreBiOpTechMemo_2026_DcCc_SSLL/
├── R/
│   ├── 00_config.R                # Global configuration, biological parameters & MCMC settings
│   ├── 01_prepare_data.R          # Data auditing, nesting grid preparation & interaction cleaning
│   ├── 02_nesting_models.R        # Loggerhead singleUQ trend & Leatherback monthly imputation/trend
│   ├── 03_historical_take.R       # Historical SSLL ANE calculation & No-SSLL-Take state-space fits
│   ├── 04_atl.R                   # Annual Conway-Maxwell-Poisson (CMP) Anticipated Take Level (ATL)
│   ├── 05_take_demographics.R     # Stan bivariate model for length (SCL) and post-release mortality
│   ├── 06_pva.R                   # 100-year stochastic Population Viability Analysis (Take vs No-Take)
│   └── 07_tech_memo_figures.R     # High-resolution figure generation and Table 1 summary export
├── inputs/                        # Raw source CSV data (excluded from Git tracking)
├── outputs/                       # Model fits, posterior summaries, audit logs, and figures
├── .gitignore                     # Version control tracking exclusions
└── README.md                      # Project documentation
```

---

## Prerequisites & Dependencies

The pipeline requires **R (v4.0 or higher)** and the following packages:

* **State-Space & Bayesian Modeling:** `jagsUI`, `rstan`
* **Statistical & Distribution Functions:** `truncnorm`, `mvtnorm`, `boot`
* **Visualization & Formatting:** `ggplot2`, `patchwork`
* **Data Manipulation:** `tidyverse` (or core R `stats` / `graphics`)

---

## Workflow Sequence

Execute the scripts sequentially from the root project directory:

1. **`R/00_config.R`**: Establishes central paths, analysis windows (1986–2025 for CC; 2001–2024 for DC), biological parameters (CF, RI, growth constants), JAGS MCMC settings, and PVA parameters.
2. **`R/01_prepare_data.R`**: Audits input CSVs, processes Yakushima loggerhead nesting (Inakahama, Maehama, Yotsusehama) and West Papua leatherback monthly nesting (Jamursba Medi / Jeen Yessa and Wermon / Jeen Syuab), and cleans PIRO shallow-set longline interaction records.
3. **`R/02_nesting_models.R`**: Runs the singleUQ JAGS state-space trend model for loggerheads (1986–2025) and monthly JAGS imputation plus trend models across Median, Low, and High branches for leatherbacks (2001–2024 seasons).
4. **`R/03_historical_take.R`**: Back-calculates historical Adult Nester Equivalents (ANE) from shallow-set interactions, reallocates ANE to nesting beaches, and reruns state-space models on "No-SSLL-Take" nesting series to produce baseline population posteriors.
5. **`R/04_atl.R`**: Fits a numerically stable Conway-Maxwell-Poisson (CMP) distribution evaluated in log space to complete annual shallow-set interaction counts to establish Anticipated Take Levels (ATL).
6. **`R/05_take_demographics.R`**: Fits a Stan bivariate normal model linking log(SCL) and logit(post-release mortality) to annual realized take to parameterize future shallow-set interaction demographics.
7. **`R/06_pva.R`**: Executes a 100-year stochastic Population Viability Analysis ($10,000 \text{ simulations} \times 100 \text{ years}$) comparing Take vs. No-Take trajectories against decline thresholds (50%, 25%, 12.5% of $N_0$), including loggerhead remigration interval ($RI \ge 1$) sensitivity analysis.
8. **`R/07_tech_memo_figures.R`**: Reads saved outputs and exports all technical memorandum figures (400 DPI PNG format) and PVA threshold probability tables without re-fitting models.

---

## Data & Methodological Rules

* **Recreated Baseline Analyses:** Baseline state-space trends, ANE calculations, and PVA scenarios recreate and extend the original methodology of Martin et al. (2020) using updated PIRO interaction datasets and nesting observations from UNIPA, STAJ, and Umigame-kan.
* **Pending Interaction Data:** Loggerhead shallow-set longline interactions for 2022 and 2023 are pending from PIRO; they are currently excluded from analysis.
