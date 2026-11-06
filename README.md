# PreBiOpTechMemo_2026_DcCc_SSLL

Analytical pipeline and codebase supporting the NOAA Fisheries 2026 Pre-Biological Opinion (Pre-BiOp) Technical Memorandum evaluating North Pacific loggerhead (*Caretta caretta*) and West Pacific leatherback (*Dermochelys coriacea*) sea turtle interactions in the Hawaiʻi-based **shallow-set longline (SSLL)** fishery.

> **Methodological Note:** The analytical workflows, state-space trend models, and population viability analyses in this repository recreate and extend the baseline methods established in the **2020 assessment** (Martin et al. 2020) using updated shallow-set longline interaction records and nesting observations through 2025.

---

## Overview

This repository contains the complete R analytical workflow for data preparation, nesting trend state-space modeling, historical shallow-set longline interaction and Adult Nester Equivalent (ANE) estimation, Conway-Maxwell-Poisson (CMP) Anticipated Take Level (ATL) fitting, Stan-based take-demographic modeling, 10-year multi-scenario Population Viability Analysis (PVA), and production figure/table generation.

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
│   ├── 05_take_demographics.R     # Stan bivariate model for length (SCL) and post-interaction mortality
│   ├── 06_pva.R                   # 10-year multi-scenario stochastic PVA (Status Quo & Counterfactuals)
│   └── 07_tech_memo_figures.R     # High-resolution figure generation (Figs 1–7) & summary table exports
├── inputs/                        # Raw source CSV data (confidential observer records; excluded from Git)
├── outputs/                       # Model fits, posterior summaries, audit logs, and publication figures
├── .gitignore                     # Version control tracking exclusions
└── README.md                      # Project documentation
```

---

## Prerequisites & Dependencies

The pipeline requires **R (v4.0 or higher)** and the following packages:

- **State-Space & Bayesian Modeling:** `jagsUI`, `rstan`
- **Statistical & Distribution Functions:** `truncnorm`, `mvtnorm`, `boot`
- **Visualization & Formatting:** `ggplot2`, `patchwork`
- **Data Manipulation:** `stats`, `utils`, `graphics` (or `tidyverse`)
- **External MCMC Dependency:** `JAGS` (Just Another Gibbs Sampler v4.3.0+)

---

## Workflow Sequence

Execute the scripts sequentially from the root project directory:

1. **`R/00_config.R`**: Establishes central paths, analysis windows (1986–2025 for loggerheads; 2001–2025 for leatherbacks), biological parameters (CF, RI, growth constants), JAGS MCMC settings, and PVA simulation horizons.
2. **`R/01_prepare_data.R`**: Audits input CSVs, processes Yakushima loggerhead nesting (Inakahama, Maehama, Yotsusehama; 2021 excluded) and West Papua leatherback monthly nesting (Jeen Yessa and Jeen Syuab), and cleans PIRO shallow-set longline interaction records (2005–2025 complete calendar years).
3. **`R/02_nesting_models.R`**: Runs the singleUQ JAGS state-space trend model for loggerheads (1986–2025) and monthly JAGS imputation plus trend models across Median, Low, and High branches for leatherbacks (2001–2025 seasons; 1,000,000 MCMC iterations).
4. **`R/03_historical_take.R`**: Back-calculates historical Adult Nester Equivalents (ANE; 8.02 for loggerheads, 39.33 for leatherbacks), reallocates ANE to nesting beaches, and fits state-space models to "No-SSLL-Take" counterfactual nesting series to establish baseline starting states.
5. **`R/04_atl.R`**: Fits Conway-Maxwell-Poisson (CMP) distributions on the log scale to complete annual interaction totals (2005–2025) to parameterize status quo ATL distributions and counterfactual uncapped CMP distributions.
6. **`R/05_take_demographics.R`**: Fits a Stan bivariate normal model linking log(SCL) and logit(post-interaction mortality) to annual realized interactions to parameterize future shallow-set take demographics.
7. **`R/06_pva.R`**: Executes 10,000 stochastic population projections over a 10-year consultation horizon comparing five scenarios: (1) No Future SSLL Take, (2) Status Quo Take, (3) No Trip Limits, (4) No Hard Caps, and (5) Both Removed. Computes depletion risk probabilities below 50%, 25%, and 12.5% of starting female nester abundance ($N_0$).
8. **`R/07_tech_memo_figures.R`**: Reads existing saved RDS/CSV outputs and exports all technical memorandum figures (Figures 1–7; 400 DPI PNG format) and PVA threshold probability summary tables (Table 1 Status Quo and Table 1B Multi-Scenario Comparison) without re-fitting models.

---

## Data & Methodological Rules

- **Recreated Baseline Analyses:** Baseline state-space trends, ANE calculations, and PVA frameworks recreate and extend the foundational methodology of the 2020 assessment (Martin et al. 2020) using updated nesting observations from UNIPA, STAJ, and Yakushima Umigame-kan.
- **Index Nesting Beach Definitions:** Consistent with Toha et al. (2025) and the 2020 assessment, Jeen Yessa comprises three contiguous monitoring beaches (Warmamedi, Batu Rumah, Wembrak) treated as a single index area, while Jeen Syuab represents a distinct index beach ~30 km to the east.
- **Leatherback Imputation Coverage:** Monthly nesting counts reflect Bayesian imputations across unmonitored months (17.0% of months imputed for Jeen Yessa [2001–2025] and 7.1% of months imputed for Jeen Syuab during its active monitoring period [2006–2025]).
- **Fishery Interaction Series:** Evaluates 21 complete calendar years (2005–2025) under ~100% observer coverage.
- **Data Confidentiality:** To protect proprietary vessel locations and comply with NOAA data privacy policies, raw observer interaction records and precise coordinates are excluded from public version control. Authorized users within NOAA Fisheries can request access via the NOAA Fisheries Pacific Islands Regional Office (PIRO) Observer Program.

---

## References

Martin, S. L., Siders, Z. A., Eguchi, T., Langseth, B., Yau, A., Baker, J., Ahrens, R. N. M., & Jones, T. T. (2020). *Assessing the population-level impacts of North Pacific loggerhead and western Pacific leatherback turtle interactions in the Hawaii-based shallow-set longline fishery* (NOAA Technical Memorandum NOAA-TM-NMFS-PIFSC-95). Pacific Islands Fisheries Science Center. https://doi.org/10.25923/ydp1-f891

---

## NOAA Open Source & Scientific Integrity Disclaimers

**Disclaimer:** This repository is a scientific product and is not an official communication of the National Oceanic and Atmospheric Administration (NOAA) or the United States Department of Commerce. All NOAA GitHub code is provided on an "as is" basis, and the user assumes responsibility for its use. NMFS and NOAA make no representation or warranty, express or implied, including the warranties of merchantability and fitness for a particular purpose, as to the completeness, accuracy, or utility of this software or its outputs.

**Commercial Endorsement:** Mention of any trade names, commercial products, or external software libraries throughout this codebase and documentation is for identification purposes only and does not imply endorsement by the National Marine Fisheries Service, NOAA, or the U.S. Department of Commerce.
