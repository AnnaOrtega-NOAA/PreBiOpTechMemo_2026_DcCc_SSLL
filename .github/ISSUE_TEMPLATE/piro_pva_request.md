---
name: PIRO PVA Request
about: Workflow checklist for tracking a PIRO marine turtle PVA assessment
title: '[PVA] <Species> - <Fishery / Action>'
labels: 'PVA, PIRO-Request'
assignees: ''
---

## Request Overview
- **PIRO Divisions:** Protected Resources Division (PRD) & Sustainable Fisheries (SF)
- **PIRO Contacts:** <!-- Add PIRO PRD / SF lead names -->
- **Target Species / Population:** <!-- North Pacific Loggerhead, Western Pacific Leatherback, etc. -->
- **Fishery / Action:** <!-- e.g., Hawaiʻi Shallow-Set Longline -->
- **Target Delivery Date:** 

## Data Gathering & Partner Coordination
*Reminder: Unaggregated observer logs and partner nesting datasets must remain on local/network drives and NEVER be committed to GitHub.*

- [ ] **Fisheries Data:** Coordinate with PIRO (PRD & SF) for updated observer interaction records
- [ ] **Loggerhead Nesting Data:** Coordinate with Sea Turtle Association of Japan (STAJ) & Yakushima Umigame-kan
- [ ] **Leatherback Nesting Data:** Coordinate with Universitas Papua (UNIPA)
- [ ] Place raw CSVs in local `inputs/` workspace (verified `.gitignore` protection)

## Analytical Pipeline Milestones
- [ ] **Data Prep & Audit:** Run `R/01_prepare_data.R` and audit output summaries
- [ ] **Nesting Models:** Run `R/02_nesting_models.R` (JAGS) for target species
- [ ] **Historical Take:** Run `R/03_historical_take.R` to compute ANE and no-historical-take posteriors
- [ ] **Anticipated Take Level:** Run `R/04_atl.R` to fit updated CMP distributions
- [ ] **Take Demographics:** Run `R/05_take_demographics.R` (Stan)
- [ ] **PVA Simulations:** Run `R/06_pva.R` for 100-year stochastic projections
- [ ] **Reporting & Figures:** Run `R/07_tech_memo_figures.R` to generate figures and tables

## Review, Clearance & Publication
- [ ] Draft Technical Memorandum (text, figures, and analytical code) submitted for **PSD Review**
- [ ] Address PSD review comments and finalize manuscript
- [ ] Submit package for formal **PIFSC Internal Center Clearance**
- [ ] Transmit final report package to **PIRO PRD & SF**
- [ ] Deposit published Tech Memo in the **NOAA Institutional Repository**
