# Examining 400 years of epidemics and pandemics

This repository contains the complete analysis code for the paper "Examining 400 years of epidemics and pandemics" by Horacio Chacon-Torrico and Christopher J. L. Murray (Institute for Health Metrics and Evaluation, University of Washington).

## Research Summary

We forecast that pandemics will cause 74.9 million deaths (95% uncertainty interval 18.0–233) and 3.66 billion DALYs between 2026 and 2100, by integrating four centuries of global mortality data (105 epidemics and pandemics from 1600–2025) with extreme-value statistical models. Our approach uses a two-component mixture model combining truncated lognormal distribution for moderate events and Generalized Pareto Distribution (GPD) for extreme tail events, with time-varying occurrence rates and rigorous out-of-sample validation.

## Model Overview

**Key Methodological Components:**
- **Two-stage mixture model**: Truncated lognormal (100k-3M population-scaled deaths) + GPD tail (>3M deaths)
- **Population scaling**: All historical deaths standardized to 2025 global population (8.2 billion)
- **Shadow mean calculation**: Finite tail expectations using Cirillo & Taleb (2020) approach for heavy-tailed distributions
- **Time-varying rates**: Poisson process with 25-year sliding windows to estimate declining pandemic frequency
- **Bootstrap uncertainty**: 10,000 draws for complete uncertainty quantification
- **Historical validation**: Training on pre-1950 data successfully predicts 1951-2025 pandemic mortality

## Requirements

**R Environment:**
- R ≥ 4.0.0
- Required packages: `tidyverse`, `mev`, `evir`, `config`, `ggplot2`, `gridExtra`, `scales`, `ggrepel`, `patchwork`, `lubridate`, `MASS`, `cowplot`, `viridis`
- Optional packages for sensitivity analysis: `future`, `future.apply`, `foreach`, `doParallel`, `parallel`

**Configuration:**
All analysis parameters are pre-configured in `config.yml` with publication-ready default values. No manual parameter adjustment needed for reproduction.

## Reproducing the Analysis

### Core Results (Table 1, Figures 1-3)
```bash
Rscript src/run_PoT.R
```
Generates main pandemic burden estimates, severity distribution plots, and rate projections.

### Out-of-Sample Validation (Figure 4)
```bash
Rscript src/run_validation.R
```
Reproduces validation showing model trained on 1600-1950 data predicting 1951-2025 outcomes.

### Sensitivity Analysis (Extended Data Figures 5-7)
```bash
Rscript src/run_sens.R
```
Tests robustness across threshold choices, time periods, and cutoff values.

### Model Diagnostics
```bash
# Bulk distribution diagnostics (Extended Data Table 2)
Rscript src/bulk_diagnostics.R

# Tail threshold selection (Extended Data Figures 2-3)
Rscript src/tail_diagnostics.R
```

## Project Structure

```
pandemic_estimation/
├── R/                      # Modular analysis functions
│   ├── analysis.R          # Core analysis orchestration
│   ├── model_fitting.R     # Mixture model implementation
│   ├── bootstrap.R         # Uncertainty quantification
│   ├── burden.R            # DALY calculations and scenarios
│   ├── plotting.R          # Visualization functions
│   ├── validation.R        # Out-of-sample testing
│   ├── sensitivity.R       # Sensitivity analysis functions
│   └── ...
├── src/                    # Execution scripts
│   ├── run_PoT.R          # Main analysis
│   ├── run_validation.R   # Validation analysis
│   ├── run_sens.R         # Sensitivity analysis
│   └── *_diagnostics.R    # Model diagnostic scripts
├── data/                   # Input datasets
│   └── epidemics_data_2025.csv  # Historical pandemic data
├── output/                 # Generated results
├── config.yml             # Analysis parameters
└── README.md
```

## Configuration Modes

The `config.yml` file includes multiple analysis configurations:

- **`default`**: Standard analysis for main results
- **`validation`**: Out-of-sample testing parameters  
- **`sensitivity`**: Multi-parameter sensitivity testing

Scripts automatically select appropriate configurations. Advanced users can modify parameters in `config.yml` or set environment variables:

```bash
# Use specific configuration
export R_CONFIG_ACTIVE=validation
Rscript src/run_PoT.R
```

## Key Parameters (Pre-configured)

- **Population reference**: 8.2 billion (2025 global population)
- **Analysis threshold**: 100,000 population-scaled deaths
- **Tail threshold**: 3 million population-scaled deaths  
- **Bootstrap draws**: 10,000 for uncertainty intervals
- **Forecast horizon**: 2026-2100 (75 years)
- **Historical data**: 1600-2025 (425 years)
- **Scenario windows**: 75-year periods for best/worst rate scenarios

## Expected Outputs

**Main Results:**
- Summary tables with point estimates and 95% uncertainty intervals
- Pandemic severity distribution plots
- Time-varying rate projections
- Scenario comparison tables
- DALY burden estimates

**Validation:**
- Out-of-sample prediction accuracy metrics
- Observed vs. predicted cumulative mortality plots

**Sensitivity Analysis:**
- Threshold sensitivity plots
- Robustness assessment tables

