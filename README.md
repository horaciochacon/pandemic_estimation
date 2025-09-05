# Examining 400 years of epidemics and pandemics

This repository contains the complete analysis code for the paper "Examining 400 years of epidemics and pandemics" by Horacio Chacon-Torrico and Christopher J. L. Murray (Institute for Health Metrics and Evaluation, University of Washington).

## Data Availability

All required data is included in the `data/` folder:
- **Pandemic/epidemic data**: Historical mortality data from 1600-2025 (`epidemics_data_2025.csv`)
- **Population forecasts**: IHME population projections at draw level (`pop_draws.csv`, `pop_forecast.csv`)
- **DALY ratios**: Disability-Adjusted Life Years to death ratios (`daly_death_ratio_draws.rds`, `dalys_draws.rds`, `deaths_draws.rds`)

This dataset supports all analysis sections including burden calculations, validation, and sensitivity analyses.

## Model Overview

**Key Methodological Components:**
- **Two-stage mixture model**: Truncated lognormal (100k-3M population-scaled deaths) + GPD tail (>3M deaths)
- **Population scaling**: All historical deaths standardized to 2025 global population (8.2 billion)
- **Shadow mean calculation**: Finite tail expectations using Cirillo & Taleb (2020) approach for heavy-tailed distributions
- **Time-varying rates**: Poisson process with 25-year sliding windows to estimate declining pandemic frequency
- **Bootstrap uncertainty**: 10,000 draws for complete uncertainty quantification
- **Historical validation**: Training on pre-1950 data successfully predicts 1951-2025 pandemic mortality

## System Requirements

**Tested Operating Systems:**
- macOS (current version)
- Linux Ubuntu

**R Environment:**
- R ≥ 4.0.0
- Required packages (install all with command below)

**Installation:**
Install all required R packages with a single command:
```r
install.packages(c(
  "tidyverse", "dplyr", "purrr", "stringr", "ggplot2", "ggrepel", "patchwork",
  "gridExtra", "scales", "data.table", "config", "evir", "MASS", "stats",
  "lubridate", "cowplot", "viridis", "mev",
  # Optional packages for sensitivity analysis (parallel processing)
  "future", "future.apply", "foreach", "doParallel", "parallel"
))
```

**Installation time:** Approximately 5-10 minutes on a typical desktop computer.

**Configuration:**
All analysis parameters are pre-configured in `config.yml` with publication-ready default values. No manual parameter adjustment needed for reproduction. For custom data analysis, see [Configuration Guide](CONFIG.md).

## Reproducing All Manuscript Results

### Complete Analysis (Recommended)
To reproduce ALL manuscript results with a single command:
```bash
Rscript src/run_all.R
```

**What this generates:**
- Main results: Table 1, Figures 1-3 (pandemic burden estimates, severity distributions, rate projections)
- Validation: Figure 4 (out-of-sample validation showing model trained on 1600-1950 predicting 1951-2025)
- Sensitivity analysis: Extended Data Figures 5-7 (robustness across thresholds and time periods)
- Model diagnostics: Extended Data Table 2, Figures 2-3 (bulk/tail distribution diagnostics)

**Expected runtime:** ~10 minutes on MacBook Pro M1 Pro, may vary on other systems.

**Output structure:**
```
output/run_all_2024-09-04_14-30-25/
├── run_pot/           # Main analysis results
├── run_validation/    # Validation results
├── run_sens/          # Sensitivity analysis
├── bulk_diagnostics/  # Bulk distribution diagnostics
├── tail_diagnostics/  # Tail threshold diagnostics
├── execution_log.csv  # Detailed run log with timing
├── session_info.txt   # R environment information
└── run_summary.md     # Human-readable execution summary
```

### Individual Analysis Components
For running specific analysis sections separately:
```bash
# Core results only
Rscript src/run_PoT.R

# Validation only  
Rscript src/run_validation.R

# Sensitivity analysis only
Rscript src/run_sens.R

# Diagnostics only
Rscript src/bulk_diagnostics.R
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

## Using Your Own Data

To analyze different pandemic datasets, modify parameters in `config.yml` or create custom configurations. See the detailed [Configuration Guide](CONFIG.md) for:

- How to specify custom data files
- Required data format and column specifications
- Adjusting analysis parameters for different datasets
- Population scaling and threshold configuration
- Output customization options

**Basic example for custom data:**
```yaml
my_analysis:
  inherit: default
  file_path: "data/my_pandemics.csv"
  pop_reference: 5e9  # Your reference population
```

Then run with:
```bash
export R_CONFIG_ACTIVE=my_analysis
Rscript src/run_all.R
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

Running `Rscript src/run_all.R` generates a comprehensive output directory with:

**Script-specific results:**
- **run_pot/**: Main analysis results with summary tables, burden estimates, severity distribution plots, and forecasts
- **run_validation/**: Out-of-sample validation results with prediction accuracy metrics and observed vs. predicted plots  
- **run_sens/**: Sensitivity analysis with threshold robustness plots and parameter stability assessments
- **bulk_diagnostics/**: Bulk distribution diagnostics with AIC/BIC comparisons and goodness-of-fit tests
- **tail_diagnostics/**: Tail threshold selection diagnostics with stability plots and recommendations

**Session documentation:**
- **execution_log.csv**: Detailed log with timing, status, and configuration for each analysis component
- **session_info.txt**: Complete R session information including package versions
- **run_summary.md**: Human-readable summary of execution results and output structure

All plots are saved in SVG format for publication quality. Tables are provided in CSV format for further analysis.

