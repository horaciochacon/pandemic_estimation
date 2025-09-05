## Centralized loader for project R scripts & packages
## Ensures deterministic, explicit sourcing order (utilities before high-level orchestration)
## Usage: source("R/load_all.R") early in run scripts.

# -----------------------------
# Package loading (deduplicated)
# -----------------------------
pkgs <- c(
  "tidyverse", "dplyr", "purrr", "stringr", "ggplot2", "ggrepel", "patchwork",
  "gridExtra", "scales", "data.table", "config", "evir", "MASS", "stats",
  "lubridate", "cowplot", "viridis"
)
for (p in pkgs) {
  suppressPackageStartupMessages({
    if (!requireNamespace(p, quietly = TRUE)) {
      message(sprintf("Package '%s' not installed; some functions may be unavailable.", p))
    } else if (!paste0("package:", p) %in% search()) {
      try(library(p, character.only = TRUE), silent = TRUE)
    }
  })
}

# -----------------------------
# Source order (low-level utils -> modeling -> analysis wrappers)
# -----------------------------
source_order <- c(
  "R/config_utils.R",
  "R/output_manager.R",
  "R/transforms.R",
  "R/parallel_utils.R",
  "R/time_rate_extremes.R",
  "R/time_trend_analysis.R",
  "R/time_forecast.R",
  "R/bootstrap.R",
  "R/model_fitting.R",
  "R/burden.R",
  "R/plotting.R",
  "R/tables.R",
  "R/yaml_export.R",
  "R/diagnostics.R",
  "R/sensitivity.R",
  "R/validation.R",
  "R/analysis.R",
  "R/mixture_plot.R"
)
for (f in source_order) {
  if (file.exists(f)) {
    source(f)
  } else {
    message(sprintf("Expected source file missing: %s", f))
  }
}

invisible(TRUE)
