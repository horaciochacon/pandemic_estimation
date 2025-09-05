####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Orchestration script for pandemic model sensitivity analysis.
##              This script coordinates the entire sensitivity analysis workflow by:
##              1. Setting up the environment with sensitivity configuration
##              2. Running sensitivity analyses across multiple dimensions:
##                 a. Time period sensitivity - how results change with data inclusion years
##                 b. Threshold sensitivity - impact of varying the tail threshold
##                 c. Cutoff sensitivity - effect of different lower cutoff values
##              3. Generating visualizations for each sensitivity dimension
##
## Requires:    - Configuration file ("config.yml") with sensitivity parameters
##
## Outputs:     - Sensitivity analysis results for each dimension
##              - Visualizations showing parameter impact on:
##                  * Expected annual burden estimates
##                  * Model parameter stability
##                  * Cumulative death projections
##              - Plots
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================

# Set active configuration for sensitivity analysis
Sys.setenv(R_CONFIG_ACTIVE = "sensitivity")
conf <- load_config()

# Read and transform the pandemic data
data <- read_transform_data(conf)

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup script-specific timestamped output directory
run_info <- setup_script_run("run_sens", conf)
output_dir <- run_info$output_dir

# Ensure required packages
for (pkg in c("future", "future.apply", "foreach", "doParallel", "parallel", "viridis")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Package", pkg, "is not available."))
  }
}

# Set parallelization strategy from config
parallel_strategy <- ifelse(
  isTRUE(conf$parallel$enabled),
  conf$parallel$strategy,
  "none"
)
workers <- conf$parallel$workers # NULL means use detectCores() - 1

# Log parallelization settings
cat(
  "Parallelization:",
  ifelse(parallel_strategy == "none", "Disabled", "Enabled"),
  ifelse(parallel_strategy != "none", paste0("using '", parallel_strategy, "' strategy"), ""),
  "\n"
)

# Performance timing - start time
start_time <- Sys.time()
cat("Starting sensitivity analysis at:", format(start_time), "\n")

# ===================================================================================
# Time period sensitivity analysis
# ===================================================================================
cat("Running time period sensitivity analysis...\n")

# Generate year range for sensitivity testing
years <- get_year_range(data, conf)
thresholds <- unlist(conf$sens$base_thresholds)

# Configure analysis parameters
sensitivity_grid <- expand.grid(
  analysis_type = "year",
  plot_type = c("cumulative"),
  var = conf$variables,
  stringsAsFactors = FALSE
)

# Execute year-based sensitivity analysis
results_year <- sensitivity_grid |>
  pmap(
    function(var, analysis_type, plot_type) {
      result <- run_parallel_sensitivity(
        data = data,
        var = var,
        analysis_type = analysis_type,
        values = years,
        plot_type = plot_type,
        conf = conf,
        strategy = parallel_strategy,
        workers = workers,
        output_dir = output_dir
      )

      # Return results for further processing
      result
    }
  )

# ===================================================================================
# Threshold sensitivity analysis
# ===================================================================================
cat("Running threshold sensitivity analysis...\n")

# Extract threshold ranges for testing
thresholds <- list(
  deaths = threshold_extract(data, "deaths", conf$base_year, conf),
  deaths_scaled = threshold_extract(data, "deaths_scaled", conf$base_year, conf)
) |>
  map(~ .x[65:90])

# Configure analysis parameters
sensitivity_grid <- expand.grid(
  analysis_type = "threshold",
  plot_type = c("cumulative"),
  var = conf$variables,
  stringsAsFactors = FALSE
)

# Execute threshold-based sensitivity analysis
results_threshold <- sensitivity_grid |>
  pmap(
    function(var, analysis_type, plot_type) {
      result <- run_parallel_sensitivity(
        data = data,
        var = var,
        analysis_type = analysis_type,
        values = thresholds[[var]],
        plot_type = plot_type,
        conf = conf,
        strategy = parallel_strategy,
        workers = workers,
        output_dir = output_dir
      )

      # Return results for further processing
      result
    }
  )

# ===================================================================================
# Lower cutoff sensitivity analysis
# ===================================================================================
cat("Running lower cutoff sensitivity analysis...\n")

# Extract cutoff ranges for testing
cutoffs <- list(
  deaths = threshold_extract(data, "deaths", conf$base_year, conf),
  deaths_scaled = threshold_extract(data, "deaths_scaled", conf$base_year, conf)
) |>
  map(~ .x[5:75])

# Configure analysis parameters
sensitivity_grid <- expand.grid(
  analysis_type = "cutoff",
  plot_type = c("cumulative"),
  var = conf$variables,
  stringsAsFactors = FALSE
)

# Execute cutoff-based sensitivity analysis
results_cutoff <- sensitivity_grid |>
  pmap(
    function(var, analysis_type, plot_type) {
      result <- run_parallel_sensitivity(
        data = data,
        var = var,
        analysis_type = analysis_type,
        values = cutoffs[[var]],
        plot_type = plot_type,
        conf = conf,
        strategy = parallel_strategy,
        workers = workers,
        output_dir = output_dir
      )

      # Return results for further processing
      result
    }
  )
