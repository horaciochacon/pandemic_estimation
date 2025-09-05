####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Orchestration script for pandemic model validation.
##              This script coordinates the model validation workflow by:
##              1. Setting up the environment with validation configuration
##              2. Reading and transforming historical pandemic data
##              3. Splitting data into training and validation periods
##              4. Training the model on historical data and validating on recent data
##              5. Generating visualization comparing predicted vs observed pandemic patterns
##
## Requires:    - Configuration file ("config.yml") with validation parameters
##
## Outputs:     - Model validation results measuring prediction accuracy
##              - Visualization comparing predicted vs. observed pandemic patterns
##              - Model performance metrics
##              - Plots
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================
# Set active configuration for validation
Sys.setenv(R_CONFIG_ACTIVE = "validation")
conf <- load_config()

# ===================================================================================
# Data preparation
# ===================================================================================
# Read and transform pandemic data for validation
data <- read_transform_data(conf)

# ===================================================================================
# Model validation execution
# ===================================================================================
# Select variable for validation (first configured variable)
variable <- conf$variables[1]

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup script-specific timestamped output directory
run_info <- setup_script_run("run_validation", conf)
output_dir <- run_info$output_dir

# Run model validation process
validation_results <- validate_model(data, variable, conf, output_dir = output_dir)

# ===================================================================================
# Results generation and output
# ===================================================================================
# Generate validation summary table
generate_validation_table(validation_results, output_dir = output_dir, conf = conf)

# Display and save validation comparison plot
if (!is.null(validation_results$plot)) {
  save_plot_if_enabled(validation_results$plot, "validation_comparison_plot", output_dir, conf)
  print(validation_results$plot)
}
