####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Master orchestration script for running all pandemic estimation analyses.
##              This script coordinates the execution of all analysis scripts:
##              1. Main peaks-over-threshold mixture model analysis
##              2. Model validation with out-of-sample testing
##              3. Comprehensive sensitivity analysis
##              4. Bulk and tail distribution component diagnostics
##
## Outputs:     - Unified output directory with script-specific subdirectories
##              - Comprehensive execution summary and logs
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================
# Set active configuration
Sys.setenv(R_CONFIG_ACTIVE = "default")
conf <- load_config()

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup run_all session with clean directory structure
run_info <- setup_run_all_session(conf)
output_dir <- run_info$output_dir

# ===================================================================================
# Execute all analysis scripts
# ===================================================================================
# Run all scripts with unified output management
results <- run_all_scripts(conf, output_dir = output_dir)

# ===================================================================================
# Generate execution summary
# ===================================================================================
# Create comprehensive summary tables and logs
generate_run_all_summary(results, output_dir = output_dir, conf = conf)