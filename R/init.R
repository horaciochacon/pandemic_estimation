# Initialization helper for pandemic_estimation project
# Source this early (e.g., from run_PoT.R) to load config and propagate options.
# Usage:
#   conf <- initialize_analysis(environment = Sys.getenv("R_CONFIG_ACTIVE", "default"))
#   # pass conf to downstream functions

if (!exists("initialize_analysis")) {
  source("R/config_utils.R")
}

conf <- initialize_analysis()
