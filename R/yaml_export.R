####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Functions for extracting model results and exporting to YAML format for use in
##              manuscript templates and automated reporting. This module provides utilities for:
##              1. Formatting numerical values with proper precision and comma separators
##              2. Extracting uncertainty intervals from model results
##              3. Converting model outputs to YAML-compatible structure
##              4. Writing formatted _variables.yml file for manuscript integration
##
## Functions:   - format_yaml_number(): Format numbers with commas for YAML output
##              - format_yaml_percent(): Format percentages with % sign
##              - format_yaml_coefficient(): Format coefficients with proper sign notation
##              - extract_ui(): Extract uncertainty interval values (median, low, upp)
##              - extract_variables_from_results(): Main extraction function for all variables
##              - write_variables_yaml(): Write formatted YAML file with proper structure
##
## Requires:    - yaml package for YAML output
##              - Model results objects with uncertainty quantification
##              - Validation results objects
##
## Outputs:     - _variables.yml file with all model results formatted as strings
##
####################################################################################################

#' Format number with commas for thousands, matching _variables.yml format
#'
#' @param value Numeric value to format
#' @param digits Number of decimal places (default: 0)
#' @return Character string with formatted number
format_yaml_number <- function(value, digits = 0) {
  if (is.null(value)) return("NA")
  
  # Handle vectors by taking the first element or mean
  if (length(value) > 1) {
    value <- value[1]  # Take first element for consistency
  }
  
  if (is.na(value)) return("NA")
  
  if (digits == 0) {
    # For whole numbers, format with commas
    formatted <- formatC(round(value), format = "f", digits = 0, big.mark = ",")
  } else {
    # For decimal numbers, round first then format
    formatted <- formatC(round(value, digits), format = "f", digits = digits, big.mark = ",")
  }
  
  return(formatted)
}

#' Format percentage as character with % sign
#'
#' @param value Numeric value to convert to percentage
#' @param digits Number of decimal places (default: 0)
#' @param already_percent Whether value is already in percentage form (default: FALSE)
#' @return Character string with formatted percentage
format_yaml_percent <- function(value, digits = 0, already_percent = FALSE) {
  if (is.null(value)) return("NA")
  
  # Handle vectors by taking the first element
  if (length(value) > 1) {
    value <- value[1]
  }
  
  if (is.na(value)) return("NA")
  
  # If already in percentage, don't multiply by 100
  percent_val <- if (already_percent) value else value * 100
  formatted <- formatC(round(percent_val, digits), format = "f", digits = digits)
  return(paste0(formatted, "%"))
}

#' Format coefficient with proper sign and precision
#'
#' @param value Numeric coefficient value
#' @param digits Number of decimal places (default: 3)
#' @return Character string with formatted coefficient using en-dash for negative values
format_yaml_coefficient <- function(value, digits = 3) {
  if (is.null(value)) return("NA")
  
  # Handle vectors by taking the first element
  if (length(value) > 1) {
    value <- value[1]
  }
  
  if (is.na(value)) return("NA")
  
  # Format the number with proper sign
  formatted <- formatC(value, format = "f", digits = digits)
  
  # The template shows –0.171, so negative values should be shown with en-dash
  if (value < 0) {
    # Replace minus sign with en-dash
    formatted <- gsub("^-", "–", formatted)
  }
  
  return(formatted)
}

#' Extract uncertainty interval values (median, low, upp)
#'
#' @param median_val Median value
#' @param low_val Lower bound of uncertainty interval
#' @param upp_val Upper bound of uncertainty interval
#' @param digits Number of decimal places for median (default: 0)
#' @param digits_low Number of decimal places for lower bound (default: same as digits)
#' @param digits_upp Number of decimal places for upper bound (default: same as digits)
#' @return List with formatted median, low, and upp values
extract_ui <- function(median_val, low_val, upp_val, digits = 0, digits_low = NULL, digits_upp = NULL) {
  # Allow different digits for low and upp if specified
  if (is.null(digits_low)) digits_low <- digits
  if (is.null(digits_upp)) digits_upp <- digits
  
  list(
    median = format_yaml_number(median_val, digits),
    low = format_yaml_number(low_val, digits_low),
    upp = format_yaml_number(upp_val, digits_upp)
  )
}

#' Extract all variables from model results and validation results
#'
#' @param results Main analysis results object
#' @param validation_results Validation analysis results object
#' @param conf Main analysis configuration
#' @param conf_validation Validation analysis configuration
#' @return List with all extracted and formatted variables
extract_variables_from_results <- function(results, validation_results, conf, conf_validation) {
  
  # Historical data - correct the variables
  total_events <- results$specs$n  # Total events in the dataset after base_year
  events_above_threshold <- results$specs$PoB  # Events above lower cutoff (Peaks over Bulk threshold)
  data_min <- format_yaml_number(conf$pop_min)
  period_start <- as.character(conf$base_year)
  period_end <- as.character(conf$time$current_year)
  forecast_start <- as.character(conf$time$current_year + 1)
  forecast_end <- as.character(conf$time$current_year + conf$time$future_years)
  
  # Out-of-sample validation periods
  oos_obs_end <- as.character(conf_validation$train_year_end)
  oos_pred_start <- as.character(conf_validation$test_year_start + 1)  # Should be 1951, not 1950
  oos_pred_end <- as.character(conf_validation$test_year_end)
  
  # Model parameters
  threshold <- format_yaml_number(conf$thresholds$lower_cutoff)
  tail_threshold <- format_yaml_number(results$specs$Threshold / 1e6, 0)  # Convert to millions, no decimals
  
  # Calculate median pandemic size from severity (convert to millions)
  # Template shows 4.14, 2.67, 11.3 (different precision for upp)
  median_pandemic_size <- extract_ui(
    results$results$boot / 1e6,
    results$results$boot_low / 1e6, 
    results$results$boot_upp / 1e6,
    digits = 2, digits_upp = 1
  )
  
  # Probability estimates
  annual_rate <- results$burden$yearly_rate
  
  # Handle time trend results
  if (is.list(annual_rate)) {
    # Time trend case - use the specific year rates
    rate_2025 <- annual_rate$rate_2026  # Use 2026 rate as closest to 2025
    rate_2100 <- annual_rate$rate_2100  # Use actual 2100 rate
    
    # Extract trend information
    beta_1_coeff <- annual_rate$beta_1_coefficient
    beta_1_se <- annual_rate$beta_1_se
    beta_1_pvalue <- annual_rate$beta_1_pvalue
    
    # Calculate CI
    beta_1_ci_low <- beta_1_coeff - 1.96 * beta_1_se
    beta_1_ci_upp <- beta_1_coeff + 1.96 * beta_1_se
    
  } else {
    # Constant rate case
    rate_2025 <- annual_rate
    rate_2100 <- annual_rate
    beta_1_coeff <- 0
    beta_1_ci_low <- 0
    beta_1_ci_upp <- 0
    beta_1_pvalue <- 1
  }
  
  # Calculate annual decrease percentage 
  # This is the annual percentage decrease based on the trend
  # From the trend coefficient: (1 - 0.9993594) * 100 = 0.06406%
  if (is.list(annual_rate) && !is.null(annual_rate$trend_coefficient)) {
    annual_decrease_pct <- (1 - annual_rate$trend_coefficient) * 100
  } else {
    annual_decrease_pct <- 0  # No trend for constant rate model
  }
  
  # Reference scenario (from burden calculations)
  # Annual deaths is the average: cumulative deaths divided by forecast years
  ref_annual_deaths <- extract_ui(
    (results$burden$cum_deaths / conf$time$future_years) / 1e6,  # Average annual deaths in millions
    (results$burden$cum_deaths_low / conf$time$future_years) / 1e6,
    (results$burden$cum_deaths_up / conf$time$future_years) / 1e6,
    digits = 0,  # Median with no digits (rounded to "1")
    digits_low = 3,  # Low with 3 digits ("0.239")
    digits_upp = 2   # Upper with 2 digits ("3.11")
  )
  
  ref_cumulative_deaths <- extract_ui(
    results$burden$cum_deaths / 1e6,  # Convert to millions
    results$burden$cum_deaths_low / 1e6,
    results$burden$cum_deaths_up / 1e6,
    digits = 1, digits_upp = 0  # Template shows "233" without decimal
  )
  
  # DALYs (must be calculated)
  if (is.null(results$dalys)) {
    stop("DALY calculations are required but not found in results. Check calc_dalys setting in config.")
  }
  total_dalys <- extract_ui(
    results$dalys$cum_dalys / 1e9,  # Convert to billions
    results$dalys$cum_dalys_low / 1e9,
    results$dalys$cum_dalys_up / 1e9,
    digits = 2
  )
  
  # Scenario analysis - check if available, if not calculate from burden/rate extremes
  if (!is.null(results$scenarios)) {
    # Scenarios are available directly
    best_deaths <- extract_ui(
      results$scenarios$best_case$cum_deaths / 1e6,
      results$scenarios$best_case$cum_deaths_low / 1e6,
      results$scenarios$best_case$cum_deaths_up / 1e6,
      digits = 1
    )
    
    worst_deaths <- extract_ui(
      results$scenarios$worst_case$cum_deaths / 1e6,
      results$scenarios$worst_case$cum_deaths_low / 1e6,
      results$scenarios$worst_case$cum_deaths_up / 1e6,
      digits = 0
    )
  } else {
    # Scenarios not in results - this happens when plots.scenario_comparison is FALSE
    # We need to either:
    # 1. Enable scenario_comparison plotting in config, OR
    # 2. Calculate these values some other way
    stop("Scenarios are required but not found. Set plots.scenario_comparison: true in config.yml")
  }
  
  # Burden gap calculations - only if scenarios exist
  if (!is.null(results$scenarios) && !is.null(results$scenarios$deaths_gap)) {
    burden_gap_deaths <- extract_ui(
      results$scenarios$deaths_gap$absolute$median / 1e6,
      results$scenarios$deaths_gap$absolute$lower_95 / 1e6,
      results$scenarios$deaths_gap$absolute$upper_95 / 1e6,
      digits = 0
    )
    
    if (!is.null(results$scenarios$dalys_gap)) {
      burden_gap_dalys <- extract_ui(
        results$scenarios$dalys_gap$absolute$median / 1e9,
        results$scenarios$dalys_gap$absolute$lower_95 / 1e9,
        results$scenarios$dalys_gap$absolute$upper_95 / 1e9,
        digits = 2
      )
    } else {
      # Calculate gap from worst - best if not provided
      stop("DALY gap calculations required but not found. Set plots.scenario_comparison: true in config.yml")
    }
    
    relative_increase <- extract_ui(
      results$scenarios$deaths_gap$relative$median,
      results$scenarios$deaths_gap$relative$lower_95,
      results$scenarios$deaths_gap$relative$upper_95,
      digits = 1
    )
    
    # Convert to percentages
    relative_increase$median <- paste0(relative_increase$median, "%")
    relative_increase$low <- paste0(relative_increase$low, "%")
    relative_increase$upp <- paste0(relative_increase$upp, "%")
  } else {
    stop("Burden gap calculations required but not found. Set plots.scenario_comparison: true in config.yml")
  }
  
  # DALY ratio (already checked above that results$dalys exists)
  daly_ratio <- extract_ui(
    results$dalys$daly_ratio,
    results$dalys$daly_ratio_low,
    results$dalys$daly_ratio_up,
    digits = 1
  )
  
  # Validation results (must exist)
  if (is.null(validation_results$predicted) || is.null(validation_results$predicted$result)) {
    stop("Validation predicted results are required but not found.")
  }
  
  # Get cumulative predicted deaths for validation period
  pred_cum <- validation_results$predicted$result$cum_deaths
  pred_low <- validation_results$predicted$result$cum_deaths_low
  pred_upp <- validation_results$predicted$result$cum_deaths_up
  
  validation_predicted <- extract_ui(pred_cum / 1e6, pred_low / 1e6, pred_upp / 1e6, digits = 1, digits_upp = 0)
  
  if (is.null(validation_results$observed)) {
    stop("Validation observed results are required but not found.")
  }
  
  # Get the maximum cumulative observation
  total_obs <- max(validation_results$observed$obs, na.rm = TRUE)
  validation_observed <- format_yaml_number(total_obs / 1e6, 1)
  
  # Return structured list of all variables
  list(
    historical = list(
      total_events = as.character(total_events),
      events_above_threshold = as.character(events_above_threshold),
      data_min = data_min,
      period_start = period_start,
      period_end = period_end,
      forecast_start = forecast_start,
      forecast_end = forecast_end,
      oos_obs_end = oos_obs_end,
      oos_pred_start = oos_pred_start,
      oos_pred_end = oos_pred_end
    ),
    
    model = list(
      threshold = threshold,
      tail_threshold = tail_threshold,
      forecast_years = as.character(conf$time$future_years),
      median_pandemic_size = median_pandemic_size
    ),
    
    probability = list(
      annual_percent = format_yaml_percent(rate_2025),
      at_2025 = format_yaml_number(rate_2025, 3),
      at_2100 = format_yaml_number(rate_2100, 3),
      mean_2026_2100 = if (is.list(annual_rate)) {
        extract_ui(
          annual_rate$mean,
          annual_rate$lower_95,
          annual_rate$upper_95,
          3
        )
      } else {
        # For constant rate, use the same rate with no uncertainty
        extract_ui(annual_rate, annual_rate, annual_rate, 3)
      },
      trend = list(
        coefficient = format_yaml_coefficient(beta_1_coeff, 3),
        ci_low = format_yaml_coefficient(beta_1_ci_low, 3),
        ci_upp = format_yaml_coefficient(beta_1_ci_upp, 3),
        p_value = format_yaml_number(beta_1_pvalue, 2),
        annual_decrease = format_yaml_percent(annual_decrease_pct, 2, already_percent = TRUE)
      )
    ),
    
    reference = list(
      annual_deaths = ref_annual_deaths,
      cumulative_deaths = ref_cumulative_deaths,
      total_dalys = total_dalys
    ),
    
    scenarios = list(
      best_plausible = list(deaths = best_deaths),
      worst_plausible = list(deaths = worst_deaths),
      burden_gap = list(
        deaths = burden_gap_deaths,
        dalys = burden_gap_dalys,
        relative_increase = relative_increase
      )
    ),
    
    daly_ratio = daly_ratio,
    
    validation = list(
      period = paste0(oos_pred_start, "–", oos_pred_end),
      predicted = validation_predicted,
      observed = validation_observed
    ),
    
    summary = list(
      main_deaths = ref_cumulative_deaths,
      main_dalys = total_dalys,
      scenario_range = paste(best_deaths$median, "to", worst_deaths$median),
      burden_gap_dalys = format_yaml_number(as.numeric(gsub(",", "", burden_gap_dalys$median)), 0)
    ),
    
    periods = list(
      projection = paste(forecast_start, "and", forecast_end),
      historical = paste0(period_start, "–", period_end)
    ),
    
    population = list(
      reference_year = period_end,
      reference_size = format_yaml_number(conf$pop_reference / 1e9, 1)
    ),
    
    statistics = list(
      tail_threshold = tail_threshold,
      poisson_window = as.character(conf$window$default_size),
      bootstrap_draws = format_yaml_number(conf$draws)
    )
  )
}

#' Write YAML file with exact formatting, comments, and spacing
#'
#' @param yaml_data List of extracted variables
#' @param output_file Path to output YAML file (default: "output/_variables.yml")
write_variables_yaml <- function(yaml_data, output_file = "output/_variables.yml") {
  
  # Create output directory if it doesn't exist
  dir.create("output", showWarnings = FALSE)
  
  # Create the YAML content as a string with exact formatting
  yaml_content <- paste0(
'# Historical data
historical:
  total_events: "', yaml_data$historical$total_events, '"
  events_above_threshold: "', yaml_data$historical$events_above_threshold, '"
  data_min: "', yaml_data$historical$data_min, '"
  period_start: "', yaml_data$historical$period_start, '"
  period_end: "', yaml_data$historical$period_end, '"
  forecast_start: "', yaml_data$historical$forecast_start, '"
  forecast_end: "', yaml_data$historical$forecast_end, '"
  oos_obs_end: "', yaml_data$historical$oos_obs_end, '"
  oos_pred_start: "', yaml_data$historical$oos_pred_start, '"
  oos_pred_end: "', yaml_data$historical$oos_pred_end, '"

# Model parameters
model:
  threshold: "', yaml_data$model$threshold, '"
  tail_threshold: "', yaml_data$model$tail_threshold, '"
  forecast_years: "', yaml_data$model$forecast_years, '"
  median_pandemic_size:
    median: "', yaml_data$model$median_pandemic_size$median, '"
    low: "', yaml_data$model$median_pandemic_size$low, '"
    upp: "', yaml_data$model$median_pandemic_size$upp, '"

# Probability estimates
probability:
  annual_percent: "', yaml_data$probability$annual_percent, '"
  at_2025: "', yaml_data$probability$at_2025, '"
  at_2100: "', yaml_data$probability$at_2100, '"
  mean_2026_2100:
    median: "', yaml_data$probability$mean_2026_2100$median, '"
    low: "', yaml_data$probability$mean_2026_2100$low, '"
    upp: "', yaml_data$probability$mean_2026_2100$upp, '"
  trend:
    coefficient: ', yaml_data$probability$trend$coefficient, '
    ci_low: ', yaml_data$probability$trend$ci_low, '
    ci_upp: ', yaml_data$probability$trend$ci_upp, '
    p_value: "', yaml_data$probability$trend$p_value, '"
    annual_decrease: "', yaml_data$probability$trend$annual_decrease, '"

# Reference scenario
reference:
  annual_deaths:
    median: "', yaml_data$reference$annual_deaths$median, '"
    low: "', yaml_data$reference$annual_deaths$low, '"
    upp: "', yaml_data$reference$annual_deaths$upp, '"
  cumulative_deaths:
    median: "', yaml_data$reference$cumulative_deaths$median, '"
    low: "', yaml_data$reference$cumulative_deaths$low, '"
    upp: "', yaml_data$reference$cumulative_deaths$upp, '"
  total_dalys:
    median: "', yaml_data$reference$total_dalys$median, '"
    low: "', yaml_data$reference$total_dalys$low, '"
    upp: "', yaml_data$reference$total_dalys$upp, '"

# Scenario analysis
scenarios:
  best_plausible:
    deaths:
      median: "', yaml_data$scenarios$best_plausible$deaths$median, '"
      low: "', yaml_data$scenarios$best_plausible$deaths$low, '"
      upp: "', yaml_data$scenarios$best_plausible$deaths$upp, '"
  worst_plausible:
    deaths:
      median: "', yaml_data$scenarios$worst_plausible$deaths$median, '"
      low: "', yaml_data$scenarios$worst_plausible$deaths$low, '"
      upp: "', yaml_data$scenarios$worst_plausible$deaths$upp, '"
  burden_gap:
    deaths:
      median: "', yaml_data$scenarios$burden_gap$deaths$median, '"
      low: "', yaml_data$scenarios$burden_gap$deaths$low, '"
      upp: "', yaml_data$scenarios$burden_gap$deaths$upp, '"
    dalys:
      median: "', yaml_data$scenarios$burden_gap$dalys$median, '"
      low: "', yaml_data$scenarios$burden_gap$dalys$low, '"
      upp: "', yaml_data$scenarios$burden_gap$dalys$upp, '"
    relative_increase:
      median: "', yaml_data$scenarios$burden_gap$relative_increase$median, '"
      low: "', yaml_data$scenarios$burden_gap$relative_increase$low, '"
      upp: "', yaml_data$scenarios$burden_gap$relative_increase$upp, '"

# DALY conversion
daly_ratio:
  median: "', yaml_data$daly_ratio$median, '"
  low: "', yaml_data$daly_ratio$low, '"
  upp: "', yaml_data$daly_ratio$upp, '"

# Validation
validation:
  period: "', yaml_data$validation$period, '"
  predicted:
    median: "', yaml_data$validation$predicted$median, '"
    low: "', yaml_data$validation$predicted$low, '"
    upp: "', yaml_data$validation$predicted$upp, '"
  observed: "', yaml_data$validation$observed, '"

# Summary
summary:
  main_deaths:
    median: "', yaml_data$summary$main_deaths$median, '"
    low: "', yaml_data$summary$main_deaths$low, '"
    upp: "', yaml_data$summary$main_deaths$upp, '"
  main_dalys:
    median: "', yaml_data$summary$main_dalys$median, '"
    low: "', yaml_data$summary$main_dalys$low, '"
    upp: "', yaml_data$summary$main_dalys$upp, '"
  scenario_range: "', yaml_data$summary$scenario_range, '"
  burden_gap_dalys: "', yaml_data$summary$burden_gap_dalys, '"

# Time periods
periods:
  projection: "', yaml_data$periods$projection, '"
  historical: "', yaml_data$periods$historical, '"

# Population
population:
  reference_year: "', yaml_data$population$reference_year, '"
  reference_size: "', yaml_data$population$reference_size, '"

# Statistical methods
statistics:
  tail_threshold: "', yaml_data$statistics$tail_threshold, '"
  poisson_window: "', yaml_data$statistics$poisson_window, '"
  bootstrap_draws: "', yaml_data$statistics$bootstrap_draws, '"
')
  
  # Write the formatted YAML content to file
  writeLines(yaml_content, output_file)
  
  cat("Successfully generated", output_file, "\n")
  cat("Variables extracted from model run and formatted for YAML output.\n")
}