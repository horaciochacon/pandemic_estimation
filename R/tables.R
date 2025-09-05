#' ===================================================================
#' LANCET FORMATTING HELPER FUNCTIONS
#' ===================================================================
#' These functions ensure compliance with The Lancet formatting requirements:
#' - Counts: 3 significant figures
#' - Rates/ratios/percentages: 1 decimal place
#' - Numbers ≥5 digits: use spaces (12 345)
#' - Numbers <5 digits: no spaces/commas (1234)
#' - Uncertainty intervals: (95% UI XXX–XXX) with en-dash
#' - First UI in paragraph includes "95% UI"
#' - Subsequent UIs just (XXX–XXX)

#' Format numbers according to The Lancet guidelines
#'
#' @param value Numeric value(s) to format - can be single value or vector
#' @param digits Number of decimal places (default: 0)
#' @param sig_figs Number of significant figures for counts (default: 3)
#' @param is_count Whether this is a count (uses sig figs) vs rate (uses decimal places)
#' @return Formatted string(s) following Lancet guidelines
format_lancet_number <- function(value, digits = 0, sig_figs = 3, is_count = FALSE) {
  # Handle NULL or all NA values
  if (is.null(value) || all(is.na(value))) {
    return("NA")
  }

  # Vectorized function for formatting individual values
  format_single <- function(val) {
    if (is.na(val)) {
      return("NA")
    }

    if (is_count) {
      # For counts: use significant figures
      formatted <- signif(val, sig_figs)
      # Convert to string and format with spaces if ≥5 digits
      str_val <- formatC(formatted, format = "f", digits = 0, big.mark = "")
    } else {
      # For rates/ratios: use decimal places
      # For very small rates, use more decimal places to avoid rounding to zero
      adjusted_digits <- digits
      if (val > 0 && val < 0.1 && digits == 1) {
        adjusted_digits <- 3 # Use 3 decimal places for rates < 0.1
      }
      str_val <- formatC(val, format = "f", digits = adjusted_digits, big.mark = "")
    }

    # Apply Lancet spacing rule: spaces for numbers ≥5 digits
    if (nchar(gsub("[^0-9]", "", str_val)) >= 5) {
      # Add spaces for thousands separator
      final_digits <- ifelse(is_count, 0, adjusted_digits)
      str_val <- formatC(as.numeric(str_val), format = "f", digits = final_digits, big.mark = " ")
    }

    return(str_val)
  }

  # Apply formatting to each value in the vector
  sapply(value, format_single, USE.NAMES = FALSE)
}

#' Format uncertainty intervals according to The Lancet guidelines
#'
#' @param mean_val Mean/median value(s) - can be single value or vector
#' @param lower_val Lower bound(s) of uncertainty interval
#' @param upper_val Upper bound(s) of uncertainty interval
#' @param digits Number of decimal places (default: 0)
#' @param sig_figs Number of significant figures for counts (default: 3)
#' @param is_count Whether this is a count vs rate/ratio
#' @param include_ui_label Whether to include "95% UI" (TRUE for first in paragraph)
#' @param scale_suffix Suffix to add (e.g., "M", "B")
#' @return Formatted uncertainty interval string(s)
format_lancet_uncertainty <- function(mean_val, lower_val, upper_val,
                                      digits = 0, sig_figs = 3, is_count = FALSE,
                                      include_ui_label = TRUE, scale_suffix = "") {
  # Check for any NA values across all vectors
  if (is.null(mean_val) || is.null(lower_val) || is.null(upper_val)) {
    return("NA")
  }

  # Handle single NA check - use any() for vectors
  has_na <- any(is.na(mean_val)) || any(is.na(lower_val)) || any(is.na(upper_val))
  if (has_na && length(mean_val) == 1) {
    return("NA")
  }

  # Vectorized function for formatting individual uncertainty intervals
  format_single_ui <- function(m_val, l_val, u_val, include_label = include_ui_label) {
    if (is.na(m_val) || is.na(l_val) || is.na(u_val)) {
      return("NA")
    }

    # Format individual values
    mean_str <- format_lancet_number(m_val, digits, sig_figs, is_count)
    lower_str <- format_lancet_number(l_val, digits, sig_figs, is_count)
    upper_str <- format_lancet_number(u_val, digits, sig_figs, is_count)

    # Add suffix if provided
    if (scale_suffix != "") {
      mean_str <- paste0(mean_str, scale_suffix)
      lower_str <- paste0(lower_str, scale_suffix)
      upper_str <- paste0(upper_str, scale_suffix)
    }

    # Create uncertainty interval with en-dash
    if (include_label) {
      ui_text <- paste0("(95% UI ", lower_str, "–", upper_str, ")")
    } else {
      ui_text <- paste0("(", lower_str, "–", upper_str, ")")
    }

    return(paste0(mean_str, " ", ui_text))
  }

  # Handle vectorized input
  if (length(mean_val) > 1 || length(lower_val) > 1 || length(upper_val) > 1) {
    # For vectors, only include UI label on first element
    include_labels <- c(include_ui_label, rep(FALSE, max(length(mean_val), length(lower_val), length(upper_val)) - 1))
    return(mapply(format_single_ui, mean_val, lower_val, upper_val, include_labels, USE.NAMES = FALSE))
  } else {
    # For single values
    return(format_single_ui(mean_val, lower_val, upper_val, include_ui_label))
  }
}

#' Format numbers with scale suffixes (M, B) according to Lancet guidelines
#'
#' @param value Numeric value
#' @param lower_val Lower uncertainty bound
#' @param upper_val Upper uncertainty bound
#' @param is_count Whether this is a count vs rate
#' @param include_ui_label Whether to include "95% UI" label
#' @return Formatted string with appropriate scale
format_lancet_scaled <- function(value, lower_val, upper_val, is_count = TRUE, include_ui_label = TRUE) {
  if (any(is.na(c(value, lower_val, upper_val)))) {
    return("NA")
  }

  # Determine scale based on magnitude
  if (value >= 1e9) {
    scale_factor <- 1e9
    suffix <- "B"
    digits <- 2 # 2 decimals for billions per Lancet
  } else if (value >= 1e6) {
    scale_factor <- 1e6
    suffix <- "M"
    digits <- 1 # 1 decimal for millions per Lancet
  } else {
    scale_factor <- 1
    suffix <- ""
    digits <- ifelse(is_count, 0, 1)
  }

  # Scale values
  scaled_value <- value / scale_factor
  scaled_lower <- lower_val / scale_factor
  scaled_upper <- upper_val / scale_factor

  return(format_lancet_uncertainty(
    scaled_value, scaled_lower, scaled_upper,
    digits = digits, is_count = FALSE,
    include_ui_label = include_ui_label, scale_suffix = suffix
  ))
}

# Internal utility (Phase 1 cleanup): generic scale + suffix decision returning vector of list(mean,low,up)
.format_scale_triplet <- function(value, low = NULL, upp = NULL) {
  if (is.null(value)) {
    return(list(value = NA, low = NA, upp = NA, suffix = ""))
  }
  if (is.list(value)) value <- value[[1]]
  if (!is.null(low) && is.list(low)) low <- low[[1]]
  if (!is.null(upp) && is.list(upp)) upp <- upp[[1]]
  if (length(value) > 1) value <- value[1]
  if (!is.null(low) && length(low) > 1) low <- low[1]
  if (!is.null(upp) && length(upp) > 1) upp <- upp[1]
  # Determine scale
  if (is.na(value)) {
    sf <- 1
    suf <- ""
  } else if (value >= 1e9) {
    sf <- 1e9
    suf <- "B"
  } else if (value >= 1e6) {
    sf <- 1e6
    suf <- "M"
  } else {
    sf <- 1
    suf <- ""
  }
  list(
    value = value / sf,
    low = if (!is.null(low) && !is.na(low)) low / sf else low,
    upp = if (!is.null(upp) && !is.na(upp)) upp / sf else upp,
    suffix = suf
  )
}

# Internal simplified wrapper for uncertainty formatting with automatic scaling
.format_scaled_uncertainty <- function(value, low = NULL, upp = NULL, digits_count = 0, include_ui_label = FALSE, is_count = TRUE) {
  trip <- .format_scale_triplet(value, low, upp)
  if (is.null(low) || is.null(upp)) {
    # point value only
    val_str <- format_lancet_number(trip$value, digits = ifelse(is_count && trip$suffix == "", 0, digits_count), is_count = is_count)
    return(paste0(val_str, trip$suffix))
  }
  # Determine digits if not provided explicitly (keep existing semantics)
  digits_use <- ifelse(trip$suffix == "B", 2, ifelse(trip$suffix == "M", 1, ifelse(is_count, 0, digits_count)))
  format_lancet_uncertainty(
    trip$value, trip$low, trip$upp,
    digits = digits_use, is_count = FALSE,
    include_ui_label = include_ui_label, scale_suffix = trip$suffix
  )
}

#' Generate a summary table of pandemic model projections
#'
#' Creates a markdown table summarizing key pandemic model projections including
#' expected pandemic size, annual deaths, and cumulative burden
#'
#' @param results Results object from run_PoT containing model outputs
#' @param period_years Number of years in the projection period (default: 75, for 2025-2100)
#' @param save_to_file Optional path to save the table to a markdown file
#' @param digits Number of decimal places for scaled values (default: 1)
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_summary_table <- function(results, period_years = NULL, save_to_file = NULL, digits = 1, conf = NULL, output_dir = NULL) {
  if (is.null(period_years)) {
    period_years <- if (!is.null(conf) && !is.null(conf$tables$projection_years)) conf$tables$projection_years else getOption("pandemic.tables.projection_years", 75)
  }
  # Validate input
  if (is.null(results$results) || is.null(results$burden)) {
    stop("Results object must contain model results and burden calculations")
  }

  # Extract required values
  single_event_mean <- results$results$boot
  single_event_low <- results$results$boot_low
  single_event_upp <- results$results$boot_upp

  yearly_deaths_mean <- results$burden$yearly_deaths
  yearly_deaths_low <- results$burden$yearly_deaths_low
  yearly_deaths_up <- results$burden$yearly_deaths_up

  cum_deaths_mean <- results$burden$cum_deaths
  cum_deaths_low <- results$burden$cum_deaths_low
  cum_deaths_up <- results$burden$cum_deaths_up

  # DALY calculations - first check if they exist
  if (!is.null(results$dalys)) {
    cum_dalys_mean <- results$dalys$cum_dalys
    cum_dalys_low <- results$dalys$cum_dalys_low
    cum_dalys_up <- results$dalys$cum_dalys_up
    has_dalys <- TRUE
  } else {
    has_dalys <- FALSE
  }

  # Format function using Lancet-compliant formatting with custom digits
  format_with_scale <- function(value, low, upp, include_ui_label = TRUE) {
    # Ensure we're working with single values, not vectors or lists
    if (is.list(value)) value <- value[[1]]
    if (is.list(low)) low <- low[[1]]
    if (is.list(upp)) upp <- upp[[1]]
    if (length(value) > 1) value <- value[1]
    if (length(low) > 1) low <- low[1]
    if (length(upp) > 1) upp <- upp[1]

    # Custom scaling with user-specified digits
    if (any(is.na(c(value, low, upp)))) {
      return("NA")
    }

    # Determine scale based on magnitude
    if (value >= 1e9) {
      scale_factor <- 1e9
      suffix <- "B"
    } else if (value >= 1e6) {
      scale_factor <- 1e6
      suffix <- "M"
    } else {
      scale_factor <- 1
      suffix <- ""
    }

    # Scale values
    scaled_value <- value / scale_factor
    scaled_lower <- low / scale_factor
    scaled_upper <- upp / scale_factor

    return(format_lancet_uncertainty(
      scaled_value, scaled_lower, scaled_upper,
      digits = digits, is_count = FALSE,
      include_ui_label = include_ui_label, scale_suffix = suffix
    ))
  }

  # Format for raw numbers with suffix using Lancet formatting with custom digits
  format_number_suffix <- function(value, low, upp, include_ui_label = FALSE) {
    .format_scaled_uncertainty(value, low, upp, digits_count = digits, include_ui_label = include_ui_label, is_count = TRUE)
  }

  # Prepare table rows
  rows <- c()

  # Single pandemic size (first UI includes label)
  rows <- c(rows, paste0(
    "|**Expected Size of a Single Pandemic**|",
    format_with_scale(single_event_mean, single_event_low, single_event_upp, include_ui_label = TRUE),
    "|"
  ))

  # Annual expected deaths (subsequent UI without label)
  rows <- c(rows, paste0(
    "|**Annual Expected Deaths**|",
    format_number_suffix(yearly_deaths_mean, yearly_deaths_low, yearly_deaths_up, include_ui_label = FALSE),
    "|"
  ))

  # Calculate average annual deaths over forecast period (2026-2100)
  avg_annual_deaths_mean <- cum_deaths_mean / period_years
  avg_annual_deaths_low <- cum_deaths_low / period_years
  avg_annual_deaths_up <- cum_deaths_up / period_years

  # Average annual deaths over forecast period (subsequent UI without label)
  rows <- c(rows, paste0(
    "|**Average Annual Deaths (2026-2100)**|",
    format_number_suffix(avg_annual_deaths_mean, avg_annual_deaths_low, avg_annual_deaths_up, include_ui_label = FALSE),
    "|"
  ))

  # Cumulative deaths (subsequent UI without label)
  rows <- c(rows, paste0(
    "|**Cumulative Deaths (2026-2100)**|",
    format_with_scale(cum_deaths_mean, cum_deaths_low, cum_deaths_up, include_ui_label = FALSE),
    "|"
  ))

  # Cumulative DALYs (if available) (subsequent UI without label)
  if (has_dalys) {
    rows <- c(rows, paste0(
      "|**Cumulative DALYs (2026-2100)**|",
      format_with_scale(cum_dalys_mean, cum_dalys_low, cum_dalys_up, include_ui_label = FALSE),
      "|"
    ))
  }

  # Create header with time period
  end_year <- 2100
  header <- c(
    paste0("# Summary of Pandemic Model Projections (2026-", end_year, ")"),
    "",
    "|Metric|Estimate (95% UI)|",
    "|---|---|"
  )

  # Combine header and rows into a single string with line breaks
  md_table <- paste(c(header, rows), collapse = "\n")

  # Print the table to console for easy copying
  cat(md_table, "\n")

  # Save to file if requested (legacy parameter)
  if (!is.null(save_to_file)) {
    dir_path <- dirname(save_to_file)
    if (!dir.exists(dir_path) && dir_path != ".") {
      dir.create(dir_path, recursive = TRUE)
    }
    writeLines(md_table, save_to_file)
    message("Summary table saved to: ", save_to_file)
  }
  
  # Save using output manager if enabled
  save_table_if_enabled(md_table, "summary_table.md", output_dir, conf)
}

#' Generate a scenario comparison table for pandemic model results
#'
#' Creates a markdown table comparing baseline, best case, and worst case scenarios
#' for pandemic projections over the specified time period
#'
#' @param results Results object from run_PoT containing scenario outputs
#' @param period_years Number of years in the projection period (default: 75, for 2025-2100)
#' @param save_to_file Optional path to save the table to a markdown file
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_scenario_table <- function(results, period_years = NULL, save_to_file = NULL, conf = NULL, output_dir = NULL) {
  if (is.null(period_years)) {
    period_years <- if (!is.null(conf) && !is.null(conf$tables$projection_years)) conf$tables$projection_years else getOption("pandemic.tables.projection_years", 75)
  }
  # Validate input
  if (is.null(results$burden) || is.null(results$scenarios)) {
    stop("Results object must contain burden calculations and scenario results")
  }

  # Extract scenario data
  baseline <- results$burden
  best_case <- results$scenarios$best_case
  worst_case <- results$scenarios$worst_case

  # Format function for values with suffix using Lancet formatting
  format_number_suffix <- function(value, low = NULL, upp = NULL, include_ui_label = FALSE) {
    .format_scaled_uncertainty(value, low, upp, digits_count = 0, include_ui_label = include_ui_label, is_count = TRUE)
  }

  # Create the rows - starting with annual deaths
  rows <- c()

  # Annual exceedance rate
  # Handle both simple numeric and list yearly_rate formats
  get_rate_value <- function(rate) {
    if (is.list(rate)) {
      # If it's a list (time trend case), get the mean value
      return(rate$mean)
    } else {
      # If it's numeric (constant rate case), use directly
      return(rate)
    }
  }

  # Helper function to extract rate uncertainty from different data structures
  get_rate_uncertainty <- function(scenario_data, rate_extremes = NULL, scenario_type = "baseline") {
    if (scenario_type == "baseline") {
      # For baseline, check if uncertainty is available in yearly_rate structure
      rate_data <- scenario_data$yearly_rate
      if (is.list(rate_data) && !is.null(rate_data$lower_95) && !is.null(rate_data$upper_95)) {
        # Time trend case with uncertainty
        return(list(
          mean = rate_data$mean,
          lower_95 = rate_data$lower_95,
          upper_95 = rate_data$upper_95
        ))
      } else {
        # No uncertainty available, return just point estimate
        return(list(
          mean = get_rate_value(rate_data),
          lower_95 = NA,
          upper_95 = NA
        ))
      }
    } else {
      # For best/worst case scenarios, use rate extremes draws
      if (!is.null(rate_extremes)) {
        # Handle both rate_extremes$draws structure and direct draws structure
        draws_data <- if (!is.null(rate_extremes$draws)) rate_extremes$draws else rate_extremes

        if (scenario_type == "best_case" && !is.null(draws_data$lowest_year_draws)) {
          draws <- draws_data$lowest_year_draws
          return(list(
            mean = mean(draws),
            lower_95 = quantile(draws, 0.025),
            upper_95 = quantile(draws, 0.975)
          ))
        } else if (scenario_type == "worst_case" && !is.null(draws_data$highest_year_draws)) {
          draws <- draws_data$highest_year_draws
          return(list(
            mean = mean(draws),
            lower_95 = quantile(draws, 0.025),
            upper_95 = quantile(draws, 0.975)
          ))
        }
      }
      # Fallback to point estimate if draws not available
      return(list(
        mean = get_rate_value(scenario_data$yearly_rate),
        lower_95 = NA,
        upper_95 = NA
      ))
    }
  }

  # Extract rate uncertainty for all scenarios
  rate_extremes_data <- if (!is.null(results$scenarios$rate_extremes_draws)) {
    results$scenarios$rate_extremes_draws
  } else {
    NULL
  }

  baseline_rate <- get_rate_uncertainty(baseline, rate_extremes_data, "baseline")
  best_case_rate <- get_rate_uncertainty(best_case, rate_extremes_data, "best_case")
  worst_case_rate <- get_rate_uncertainty(worst_case, rate_extremes_data, "worst_case")

  # Format rate values with uncertainty intervals
  baseline_rate_formatted <- if (!is.na(baseline_rate$lower_95) && !is.na(baseline_rate$upper_95)) {
    format_lancet_uncertainty(baseline_rate$mean, baseline_rate$lower_95, baseline_rate$upper_95,
      digits = 3, is_count = FALSE, include_ui_label = TRUE
    )
  } else {
    format_lancet_number(baseline_rate$mean, digits = 3, is_count = FALSE)
  }

  best_case_rate_formatted <- if (!is.na(best_case_rate$lower_95) && !is.na(best_case_rate$upper_95)) {
    format_lancet_uncertainty(best_case_rate$mean, best_case_rate$lower_95, best_case_rate$upper_95,
      digits = 3, is_count = FALSE, include_ui_label = FALSE
    )
  } else {
    format_lancet_number(best_case_rate$mean, digits = 3, is_count = FALSE)
  }

  worst_case_rate_formatted <- if (!is.na(worst_case_rate$lower_95) && !is.na(worst_case_rate$upper_95)) {
    format_lancet_uncertainty(worst_case_rate$mean, worst_case_rate$lower_95, worst_case_rate$upper_95,
      digits = 3, is_count = FALSE, include_ui_label = FALSE
    )
  } else {
    format_lancet_number(worst_case_rate$mean, digits = 3, is_count = FALSE)
  }

  rows <- c(rows, paste0(
    "|**Annual Exceedance Rate**|",
    baseline_rate_formatted, "|",
    best_case_rate_formatted, "|",
    worst_case_rate_formatted, "|"
  ))

  # Annual expected deaths (first UI includes label)
  rows <- c(rows, paste0(
    "|**Annual Deaths**|",
    format_number_suffix(baseline$yearly_deaths, baseline$yearly_deaths_low, baseline$yearly_deaths_up, include_ui_label = TRUE), "|",
    format_number_suffix(best_case$yearly_deaths, best_case$yearly_deaths_low, best_case$yearly_deaths_up, include_ui_label = FALSE), "|",
    format_number_suffix(worst_case$yearly_deaths, worst_case$yearly_deaths_low, worst_case$yearly_deaths_up, include_ui_label = FALSE), "|"
  ))

  # Cumulative deaths (subsequent UIs without label)
  rows <- c(rows, paste0(
    "|**Cumulative Deaths**|",
    format_number_suffix(baseline$cum_deaths, baseline$cum_deaths_low, baseline$cum_deaths_up, include_ui_label = FALSE), "|",
    format_number_suffix(best_case$cum_deaths, best_case$cum_deaths_low, best_case$cum_deaths_up, include_ui_label = FALSE), "|",
    format_number_suffix(worst_case$cum_deaths, worst_case$cum_deaths_low, worst_case$cum_deaths_up, include_ui_label = FALSE), "|"
  ))

  # Add DALY data if available
  if (!is.null(results$dalys)) {
    rows <- c(rows, paste0(
      "|**Cumulative DALYs**|",
      format_number_suffix(results$dalys$cum_dalys, results$dalys$cum_dalys_low, results$dalys$cum_dalys_up, include_ui_label = FALSE), "|",
      format_number_suffix(results$dalys$cum_dalys_best, results$dalys$cum_dalys_best_low, results$dalys$cum_dalys_best_up, include_ui_label = FALSE), "|",
      format_number_suffix(results$dalys$cum_dalys_worst, results$dalys$cum_dalys_worst_low, results$dalys$cum_dalys_worst_up, include_ui_label = FALSE), "|"
    ))
  }

  # Create header
  end_year <- 2100
  header <- c(
    paste0("# Pandemic Scenario Comparison (2026-", end_year, ")"),
    "",
    "|Metric|Baseline Projection (95% UI)|Best Case Scenario (95% UI)|Worst Case Scenario (95% UI)|",
    "|---|---|---|---|"
  )

  # Combine header and rows into a single string with line breaks
  md_table <- paste(c(header, rows), collapse = "\n")

  # Print the table to console for easy copying
  cat(md_table, "\n")

  # Save to file if requested (legacy parameter)
  if (!is.null(save_to_file)) {
    dir_path <- dirname(save_to_file)
    if (!dir.exists(dir_path) && dir_path != ".") {
      dir.create(dir_path, recursive = TRUE)
    }
    writeLines(md_table, save_to_file)
    message("Scenario table saved to: ", save_to_file)
  }
  
  # Save using output manager if enabled
  save_table_if_enabled(md_table, "scenario_table.md", output_dir, conf)
}


#' Generate a validation summary table for pandemic model results
#'
#' Creates a markdown table summarizing model validation results including
#' year-specific observed vs predicted comparisons, model performance metrics,
#' and overall validation statistics
#'
#' @param validation_results Results object from validate_model containing validation outputs
#' @param save_to_file Optional path to save the table to a markdown file
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_validation_table <- function(validation_results, save_to_file = NULL, output_dir = NULL, conf = NULL) {
  # Validate input
  if (is.null(validation_results$model) || is.null(validation_results$observed) || is.null(validation_results$predicted)) {
    stop("Validation results object must contain model, observed, and predicted components")
  }

  # Extract components
  model_results <- validation_results$model
  observed_data <- validation_results$observed
  predicted_data <- validation_results$predicted

  # Format function for values with suffix using Lancet formatting
  format_number_suffix <- function(value, low = NULL, upp = NULL, include_ui_label = FALSE) {
    .format_scaled_uncertainty(value, low, upp, digits_count = 0, include_ui_label = include_ui_label, is_count = TRUE)
  }

  # Create model information section
  model_info <- c(
    "## Model Training Information",
    "",
    "|Parameter|Value|",
    "|---|---|",
    "|Training Period|1600-1950|",
    "|Validation Period|1950-2025|",
    sprintf("|Training Sample Size|%s events|", format_lancet_number(model_results$specs$n, digits = 0, is_count = TRUE)),
    sprintf("|Population-scaled Threshold|%s deaths|", format_number_suffix(model_results$specs$Threshold)),
    sprintf("|Lower Cutoff|%s deaths|", format_number_suffix(model_results$specs$Cutoff)),
    sprintf(
      "|Shape Parameter (ξ)|%s|",
      if (!is.null(model_results$params$xi_ci_lower) && !is.null(model_results$params$xi_ci_upper)) {
        # Use pre-calculated confidence intervals from analysis
        format_lancet_uncertainty(
          model_results$params$xi,
          model_results$params$xi_ci_lower,
          model_results$params$xi_ci_upper,
          digits = 3, is_count = FALSE, include_ui_label = TRUE
        )
      } else {
        # Fallback to point estimate with SE
        paste0(
          format_lancet_number(model_results$params$xi, digits = 3, is_count = FALSE),
          " (SE: ", format_lancet_number(model_results$params$xi_se, digits = 3, is_count = FALSE), ")"
        )
      }
    ),
    sprintf(
      "|Scale Parameter (β)|%s|",
      if (!is.null(model_results$params$beta_ci_lower) && !is.null(model_results$params$beta_ci_upper)) {
        # Use pre-calculated confidence intervals from analysis
        format_lancet_uncertainty(
          model_results$params$beta,
          model_results$params$beta_ci_lower,
          model_results$params$beta_ci_upper,
          digits = 0, is_count = TRUE, include_ui_label = FALSE
        )
      } else {
        # Fallback to point estimate only
        format_number_suffix(model_results$params$beta)
      }
    ),
    ""
  )

  # Get predicted values for observation years
  # For validation, we need to interpolate or calculate predicted cumulative values
  obs_years <- observed_data$year
  predicted_cum <- predicted_data$result$cum_deaths
  predicted_cum_low <- predicted_data$result$cum_deaths_low
  predicted_cum_up <- predicted_data$result$cum_deaths_up

  # Create year-specific comparison table
  comparison_header <- c(
    "## Year-Specific Observed vs Predicted Comparison",
    "",
    "|Year|Observed Cumulative Deaths|Predicted Cumulative Deaths (95% UI)|",
    "|---|---|---|"
  )

  # Create comparison rows
  comparison_rows <- c()
  for (i in seq_len(nrow(observed_data))) {
    year <- observed_data$year[i]
    obs_val <- observed_data$obs[i]

    # Calculate predicted cumulative up to this year (simple approximation)
    # This assumes linear accumulation over the validation period
    years_elapsed <- year - 1950
    total_years <- 2025 - 1950
    pred_fraction <- years_elapsed / total_years

    pred_val <- predicted_cum * pred_fraction
    pred_low <- predicted_cum_low * pred_fraction
    pred_up <- predicted_cum_up * pred_fraction

    obs_formatted <- format_number_suffix(obs_val)
    pred_formatted <- format_number_suffix(pred_val, pred_low, pred_up, include_ui_label = (i == 1))

    comparison_rows <- c(comparison_rows, paste0("|", year, "|", obs_formatted, "|", pred_formatted, "|"))
  }

  # Create overall performance summary
  total_observed <- max(observed_data$obs)
  performance_summary <- c(
    "",
    "## Overall Validation Summary",
    "",
    "|Metric|Value|",
    "|---|---|",
    sprintf("|Total Observed Deaths (1950-2025)|%s|", format_number_suffix(total_observed)),
    sprintf("|Total Predicted Deaths (95%% UI)|%s|", format_number_suffix(predicted_cum, predicted_cum_low, predicted_cum_up, include_ui_label = FALSE)),
    sprintf(
      "|Annual Predicted Rate|%s events/year|",
      if (is.list(predicted_data$result$yearly_rate)) {
        # Extract uncertainty intervals if available
        rate_mean <- predicted_data$result$yearly_rate$mean
        rate_low <- if (!is.null(predicted_data$result$yearly_rate$lower_95)) {
          predicted_data$result$yearly_rate$lower_95
        } else {
          NULL
        }
        rate_upp <- if (!is.null(predicted_data$result$yearly_rate$upper_95)) {
          predicted_data$result$yearly_rate$upper_95
        } else {
          NULL
        }

        # Format with uncertainty if available
        if (!is.null(rate_low) && !is.null(rate_upp)) {
          format_lancet_uncertainty(
            rate_mean, rate_low, rate_upp,
            digits = 3, is_count = FALSE, include_ui_label = TRUE
          )
        } else {
          format_lancet_number(rate_mean, digits = 3, is_count = FALSE)
        }
      } else {
        format_lancet_number(predicted_data$result$yearly_rate, digits = 3, is_count = FALSE)
      }
    ),
    sprintf("|Expected Single Event Size|%s deaths|", format_number_suffix(model_results$results$boot, model_results$results$boot_low, model_results$results$boot_upp, include_ui_label = FALSE)),
    ""
  )

  # Create header
  header <- c(
    "# Pandemic Model Validation Results",
    ""
  )

  # Combine all sections
  full_table <- c(
    header,
    model_info,
    comparison_header,
    comparison_rows,
    performance_summary
  )

  # Combine into a single string with line breaks
  md_table <- paste(full_table, collapse = "\n")

  # Print the table to console for easy copying
  cat(md_table, "\n")

  # Save to file if requested
  if (!is.null(save_to_file)) {
    dir_path <- dirname(save_to_file)
    if (!dir.exists(dir_path) && dir_path != ".") {
      dir.create(dir_path, recursive = TRUE)
    }
    writeLines(md_table, save_to_file)
    message("Validation table saved to: ", save_to_file)
  }
  
  # Save using output manager if enabled
  save_table_if_enabled(md_table, "validation_table.md", output_dir, conf)

  return(md_table)
}

#' Generate a health burden gap table between scenarios
#'
#' Creates a markdown table summarizing the health burden gap between worst and best
#' plausible scenarios, including both deaths and DALYs with uncertainty intervals
#'
#' @param results Results object from run_PoT containing scenario gap calculations
#' @param save_to_file Optional path to save the table to a markdown file
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_gap_table <- function(results, save_to_file = NULL, output_dir = NULL, conf = NULL) {
  # Validate input
  if (is.null(results$scenarios)) {
    stop("Results object must contain scenario calculations")
  }

  # Check if gap calculations exist
  if (is.null(results$scenarios$deaths_gap)) {
    stop("Gap calculations not found. Ensure scenarios were calculated with the latest version of run_analysis()")
  }

  # Format function for values with suffix and uncertainty using Lancet formatting
  format_gap_value <- function(median_val, low_val, upp_val, is_percentage = FALSE, include_ui_label = FALSE) {
    if (is_percentage) {
      # Format as percentage with 1 decimal place using Lancet formatting
      return(format_lancet_uncertainty(
        median_val, low_val, upp_val,
        digits = 1, is_count = FALSE,
        include_ui_label = include_ui_label, scale_suffix = "%"
      ))
    } else {
      # Use Lancet scaled formatting for absolute values
      return(format_lancet_scaled(median_val, low_val, upp_val, is_count = TRUE, include_ui_label = include_ui_label))
    }
  }

  # Extract deaths gap data
  deaths_gap <- results$scenarios$deaths_gap
  deaths_absolute <- format_gap_value(
    deaths_gap$absolute$median,
    deaths_gap$absolute$lower_95,
    deaths_gap$absolute$upper_95,
    is_percentage = FALSE,
    include_ui_label = TRUE
  )
  deaths_relative <- format_gap_value(
    deaths_gap$relative$median,
    deaths_gap$relative$lower_95,
    deaths_gap$relative$upper_95,
    is_percentage = TRUE,
    include_ui_label = FALSE
  )

  # Start building rows
  rows <- c()

  # Deaths gap row
  rows <- c(rows, paste0(
    "|**Deaths**|",
    deaths_absolute, "|",
    deaths_relative, "|"
  ))

  # Check if DALY gap exists and add if available
  if (!is.null(results$scenarios$dalys_gap)) {
    dalys_gap <- results$scenarios$dalys_gap
    dalys_absolute <- format_gap_value(
      dalys_gap$absolute$median,
      dalys_gap$absolute$lower_95,
      dalys_gap$absolute$upper_95,
      is_percentage = FALSE,
      include_ui_label = FALSE
    )
    dalys_relative <- format_gap_value(
      dalys_gap$relative$median,
      dalys_gap$relative$lower_95,
      dalys_gap$relative$upper_95,
      is_percentage = TRUE,
      include_ui_label = FALSE
    )

    rows <- c(rows, paste0(
      "|**DALYs**|",
      dalys_absolute, "|",
      dalys_relative, "|"
    ))
  }

  # Create header
  header <- c(
    "# Health Burden Gap: Worst vs Best Plausible Scenarios (2026-2100)",
    "",
    "|Health Metric|Absolute Gap (95% UI)|Relative Gap (95% UI)|",
    "|---|---|---|"
  )

  # Add note about interpretation
  note <- c(
    "",
    "_Note: Gap represents the additional burden in the worst-case scenario compared to the best-case scenario._"
  )

  # Combine all parts
  md_table <- paste(c(header, rows, note), collapse = "\n")

  # Print the table to console
  cat(md_table, "\n")

  # Save to file if requested (legacy parameter)
  if (!is.null(save_to_file)) {
    dir_path <- dirname(save_to_file)
    if (!dir.exists(dir_path) && dir_path != ".") {
      dir.create(dir_path, recursive = TRUE)
    }
    writeLines(md_table, save_to_file)
    message("Gap table saved to: ", save_to_file)
  }
  
  # Save using output manager if enabled
  save_table_if_enabled(md_table, "gap_table.md", output_dir, conf)
}

#' Generate a parameter summary table for pandemic model results
#'
#' Creates a markdown table summarizing key model parameters including shape and scale
#' parameters for the GPD tail distribution, DALY-to-death ratio, and linear trend
#' coefficients with their uncertainty intervals and statistical significance.
#'
#' @param results Results object from run_analysis containing model parameters and estimates
#' @param save_to_file Optional path to save the table to a markdown file
#' @return A string containing the markdown table with proper line breaks
#' @export
generate_parameter_summary <- function(results, save_to_file = NULL, output_dir = NULL, conf = NULL) {
  # Validate input
  if (is.null(results$params)) {
    stop("Results object must contain model parameters")
  }

  # Initialize rows container
  rows <- c()

  #--------------------------------------------------
  # GPD TAIL PARAMETERS
  #--------------------------------------------------
  if (!is.null(results$params$xi) && !is.null(results$params$beta)) {
    # Shape parameter (ξ) - prefer bootstrap intervals if available
    xi_estimate <- results$params$xi
    xi_formatted <- NULL
    
    # Check for bootstrap-derived intervals first
    if (!is.null(results$params$xi_boot_median) && !is.null(results$params$xi_boot_lower) && !is.null(results$params$xi_boot_upper)) {
      xi_formatted <- format_lancet_uncertainty(
        results$params$xi_boot_median, results$params$xi_boot_lower, results$params$xi_boot_upper,
        digits = 3, is_count = FALSE, include_ui_label = TRUE
      )
    } else if (!is.null(results$params$xi_se)) {
      # Fallback to SE-based intervals
      xi_se <- results$params$xi_se
      xi_lower <- xi_estimate - 1.96 * xi_se
      xi_upper <- xi_estimate + 1.96 * xi_se
      xi_formatted <- format_lancet_uncertainty(
        xi_estimate, xi_lower, xi_upper,
        digits = 3, is_count = FALSE, include_ui_label = TRUE
      )
    } else {
      xi_formatted <- format_lancet_number(xi_estimate, digits = 3, is_count = FALSE)
    }

    rows <- c(rows, paste0(
      "|**Shape Parameter (ξ)**|",
      xi_formatted, "|",
      "GPD tail distribution|"
    ))

    # Scale parameter (β) - prefer bootstrap intervals if available
    beta_estimate <- results$params$beta
    beta_formatted <- NULL
    
    # Check for bootstrap-derived intervals first
    if (!is.null(results$params$beta_boot_median) && !is.null(results$params$beta_boot_lower) && !is.null(results$params$beta_boot_upper)) {
      beta_formatted <- format_lancet_uncertainty(
        results$params$beta_boot_median, results$params$beta_boot_lower, results$params$beta_boot_upper,
        digits = 0, is_count = TRUE, include_ui_label = FALSE
      )
    } else if (!is.null(results$params$beta_se)) {
      # Fallback to SE-based intervals
      beta_se <- results$params$beta_se
      beta_lower <- max(0, beta_estimate - 1.96 * beta_se)  # Ensure positive
      beta_upper <- beta_estimate + 1.96 * beta_se
      beta_formatted <- format_lancet_uncertainty(
        beta_estimate, beta_lower, beta_upper,
        digits = 0, is_count = TRUE, include_ui_label = FALSE
      )
    } else {
      beta_formatted <- format_lancet_number(beta_estimate, digits = 0, is_count = TRUE)
    }

    rows <- c(rows, paste0(
      "|**Scale Parameter (β)**|",
      beta_formatted, "|",
      "GPD tail distribution|"
    ))
  }

  #--------------------------------------------------
  # BULK DISTRIBUTION PARAMETERS
  #--------------------------------------------------
  # Add bulk distribution parameters if available from bootstrap
  if (!is.null(results$param_summary)) {
    bulk_params <- results$param_summary[results$param_summary$component == "bulk", ]
    
    if (nrow(bulk_params) > 0) {
      # Get the distribution type from the first bulk parameter
      dist_type <- if (!is.null(results$params$bulk_distribution)) {
        results$params$bulk_distribution
      } else {
        "Bulk distribution"
      }
      
      for (i in seq_len(nrow(bulk_params))) {
        param_name <- bulk_params$parameter[i]
        
        # Format parameter name for display
        display_name <- switch(param_name,
          "meanlog" = "Log Mean (μ)",
          "sdlog" = "Log SD (σ)", 
          "shape" = "Shape (α)",
          "scale" = "Scale (β)",
          "rate" = "Rate (β)",
          "bulk_prob" = "Bulk Probability",
          paste0(tools::toTitleCase(param_name))
        )
        
        # Format the parameter estimate with uncertainty
        param_formatted <- format_lancet_uncertainty(
          bulk_params$median[i], bulk_params$lower_95[i], bulk_params$upper_95[i],
          digits = if (param_name == "bulk_prob") 3 else 2,
          is_count = FALSE,
          include_ui_label = FALSE
        )
        
        rows <- c(rows, paste0(
          "|**", display_name, "**|",
          param_formatted, "|",
          paste0(tools::toTitleCase(dist_type), " bulk distribution|")
        ))
      }
    }
  }

  #--------------------------------------------------
  # DALY-TO-DEATH RATIO
  #--------------------------------------------------
  if (!is.null(results$dalys) && !is.null(results$dalys$daly_ratio)) {
    daly_ratio_formatted <- format_lancet_uncertainty(
      results$dalys$daly_ratio,
      results$dalys$daly_ratio_low,
      results$dalys$daly_ratio_up,
      digits = 1, is_count = FALSE, include_ui_label = FALSE
    )

    rows <- c(rows, paste0(
      "|**DALY-to-Death Ratio**|",
      daly_ratio_formatted, "|",
      "Burden conversion factor|"
    ))
  }

  #--------------------------------------------------
  # LINEAR TREND COEFFICIENT WITH P-VALUE
  #--------------------------------------------------
  # Check if time trend analysis was performed and extract trend coefficient
  if (!is.null(results$burden) && is.list(results$burden$yearly_rate)) {
    # Time trend was used - extract trend information
    yearly_rate_info <- results$burden$yearly_rate


    # Check if we have coefficients
    model_type <- if (!is.null(yearly_rate_info$model_type)) yearly_rate_info$model_type else "Unknown"

    # Extract beta_0 (intercept) if available
    if (!is.null(yearly_rate_info$beta_0_coefficient) &&
      !is.na(yearly_rate_info$beta_0_coefficient)) {
      beta_0_estimate <- yearly_rate_info$beta_0_coefficient
      beta_0_se <- yearly_rate_info$beta_0_se
      beta_0_pvalue <- yearly_rate_info$beta_0_pvalue

      # Format p-value according to standard conventions
      if (is.na(beta_0_pvalue)) {
        p_value_formatted <- ""
      } else if (beta_0_pvalue < 0.001) {
        p_value_formatted <- ", p<0.001"
      } else if (beta_0_pvalue < 0.01) {
        p_value_formatted <- paste0(", p=", format_lancet_number(beta_0_pvalue, digits = 3, is_count = FALSE))
      } else {
        p_value_formatted <- paste0(", p=", format_lancet_number(beta_0_pvalue, digits = 2, is_count = FALSE))
      }

      # Calculate 95% CI for beta_0 if SE is available
      if (!is.na(beta_0_se)) {
        beta_0_lower <- beta_0_estimate - 1.96 * beta_0_se
        beta_0_upper <- beta_0_estimate + 1.96 * beta_0_se

        # Format beta_0 with uncertainty
        beta_0_formatted <- format_lancet_uncertainty(
          beta_0_estimate, beta_0_lower, beta_0_upper,
          digits = 3, is_count = FALSE, include_ui_label = FALSE
        )
      } else {
        # Format without uncertainty if SE not available
        beta_0_formatted <- format_lancet_number(beta_0_estimate, digits = 3, is_count = FALSE)
      }

      rows <- c(rows, paste0(
        "|**Intercept (β₀)**|",
        beta_0_formatted, p_value_formatted, "|",
        "Log-scale intercept (", model_type, ")|"
      ))
    }

    # Check if we have the beta_1 coefficient and p-value
    if (!is.null(yearly_rate_info$beta_1_coefficient) &&
      !is.na(yearly_rate_info$beta_1_coefficient)) {
      beta_1_estimate <- yearly_rate_info$beta_1_coefficient
      beta_1_se <- yearly_rate_info$beta_1_se
      beta_1_pvalue <- yearly_rate_info$beta_1_pvalue

      # Format p-value according to standard conventions
      if (is.na(beta_1_pvalue)) {
        p_value_formatted <- ""
      } else if (beta_1_pvalue < 0.001) {
        p_value_formatted <- ", p<0.001"
      } else if (beta_1_pvalue < 0.01) {
        p_value_formatted <- paste0(", p=", format_lancet_number(beta_1_pvalue, digits = 3, is_count = FALSE))
      } else {
        p_value_formatted <- paste0(", p=", format_lancet_number(beta_1_pvalue, digits = 2, is_count = FALSE))
      }

      # Calculate 95% CI for beta_1 if SE is available
      if (!is.na(beta_1_se)) {
        beta_1_lower <- beta_1_estimate - 1.96 * beta_1_se
        beta_1_upper <- beta_1_estimate + 1.96 * beta_1_se

        # Format beta_1 with uncertainty
        beta_1_formatted <- format_lancet_uncertainty(
          beta_1_estimate, beta_1_lower, beta_1_upper,
          digits = 3, is_count = FALSE, include_ui_label = FALSE
        )
      } else {
        # Format without uncertainty if SE not available
        beta_1_formatted <- format_lancet_number(beta_1_estimate, digits = 3, is_count = FALSE)
      }

      rows <- c(rows, paste0(
        "|**Linear Trend Coefficient (β₁)**|",
        beta_1_formatted, p_value_formatted, "|",
        "Log-scale coefficient (", model_type, ")|"
      ))
    } else if (!is.null(yearly_rate_info$model_type) &&
      grepl("Constant", yearly_rate_info$model_type)) {
      # Constant rate model - no trend
      rows <- c(rows, paste0(
        "|**Linear Trend Coefficient (β₁)**|",
        "0.000 (reference)|",
        "No temporal trend (Constant Rate)|"
      ))
    }

    # Also add the annualized rate change if available
    if (!is.null(yearly_rate_info$trend_coefficient)) {
      trend_coeff <- yearly_rate_info$trend_coefficient

      if (!is.na(trend_coeff) && trend_coeff != 1.0) {
        # Convert to annualized rate change (as percentage)
        annual_change_pct <- (trend_coeff - 1) * 100

        # Check for uncertainty bounds
        if (!is.null(yearly_rate_info$trend_coefficient_lower) && !is.null(yearly_rate_info$trend_coefficient_upper)) {
          # Convert uncertainty bounds to annual rate change percentages
          annual_change_pct_lower <- (yearly_rate_info$trend_coefficient_lower - 1) * 100
          annual_change_pct_upper <- (yearly_rate_info$trend_coefficient_upper - 1) * 100
          
          # Format with uncertainty interval
          trend_formatted <- paste0(
            format_lancet_number(annual_change_pct, digits = 2, is_count = FALSE),
            "% (95% UI ", 
            format_lancet_number(annual_change_pct_lower, digits = 2, is_count = FALSE),
            "–",
            format_lancet_number(annual_change_pct_upper, digits = 2, is_count = FALSE),
            "%)"
          )
        } else {
          # Format without uncertainty interval (legacy format)
          trend_formatted <- paste0(
            format_lancet_number(annual_change_pct, digits = 2, is_count = FALSE), "%"
          )
        }

        rows <- c(rows, paste0(
          "|**Annual Rate Change**|",
          trend_formatted, "|",
          "Yearly change in exceedance rate|"
        ))
      }
    }
  } else {
    # Time trend analysis not performed
    rows <- c(rows, paste0(
      "|**Linear Trend Coefficient (β₁)**|",
      "Not available|",
      "Time trend analysis not performed|"
    ))
  }

  #--------------------------------------------------
  # CREATE TABLE
  #--------------------------------------------------
  # Create header
  header <- c(
    "# Model Parameter Summary",
    "",
    "|Parameter|Estimate (95% UI)|Description|",
    "|---|---|---|"
  )

  # Add note about methodology
  note <- c(
    "",
    "_Note: Uncertainty intervals calculated using bootstrap quantiles (2.5%-97.5%) when available, otherwise normal approximation (parameter ± 1.96 × SE)._"
  )

  # Combine all parts
  md_table <- paste(c(header, rows, note), collapse = "\n")

  # Print the table to console
  cat(md_table, "\n")

  # Save to file if requested (legacy parameter)
  if (!is.null(save_to_file)) {
    dir_path <- dirname(save_to_file)
    if (!dir.exists(dir_path) && dir_path != ".") {
      dir.create(dir_path, recursive = TRUE)
    }
    writeLines(md_table, save_to_file)
    message("Parameter summary table saved to: ", save_to_file)
  }
  
  # Save using output manager if enabled
  save_table_if_enabled(md_table, "parameter_summary.md", output_dir, conf)
}
