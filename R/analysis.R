#' Run Complete Pandemic Analysis on Dataset
#'
#' Performs a complete extreme value analysis on pandemic data using a two-stage
#' mixture model approach. This function orchestrates the entire analytical pipeline,
#' including data transformation, threshold selection, model fitting, bootstrap
#' simulation, burden calculation, and DALY estimation.
#'
#' @param data Data frame containing pandemic data with at least the columns:
#'   - deaths: Number of deaths in the pandemic event
#'   - start_year: Year when the pandemic started
#'   - pop_affected: Population affected by the pandemic (optional)
#' @param variable String, name of variable to analyze (typically "deaths")
#' @param conf Configuration object from config package containing analysis parameters:
#'   - pop_reference: Reference population for standardization
#'   - pop_min: Minimum population threshold for inclusion
#'   - base_year: Default starting year for analysis
#'   - thresholds: Various cutoff parameters for death counts
#'   - bootstrap: Configuration for uncertainty quantification
#'   - bulk_distribution: Distribution type for bulk modeling (lognormal, weibull, gamma)
#'   - use_constrained_model: Enable constrained survival mixture model
#'   - plots: Plot generation settings and options
#'   - calc_burden: Enable burden calculation
#'   - calc_dalys: Enable DALY estimation
#'   - calc_scenarios: Enable best/worst case scenario projections
#'   - time: Parameters for time trend modeling
#' @param threshold Optional manual threshold override for peaks-over-threshold analysis
#' @param year Optional year filter to restrict analysis to events after this year
#' @param lower_cutoff Optional override for the lower cutoff threshold
#' @param return_severity_draws Logical, if TRUE, the return object will include year-by-year severity draws matrix (rows: years 2026-2100, columns: draws) (default: FALSE)
#' @param return_cumulative_draws Logical, if TRUE, the return object will include cumulative_draws object (default: FALSE)
#' @return A list containing multiple analysis components:
#'   - specs: Data frame with model specifications and parameters
#'   - params: List of fitted model parameters (distribution types, shape, scale, etc.)
#'   - results: Expected values and point estimates with bootstrap confidence intervals
#'   - burden: Population-scaled death burden estimates (yearly and cumulative)
#'   - dalys: Global DALY estimates with uncertainty bounds
#'   - scenarios: Results from best/worst case scenario analysis (if calculated)
#'   - severity_draws: Year-by-year severity draws matrix with rows as years (2026-2100) and columns as draws (if return_severity_draws=TRUE)
#'   - cumulative_draws: Draw values for cumulative death estimates (if return_cumulative_draws=TRUE)
run_analysis <- function(
    data,
    variable,
    conf = config::get(),
    threshold = NULL,
    year = NULL,
    lower_cutoff = NULL,
    return_severity_draws = FALSE,
    return_cumulative_draws = FALSE) {
  #--------------------------------------------------
  # 1. DATA PREPARATION AND FILTERING
  #--------------------------------------------------
  # Set master random seed for reproducibility
  if (!is.null(conf$seed)) {
    set.seed(conf$seed)
    cat("Setting master random seed:", conf$seed, "\n")
  }

  # Keep copy of original data
  original_data <- data

  # Set base year
  base_year <- ifelse(is.null(year), conf$base_year, year)

  # Filter data
  data <- data |>
    filter(
      start_year >= base_year,
      deaths > conf$pop_min
    )

  # Generate descriptive plots if enabled
  if (conf$plots$descriptive) {
    print(create_descriptive_plots(data, variable))
  }

  # Generate time-scaled deaths plot if enabled
  if (conf$plots$time_scaled_deaths) {
    print(plot_time_scaled_deaths(data, conf$thresholds$lower_cutoff,
      reference_year = conf$time$current_year,
      include_labels = conf$plots$time_scaled_deaths_labels
    ))
  }

  #--------------------------------------------------
  # 2. DATA TRANSFORMATION AND THRESHOLD SELECTION
  #--------------------------------------------------
  # Apply dual transformation to standardize pandemic size relative to population
  data$dual <- dual_transform(data, variable, h = conf$pop_reference, l = conf$pop_min)

  # Determine threshold
  u_best <- determine_threshold(data, variable, threshold, conf)

  # Determine lower cutoff
  variable_specific_cutoff_name <- paste0("lower_cutoff_", variable)

  # Use the appropriate lower cutoff based on the variable
  if (!is.null(lower_cutoff)) {
    # User-specified cutoff takes precedence
    conf$thresholds$lower_cutoff <- lower_cutoff
  } else if (!is.null(conf$thresholds[[variable_specific_cutoff_name]])) {
    # Variable-specific cutoff from config
    conf$thresholds$lower_cutoff <- conf$thresholds[[variable_specific_cutoff_name]]
  }
  # Else use the default lower_cutoff already in conf

  #--------------------------------------------------
  # 3. INITIAL SUMMARY STATISTICS
  #--------------------------------------------------
  # Create initial summary of data characteristics and thresholds
  specs <- tibble(
    Variable = variable,
    Year = base_year,
    Threshold = inv_dual_transform(u_best, h = conf$pop_reference, l = conf$pop_min),
    Cutoff = conf$thresholds$lower_cutoff,
    n = nrow(data),
    PoT = sum(data$dual > u_best),
    PcoT = mean(data$dual > u_best),
    PoB = sum(data$dual > conf$thresholds$lower_cutoff)
  )

  #--------------------------------------------------
  # 4. INITIALIZE ANALYSIS OBJECTS
  #--------------------------------------------------
  # Initialize return objects and containers for results
  model_results <- NULL # Will store fitted model components
  boot_results <- NULL # Will store bootstrap results
  burden <- NULL # Will store burden calculation results
  params <- list() # Will store model parameters
  results <- list() # Will store key result metrics
  forecast <- list() # Will store forecasted values

  # Ensure scenario-related objects always exist (validation may disable burden/scenarios)
  scenario_results <- NULL
  worst_case <- NULL
  best_case <- NULL

  #--------------------------------------------------
  # 5. MIXTURE MODEL FITTING
  #--------------------------------------------------
  # Fit two-stage mixture model with flexible bulk and GPD tail
  model_results <- fit_mixture_model(
    data, variable, conf, u_best, conf$bulk_distribution,
    constrain_survival = conf$use_constrained_model,
    summary = specs
  )

  if (!is.null(model_results)) {
    #--------------------------------------------------
    # 6. EXTRACT MODEL PARAMETERS
    #--------------------------------------------------
    # Extract key parameters from fitted model
    params <- list(
      model_type = ifelse(is.null(model_results$bulk), "gpd_only", "mixture"),
      xi = model_results$tail$xi,
      xi_se = model_results$tail$xi_se,
      beta = model_results$tail$beta,
      beta_se = model_results$tail$beta_se
    )

    # Add confidence intervals using confint() if model fit is available
    if (!is.null(model_results$tail$fit)) {
      tryCatch(
        {
          ci <- confint(model_results$tail$fit)
          if (!is.null(ci) && nrow(ci) >= 2) {
            # Extract confidence intervals for both parameters
            params$xi_ci_lower <- ci[1, 1] # Shape parameter lower bound
            params$xi_ci_upper <- ci[1, 2] # Shape parameter upper bound
            params$beta_ci_lower <- ci[2, 1] # Scale parameter lower bound
            params$beta_ci_upper <- ci[2, 2] # Scale parameter upper bound
          }
        },
        error = function(e) {
          # If confint fails, calculate using normal approximation with constraints
          params$xi_ci_lower <<- model_results$tail$xi - 1.96 * model_results$tail$xi_se
          params$xi_ci_upper <<- model_results$tail$xi + 1.96 * model_results$tail$xi_se
          # Ensure beta confidence intervals are positive
          params$beta_ci_lower <<- max(0, model_results$tail$beta - 1.96 * model_results$tail$beta_se)
          params$beta_ci_upper <<- model_results$tail$beta + 1.96 * model_results$tail$beta_se
        }
      )
    }

    # Add bulk parameters if it's a mixture model
    if (!is.null(model_results$bulk)) {
      params$bulk_distribution <- model_results$bulk$dist_type
      params <- c(params, model_results$bulk$parameters)
      params$constrained_survival <- conf$use_constrained_model
    }

    #--------------------------------------------------
    # 7. VISUALIZATION
    #--------------------------------------------------
    # Generate plots if enabled in configuration
    if (conf$plots$tail) {
      tail_plot(
        data = data,
        variable = variable,
        fit = model_results$tail$fit,
        xi = model_results$tail$xi,
        u_best = u_best,
        beta = model_results$tail$beta,
        pop_ref = conf$pop_reference,
        shadow = model_results$mean,
        tail_limit = conf$plots$tail_limit
      )
    }

    # Note: tail_mixture plot moved to after bootstrap results are available

    #--------------------------------------------------
    # 8. UNCERTAINTY QUANTIFICATION
    #--------------------------------------------------
    # Initialize results list with shadow mean (expected value)
    results$shadow <- model_results$mean
    
    # Add point estimates for individual component means
    if (!is.null(model_results$bulk)) {
      results$bulk_mean <- model_results$bulk$conditional_mean
    } else {
      results$bulk_mean <- 0  # No bulk component in GPD-only model
    }
    results$tail_mean <- model_results$tail$shadow

    # Run bootstrap simulation if enabled for uncertainty estimation
    if (conf$bootstrap$enabled) {
      boot_results <- bootstrap_mixture_mean(
        data = data,
        variable = variable,
        conf = conf,
        mixture_pe = model_results$mean,
        u_best = u_best,
        summary = specs
      )

      # Add bootstrap results to results list
      results$boot <- boot_results$estimate
      results$boot_low <- boot_results$ci_lower
      results$boot_upp <- boot_results$ci_upper
      
      # Add bootstrap results for individual component means
      results$bulk_mean_boot <- boot_results$bulk_mean_estimate
      results$bulk_mean_boot_low <- boot_results$bulk_mean_ci_lower
      results$bulk_mean_boot_upp <- boot_results$bulk_mean_ci_upper
      results$tail_mean_boot <- boot_results$tail_mean_estimate
      results$tail_mean_boot_low <- boot_results$tail_mean_ci_lower
      results$tail_mean_boot_upp <- boot_results$tail_mean_ci_upper
      
      # Process parameter draws and create summaries (store separately to avoid tibble conversion issues)
      bootstrap_param_draws <- NULL
      bootstrap_param_summary <- NULL
      if (!is.null(boot_results$param_draws)) {
        bootstrap_param_summary <- calculate_parameter_summaries(boot_results$param_draws)
        bootstrap_param_draws <- boot_results$param_draws
        
        # Update params with bootstrap-derived uncertainty intervals
        params <- merge_bootstrap_params(params, bootstrap_param_summary)
      }
    }

    # Generate tail mixture plot after bootstrap results are available
    if (conf$plots$tail_mixture) {
      # Use bootstrapped median if available, otherwise use point estimate
      expected_severity_value <- if (conf$bootstrap$enabled && exists("boot_results")) {
        boot_results$estimate
      } else {
        NULL # Will fall back to mixture_fit$mean in the function
      }

      truncated_bulk_tail_plot(
        data = data,
        variable = variable,
        mixture_fit = model_results,
        pop_ref = conf$pop_reference,
        log10 = conf$plots$log_scale,
        y_title = conf$plots$tail_mixture_y_title,
        expected_severity = expected_severity_value,
        param_draws = bootstrap_param_draws,
        show_uncertainty = conf$bootstrap$enabled && !is.null(bootstrap_param_draws),
        shade_regions = FALSE,
        point_size = if (!is.null(conf$plots$point_size)) conf$plots$point_size else 3,
        severity_marker_height = if (!is.null(conf$plots$severity_marker_height)) conf$plots$severity_marker_height else 0.015
      )
    }

    #--------------------------------------------------
    # 9. BURDEN CALCULATION
    #--------------------------------------------------
    # Initialize burden variables for death forecasts
    burden <- NULL
    burden_results <- tibble()

    # Calculate population-scaled burden if enabled in config
    if (isTRUE(conf$calc_burden)) {
      burden <- calc_burden(
        data = original_data,
        var = variable,
        thresh = conf$thresholds$lower_cutoff,
        start = base_year,
        end = conf$time$current_year,
        future_yrs = conf$time$future_years,
        plot_forecast = conf$plots$forecast,
        shadow_mean = results$boot,
        shadow_low = results$boot_low,
        shadow_up = results$boot_upp,
        resp = conf$respiratory_only,
        use_time_trend = conf$time$use_time_trend,
        severity_draws = boot_results$severity_draws,
        conf = conf,
        return_yearly_deaths_draws = return_severity_draws,
        window_sizes = conf$window$default_size
      )

      if (conf$plots$forecast) {
        print(burden$plot)
      }

      if (conf$plots$time_trend) {
        print(burden$rate_plot)
      }

      # Add yearly rate to results
      results$yearly_rate <- burden$result$yearly_rate

      # Add burden results to specs
      # Handle both simple numeric and detailed list yearly_rate formats
      if (is.list(burden$result$yearly_rate[[1]])) {
        # Time trend case - extract mean rate for backward compatibility
        burden_results <- tibble(
          yearly_rate = burden$result$yearly_rate[[1]]$mean,
          yearly_rate_detailed = burden$result$yearly_rate,
          yearly_deaths = burden$result$yearly_deaths,
          yearly_deaths_low = burden$result$yearly_deaths_low,
          yearly_deaths_up = burden$result$yearly_deaths_up,
          cum_deaths = burden$result$cum_deaths,
          cum_deaths_low = burden$result$cum_deaths_low,
          cum_deaths_up = burden$result$cum_deaths_up
        )
      } else {
        # Constant rate case - keep simple structure
        burden_results <- tibble(
          yearly_rate = burden$result$yearly_rate,
          yearly_deaths = burden$result$yearly_deaths,
          yearly_deaths_low = burden$result$yearly_deaths_low,
          yearly_deaths_up = burden$result$yearly_deaths_up,
          cum_deaths = burden$result$cum_deaths,
          cum_deaths_low = burden$result$cum_deaths_low,
          cum_deaths_up = burden$result$cum_deaths_up
        )
      }

      #--------------------------------------------------
      # 10. SCENARIO ANALYSIS
      #--------------------------------------------------
      # This section identifies historical periods with extreme pandemic rates
      # and uses them to create best- and worst-case future projections
      scenario_results <- NULL
      worst_case <- NULL
      best_case <- NULL
      if (isTRUE(conf$calc_scenarios)) {
        # Use sliding window analysis to find periods with extreme pandemic rates
        window_size <- conf$scenario_window
        need_plot <- isTRUE(conf$plots$rate_extremes)

        rate_extremes <- find_rate_extremes(
          data = original_data,
          variable = variable,
          threshold = conf$thresholds$lower_cutoff,
          start_year = base_year,
          end_year = conf$time$current_year,
          window_size = window_size,
          return_all = need_plot,
          generate_draws = TRUE,
          n_draws = conf$draws
        )

        # Generate visualization if enabled
        if (need_plot) {
          rate_extremes_plot <- plot_rate_extremes(
            extreme_periods = rate_extremes,
            threshold_value = conf$thresholds$lower_cutoff
          )
          print(rate_extremes_plot)
        }

        if (conf$plots$scenario_rate_draws) {
          rate_extremes$draws$draw_plot
        }

        # Extract rate data from most extreme historical periods
        highest_rate <- rate_extremes$draws$highest_year_draws
        lowest_rate <- rate_extremes$draws$lowest_year_draws

        # Create period labels for reporting
        period_len <- if (!is.null(conf$tables$projection_years)) conf$tables$projection_years else (conf$time$future_years %||% getOption("pandemic.tables.projection_years", 75))
        highest_period <- paste0(
          rate_extremes$draws$highest_draw_year - period_len, "-",
          rate_extremes$draws$highest_draw_year
        )
        lowest_period <- paste0(
          rate_extremes$draws$lowest_draw_year - period_len, "-",
          rate_extremes$draws$lowest_draw_year
        )

        # Calculate worst-case scenario using highest historical pandemic rate
        # This represents a future where pandemics occur as frequently as in the worst historical period
        worst_case <- calc_burden(
          data = original_data,
          var = variable,
          thresh = conf$thresholds$lower_cutoff,
          start = base_year,
          end = conf$time$current_year,
          future_yrs = conf$time$future_years,
          plot_forecast = FALSE,
          shadow_mean = results$boot,
          shadow_low = results$boot_low,
          shadow_up = results$boot_upp,
          resp = conf$respiratory_only,
          use_time_trend = FALSE, # Force constant rate
          fixed_rate = highest_rate, # Override rate with highest observed
          severity_draws = boot_results$severity_draws,
          conf = conf,
          return_yearly_deaths_draws = return_severity_draws
        )

        # Calculate best-case scenario using lowest historical pandemic rate
        # This represents a future where pandemics occur as infrequently as in the best historical period
        best_case <- calc_burden(
          data = original_data,
          var = variable,
          thresh = conf$thresholds$lower_cutoff,
          start = base_year,
          end = conf$time$current_year,
          future_yrs = conf$time$future_years,
          plot_forecast = FALSE,
          shadow_mean = results$boot,
          shadow_low = results$boot_low,
          shadow_up = results$boot_upp,
          resp = conf$respiratory_only,
          use_time_trend = FALSE, # Force constant rate
          fixed_rate = lowest_rate, # Override with lowest observed rate
          severity_draws = boot_results$severity_draws,
          conf = conf,
          return_yearly_deaths_draws = return_severity_draws
        )

        # Compile scenario results
        scenario_results <- list(
          rate_extremes = rate_extremes$summary,
          rate_extremes_draws = rate_extremes$draws, # Store full draws for uncertainty quantification
          worst_case = worst_case$result,
          best_case = best_case$result
        )

        # Calculate health burden gap between scenarios at draw level
        if (!is.null(worst_case) && !is.null(best_case)) {
          # Get cumulative draws for both scenarios
          cumulative_worst_draws <- as.vector(colSums(worst_case$yearly_deaths_adjusted_draws))
          cumulative_best_draws <- as.vector(colSums(best_case$yearly_deaths_adjusted_draws))

          # Calculate gap for each draw (worst - best)
          gap_draws <- cumulative_worst_draws - cumulative_best_draws

          # Calculate summary statistics for absolute gap
          gap_median <- median(gap_draws)
          gap_lower <- quantile(gap_draws, 0.025)
          gap_upper <- quantile(gap_draws, 0.975)

          # Calculate relative gap (as percentage of worst case)
          relative_gap_draws <- gap_draws / cumulative_worst_draws * 100

          # Add gap calculations to scenario results
          scenario_results$deaths_gap <- list(
            absolute = list(
              median = gap_median,
              lower_95 = gap_lower,
              upper_95 = gap_upper
            ),
            relative = list(
              median = median(relative_gap_draws),
              lower_95 = quantile(relative_gap_draws, 0.025),
              upper_95 = quantile(relative_gap_draws, 0.975)
            ),
            draws = gap_draws # Keep draws for DALY calculation later
          )
        }
      }
    }

    # Populate forecast list
    forecast <- burden
  }

  #--------------------------------------------------
  # 11. DALY CALCULATIONS
  #--------------------------------------------------
  # Initialize DALY result containers
  daly_results <- NULL
  regional_daly_results <- NULL

  # Check if DALY calculations are enabled in config
  if (isTRUE(conf$calc_burden) && isTRUE(conf$calc_dalys) && isTRUE(conf$calc_scenarios) && !is.null(burden_results)) {
    # Calculate DALY to death ratio
    daly_ratio_draws <- calc_daly_ratio_draws()
    cumulative_draws <- as.vector(colSums(burden$yearly_deaths_adjusted_draws))
    cumulative_worst_draws <- as.vector(colSums(worst_case$yearly_deaths_adjusted_draws))
    cumulative_best_draws <- as.vector(colSums(best_case$yearly_deaths_adjusted_draws))

    # Resample daly-to-death ratio draws (with reproducible seed)
    if (!is.null(conf$seed)) {
      set.seed(conf$seed + 2) # Offset to avoid interfering with other random operations
    }
    daly_ratio_draws <- sample(daly_ratio_draws, conf$draws, replace = TRUE)

    # Multiply DALYs while propagating uncertainty
    cum_dalys_base_draws <- cumulative_draws * daly_ratio_draws
    cum_dalys_worst_draws <- cumulative_worst_draws * daly_ratio_draws
    cum_dalys_best_draws <- cumulative_best_draws * daly_ratio_draws

    # Calculate DALYs based on death estimates using only the central ratio value
    daly_results <- tibble(
      cum_dalys = median(cum_dalys_base_draws),
      cum_dalys_low = quantile(cum_dalys_base_draws, prob = 0.025),
      cum_dalys_up = quantile(cum_dalys_base_draws, prob = 0.975),
      daly_ratio = median(daly_ratio_draws),
      daly_ratio_low = quantile(daly_ratio_draws, prob = 0.025),
      daly_ratio_up = quantile(daly_ratio_draws, prob = 0.975),
      cum_dalys_best = median(cum_dalys_best_draws),
      cum_dalys_best_low = quantile(cum_dalys_best_draws, prob = 0.025),
      cum_dalys_best_up = quantile(cum_dalys_best_draws, prob = 0.975),
      cum_dalys_worst = median(cum_dalys_worst_draws),
      cum_dalys_worst_low = quantile(cum_dalys_worst_draws, prob = 0.025),
      cum_dalys_worst_up = quantile(cum_dalys_worst_draws, prob = 0.975)
    )

    # Calculate DALY gap between scenarios at draw level
    if (!is.null(scenario_results) && !is.null(cum_dalys_worst_draws) && !is.null(cum_dalys_best_draws)) {
      # Calculate gap for each draw (worst - best)
      dalys_gap_draws <- cum_dalys_worst_draws - cum_dalys_best_draws

      # Calculate summary statistics for absolute gap
      dalys_gap_median <- median(dalys_gap_draws)
      dalys_gap_lower <- quantile(dalys_gap_draws, 0.025)
      dalys_gap_upper <- quantile(dalys_gap_draws, 0.975)

      # Calculate relative gap (as percentage of worst case)
      dalys_relative_gap_draws <- dalys_gap_draws / cum_dalys_worst_draws * 100

      # Add DALY gap calculations to scenario results
      scenario_results$dalys_gap <- list(
        absolute = list(
          median = dalys_gap_median,
          lower_95 = dalys_gap_lower,
          upper_95 = dalys_gap_upper
        ),
        relative = list(
          median = median(dalys_relative_gap_draws),
          lower_95 = quantile(dalys_relative_gap_draws, 0.025),
          upper_95 = quantile(dalys_relative_gap_draws, 0.975)
        )
      )
    }
  }



  #--------------------------------------------------
  # 12. RESULT COMPILATION AND RETURN
  #--------------------------------------------------
  # Create and return the structured output
  result_list <- list(
    specs = specs,
    params = as_tibble(params),
    results = as_tibble(results),
    burden = burden_results,
    dalys = daly_results,
    scenarios = scenario_results
  )
  
  # Add parameter draws and summaries if available
  if (!is.null(bootstrap_param_draws)) {
    result_list$param_draws <- bootstrap_param_draws
  }
  if (!is.null(bootstrap_param_summary)) {
    result_list$param_summary <- bootstrap_param_summary
  }

  # Add severity draws if requested
  if (return_severity_draws) {
    if (!is.null(burden) && !is.null(burden$yearly_deaths_draws)) {
      # Preferred: use yearly death draws produced during burden calculation
      result_list$severity_draws <- burden$yearly_deaths_draws
      if (!is.null(worst_case)) {
        result_list$severity_draws_worst  <- worst_case$yearly_deaths_draws 
      }
      if (!is.null(best_case)) {
        result_list$severity_draws_best  <- best_case$yearly_deaths_draws
      }
      
      # Add population-adjusted draws for FHS upload (these match the scenario table values)
      if (!is.null(burden$yearly_deaths_adjusted_draws)) {
        result_list$yearly_deaths_adjusted_draws <- burden$yearly_deaths_adjusted_draws
      }
      if (!is.null(worst_case) && !is.null(worst_case$yearly_deaths_adjusted_draws)) {
        result_list$yearly_deaths_adjusted_draws_worst <- worst_case$yearly_deaths_adjusted_draws
      }
      if (!is.null(best_case) && !is.null(best_case$yearly_deaths_adjusted_draws)) {
        result_list$yearly_deaths_adjusted_draws_best <- best_case$yearly_deaths_adjusted_draws
      }
    } else if (!is.null(boot_results) && !is.null(boot_results$severity_draws)) {
      # Fallback: when burden is disabled (e.g., validation), expose bootstrap severity draws
      result_list$severity_draws <- boot_results$severity_draws
    }
  }

  # Add cumulative draws if requested
  if (return_cumulative_draws && exists("cumulative_draws")) {
    result_list$cumulative_draws <- cumulative_draws
  }

  return(result_list)
}

#' Calculate Parameter Summaries from Bootstrap Draws
#'
#' @param param_draws List containing parameter draws from bootstrap
#' @return Tibble with parameter summaries including median and 95% UI
calculate_parameter_summaries <- function(param_draws) {
  # Initialize results
  summaries <- tibble()
  
  # Tail parameters
  tail_params <- c("xi", "beta", "sigma", "tail_prob")
  for (param in tail_params) {
    if (param %in% names(param_draws) && length(param_draws[[param]]) > 0) {
      values <- param_draws[[param]][!is.na(param_draws[[param]])]
      if (length(values) > 0) {
        summaries <- bind_rows(summaries, tibble(
          parameter = param,
          component = "tail",
          median = median(values, na.rm = TRUE),
          lower_95 = quantile(values, 0.025, na.rm = TRUE),
          upper_95 = quantile(values, 0.975, na.rm = TRUE),
          n_valid = length(values)
        ))
      }
    }
  }
  
  # Bulk probability 
  if ("bulk_prob" %in% names(param_draws) && length(param_draws$bulk_prob) > 0) {
    values <- param_draws$bulk_prob[!is.na(param_draws$bulk_prob)]
    if (length(values) > 0) {
      summaries <- bind_rows(summaries, tibble(
        parameter = "bulk_prob",
        component = "bulk",
        median = median(values, na.rm = TRUE),
        lower_95 = quantile(values, 0.025, na.rm = TRUE),
        upper_95 = quantile(values, 0.975, na.rm = TRUE),
        n_valid = length(values)
      ))
    }
  }
  
  # Bulk distribution parameters
  if ("bulk_params" %in% names(param_draws) && length(param_draws$bulk_params) > 0) {
    for (bulk_param in names(param_draws$bulk_params)) {
      values <- param_draws$bulk_params[[bulk_param]]
      values <- values[!is.na(values) & values != 0]  # Remove NA and zero values
      if (length(values) > 0) {
        summaries <- bind_rows(summaries, tibble(
          parameter = bulk_param,
          component = "bulk",
          median = median(values, na.rm = TRUE),
          lower_95 = quantile(values, 0.025, na.rm = TRUE),
          upper_95 = quantile(values, 0.975, na.rm = TRUE),
          n_valid = length(values)
        ))
      }
    }
  }
  
  return(summaries)
}

#' Merge Bootstrap Parameter Summaries with Point Estimates
#'
#' @param params Original parameter list from point estimation
#' @param param_summary Parameter summaries from bootstrap draws
#' @return Updated parameter list with bootstrap-derived intervals
merge_bootstrap_params <- function(params, param_summary) {
  # Update tail parameters with bootstrap intervals
  tail_param_map <- list(
    "xi" = c("xi_boot_median", "xi_boot_lower", "xi_boot_upper"),
    "beta" = c("beta_boot_median", "beta_boot_lower", "beta_boot_upper"),
    "sigma" = c("sigma_boot_median", "sigma_boot_lower", "sigma_boot_upper"),
    "tail_prob" = c("tail_prob_boot_median", "tail_prob_boot_lower", "tail_prob_boot_upper")
  )
  
  for (param in names(tail_param_map)) {
    param_row <- param_summary[param_summary$parameter == param & param_summary$component == "tail", ]
    if (nrow(param_row) > 0) {
      param_names <- tail_param_map[[param]]
      params[[param_names[1]]] <- param_row$median
      params[[param_names[2]]] <- param_row$lower_95
      params[[param_names[3]]] <- param_row$upper_95
    }
  }
  
  # Add bulk parameter summaries
  bulk_params <- param_summary[param_summary$component == "bulk", ]
  for (i in seq_len(nrow(bulk_params))) {
    param_name <- bulk_params$parameter[i]
    base_name <- paste0(param_name, "_boot")
    params[[paste0(base_name, "_median")]] <- bulk_params$median[i]
    params[[paste0(base_name, "_lower")]] <- bulk_params$lower_95[i]
    params[[paste0(base_name, "_upper")]] <- bulk_params$upper_95[i]
  }
  
  return(params)
}
