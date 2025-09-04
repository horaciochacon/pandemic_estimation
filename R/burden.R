#' Calculate global DALY to Death ratio using GBD data
#'
#' @param gbd_data_path Path to the GBD data file with DALY and death values
#'
#' @return A tibble with daly_ratio, daly_ratio_lower, and daly_ratio_upper
#' @export
calc_daly_ratio <- function(gbd_data_path = "data/gbd2023_global_DALY_death_ratio.csv") {
  # Read GBD data
  gbd_data <- read.csv(gbd_data_path)

  # Calculate total deaths
  deaths <- gbd_data |>
    dplyr::filter(measure == "Deaths") |>
    dplyr::pull(val) |>
    sum()

  # Calculate total DALYs
  dalys <- gbd_data |>
    dplyr::filter(measure == "DALYs") |>
    dplyr::pull(val) |>
    sum()

  # Calculate DALY to death ratio
  daly_ratio <- dalys / deaths

  # Include uncertainty bounds
  deaths_lower <- gbd_data |>
    dplyr::filter(measure == "Deaths") |>
    dplyr::pull(lower) |>
    sum()

  dalys_lower <- gbd_data |>
    dplyr::filter(measure == "DALYs") |>
    dplyr::pull(lower) |>
    sum()

  deaths_upper <- gbd_data |>
    dplyr::filter(measure == "Deaths") |>
    dplyr::pull(upper) |>
    sum()

  dalys_upper <- gbd_data |>
    dplyr::filter(measure == "DALYs") |>
    dplyr::pull(upper) |>
    sum()

  daly_ratio_lower <- dalys_lower / deaths_lower
  daly_ratio_upper <- dalys_upper / deaths_upper

  # Return as tibble
  daly_ratio_df <- tibble(daly_ratio, daly_ratio_lower, daly_ratio_upper)
  return(daly_ratio_df)
}

#' Get DALY to Death Ratio Draws for Uncertainty Propagation
#'
#' Retrieves pre-generated draws of the DALY-to-death ratio for uncertainty
#' quantification in burden estimation. These draws capture the uncertainty in
#' the relationship between deaths and disability-adjusted life years.
#'
#' @details
#' The function loads a set of pre-generated draws that represent the uncertainty
#' in the DALY-to-death ratio. Each draw is a possible value of the ratio based on
#' sampling from the joint distribution of DALYs and deaths. These draws are used
#' for Monte Carlo simulation in the burden calculation pipeline to properly
#' propagate uncertainty from multiple sources.
#'
#' @param path Path to the RDS file containing the DALY-to-death ratio draws
#' @return A vector of DALY-to-death ratio draws
#' @export
calc_daly_ratio_draws <- function(path = "data/daly_death_ratio_draws.rds") {
  ratio_draws <- readRDS(path)
  return(ratio_draws)
}


#' Calculate Expected Annual and Total Burden with Confidence Intervals
#'
#' @param data A data frame containing the event data.
#' @param var The variable of interest for which the burden is calculated.
#' @param thresh The threshold value for the variable of interest.
#' @param start The start year of the time period for analysis.
#' @param end The end year of the time period for analysis.
#' @param future_yrs The number of future years for which to forecast the burden.
#' @param plot_forecast Logical, if TRUE, a plot of the forecasted burden is returned.
#' @param shadow_mean The mean shadow value for the burden calculation.
#' @param shadow_low The lower bound shadow value for the burden calculation.
#' @param shadow_up The upper bound shadow value for the burden calculation.
#' @param resp Logical, if TRUE, only respiratory events are considered.
#' @param use_time_trend Logical, if TRUE, use time-variable rate for forecasting.
#' @param fixed_rate Optional fixed rate to use instead of calculating from data (for scenario analysis).
#' @param severity_draws Matrix of severity parameter draws for uncertainty propagation.
#' @param conf Configuration parameters from config.yml.
#' @param return_rate_draws Logical, if TRUE, returns the rate draws matrix (default: FALSE).
#' @param return_severity_draws Logical, if TRUE, returns the severity draws (default: FALSE).
#' @param return_yearly_deaths_draws Logical, if TRUE, returns the unadjusted yearly deaths draws (default: FALSE).
#' @param return_yearly_deaths_adjusted_draws Logical, if TRUE, returns the population-adjusted yearly deaths draws (default: TRUE).
#'
#' @return A list containing:
#'   \item{result}{Tibble with burden summary statistics}
#'   \item{plot}{Forecast plot if plot_forecast=TRUE}
#'   \item{rate_plot}{Rate trend plot if use_time_trend=TRUE and plot_forecast=TRUE}
#'   \item{rate_data}{Rate data if use_time_trend=TRUE}
#'   \item{yearly_deaths_adjusted_draws}{Matrix of population-adjusted yearly deaths draws if return_yearly_deaths_adjusted_draws=TRUE}
#'   \item{yearly_deaths_draws}{Matrix of unadjusted yearly deaths draws if return_yearly_deaths_draws=TRUE}
#'   \item{rate_draws}{Matrix of rate draws if return_rate_draws=TRUE}
#'   \item{severity_draws}{Matrix of severity parameter draws if return_severity_draws=TRUE}
#' @export
calc_burden <- function(
    data, var, thresh, start, end, future_yrs, plot_forecast = FALSE,
    shadow_mean, shadow_low, shadow_up, resp = FALSE, use_time_trend = FALSE,
    fixed_rate = NULL, severity_draws, conf,
    return_rate_draws = FALSE,
    return_severity_draws = FALSE,
    return_yearly_deaths_draws = FALSE,
    return_yearly_deaths_adjusted_draws = TRUE,
    return_forecast = FALSE,
    validation = FALSE,
    window_sizes = NULL,
    ...) {
  # Filter data for the time period
  data_period <- data %>%
    dplyr::filter(start_year >= start, end_year <= end)

  # Optionally filter for respiratory-only events
  if (resp) {
    data_period <- data_period %>% dplyr::filter(type == "Respiratory")
  }

  # Total time span in years
  time_span <- end - start + 1

  # Draws
  n_draws <- conf$draws

  # Get population forecast data
  pop_forecast <- get_pop_forecast() %>% dplyr::select(1, 2)
  pop_draws <- get_pop_draws()

  # Annual rate calculation - either static or time-variable approach
  if (!use_time_trend) {
    # Traditional approach - constant rate

    # Prepare for draws outputs
    yearly_deaths_draws <- NULL
    yearly_deaths_adjusted_draws <- NULL
    rate_draws <- NULL
    n_draws <- conf$draws

    # Process fixed rate for scenario analysis
    if (!is.null(fixed_rate)) {
      # Check if fixed_rate is a vector of draws or a single value
      if (is.vector(fixed_rate) && length(fixed_rate) > 1) {
        # Handle fixed_rate as draws for uncertainty propagation
        # Ensure rate_draws matches the number of severity draws (after potential xi filtering)
        actual_n_draws <- length(severity_draws)
        if (length(fixed_rate) > actual_n_draws) {
          # Resample fixed_rate to match severity_draws length
          fixed_rate_resampled <- sample(fixed_rate, actual_n_draws, replace = TRUE)
        } else if (length(fixed_rate) < actual_n_draws) {
          # This shouldn't happen, but handle it just in case
          fixed_rate_resampled <- sample(fixed_rate, actual_n_draws, replace = TRUE)
        } else {
          fixed_rate_resampled <- fixed_rate
        }

        rate_draws <- matrix(rep(fixed_rate_resampled, future_yrs), nrow = future_yrs, byrow = TRUE)

        # Calculate mean rate for summary
        lambda <- mean(fixed_rate)

        # Multiply rate and severity draws to propagate uncertainty
        # Similar to time trend approach in lines 238-242
        yearly_deaths_draws <- t(t(rate_draws) * severity_draws)

        # Resample population draws to match the actual number of severity draws
        # (which may be less than n_draws after xi filtering)
        idx_pop <- sample(seq_len(ncol(pop_draws) - 1), actual_n_draws, replace = TRUE)
        pop_draws_filtered <- pop_draws |>
          filter(year > end) |>
          dplyr::select(-year) |>
          as.matrix()
        pop_draws_filtered <- pop_draws_filtered[, idx_pop]

        # Calculate population-adjusted yearly deaths with draws
        yearly_deaths_adjusted_draws <- yearly_deaths_draws * pop_draws_filtered / conf$pop_reference

        # Calculate summary statistics for yearly deaths
        yearly_deaths_summary <- apply(
          yearly_deaths_adjusted_draws,
          MARGIN = 1,
          function(x) {
            # Using median for central tendency as deaths may be skewed
            mean_val <- median(x)
            quant <- quantile(x, prob = c(0.025, 0.975))
            c(mean_val, quant)
          }
        ) |>
          t() |>
          as_tibble(.name_repair = ~ c("mean", "lower", "upper")) |>
          mutate(year = (end + 1):(end + future_yrs), .before = mean)

        # Create forecast data frame with uncertainty from draws
        data_forecast <- tibble(
          year = (end + 1):(end + future_yrs),
          cum_years = 1:future_yrs,
          rate = lambda,
          pred_deaths = yearly_deaths_summary$mean,
          pred_deaths_low = yearly_deaths_summary$lower,
          pred_deaths_upp = yearly_deaths_summary$upper
        ) %>%
          mutate(
            cum_deaths = cumsum(pred_deaths),
            cum_deaths_low = cumsum(pred_deaths_low),
            cum_deaths_up = cumsum(pred_deaths_upp)
          )

        # Calculate conventional yearly values for compatibility
        yearly_deaths <- lambda * shadow_mean
        yearly_deaths_low <- lambda * shadow_low
        yearly_deaths_up <- lambda * shadow_up
      } else {
        # Traditional single-value fixed rate
        lambda <- fixed_rate

        # Expected annual burden (mean, lower, upper)
        yearly_deaths <- lambda * shadow_mean
        yearly_deaths_low <- lambda * shadow_low
        yearly_deaths_up <- lambda * shadow_up

        # Create forecast data frame with constant rate
        data_forecast <- tibble(
          year = (end + 1):(end + future_yrs),
          cum_years = 1:future_yrs,
          rate = lambda
        ) %>%
          left_join(pop_forecast, by = c("year" = "year")) %>%
          mutate(
            pred_deaths = yearly_deaths * val / conf$pop_reference,
            pred_deaths_low = yearly_deaths_low * val / conf$pop_reference,
            pred_deaths_upp = yearly_deaths_up * val / conf$pop_reference
          ) %>%
          mutate(
            cum_deaths = cumsum(pred_deaths),
            cum_deaths_low = cumsum(pred_deaths_low),
            cum_deaths_up = cumsum(pred_deaths_upp)
          )
      }
    } else {
      # Calculate rate from data - Number of exceedances over threshold
      events_over_thresh <- sum(threshold_extract(data_period, var, start, conf) > thresh)

      # Annual rate of exceedances (lambda)
      lambda <- events_over_thresh / time_span

      # Expected annual burden (mean, lower, upper)
      yearly_deaths <- lambda * shadow_mean
      yearly_deaths_low <- lambda * shadow_low
      yearly_deaths_up <- lambda * shadow_up

      # Create forecast data frame with constant rate
      data_forecast <- tibble(
        year = (end + 1):(end + future_yrs),
        cum_years = 1:future_yrs,
        rate = lambda
      ) %>%
        left_join(pop_forecast, by = c("year" = "year")) %>%
        mutate(
          pred_deaths = yearly_deaths * val / conf$pop_reference,
          pred_deaths_low = yearly_deaths_low * val / conf$pop_reference,
          pred_deaths_upp = yearly_deaths_up * val / conf$pop_reference
        ) %>%
        mutate(
          cum_deaths = cumsum(pred_deaths),
          cum_deaths_low = cumsum(pred_deaths_low),
          cum_deaths_up = cumsum(pred_deaths_upp)
        )
    }

    # Return results
    result <- tibble(
      yearly_rate = lambda,
      yearly_deaths,
      yearly_deaths_low,
      yearly_deaths_up,
      cum_deaths = data_forecast[future_yrs, ]$cum_deaths,
      cum_deaths_low = data_forecast[future_yrs, ]$cum_deaths_low,
      cum_deaths_up = data_forecast[future_yrs, ]$cum_deaths_up
    )
  } else {
    # Time-variable approach - use analyze_pandemic_trend function
    set.seed(42)
    trend_analysis <- analyze_pandemic_trend(
      data = data_period,
      variable = var,
      lower_cutoff = thresh,
      start_year = start,
      end_year = end,
      window_sizes = window_sizes,
      future_years = future_yrs,
      with_plots = plot_forecast,
      model_type = conf$time$model_type,
      n_draws = n_draws,
      generate_draws = TRUE,
      conf = conf
    )

    # Extract draws of prediction rates from time trend analysis
    rate_draws <- trend_analysis$draws$annual_rates

    # Get year indicator from time trend analysis
    years <- trend_analysis$draws$years

    # Filter only prediction period (years >= current year)
    current_year <- conf$time$current_year
    rate_draws <- rate_draws[years > current_year, ]

    # Get actual number of draws after potential xi filtering
    actual_n_draws <- length(severity_draws)

    # Ensure rate_draws matches the number of severity draws
    if (ncol(rate_draws) != actual_n_draws) {
      # Resample columns to match severity_draws length
      idx_rate <- sample(seq_len(ncol(rate_draws)), actual_n_draws, replace = TRUE)
      rate_draws <- rate_draws[, idx_rate, drop = FALSE]
    }

    # Multiply rate and severity draws, propagating uncertainty
    # This creates a matrix where each column is a draw of yearly deaths
    yearly_deaths_draws <- t(t(rate_draws) * severity_draws)

    if (validation) {
      pop_forecast <- pop_forecast |>
        filter(year > end, year <= end + future_yrs) |>
        pull(val)
      pop_draws_filtered <- matrix(pop_forecast, nrow = length(pop_forecast), ncol = actual_n_draws)
    } else if (!validation) {
      # Resample population draws to match the actual number of severity draws
      idx_pop <- sample(seq_len(ncol(pop_draws) - 1), actual_n_draws, replace = TRUE)
      pop_draws_filtered <- pop_draws |>
        filter(year > current_year) |>
        dplyr::select(-year) |>
        as.matrix()
      pop_draws_filtered <- pop_draws_filtered[, idx_pop]
    }

    # Estimate the draw-level population adjusted future yearly deaths
    # This accounts for changing global population over the forecast period
    yearly_deaths_adjusted_draws <- yearly_deaths_draws * pop_draws_filtered / conf$pop_reference

    # Calculate summary statistics for forecast rates
    forecast_rates <- apply(
      rate_draws,
      MARGIN = 1,
      function(x) {
        mean_val <- mean(x)
        quant <- quantile(x, prob = c(0.025, 0.975))
        c(mean_val, quant)
      }
    ) |>
      t() |>
      as_tibble(.name_repair = ~ c("fitted_rate", "lower_95", "upper_95")) |>
      mutate(year = (current_year + 1):(current_year + future_yrs), .before = fitted_rate)

    # Calculate summary statistics for yearly deaths
    yearly_deaths_summary <- apply(
      yearly_deaths_adjusted_draws,
      MARGIN = 1,
      function(x) {
        # Using median for central tendency as deaths may be skewed
        mean_val <- median(x)
        quant <- quantile(x, prob = c(0.025, 0.975))
        c(mean_val, quant)
      }
    ) |>
      t() |>
      as_tibble(.name_repair = ~ c("mean", "lower", "upper")) |>
      mutate(year = current_year:(current_year + future_yrs - 1), .before = mean)


    # Calculate year-by-year burden using time-variable rates
    forecast_burden <- forecast_rates %>%
      mutate(
        yearly_deaths = yearly_deaths_summary$mean,
        yearly_deaths_low = yearly_deaths_summary$lower,
        yearly_deaths_up = yearly_deaths_summary$upper
      ) %>%
      mutate(
        # Scale by population
        pred_deaths = yearly_deaths_summary$mean,
        pred_deaths_low = yearly_deaths_summary$lower,
        pred_deaths_upp = yearly_deaths_summary$upper
      )

    # Calculate cumulative values
    data_forecast <- forecast_burden %>%
      mutate(
        cum_deaths = cumsum(pred_deaths),
        cum_deaths_low = cumsum(pred_deaths_low),
        cum_deaths_up = cumsum(pred_deaths_upp),
        cum_years = seq_len(n())
      )

    # Average rate for summary statistics
    mean_rate <- mean(forecast_rates$fitted_rate)
    mean_rate_low <- mean(forecast_rates$lower_95)
    mean_rate_up <- mean(forecast_rates$upper_95)

    # Extract specific year rates
    # First year of forecast (typically 2026)
    first_year_rate <- forecast_rates$fitted_rate[1]
    first_year <- forecast_rates$year[1]
    
    # Extract uncertainty bounds for first year
    first_year_rate_lower <- forecast_rates$lower_95[1]
    first_year_rate_upper <- forecast_rates$upper_95[1]

    # Last year of forecast (typically 2100)
    last_year_rate <- forecast_rates$fitted_rate[nrow(forecast_rates)]
    last_year <- forecast_rates$year[nrow(forecast_rates)]
    
    # Extract uncertainty bounds for last year
    last_year_rate_lower <- forecast_rates$lower_95[nrow(forecast_rates)]
    last_year_rate_upper <- forecast_rates$upper_95[nrow(forecast_rates)]

    # Extract model information from trend_analysis
    best_model_name <- trend_analysis$model_comparison$Model[1]

    # Extract statistical model coefficients and p-values
    best_model <- NULL
    beta_0_coefficient <- NA
    beta_0_se <- NA
    beta_0_pvalue <- NA
    beta_1_coefficient <- NA
    beta_1_se <- NA
    beta_1_pvalue <- NA

    # Get the best model object directly from trend_analysis
    if (!is.null(trend_analysis$best_model)) {
      best_model <- trend_analysis$best_model
    }

    # Extract coefficient information if we have a model
    if (!is.null(best_model) && (inherits(best_model, "glm") || inherits(best_model, "negbin"))) {
      model_summary <- summary(best_model)

      # Extract coefficients
      coef_names <- rownames(model_summary$coefficients)
      
      # Look for intercept coefficient (beta_0)
      intercept_idx <- which(coef_names == "(Intercept)")
      if (length(intercept_idx) > 0) {
        beta_0_coefficient <- model_summary$coefficients[intercept_idx, "Estimate"]
        beta_0_se <- model_summary$coefficients[intercept_idx, "Std. Error"]
        beta_0_pvalue <- model_summary$coefficients[intercept_idx, "Pr(>|z|)"]
      }
      
      # Look for linear trend coefficient (beta_1)
      # This could be named "year_scaled", "year_centered", or similar
      trend_coef_idx <- which(grepl("year", coef_names, ignore.case = TRUE) & !grepl("squared|quad", coef_names, ignore.case = TRUE))

      if (length(trend_coef_idx) > 0) {
        # Take the first linear year coefficient
        idx <- trend_coef_idx[1]
        beta_1_coefficient <- model_summary$coefficients[idx, "Estimate"]
        beta_1_se <- model_summary$coefficients[idx, "Std. Error"]
        beta_1_pvalue <- model_summary$coefficients[idx, "Pr(>|z|)"]
      }
    }

    # Calculate trend coefficient as the ratio of last to first year rates
    # This gives us a measure of how much the rate changes over the forecast period
    trend_coefficient <- NA
    trend_coefficient_lower <- NA
    trend_coefficient_upper <- NA
    
    # First try to calculate from β₁ coefficient (parameter uncertainty only)
    # This avoids using prediction intervals that include random variation
    if (!is.na(beta_1_coefficient) && !is.na(beta_1_se) && beta_1_se > 0) {
      # The model uses year_scaled which goes from 0 to 1 over the historical period
      # Try to get actual historical span from config or use reasonable default
      historical_years_span <- if (!is.null(conf$time$current_year) && !is.null(conf$base_year)) {
        conf$time$current_year - conf$base_year
      } else {
        425  # Default assumption: 1600-2025
      }
      
      # Calculate annual rate change from β₁ coefficient directly
      # β₁ represents change in log(rate) as year_scaled increases by 1 (full historical span)
      # Annual change in log space = β₁ / historical_years_span
      annual_log_change <- beta_1_coefficient / historical_years_span
      trend_coefficient <- exp(annual_log_change)
      
      # Calculate confidence bounds using β₁ standard error (not prediction intervals)
      beta_1_lower <- beta_1_coefficient - 1.96 * beta_1_se
      beta_1_upper <- beta_1_coefficient + 1.96 * beta_1_se
      
      annual_log_change_lower <- beta_1_lower / historical_years_span
      annual_log_change_upper <- beta_1_upper / historical_years_span
      
      trend_coefficient_lower <- exp(annual_log_change_lower)
      trend_coefficient_upper <- exp(annual_log_change_upper)
      
    } else if (first_year_rate > 0) {
      # Fallback: Use forecast rate ratios (includes prediction uncertainty)
      # This method includes random variation and gives wider intervals
      total_change_ratio <- last_year_rate / first_year_rate

      # Annualized trend coefficient (geometric mean of growth)
      years_span <- last_year - first_year
      if (years_span > 0) {
        trend_coefficient <- total_change_ratio^(1 / years_span)
        
        # Calculate uncertainty bounds using prediction intervals
        # Note: This includes random variation, making intervals wider
        if (first_year_rate_upper > 0 && first_year_rate_lower > 0) {
          total_change_ratio_lower <- last_year_rate_lower / first_year_rate_upper
          total_change_ratio_upper <- last_year_rate_upper / first_year_rate_lower
          trend_coefficient_lower <- total_change_ratio_lower^(1 / years_span)
          trend_coefficient_upper <- total_change_ratio_upper^(1 / years_span)
        }
      } else {
        trend_coefficient <- 1.0
        trend_coefficient_lower <- 1.0
        trend_coefficient_upper <- 1.0
      }
    }

    # For constant rate models, the coefficient should be exactly 1
    if (best_model_name == "Constant Rate") {
      trend_coefficient <- 1.0
      trend_coefficient_lower <- 1.0
      trend_coefficient_upper <- 1.0
      beta_1_coefficient <- 0.0 # No trend for constant rate
      beta_1_se <- 0.0
      beta_1_pvalue <- 1.0 # Not significant
    }

    # Calculate rate uncertainty bounds as mean of forecast period uncertainty
    # This provides uncertainty intervals for the mean annual exceedance rate
    mean_rate_lower_95 <- mean(forecast_rates$lower_95)
    mean_rate_upper_95 <- mean(forecast_rates$upper_95)

    # Create detailed yearly_rate information
    yearly_rate_detailed <- list(
      mean = mean_rate,
      lower_95 = mean_rate_lower_95,
      upper_95 = mean_rate_upper_95,
      rate_2026 = if (first_year == 2026) first_year_rate else NA,
      rate_2100 = if (last_year == 2100) last_year_rate else NA,
      first_year_rate = first_year_rate,
      first_year = first_year,
      last_year_rate = last_year_rate,
      last_year = last_year,
      trend_coefficient = trend_coefficient,
      trend_coefficient_lower = trend_coefficient_lower,
      trend_coefficient_upper = trend_coefficient_upper,
      model_type = best_model_name,
      beta_0_coefficient = beta_0_coefficient,
      beta_0_se = beta_0_se,
      beta_0_pvalue = beta_0_pvalue,
      beta_1_coefficient = beta_1_coefficient,
      beta_1_se = beta_1_se,
      beta_1_pvalue = beta_1_pvalue
    )

    # Return results
    result <- tibble(
      yearly_rate = yearly_rate_detailed,
      yearly_deaths = mean_rate * median(severity_draws),
      yearly_deaths_low = mean_rate * quantile(severity_draws, prob = 0.025),
      yearly_deaths_up = mean_rate * quantile(severity_draws, prob = 0.975),
      cum_deaths = data_forecast[future_yrs, ]$cum_deaths,
      cum_deaths_low = data_forecast[future_yrs, ]$cum_deaths_low,
      cum_deaths_up = data_forecast[future_yrs, ]$cum_deaths_up,
      dynamic_rates = TRUE
    )
  }

  if (plot_forecast) {
    # Main forecast plot
    p <- ggplot(data_forecast, aes(x = year)) +
      geom_line(aes(y = cum_deaths, colour = "Predicted"), size = 1) +
      geom_ribbon(
        aes(
          ymin = cum_deaths_low, ymax = cum_deaths_up,
          fill = "Predicted"
        ),
        alpha = 0.2
      ) +
      labs(
        title = paste0("Expected Cumulative Deaths from ", end + 1, " to ", end + future_yrs),
        subtitle = ifelse(use_time_trend, "Using time-variable rates", "Using constant rate"),
        x = "Year",
        y = "Cumulative Deaths (in million)"
      ) +
      scale_y_log10(
        labels = scales::label_number(scale_cut = scales::cut_short_scale()),
        breaks = c(1e5, 1e6, 1e7, 1e8),
        limits = c(5e4, 4e8)
      ) +
      # Apply professional theme styling
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)),
        plot.caption = element_text(size = 9, color = "gray30", hjust = 1, margin = margin(t = 10)),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.title.x = element_text(margin = margin(t = 10), size = 11),
        axis.title.y = element_text(margin = margin(r = 10), size = 11),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(20, 20, 20, 20)
      )


    # Build return object with plot
    if (use_time_trend || (!is.null(fixed_rate) && is.vector(fixed_rate) && length(fixed_rate) > 1)) {
      # Create return list with basic elements
      return_list <- list(
        result = result,
        plot = p
      )

      # Add time-trend specific elements
      if (use_time_trend) {
        return_list$rate_plot <- trend_analysis$plots$draw_plot
        return_list$rate_data <- trend_analysis$exceedance_rate_df
      }

      # Add optional draw matrices based on function parameters
      if (return_yearly_deaths_adjusted_draws && !is.null(yearly_deaths_adjusted_draws)) {
        return_list$yearly_deaths_adjusted_draws <- yearly_deaths_adjusted_draws
      }
      if (return_yearly_deaths_draws && !is.null(yearly_deaths_draws)) {
        return_list$yearly_deaths_draws <- yearly_deaths_draws
      }
      if (return_rate_draws && !is.null(rate_draws)) {
        return_list$rate_draws <- rate_draws
      }
      if (return_severity_draws) {
        return_list$severity_draws <- severity_draws
      }
      if (return_forecast) {
        return_list$forecast <- data_forecast
      }



      return(return_list)
    } else {
      return(list(result = result, plot = p))
    }
  } else {
    # Build return object without plot
    if (use_time_trend || (!is.null(fixed_rate) && is.vector(fixed_rate) && length(fixed_rate) > 1)) {
      # Create return list with basic elements
      return_list <- list(
        result = result
      )

      # Add rate_data if time trend was used
      if (use_time_trend) {
        return_list$rate_data <- trend_analysis$exceedance_rate_df
      }

      # Add optional draw matrices based on function parameters
      if (return_yearly_deaths_adjusted_draws && !is.null(yearly_deaths_adjusted_draws)) {
        return_list$yearly_deaths_adjusted_draws <- yearly_deaths_adjusted_draws
      }
      if (return_yearly_deaths_draws && !is.null(yearly_deaths_draws)) {
        return_list$yearly_deaths_draws <- yearly_deaths_draws
      }
      if (return_rate_draws && !is.null(rate_draws)) {
        return_list$rate_draws <- rate_draws
      }
      if (return_severity_draws) {
        return_list$severity_draws <- severity_draws
      }
      if (return_forecast) {
        return_list$forecast <- data_forecast
      }

      return(return_list)
    } else {
      return(list(result = result))
    }
  }
}
