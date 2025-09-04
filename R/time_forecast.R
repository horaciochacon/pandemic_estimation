#' Forecast utilities (extracted from time_functions.R)
#'
#' Functions supporting forecasting of exceedance rates.
#' Logic preserved verbatim; moved during Phase 5 modularization.
# (Removed redundant library() calls; packages loaded via R/load_all.R)

#' Prepare forecast data with model variables
prepare_forecast_data <- function(forecast_years_seq, window_counts, window_size) {
  forecast_data <- tibble::tibble(
    window_end = forecast_years_seq,
    window_start = window_end - window_size + 1
  )
  last_historical <- window_counts[nrow(window_counts), ]
  if ("year_centered" %in% names(window_counts)) {
    year_increment <- window_counts$year_centered[nrow(window_counts)] -
      window_counts$year_centered[nrow(window_counts) - 1]
    forecast_data$year_centered <- last_historical$year_centered +
      year_increment * (seq_len(nrow(forecast_data)))
    if ("year_scaled" %in% names(window_counts)) {
      max_centered <- max(window_counts$year_centered)
      forecast_data$year_scaled <- forecast_data$year_centered / max_centered
      if ("year_squared" %in% names(window_counts)) {
        forecast_data$year_squared <- forecast_data$year_scaled^2
      }
    }
  }
  forecast_data
}

#' Create prediction intervals using simulation
add_prediction_intervals <- function(forecast_data, model, window_size, n_sims) {
  pred_link <- predict(model, newdata = forecast_data, type = "link", se.fit = TRUE)
  forecast_data$fit <- pred_link$fit
  forecast_data$se.fit <- pred_link$se.fit
  sims <- matrix(0, nrow = nrow(forecast_data), ncol = n_sims)
  for (i in seq_len(n_sims)) {
    sim_values <- rnorm(nrow(forecast_data), forecast_data$fit, forecast_data$se.fit)
    sim_values <- exp(sim_values)
    if (inherits(model, "negbin")) {
      theta <- model$theta
      sim_values <- rnbinom(length(sim_values), size = theta, mu = sim_values)
    } else {
      sim_values <- rpois(length(sim_values), lambda = sim_values)
    }
    sims[, i] <- sim_values / window_size
  }
  forecast_data$lower_95 <- apply(sims, 1, quantile, probs = 0.025, na.rm = TRUE)
  forecast_data$upper_95 <- apply(sims, 1, quantile, probs = 0.975, na.rm = TRUE)
  list(forecast_data = forecast_data, sims = sims)
}

#' Process historical data for visualization
prepare_visualization_data <- function(model, window_counts, window_size, forecast_data) {
  vars_to_select <- intersect(
    c("window_start", "window_end", "year_centered", "year_scaled", "year_squared"),
    names(window_counts)
  )
  historical_pred_data <- window_counts %>%
    dplyr::select(dplyr::all_of(vars_to_select))
  historical_preds <- predict(model, newdata = historical_pred_data, type = "response")
  window_counts$model_predicted_rate <- historical_preds / window_size
  combined_data <- dplyr::bind_rows(
    window_counts %>%
      dplyr::select(window_end, annual_rate = rate) %>%
      dplyr::mutate(type = "Historical"),
    forecast_data %>%
      dplyr::select(window_end, annual_rate) %>%
      dplyr::mutate(type = "Forecast")
  )
  historical_fit <- window_counts %>%
    dplyr::select(window_end, fitted_rate = model_predicted_rate)
  list(
    window_counts = window_counts,
    combined_data = combined_data,
    historical_fit = historical_fit
  )
}

#' Create continuous prediction line across all time periods
create_continuous_line <- function(history_fit, forecast_data, window_counts,
                                   gap_years, model, window_size) {
  all_years <- seq(
    min(window_counts$window_end),
    max(forecast_data$window_end)
  )
  year_df <- tibble::tibble(window_end = all_years) %>%
    dplyr::left_join(history_fit, by = "window_end") %>%
    dplyr::left_join(
      forecast_data %>% dplyr::select(window_end, forecast_rate = annual_rate),
      by = "window_end"
    )
  if (length(gap_years) > 0) {
    gap_data <- prepare_forecast_data(gap_years, window_counts, window_size)
    gap_preds <- predict(model, newdata = gap_data, type = "response")
    gap_data$gap_rate <- gap_preds / window_size
    year_df <- year_df %>%
      dplyr::left_join(
        gap_data %>% dplyr::select(window_end, gap_rate),
        by = "window_end"
      )
  }
  if (!"gap_rate" %in% names(year_df)) {
    year_df$gap_rate <- NA_real_
  }
  year_df <- year_df %>%
    dplyr::mutate(
      fitted_rate = dplyr::case_when(
        !is.na(fitted_rate) ~ fitted_rate,
        !is.na(gap_rate) ~ gap_rate,
        !is.na(forecast_rate) ~ forecast_rate,
        TRUE ~ NA_real_
      )
    )
  year_df
}

#' Create forecast plot
create_forecast_plot <- function(combined_data, continuous_line, forecast_data,
                                 last_year, current_year, forecast_years, prediction_interval,
                                 window_size) {
  plot <- ggplot() +
    geom_point(
      data = dplyr::filter(combined_data, type == "Historical"),
      aes(x = window_end, y = annual_rate, shape = "Observed Values"),
      alpha = 0.5, color = "black"
    )
  plot <- plot +
    geom_line(
      data = continuous_line %>% dplyr::filter(!is.na(fitted_rate)),
      aes(
        x = window_end, y = fitted_rate,
        color = dplyr::case_when(
          window_end <= last_year ~ "Fitted Model (Historical)",
          window_end <= current_year ~ "Fitted Model (Gap to Present)",
          TRUE ~ "Forecast"
        )
      ),
      linewidth = 1
    )
  if (prediction_interval && "lower_95" %in% names(forecast_data)) {
    plot <- plot +
      geom_ribbon(
        data = forecast_data,
        aes(x = window_end, ymin = lower_95, ymax = upper_95, fill = "95% Prediction Interval"),
        alpha = 0.2
      )
  }
  max_rate <- max(continuous_line$fitted_rate, na.rm = TRUE)
  plot <- plot +
    labs(
      title = "Historical and Projected Pandemic Exceedance Rate",
      subtitle = paste0("Forecast for ", forecast_years, " years (rates show exceedances in preceding ", window_size, " years)"),
      x = "Year (end of backwards-looking window)",
      y = "Annual Exceedance Rate"
    ) +
    geom_vline(xintercept = last_year, linetype = "dashed") +
    geom_vline(xintercept = current_year, linetype = "dashed", color = "darkgreen") +
    annotate("text",
      x = last_year - 2, y = max_rate * 0.95,
      label = paste0("Last data: ", last_year), color = "black", hjust = 1, size = 3.5
    ) +
    annotate("text",
      x = current_year + 2, y = max_rate * 0.95,
      label = paste0("Current: ", current_year), color = "darkgreen", hjust = 0, size = 3.5
    ) +
    annotate("text",
      x = (last_year + current_year) / 2,
      y = max_rate * 0.8,
      label = "Gap Period (Zero Events)",
      color = "darkgreen", hjust = 0.5, size = 3, fontface = "italic"
    ) +
    scale_color_manual(
      values = c(
        "Fitted Model (Historical)" = "blue",
        "Fitted Model (Gap to Present)" = "darkgreen",
        "Forecast" = "red"
      ),
      name = ""
    ) +
    scale_shape_manual(values = c("Observed Values" = 16), name = "") +
    scale_fill_manual(values = c("95% Prediction Interval" = "red"), name = "") +
    guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  plot
}

#' Calculate forecast statistics
calculate_forecast_stats <- function(forecast_data, conf = NULL) {
  summary <- forecast_data %>%
    dplyr::summarise(
      mean_annual_rate = mean(annual_rate),
      median_annual_rate = median(annual_rate),
      min_annual_rate = min(annual_rate),
      max_annual_rate = max(annual_rate)
    )
  mean_rate <- summary$mean_annual_rate
  wt_draws <- if (!is.null(conf) && !is.null(conf$draws)) conf$draws else 10000
  waiting_times <- rexp(wt_draws, rate = mean_rate)
  waiting_quantiles <- quantile(waiting_times, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  waiting_summary <- tibble::tibble(
    statistic = c("2.5%", "25%", "50% (median)", "75%", "97.5%"),
    waiting_years = as.numeric(waiting_quantiles)
  )
  list(
    forecast_summary = summary,
    waiting_summary = waiting_summary
  )
}

#' Create forecast of future exceedance rates using best model
forecast_exceedance_rates <- function(model_results,
                                      forecast_years = if (!is.null(model_results$conf$time$future_years)) model_results$conf$time$future_years else 75,
                                      window_size = 50,
                                      prediction_interval = TRUE,
                                      n_sims = if (!is.null(model_results$conf$draws)) model_results$conf$draws else 1000) {
  model <- model_results$best_model
  window_counts <- model_results$window_counts
  last_year <- max(window_counts$window_end)
  current_year <- ifelse(
    !is.null(model_results$conf) && !is.null(model_results$conf$time$current_year),
    model_results$conf$time$current_year,
    as.numeric(format(Sys.Date(), "%Y"))
  )
  forecast_start <- current_year + 1
  forecast_end <- forecast_start + forecast_years - 1
  forecast_seq <- seq(forecast_start, forecast_end)
  forecast_data <- prepare_forecast_data(forecast_seq, window_counts, window_size)
  preds <- predict(model, newdata = forecast_data, type = "response")
  forecast_data$predicted <- preds
  forecast_data$annual_rate <- preds / window_size
  sim_results <- NULL
  if (prediction_interval) {
    interval_results <- add_prediction_intervals(
      forecast_data, model, window_size, n_sims
    )
    forecast_data <- interval_results$forecast_data
    sim_results <- interval_results$sims
  }
  viz_data <- prepare_visualization_data(
    model, window_counts, window_size, forecast_data
  )
  gap_years <- NULL
  if (last_year < current_year) {
    gap_years <- seq(last_year + 1, current_year)
  }
  continuous_line <- create_continuous_line(
    viz_data$historical_fit, forecast_data, window_counts,
    gap_years, model, window_size
  )
  forecast_plot <- create_forecast_plot(
    viz_data$combined_data, continuous_line, forecast_data,
    last_year, current_year, forecast_years, prediction_interval,
    window_size = window_size
  )
  stats <- calculate_forecast_stats(forecast_data)
  list(
    forecast_data = forecast_data,
    combined_data = viz_data$combined_data,
    plot = forecast_plot,
    summary = stats$forecast_summary,
    waiting_summary = stats$waiting_summary,
    simulations = sim_results
  )
}
