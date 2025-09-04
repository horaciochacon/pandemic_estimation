#' Time trend analysis functions (extracted from time_functions.R)
#'
#' Contains functions for analyzing pandemic exceedance rate trends and comparing temporal models.
#' Logic preserved verbatim during Phase 5 modularization.
# (Removed redundant library() calls; packages loaded via R/load_all.R)

#' Analyze pandemic exceedance rates with time trend
#'
#' Complete time series analysis framework for pandemic exceedance rates with
#' multiple model types and uncertainty quantification through Monte Carlo simulation.
#'
#' @inheritParams analyze_exceedance_rates
#' @param lower_cutoff Threshold value to define exceedance
#' @param end_year Ending year for the analysis (default: current year)
#' @param future_years Number of years to forecast (default: 75)
#' @param window_sizes Vector of window sizes to try (default from config)
#' @param n_sims Number of simulations for prediction intervals (default: 1000)
#' @param with_plots Whether to generate plots (default: TRUE)
#' @param model_type Type of temporal model to use ("linear", "quadratic",
#' "constant", "nb_linear", "nb_quadratic", or "auto")
#' @param event_only_windows If TRUE, uses non-overlapping windows of the specified size,
#' creating exactly one observation per time interval; dramatically reduces computation (default: FALSE)
#' @param conf Configuration object containing time.current_year for forecasting reference
#' @param generate_draws Whether to generate prediction draws for uncertainty propagation (default: FALSE)
#' @param n_draws Number of draws to generate for prediction intervals (default: 1000)
#'
#' @return See original documentation for structure
analyze_pandemic_trend <- function(
    data,
    variable = "deaths_scaled",
    lower_cutoff,
    start_year = NULL,
    end_year = NULL,
    future_years = if (!is.null(conf) && !is.null(conf$time$future_years)) conf$time$future_years else 75,
    window_sizes = if (!is.null(conf) && !is.null(conf$window$sizes)) conf$window$sizes else c(10, 25, 50, 75, 100),
    n_sims = if (!is.null(conf) && !is.null(conf$draws)) conf$draws else 1000,
    with_plots = TRUE,
    model_type = "linear",
    event_only_windows = FALSE,
    conf = NULL,
    generate_draws = FALSE,
    n_draws = if (!is.null(conf) && !is.null(conf$draws)) conf$draws else 1000) {
  if (is.null(start_year)) {
    start_year <- min(data$start_year)
  }
  if (is.null(end_year)) {
    end_year <- max(data$start_year, as.numeric(format(Sys.Date(), "%Y")))
  }
  filtered_data <- data %>%
    dplyr::filter(start_year >= start_year) %>%
    dplyr::arrange(start_year)
  exceedance_results <- analyze_exceedance_rates(
    data = filtered_data,
    variable = variable,
    threshold = lower_cutoff,
    base_year = start_year,
    window_sizes = window_sizes,
    plot_diagnostics = with_plots,
    event_only_windows = event_only_windows,
    conf = conf
  )
  model_comparison <- compare_temporal_models(
    data = filtered_data,
    variable = variable,
    threshold = lower_cutoff,
    window_size = exceedance_results$optimal_window,
    base_year = start_year,
    event_only_windows = event_only_windows,
    conf = conf
  )
  if (model_type != "auto") {
    valid_models <- names(model_comparison$models)
    model_map <- c(
      "constant" = "Constant Rate",
      "linear" = "Linear Trend",
      "quadratic" = "Quadratic Trend",
      "nb_linear" = "NB Linear",
      "nb_quadratic" = "NB Quadratic"
    )
    requested_model <- model_map[model_type]
    if (!requested_model %in% valid_models) {
      warning(paste0(
        "Requested model '", model_type, "' (", requested_model, ") not available. ",
        "Using best model: ", model_comparison$comparison$Model[1]
      ))
    } else {
      user_model <- model_comparison$models[[requested_model]]
      model_comparison$best_model <- user_model
      custom_order <- c(requested_model, setdiff(valid_models, requested_model))
      model_idx <- match(custom_order, model_comparison$comparison$Model)
      model_comparison$comparison <- model_comparison$comparison[model_idx, ]
      cat("Using user-specified model:", requested_model, "\n")
    }
  } else {
    cat("Using auto-selected best model:", model_comparison$comparison$Model[1], "\n")
  }
  forecast_results <- forecast_exceedance_rates(
    model_results = model_comparison,
    forecast_years = future_years,
    window_size = exceedance_results$optimal_window,
    prediction_interval = TRUE,
    n_sims = n_sims
  )
  historical_detailed <- exceedance_results$optimal_results$window_counts %>%
    dplyr::select(
      year = window_end,
      window_start,
      observed_rate = rate,
      exceedances,
      fitted_rate = annual_rate,
      predicted,
      ci_lower_rate,
      ci_upper_rate
    ) %>%
    dplyr::mutate(
      type = "historical",
      year_idx = year - min(year) + 1
    )
  historical_rates <- historical_detailed %>%
    dplyr::select(year, observed_rate, fitted_rate, type, year_idx)
  forecast_rates <- forecast_results$forecast_data %>%
    dplyr::select(year = window_end, fitted_rate = annual_rate, lower_95, upper_95) %>%
    dplyr::mutate(
      type = "forecast",
      year_idx = max(historical_rates$year_idx) + dplyr::row_number(),
      observed_rate = NA_real_
    )
  exceedance_rate_df <- dplyr::bind_rows(
    historical_detailed %>%
      dplyr::select(
        year, type, year_idx,
        observed_rate,
        fitted_rate,
        exceedances,
        lower_95 = ci_lower_rate,
        upper_95 = ci_upper_rate
      ),
    forecast_rates
  ) %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(
      lower_95 = dplyr::if_else(is.na(lower_95), fitted_rate * 0.5, lower_95),
      upper_95 = dplyr::if_else(is.na(upper_95), fitted_rate * 1.5, upper_95)
    )
  if (with_plots) {
    main_plot <- forecast_results$plot +
      labs(
        title = "Pandemic Exceedance Rate Projection",
        subtitle = paste0(
          "Events exceeding ", format(lower_cutoff, scientific = TRUE),
          " in variable '", variable, "'"
        )
      )
    model_plot <- model_comparison$model_plot +
      labs(
        title = "Model Selection for Rate Projection",
        subtitle = paste0("Best model: ", model_comparison$comparison$Model[1])
      )
    plots <- list(
      main_plot = main_plot,
      model_plot = model_plot,
      dashboard = main_plot
    )
  } else {
    plots <- NULL
  }
  draws_results <- NULL
  if (generate_draws) {
    cat("Generating prediction draws for uncertainty quantification...\n")
    best_model <- model_comparison$best_model
    window_counts <- model_comparison$window_counts
    window_size <- exceedance_results$optimal_window
    last_observed_year <- max(window_counts$window_end)
    current_year <- ifelse(is.null(conf$time$current_year),
      as.numeric(format(Sys.Date(), "%Y")),
      conf$time$current_year
    )
    forecast_end <- current_year + future_years
    years_seq <- min(window_counts$window_end):forecast_end
    years_df <- data.frame(window_end = years_seq)
    years_df$window_start <- years_df$window_end - window_size + 1
    years_df$year_centered <- years_df$window_end - min(window_counts$window_end)
    years_df$year_scaled <- years_df$year_centered / max(window_counts$year_centered)
    years_df$year_squared <- years_df$year_scaled^2
    years_df$offset <- window_size
    if (!is.null(conf$seed)) {
      set.seed(conf$seed + 3)
    } else if (!is.null(conf$seeds$secondary)) {
      set.seed(conf$seeds$secondary)
    } else {
      set.seed(42)
    }
    param_draws <- MASS::mvrnorm(
      n = n_draws,
      mu = coef(best_model),
      Sigma = vcov(best_model)
    )
    X <- model.matrix(formula(best_model)[-2], years_df)
    eta_draws <- X %*% t(param_draws)
    if (inherits(best_model, "negbin")) {
      theta <- best_model$theta
      annual_rates <- matrix(0, nrow = length(years_seq), ncol = n_draws)
      for (i in seq_len(n_draws)) {
        lambda <- exp(eta_draws[, i])
        if (!is.null(conf$seed)) {
          set.seed(conf$seed + 4 + i)
        }
        counts <- rnbinom(length(lambda), size = theta, mu = lambda)
        annual_rates[, i] <- counts / window_size
      }
    } else {
      lambda_draws <- exp(eta_draws)
      annual_rates <- matrix(0, nrow = length(years_seq), ncol = n_draws)
      for (i in seq_len(n_draws)) {
        counts <- rpois(length(years_seq), lambda = lambda_draws[, i])
        annual_rates[, i] <- counts / window_size
      }
    }
    draws_summary <- data.frame(
      year = years_seq,
      mean = rowMeans(annual_rates),
      median = apply(annual_rates, 1, median),
      lower_95 = apply(annual_rates, 1, quantile, 0.025),
      upper_95 = apply(annual_rates, 1, quantile, 0.975),
      lower_80 = apply(annual_rates, 1, quantile, 0.1),
      upper_80 = apply(annual_rates, 1, quantile, 0.9)
    )
    draws_summary$type <- ifelse(
      draws_summary$year <= current_year,
      "historical",
      "forecast"
    )
    if (with_plots) {
      n_sample_draws <- min(100, n_draws)
      sample_indices <- sample(seq_len(n_draws), n_sample_draws)
      traj_data <- data.frame(year = rep(years_seq, n_sample_draws))
      traj_data$draw_id <- rep(1:n_sample_draws, each = length(years_seq))
      traj_data$rate <- as.vector(annual_rates[, sample_indices])
      lancet_colors <- c(
        "Primary" = "#00468B",
        "Highlight" = "#ED0000",
        "Secondary" = "#42B540",
        "Tertiary" = "#0099B4",
        "Quaternary" = "#925E9F"
      )
      formatted_cutoff <- format(round(lower_cutoff), big.mark = " ")
      year_min <- floor(min(years_seq) / 50) * 50
      year_max <- ceiling(max(years_seq) / 50) * 50
      year_breaks <- seq(year_min, year_max, by = 50)
      show_labels <- if (is.null(conf$plots$draw_plot_labels)) FALSE else conf$plots$draw_plot_labels
      draw_plot <- ggplot() +
        geom_point(
          data = dplyr::filter(exceedance_rate_df, type == "historical"),
          aes(x = year, y = observed_rate),
          shape = 16,
          size = 2,
          fill = lancet_colors["Highlight"],
          alpha = 0.2
        ) +
        geom_ribbon(
          data = draws_summary,
          aes(x = year, ymin = lower_95, ymax = upper_95),
          fill = lancet_colors["Primary"],
          alpha = 0.15
        ) +
        geom_line(
          data = traj_data,
          aes(x = year, y = rate, group = draw_id),
          alpha = 0.03,
          color = lancet_colors["Primary"],
          linewidth = 0.05
        ) +
        geom_line(
          data = draws_summary,
          aes(x = year, y = mean),
          color = lancet_colors["Primary"],
          linewidth = 1
        ) +
        geom_vline(
          xintercept = current_year,
          linetype = "dashed",
          linewidth = 0.5,
          color = "grey50",
          alpha = 0.8
        ) +
        labs(
          title = if (show_labels) "Pandemic Exceedance Rate with Uncertainty" else NULL,
          subtitle = if (show_labels) {
            paste0(
              "Events exceeding ", formatted_cutoff, " deaths with ", window_size, "-year windows"
            )
          } else {
            NULL
          },
          x = "Year",
          y = paste0("Annual rate of >", round(lower_cutoff / 1000), "k\npopulation-scaled death pandemics"),
          caption = if (show_labels) "Shaded region shows 95% uncertainty interval" else NULL
        ) +
        scale_x_continuous(breaks = year_breaks) +
        theme_minimal() +
        theme(
          plot.title = element_text(
            family = "Times New Roman", size = 10, face = "bold",
            hjust = 0.5, margin = margin(b = 10)
          ),
          plot.subtitle = element_text(
            family = "Arial", size = 9, hjust = 0.5,
            margin = margin(b = 20)
          ),
          plot.caption = element_text(
            family = "Arial", size = 8, color = "gray30",
            hjust = 0, margin = margin(t = 10)
          ),
          legend.position = "bottom",
          legend.title = element_text(family = "Arial", size = 9, face = "bold"),
          legend.text = element_text(family = "Arial", size = 8),
          legend.key.size = unit(1.2, "lines"),
          legend.key = element_rect(fill = "white", color = NA),
          axis.title.x = element_text(family = "Arial", size = 11, margin = margin(t = 10)),
          axis.title.y = element_text(family = "Arial", size = 11, margin = margin(r = 10)),
          axis.text = element_text(family = "Arial", size = 8),
          axis.ticks = element_line(linewidth = 0.5),
          axis.line = element_line(color = "black", linewidth = 0.5),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(20, 20, 20, 20)
        )
      plots$draw_plot <- draw_plot
    }
    draws_results <- list(
      annual_rates = annual_rates,
      years = years_seq,
      param_draws = param_draws,
      summary = draws_summary,
      observed_years = exceedance_rate_df %>%
        dplyr::filter(type == "historical") %>%
        dplyr::mutate(observed = TRUE) %>%
        dplyr::select(year, observed_rate, observed)
    )
  }
  return(list(
    exceedance_rate_df = exceedance_rate_df,
    plots = plots,
    optimal_window = exceedance_results$optimal_window,
    model_comparison = model_comparison$comparison,
    best_model = model_comparison$best_model,
    forecast_summary = forecast_results$summary,
    waiting_time = forecast_results$waiting_summary,
    sampling_approach = ifelse(event_only_windows, "sparse", "yearly"),
    window_counts = model_comparison$window_counts,
    draws = draws_results
  ))
}

#' Analyze exceedance rates over time
#'
#' @return List containing analysis results and diagnostics
analyze_exceedance_rates <- function(data,
                                     variable = "deaths_scaled",
                                     threshold,
                                     base_year = 1700,
                                     window_sizes = c(10, 20, 30, 50, 75, 100),
                                     plot_diagnostics = TRUE,
                                     event_only_windows = FALSE,
                                     conf = NULL) {
  filtered_data <- data %>%
    dplyr::filter(start_year >= base_year) %>%
    dplyr::arrange(start_year)
  min_year <- min(filtered_data$start_year)
  max_year <- max(filtered_data$start_year)
  current_year <- ifelse(is.null(conf$time$current_year),
    as.numeric(format(Sys.Date(), "%Y")),
    conf$time$current_year
  )
  effective_max_year <- max(max_year, current_year)
  results <- list()
  results$data <- filtered_data
  results$parameters <- list(
    variable = variable,
    threshold = threshold,
    base_year = base_year,
    window_sizes = window_sizes
  )
  window_results <- list()
  for (window_size in window_sizes) {
    if (event_only_windows) {
      total_range <- effective_max_year - min_year
      num_windows <- floor(total_range / window_size)
      window_ends <- min_year + (1:num_windows) * window_size
      if (effective_max_year > max(window_ends)) {
        window_ends <- c(window_ends, effective_max_year)
      }
      cat(paste0(
        "Window size ", window_size, ": Creating ", length(window_ends),
        " non-overlapping windows (vs ", length(seq(min_year + window_size - 1, effective_max_year, by = 1)),
        " in regular approach).\n"
      ))
    } else {
      window_ends <- seq(min_year + window_size - 1, effective_max_year, by = 1)
    }
    window_starts <- window_ends - window_size + 1
    window_counts <- tibble::tibble(
      window_start = window_starts,
      window_end = window_ends,
      exceedances = 0,
      events = 0,
      rate = 0
    )
    for (i in seq_len(nrow(window_counts))) {
      window_start <- window_counts$window_start[i]
      window_end <- window_counts$window_end[i]
      window_events <- filtered_data %>%
        dplyr::filter(start_year >= window_start, start_year <= window_end)
      total_events <- nrow(window_events)
      exceedances <- sum(window_events[[variable]] > threshold)
      window_counts$exceedances[i] <- exceedances
      window_counts$events[i] <- total_events
      window_counts$rate[i] <- exceedances / window_size
    }
    window_stats <- tibble::tibble(
      window_size = window_size,
      total_exceedances = sum(window_counts$exceedances),
      total_events = sum(window_counts$events),
      mean_exceedances = mean(window_counts$exceedances),
      sd_exceedances = sd(window_counts$exceedances),
      var_exceedances = var(window_counts$exceedances),
      dispersion = var(window_counts$exceedances) / mean(window_counts$exceedances),
      zero_proportion = mean(window_counts$exceedances == 0),
      mean_rate = mean(window_counts$rate)
    )
    window_counts$year_center <- window_counts$window_end - min(window_counts$window_end)
    model <- glm(exceedances ~ year_center, family = poisson(link = "log"), data = window_counts)
    window_counts$predicted <- predict(model, type = "response")
    window_counts$annual_rate <- window_counts$predicted / window_size
    pred_ci <- predict(model,
      newdata = window_counts,
      type = "link", se.fit = TRUE
    )
    ci_data <- tibble::tibble(
      fit = pred_ci$fit,
      se = pred_ci$se.fit
    ) %>%
      dplyr::mutate(
        ci_lower = exp(fit - 1.96 * se),
        ci_upper = exp(fit + 1.96 * se),
        predicted_rate = exp(fit) / window_size,
        ci_lower_rate = ci_lower / window_size,
        ci_upper_rate = ci_upper / window_size
      )
    window_counts <- window_counts %>%
      dplyr::bind_cols(
        ci_data %>%
          dplyr::select(ci_lower, ci_upper, predicted_rate, ci_lower_rate, ci_upper_rate)
      )
    window_results[[as.character(window_size)]] <- list(
      window_counts = window_counts,
      window_stats = window_stats,
      model = model
    )
  }
  window_comparisons <- dplyr::bind_rows(
    lapply(window_results, function(x) x$window_stats)
  )
  if (plot_diagnostics) {
    dispersion_plot <- ggplot(window_comparisons, aes(x = window_size, y = dispersion)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Dispersion Index by Window Size",
        x = "Window Size (years)",
        y = "Dispersion (Var/Mean)"
      ) +
      theme_minimal()
    zero_prop_plot <- ggplot(window_comparisons, aes(x = window_size, y = zero_proportion)) +
      geom_line() +
      geom_point() +
      labs(
        title = "Proportion of Windows with Zero Exceedances",
        x = "Window Size (years)",
        y = "Proportion of Zeros"
      ) +
      theme_minimal()
    mean_var_plot <- ggplot(window_comparisons, aes(x = mean_exceedances, y = var_exceedances)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Mean-Variance Relationship",
        x = "Mean Exceedances",
        y = "Variance of Exceedances"
      ) +
      theme_minimal()
    time_plots <- list()
    for (i in seq_len(length(window_sizes))) {
      window_size <- window_sizes[i]
      window_data <- window_results[[as.character(window_size)]]$window_counts
      p <- ggplot(window_data, aes(x = window_end)) +
        geom_point(aes(y = rate), alpha = 0.5) +
        geom_line(aes(y = predicted_rate), color = "blue") +
        geom_ribbon(aes(ymin = ci_lower_rate, ymax = ci_upper_rate), alpha = 0.2, fill = "blue") +
        labs(
          title = paste0("Window Size: ", window_size, " years"),
          x = "Year (end of backwards-looking window)",
          y = "Annual Exceedance Rate"
        ) +
        theme_minimal()
      time_plots[[i]] <- p
    }
    results$diagnostics <- list(
      dispersion_plot = dispersion_plot,
      zero_prop_plot = zero_prop_plot,
      mean_var_plot = mean_var_plot,
      time_plots = time_plots
    )
  }
  dispersion_diff <- abs(window_comparisons$dispersion - 1)
  optimal_idx <- which.min(dispersion_diff)
  optimal_window <- window_comparisons$window_size[optimal_idx]
  results$optimal_window <- optimal_window
  results$optimal_results <- window_results[[as.character(optimal_window)]]
  optimal_data <- results$optimal_results$window_counts
  main_plot <- ggplot(optimal_data, aes(x = window_end)) +
    geom_line(aes(y = predicted_rate), color = "blue", linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lower_rate, ymax = ci_upper_rate), alpha = 0.2, fill = "blue") +
    geom_point(aes(y = rate), alpha = 0.6) +
    labs(
      title = paste0("Pandemic Exceedance Rate Over Time (Window: ", optimal_window, " years)"),
      subtitle = paste0(
        "Events exceeding ", format(threshold, scientific = TRUE, digits = 2),
        " in variable '", variable, "' (rates based on previous ", optimal_window, " years)"
      ),
      x = "Year (end of backwards-looking window)",
      y = "Annual Exceedance Rate"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  results$main_plot <- main_plot
  return(results)
}

#' Fit and compare multiple temporal models
#'
#' @return List of model comparison results
compare_temporal_models <- function(data,
                                    variable = "deaths_scaled",
                                    threshold,
                                    window_size = 50,
                                    base_year = 1700,
                                    event_only_windows = FALSE,
                                    conf = NULL) {
  filtered_data <- data %>%
    dplyr::filter(start_year >= base_year) %>%
    dplyr::arrange(start_year)
  min_year <- min(filtered_data$start_year)
  max_year <- max(filtered_data$start_year)
  current_year <- ifelse(is.null(conf$time$current_year),
    as.numeric(format(Sys.Date(), "%Y")),
    conf$time$current_year
  )
  effective_max_year <- max(max_year, current_year)
  if (event_only_windows) {
    total_range <- effective_max_year - min_year
    num_windows <- floor(total_range / window_size)
    window_ends <- min_year + (1:num_windows) * window_size
    if (effective_max_year > max(window_ends)) {
      window_ends <- c(window_ends, effective_max_year)
    }
    cat(paste0(
      "Creating ", length(window_ends),
      " non-overlapping windows for model comparison (vs ",
      length(seq(min_year + window_size - 1, effective_max_year, by = 1)),
      " in regular approach).\n"
    ))
  } else {
    window_ends <- seq(min_year + window_size - 1, effective_max_year, by = 1)
  }
  window_starts <- window_ends - window_size + 1
  window_counts <- tibble::tibble(
    window_start = window_starts,
    window_end = window_ends,
    exceedances = 0,
    events = 0,
    rate = 0
  )
  for (i in seq_len(nrow(window_counts))) {
    window_start <- window_counts$window_start[i]
    window_end <- window_counts$window_end[i]
    window_events <- filtered_data %>%
      dplyr::filter(start_year >= window_start, start_year <= window_end)
    total_events <- nrow(window_events)
    exceedances <- sum(window_events[[variable]] > threshold)
    window_counts$exceedances[i] <- exceedances
    window_counts$events[i] <- total_events
    window_counts$rate[i] <- exceedances / window_size
  }
  window_counts <- window_counts %>%
    dplyr::mutate(
      year_centered = window_end - min(window_end),
      year_scaled = year_centered / max(year_centered),
      year_squared = year_scaled^2
    )
  null_model <- glm(exceedances ~ 1,
    family = poisson(link = "log"),
    data = window_counts
  )
  linear_model <- glm(exceedances ~ year_scaled,
    family = poisson(link = "log"),
    data = window_counts
  )
  quadratic_model <- glm(exceedances ~ year_scaled + year_squared,
    family = poisson(link = "log"),
    data = window_counts
  )
  nb_model <- try(
    MASS::glm.nb(exceedances ~ year_scaled, data = window_counts),
    silent = TRUE
  )
  nb_quad_model <- try(
    MASS::glm.nb(exceedances ~ year_scaled + year_squared, data = window_counts),
    silent = TRUE
  )
  models <- list(
    "Constant Rate" = null_model,
    "Linear Trend" = linear_model,
    "Quadratic Trend" = quadratic_model
  )
  if (!inherits(nb_model, "try-error")) {
    models[["NB Linear"]] <- nb_model
  }
  if (!inherits(nb_quad_model, "try-error")) {
    models[["NB Quadratic"]] <- nb_quad_model
  }
  model_comparison <- tibble::tibble(
    Model = names(models),
    AIC = sapply(models, AIC),
    BIC = sapply(models, BIC),
    LogLik = sapply(models, function(m) as.numeric(logLik(m)))
  ) %>%
    dplyr::arrange(AIC)
  best_model_name <- model_comparison$Model[1]
  best_model <- models[[best_model_name]]
  prediction_data <- window_counts %>%
    dplyr::select(window_start, window_end, year_scaled, year_squared)
  for (model_name in names(models)) {
    model <- models[[model_name]]
    pred <- predict(model, newdata = prediction_data, type = "response")
    pred_rate <- pred / window_size
    prediction_data[[paste0(model_name, "_pred")]] <- pred_rate
  }
  prediction_data$actual_rate <- window_counts$rate
  pred_cols <- grep("_pred$", names(prediction_data), value = TRUE)
  models_data <- prediction_data %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(pred_cols),
      names_to = "model_name",
      values_to = "predicted_rate"
    ) %>%
    dplyr::mutate(
      model_name = gsub("_pred$", "", model_name),
      is_best = model_name == best_model_name
    )
  model_plot <- ggplot() +
    geom_point(
      data = prediction_data, aes(x = window_end, y = actual_rate),
      alpha = 0.4, color = "black"
    ) +
    geom_line(
      data = models_data %>%
        dplyr::left_join(prediction_data %>% dplyr::select(window_start, window_end), by = "window_start"),
      aes(
        x = window_end, y = predicted_rate, color = model_name,
        size = is_best, linetype = is_best
      )
    ) +
    scale_size_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide = "none") +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"), guide = "none") +
    scale_color_manual(
      values = c(
        "Constant Rate" = "red",
        "Linear Trend" = "blue",
        "Quadratic Trend" = "green",
        "NB Linear" = "purple",
        "NB Quadratic" = "orange"
      ),
      name = "Model Type"
    ) +
    labs(
      title = "Model Comparison for Pandemic Exceedance Rate",
      subtitle = paste0("Window Size: ", window_size, " years (looking back from each year)"),
      x = "Year (end of backwards-looking window)",
      y = "Annual Exceedance Rate"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  return(list(
    models = models,
    best_model = best_model,
    comparison = model_comparison,
    prediction_data = prediction_data,
    model_plot = model_plot,
    window_counts = window_counts,
    conf = conf
  ))
}
