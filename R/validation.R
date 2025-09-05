#' Run Model Validation
#'
#' Performs out-of-sample validation by training on pre-1950 data and validating on post-1950
#'
#' @param data Complete dataset
#' @param variable Variable to analyze
#' @param conf Configuration object
#' @return List containing validation results and plots
validate_model <- function(data, variable, conf, output_dir = NULL) {
  # Split data into training and test sets
  train_data <- data %>%
    filter(
      start_year >= conf$train_year_start,
      start_year <= conf$train_year_end
    )

  # Run analysis on training data
  # Note: calc_burden should be disabled in the config
  train_results <- run_analysis(
    train_data,
    variable,
    conf = conf, # Pass validation config with calc_burden=false
    year = conf$train_year_start,
    return_severity_draws = TRUE
  )

  # Calculate validation period
  validation_years <- conf$test_year_end - conf$test_year_start

  # Calculate expected burden separately
  expected_burden <- calc_burden(
    data = train_data,
    var = variable,
    thresh = conf$thresholds[[paste0("lower_cutoff_", variable)]],
    start = conf$train_year_start,
    end = conf$train_year_end,
    future_yrs = validation_years,
    shadow_mean = train_results$results$boot,
    shadow_low = train_results$results$boot_low,
    shadow_up = train_results$results$boot_upp,
    resp = conf$respiratory_only,
    use_time_trend = conf$time$use_time_trend,
    conf = conf,
    plot_forecast = TRUE,
    window_sizes = conf$window_size,
    severity_draws = train_results$severity_draws,
    return_yearly_deaths_adjusted_draws = TRUE,
    return_forecast = TRUE,
    validation = TRUE,
    output_dir = output_dir
  )

  # Store the forecast plot
  forecast_plot <- expected_burden$plot

  # Process observed data
  observed_cumulative <- process_observed_data(
    data,
    variable,
    start = conf$test_year_start,
    end = conf$test_year_end,
    obs_years = conf$observation_years
  )

  # Create validation plot
  validation_plot <- create_validation_plot(
    forecast_plot,
    observed_cumulative,
    conf,
    variable
  )


  return(list(
    model = train_results,
    plot = validation_plot,
    observed = observed_cumulative,
    predicted = expected_burden,
    forecast = expected_burden$forecast |>
      filter(year %in% observed_cumulative$year) |>
      dplyr::select(year, cum_deaths:cum_deaths_up)
  ))
}

#' Process Observed Data for Validation
#'
#' @param data Complete dataset
#' @param variable Variable name to analyze
#' @param start_year Beginning of test period
#' @param end_year End of test period
#' @param obs_years Years to extract for comparison
#' @return Dataframe with cumulative observed values
process_observed_data <- function(data, variable, start, end, obs_years) {
  # Filter to test period
  observed_data <- data %>%
    filter(start_year >= start)

  pop <- get_pop_forecast() |>
    filter(year <= 2025) |>
    dplyr::select(year, pop = val)

  # Process observed data to get cumulative deaths by year
  observed_expanded <- observed_data %>%
    dplyr::select(name, start_year, end_year, deaths) |>
    mutate(length = end_year - start_year + 1) |>
    # Expand each event to yearly rows
    mutate(year = map2(start_year, end_year, seq)) %>%
    unnest(year) %>%
    # Calculate deaths per year (assuming uniform distribution within event)
    mutate(deaths_per_year = deaths / length) %>%
    dplyr::select(year, name, deaths_per_year) %>%
    # Transform to wide format with events as columns
    pivot_wider(names_from = name, values_from = deaths_per_year, values_fill = 0) |>
    arrange(year)

  # Create full yearly dataset and fill in missing years
  observed_yearly <- tibble(year = start:end) %>%
    left_join(observed_expanded, by = "year") %>%
    replace(is.na(.), 0)

  # Calculate cumulative deaths by year
  observed_cumulative <- observed_yearly %>%
    pivot_longer(-year) %>%
    group_by(year) %>%
    summarize(deaths = sum(value)) %>%
    mutate(obs = cumsum(deaths)) %>%
    dplyr::select(year, obs) |>
    # Filter to specific observation years for comparison
    filter(year %in% obs_years)

  return(observed_cumulative)
}

#' Create Enhanced Validation Plot
#'
#' @param forecast_plot Base forecast plot from burden calculation
#' @param observed_data Observed cumulative data
#' @param conf Configuration object
#' @param variable Variable name used in analysis
#' @return ggplot object with an enhanced validation plot
create_validation_plot <- function(forecast_plot, observed_data, conf, variable) {
  # Define Lancet color palette for consistency with Lancet guidelines
  lancet_colors <- c(
    "Predicted" = "#00468B", # Primary Lancet blue for predictions
    "Observed" = "#ED0000" # Lancet highlight red for observations
  )

  # Use custom colors if defined in config
  if (!is.null(conf$plot_colors$predicted)) {
    lancet_colors["Predicted"] <- conf$plot_colors$predicted
  }
  if (!is.null(conf$plot_colors$observed)) {
    lancet_colors["Observed"] <- conf$plot_colors$observed
  }

  # Create enhanced validation plot
  validation_plot <- forecast_plot +
    # Add observed data points with Lancet styling and more transparency
    geom_point(
      data = observed_data,
      aes(x = year, y = obs, colour = "Observed"),
      shape = 21, # Filled circle with border per Lancet style
      size = 3,
      stroke = 0.5, # Keep line weight ≤ 0.5 pt per Lancet guidelines
      fill = lancet_colors["Observed"],
      alpha = 0.7 # More transparent
    ) +
    # Add connecting line for observed points
    geom_line(
      data = observed_data,
      aes(x = year, y = obs, colour = "Observed"),
      linetype = "dashed",
      linewidth = 0.5, # Keep line weight ≤ 0.5 pt per Lancet guidelines
      alpha = 0.7 # More transparent
    ) +
    # Lancet color scales
    scale_color_manual(
      values = lancet_colors,
      name = "Data Source",
      labels = c(
        "Predicted" = "Model prediction",
        "Observed" = "Historical data"
      )
    ) +
    scale_fill_manual(
      values = alpha(lancet_colors["Predicted"], 0.3), # More transparent fill for uncertainty interval
      name = "Uncertainty",
      labels = c("Predicted" = "95% uncertainty interval")
    ) +
    # Conditionally add labels based on config
    labs(
      title = if (!is.null(conf$plots$validation_plot_labels) && conf$plots$validation_plot_labels) "Model Validation Results" else NULL,
      subtitle = if (!is.null(conf$plots$validation_plot_labels) && conf$plots$validation_plot_labels) {
        sprintf(
          "Training: %d–%d, Validation: %d–%d",
          conf$train_year_start,
          conf$train_year_end,
          conf$test_year_start + 1,
          conf$test_year_end
        )
      } else {
        NULL
      },
      caption = if (!is.null(conf$plots$validation_plot_labels) && conf$plots$validation_plot_labels) {
        sprintf(
          "Projection based on %d bootstrap samples",
          conf$draws
        )
      } else {
        NULL
      },
      x = "Year",
      y = "Cumulative pandemic-related deaths",
      fill = NULL,
      color = "Source"
    ) +
    scale_y_log10(
      labels = scales::label_number(),
      breaks = scales::breaks_log(5), limits = c(5e4, NA)
    ) +
    # Set x-axis breaks to 25-year intervals starting from 1950
    scale_x_continuous(
      breaks = seq(1950, 2025, by = 25)
    ) +
    # Apply Lancet theme styling
    theme_minimal() +
    theme(
      # Typography per Lancet guidelines
      plot.title = element_text(family = "Times New Roman", size = 10, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(family = "Arial", size = 9, hjust = 0.5, margin = margin(b = 15)),
      plot.caption = element_text(family = "Arial", size = 8, color = "gray30", hjust = 0, margin = margin(t = 10)),

      # Legend styling per Lancet guidelines
      legend.position = "bottom",
      legend.title = element_text(family = "Arial", size = 9, face = "bold"),
      legend.text = element_text(family = "Arial", size = 8),
      legend.key.size = unit(1.2, "lines"),
      legend.key = element_rect(fill = "white", color = NA),

      # Axis styling per Lancet guidelines
      axis.title.x = element_text(family = "Arial", size = 9, margin = margin(t = 10)),
      axis.title.y = element_text(family = "Arial", size = 9, margin = margin(r = 10)),
      axis.text = element_text(family = "Arial", size = 8),
      axis.ticks = element_line(linewidth = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.5),

      # Grid lines - removed per Lancet guidelines
      panel.grid = element_blank(),

      # Background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    )

  return(validation_plot)
}
