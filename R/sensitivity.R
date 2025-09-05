#' Get Year Range from Data
#'
#' @param data The input dataset containing start_year column
#' @param start_idx Starting index for the range (default: 15)
#' @param end_idx Ending index or number of years to include
#' @return Vector of ordered years within the specified range
#' @export
get_year_range <- function(data, conf) {
  data %>%
    distinct(start_year) %>%
    arrange(start_year) %>%
    slice(conf$sens$year$start_idx:conf$sens$year$end_idx) %>%
    pull(start_year)
}

#' Generic sensitivity analysis function
#' @param data Analysis dataset
#' @param var Variable to analyze ("deaths" or "deaths_scaled")
#' @param analysis_type Type of sensitivity analysis ("year", "threshold", "parameter")
#' @param values Vector of values to analyze over (years, thresholds, or parameter values)
#' @param plot_type Type of plot ("shadow", "yearly", "cumulative", "parameter")
#' @param conf Configuration list
#' @return List containing results and plot
run_sensitivity <- function(data,
                            var,
                            analysis_type = "year",
                            values,
                            plot_type = "shadow",
                            conf,
                            output_dir = NULL) {
  # Explicitly define valid options
  valid_analysis_types <- c("year", "threshold", "par", "cutoff")
  valid_plot_types <- c("shadow", "yearly", "cumulative", "par")

  # Validate inputs
  analysis_type <- match.arg(analysis_type, valid_analysis_types)
  plot_type <- match.arg(plot_type, valid_plot_types)

  # Analysis function mapping
  analysis_fns <- list(
    year = function(x) run_analysis(data, var, year = x, conf = conf),
    threshold = function(x) run_analysis(data, var, threshold = x, conf = conf),
    par = function(x) run_analysis(data, var, parameter = x, conf = conf),
    cutoff = function(x) {
      # When varying lower cutoff, preserve original spacing between lower cutoff and tail threshold.
      # Use configured deaths_scaled threshold if available; fall back to historical 2.8e6 offset.
      base_tail <- conf$thresholds$deaths_scaled %||% 2.8e6
      # Maintain relative gap: set threshold to max(x, base lower cutoff) + (base_tail - conf$thresholds$lower_cutoff)
      gap <- if (!is.null(conf$thresholds$lower_cutoff)) base_tail - conf$thresholds$lower_cutoff else (base_tail - 1e5)
      adj_threshold <- x + gap
      run_analysis(
        data, var,
        lower_cutoff = x,
        threshold = adj_threshold,
        conf = conf
      )
    }
  )

  # Process results
  raw_results <- map(values, analysis_fns[[analysis_type]])

  # Extract and combine results from the nested structure
  results <- map_df(raw_results, function(result) {
    # Create a base tibble with the analysis parameter
    param_name <- case_when(
      analysis_type == "year" ~ "year",
      analysis_type == "threshold" ~ "threshold",
      analysis_type == "cutoff" ~ "cutoff",
      TRUE ~ "value"
    )

    base_tibble <- tibble(
      !!param_name := case_when(
        analysis_type == "year" ~ result$specs$Year,
        analysis_type == "threshold" ~ result$specs$Threshold,
        analysis_type == "cutoff" ~ result$specs$Cutoff,
        TRUE ~ NA_real_
      ),
      # Add commonly used fields for sensitivity plots with appropriate mappings
      shadow = ifelse(!is.null(result$results$shadow), result$results$shadow, NA_real_),
      boot = ifelse(!is.null(result$results$boot), result$results$boot, NA_real_),
      boot_low = ifelse(!is.null(result$results$boot_low), result$results$boot_low, NA_real_),
      boot_upp = ifelse(!is.null(result$results$boot_upp), result$results$boot_upp, NA_real_),
      Xi = ifelse(!is.null(result$params$xi), result$params$xi, NA_real_),
      Beta = ifelse(!is.null(result$params$beta), result$params$beta, NA_real_)
    )

    # Add burden fields if available
    if (!is.null(result$burden)) {
      burden_fields <- tibble(
        yearly_deaths = result$burden$yearly_deaths,
        yearly_deaths_low = result$burden$yearly_deaths_low,
        yearly_deaths_up = result$burden$yearly_deaths_up,
        cum_deaths = result$burden$cum_deaths,
        cum_deaths_low = result$burden$cum_deaths_low,
        cum_deaths_up = result$burden$cum_deaths_up,
        exp_annual_burden = result$burden$yearly_deaths
      )
      base_tibble <- bind_cols(base_tibble, burden_fields)
    }

    base_tibble
  }) %>%
    # Scale values
    mutate(across(
      c(
        shadow, boot, boot_low, boot_upp,
        yearly_deaths, yearly_deaths_low, yearly_deaths_up,
        cum_deaths, cum_deaths_low, cum_deaths_up,
        exp_annual_burden, Beta
      ),
      ~ ifelse(!is.na(.x), .x / 1e6, .x)
    ))

  # Plot function mapping
  plot_fns <- list(
    shadow = function(res, var, analysis_type) plot_shadow_sensitivity(res, var, analysis_type, output_dir, conf),
    yearly = function(res, var, analysis_type) plot_yearly_deaths(res, var, analysis_type, output_dir, conf),
    cumulative = function(res, var, analysis_type) plot_cumulative_deaths(res, var, analysis_type, config = conf, output_dir = output_dir, conf = conf),
    par = function(res, var, analysis_type) plot_par_sensitivity(res, var, analysis_type, output_dir, conf)
  )

  # Generate plot
  plot_fns[[plot_type]](results, var, analysis_type)

  list(results = results)
}
