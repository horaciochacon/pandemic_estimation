#' Initialize parallel processing
#'
#' Sets up parallel backend for computation with automatic detection of optimal worker count.
#' @param workers Number of parallel workers to use (default: NULL, which uses available cores minus 1)
#' @param strategy Character, parallel strategy to use: 'future', 'foreach', or 'none'
#' @return Invisible logical indicating if parallel backend was successfully initialized
#' @export
setup_parallel <- function(workers = NULL, strategy = "future") {
  # Check and validate strategy
  strategy <- match.arg(strategy, c("future", "foreach", "none"))

  if (strategy == "none") {
    message("Using sequential processing (no parallelization)")
    return(invisible(FALSE))
  }

  # Determine number of workers
  if (is.null(workers)) {
    # Default to number of available cores minus 1 (leave one for system)
    workers <- max(1, parallel::detectCores() - 1)
  }

  message(paste0("Setting up parallel processing with ", workers, " workers using ", strategy, " strategy"))

  if (strategy == "future") {
    # Future-based parallelization
    if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
      warning("Required packages 'future' and 'future.apply' not available. Using sequential processing.")
      return(invisible(FALSE))
    }

    # Set up future plan with specified workers
    future::plan(future::multisession, workers = workers)
    return(invisible(TRUE))
  } else if (strategy == "foreach") {
    # Foreach-based parallelization
    if (!requireNamespace("foreach", quietly = TRUE) ||
      !requireNamespace("doParallel", quietly = TRUE)) {
      warning("Required packages 'foreach' and 'doParallel' not available. Using sequential processing.")
      return(invisible(FALSE))
    }

    # Register parallel backend for foreach
    cl <- parallel::makeCluster(workers)
    doParallel::registerDoParallel(cl)

    # Store cluster object for cleanup
    assign("SENSITIVITY_CLUSTER", cl, envir = .GlobalEnv)
    return(invisible(TRUE))
  }
}

#' Stop parallel processing
#'
#' Shuts down parallel backend and cleans up resources.
#' @return Invisible NULL
#' @export
stop_parallel <- function() {
  # Check if future is in use
  if (requireNamespace("future", quietly = TRUE)) {
    # Reset future plan to sequential
    future::plan(future::sequential)
  }

  # Check if foreach cluster exists and stop it
  if (exists("SENSITIVITY_CLUSTER", envir = .GlobalEnv)) {
    parallel::stopCluster(get("SENSITIVITY_CLUSTER", envir = .GlobalEnv))
    rm(SENSITIVITY_CLUSTER, envir = .GlobalEnv)
  }

  return(invisible(NULL))
}

#' Run parallel sensitivity analysis
#'
#' Executes a sensitivity analysis across multiple parameter values in parallel
#' @param data Analysis dataset
#' @param var Variable to analyze ("deaths" or "deaths_scaled")
#' @param analysis_type Type of sensitivity analysis ("year", "threshold", "cutoff")
#' @param values Vector of values to analyze over (years, thresholds, or cutoff values)
#' @param plot_type Type of plot to generate
#' @param conf Configuration list
#' @param strategy Parallel strategy to use: "future", "foreach", or "none"
#' @param workers Number of workers for parallel processing (default: NULL uses available cores - 1)
#' @return A combined result object from the parallel runs
#' @export
run_parallel_sensitivity <- function(
    data,
    var,
    analysis_type = "year",
    values,
    plot_type = "shadow",
    conf,
    strategy = "future",
    workers = NULL) {
  # ---------------------------------------------------------------------------
  # Internal helpers (non-exported) to reduce repeated code blocks
  # ---------------------------------------------------------------------------
  get_analysis_fn <- function(analysis_type, data, var, conf) {
    switch(analysis_type,
      year = function(x) run_analysis(data, var, year = x, conf = conf),
      threshold = function(x) run_analysis(data, var, threshold = x, conf = conf),
      par = function(x) run_analysis(data, var, parameter = x, conf = conf),
      cutoff = function(x) {
        base_tail <- conf$thresholds$deaths_scaled %||% 2.8e6
        gap <- if (!is.null(conf$thresholds$lower_cutoff)) base_tail - conf$thresholds$lower_cutoff else (base_tail - 1e5)
        adj_threshold <- x + gap
        run_analysis(
          data, var,
          lower_cutoff = x,
          threshold = adj_threshold,
          conf = conf
        )
      },
      stop("Unsupported analysis_type: ", analysis_type)
    )
  }

  process_result <- function(result, analysis_type) {
    param_name <- switch(analysis_type,
      year = "year",
      threshold = "threshold",
      cutoff = "cutoff",
      "value"
    )

    base_tibble <- tibble::tibble(
      !!param_name := switch(analysis_type,
        year = result$specs$Year,
        threshold = result$specs$Threshold,
        cutoff = result$specs$Cutoff,
        NA_real_
      ),
      shadow = ifelse(!is.null(result$results$shadow), result$results$shadow, NA_real_),
      boot = ifelse(!is.null(result$results$boot), result$results$boot, NA_real_),
      boot_low = ifelse(!is.null(result$results$boot_low), result$results$boot_low, NA_real_),
      boot_upp = ifelse(!is.null(result$results$boot_upp), result$results$boot_upp, NA_real_),
      Xi = ifelse(!is.null(result$params$xi), result$params$xi, NA_real_),
      Beta = ifelse(!is.null(result$params$beta), result$params$beta, NA_real_)
    )

    if (!is.null(result$burden)) {
      burden_fields <- tibble::tibble(
        yearly_deaths = result$burden$yearly_deaths,
        yearly_deaths_low = result$burden$yearly_deaths_low,
        yearly_deaths_up = result$burden$yearly_deaths_up,
        cum_deaths = result$burden$cum_deaths,
        cum_deaths_low = result$burden$cum_deaths_low,
        cum_deaths_up = result$burden$cum_deaths_up,
        exp_annual_burden = result$burden$yearly_deaths
      )
      base_tibble <- dplyr::bind_cols(base_tibble, burden_fields)
    }

    base_tibble
  }

  # Wrapper used by all strategies
  analysis_fn <- get_analysis_fn(analysis_type, data, var, conf)
  run_single <- function(value) {
    res <- analysis_fn(value)
    process_result(res, analysis_type)
  }
  # Set up parallel backend
  is_parallel <- setup_parallel(workers, strategy)


  # Run analysis in parallel using the chosen strategy
  if (strategy == "future" && is_parallel) {
    # Use future.apply for parallelization
    raw_results <- future.apply::future_lapply(
      values,
      FUN = function(value) run_single(value),
      future.seed = TRUE,
      future.packages = c("dplyr", "tibble")
    )
  } else if (strategy == "foreach" && is_parallel) {
    # Use foreach for parallelization
    raw_results <- foreach::`%dopar%`(
      foreach::foreach(
        value = values,
        .packages = c("dplyr", "tibble")
      ), run_single(value)
    )
  } else {
    # Fall back to sequential processing
    raw_results <- lapply(values, run_single)
  }

  # Combine results
  results <- dplyr::bind_rows(raw_results) %>%
    # Scale values
    dplyr::mutate(dplyr::across(
      c(
        shadow, boot, boot_low, boot_upp,
        yearly_deaths, yearly_deaths_low, yearly_deaths_up,
        cum_deaths, cum_deaths_low, cum_deaths_up,
        exp_annual_burden, Beta
      ),
      ~ ifelse(!is.na(.x), .x / 1e6, .x)
    ))

  # Clean up parallel resources
  if (is_parallel) {
    stop_parallel()
  }

  # Plot function mapping
  plot_fns <- list(
    shadow = function(res, var, analysis_type) plot_shadow_sensitivity(res, var, analysis_type),
    yearly = function(res, var, analysis_type) plot_yearly_deaths(res, var, analysis_type),
    cumulative = function(res, var, analysis_type) plot_cumulative_deaths(res, var, analysis_type, conf),
    par = function(res, var, analysis_type) plot_par_sensitivity(res, var, analysis_type)
  )

  # Generate plot
  plot <- plot_fns[[plot_type]](results, var, analysis_type)

  list(results = results, plot = plot)
}
