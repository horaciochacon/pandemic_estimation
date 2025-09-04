# Required Packages --------------------------------------------------------------
library(tidyverse)
library(evir)
library(gridExtra)
library(scales)
library(purrr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)
library(config)
library(data.table)

#' Load Configuration Settings
#'
#' @param config_file Path to config file
#' @param environment Configuration environment to use
#' @return List containing configuration settings
#' @export
load_config <- function(
    config_file = "config.yml",
    environment = Sys.getenv("R_CONFIG_ACTIVE", "default")) {
  # Load the specified configuration
  conf <- config::get(file = config_file, config = environment)

  # Apply population-scaled threshold calculations for validation configuration
  if (environment == "validation") {
    # Load default config for reference values
    default_conf <- config::get(file = config_file, config = "default")

    # Calculate scaled thresholds based on population ratio
    conf <- calculate_population_scaled_thresholds(default_conf, conf)
  }

  return(conf)
}

#' Read and Transform Analysis Data
#'
#' Reads the data file and applies necessary transformations according to configuration
#'
#' @param conf List containing configuration settings
#' @return A tibble with transformed data
#' @export
read_transform_data <- function(conf) {
  read_csv(conf$file_path) %>%
    mutate(deaths = deaths * conf$deaths_unit) %>%
    mutate(
      length = end_year - start_year + 1,
      deaths_scaled = deaths * (conf$pop_reference / pop),
      deaths_scaled_year = (deaths * (conf$pop_reference / pop)) / length
    )
}

#' Get IHME Population Forecast Data
#'
#' Retrieves and processes IHME's population forecast global values.
#'
#' @return A tibble containing processed population forecast data
#' @export
get_pop_forecast <- function() {
  readr::read_csv("data/pop_forecast.csv", show_col_types = FALSE) %>%
    dplyr::select(
      year = .data$year_id,
      .data$sex,
      .data$val,
      .data$upper,
      .data$lower
    ) %>%
    dplyr::group_by(.data$year) %>%
    dplyr::summarize(
      val = sum(.data$val),
      upper = sum(.data$upper),
      lower = sum(.data$lower),
      .groups = "drop"
    )
}


#' Get IHME Population Forecast Draws
#'
#' Retrieves and processes IHME's population forecast global values with draws.
#' These draws represent uncertainty in future population projections and are
#' used for propagating uncertainty in burden calculations.
#'
#' @return A tibble containing processed population forecast draws, with years as rows
#'   and draws as columns (draw_1, draw_2, etc.)
#' @export
get_pop_draws <- function() {
  readr::read_csv("data/pop_draws.csv", show_col_types = FALSE) %>%
    dplyr::select(
      year = year_id,
      draw,
      value
    ) |>
    pivot_wider(names_from = draw, values_from = value, names_prefix = "draw_")
}

#' Get IHME Regional Population Forecast Draws
#'
#' Retrieves and processes IHME's population forecast data at the super-region level with draws.
#' This function maps location IDs to human-readable region names and formats the data
#' for use in regional burden calculations with uncertainty propagation.
#'
#' @details
#' The function loads regional population draws from the data file and joins them with
#' super-region definitions to provide meaningful region names. The draws represent
#' uncertainty in regional population projections over time and are used to generate
#' confidence intervals for regional burden estimates.
#'
#' The super-regions used are:
#' - Southeast Asia, East Asia, and Oceania
#' - Central Europe, Eastern Europe, and Central Asia
#' - High-income
#' - Latin America and Caribbean
#' - North Africa and Middle East
#' - South Asia
#' - Sub-Saharan Africa
#'
#' @return A tibble containing processed regional population forecast draws, with
#'   location ID, location name, year, and multiple draw columns (draw_1, draw_2, etc.)
#' @export
get_regional_pop_draws <- function() {
  super_region <- get_super_region_mapping()

  readr::read_csv("data/pop_regional_draws.csv", show_col_types = FALSE) %>%
    left_join(super_region) |>
    dplyr::select(
      year = year_id,
      location_id,
      location_name,
      draw,
      value
    ) |>
    pivot_wider(names_from = draw, values_from = value, names_prefix = "draw_")
}

#' Get super-region mapping table (location_id, location_name)
#'
#' @return Tibble with super-region location IDs and names
#' @export
get_super_region_mapping <- function() {
  tibble(
    location_id = c(4, 31, 64, 103, 137, 158, 166),
    location_name = c(
      "Southeast Asia, East Asia, and Oceania",
      "Central Europe, Eastern Europe, and Central Asia",
      "High-income",
      "Latin America and Caribbean",
      "North Africa and Middle East",
      "South Asia",
      "Sub-Saharan Africa"
    )
  )
}

#' Calculate Population-Scaled Thresholds
#'
#' Automatically calculates threshold values for configurations that need to scale
#' thresholds based on population reference ratios. This is particularly useful
#' for validation configurations that should inherit default thresholds but scale
#' them according to their different population reference values.
#'
#' @param default_config Configuration object containing default threshold values
#' @param target_config Configuration object that needs scaled thresholds
#' @param threshold_names Character vector of threshold names to scale (default: c("deaths_scaled"))
#' @return Modified target_config with calculated scaled thresholds
#' @export
calculate_population_scaled_thresholds <- function(
    default_config,
    target_config,
    threshold_names = c("deaths_scaled")) {
  # Calculate population scaling factor
  population_ratio <- target_config$pop_reference / default_config$pop_reference

  # Scale each specified threshold
  for (threshold_name in threshold_names) {
    if (!is.null(default_config$thresholds[[threshold_name]])) {
      target_config$thresholds[[threshold_name]] <-
        default_config$thresholds[[threshold_name]] * population_ratio
    }
  }

  return(target_config)
}

#' Initialize analysis environment
#'
#' Sets global options, seeds, and returns loaded configuration.
#' Centralizes configuration propagation so downstream code does not call getOption directly.
#'
#' @param config_file Path to config.yml
#' @param environment Config environment (default, validation, sensitivity, diagnostics)
#' @param set_seed Logical; if TRUE applies base seed
#' @return Configuration list (invisibly) for assignment
#' @export
initialize_analysis <- function(config_file = "config.yml", environment = Sys.getenv("R_CONFIG_ACTIVE", "default"), set_seed = TRUE) {
  conf <- load_config(config_file = config_file, environment = environment)
  # Propagate frequently used scalar values to options namespace (single source of truth: config)
  opts <- list(
    "pandemic.tables.projection_years" = conf$tables$projection_years %||% conf$time$future_years,
    "pandemic.draws.n" = conf$draws,
    "pandemic.seed" = conf$seed,
    "pandemic.seeds.secondary" = conf$seeds$secondary %||% NA,
    "pandemic.plots.log_breaks_major" = conf$plots$log_breaks_major %||% c(1, 10, 100, 1000)
  )
  do.call(options, opts)
  if (set_seed && !is.null(conf$seed)) {
    set.seed(conf$seed)
  }
  invisible(conf)
}
