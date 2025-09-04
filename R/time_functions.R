#' Time-based event rate analysis functions
#'
#' Functions to analyze how pandemic event rates vary over time
#' using various time windows and diagnostic methods.

library(tidyverse)
library(ggplot2)
library(lubridate)
library(patchwork)
library(MASS)
library(stats)

## NOTE (Phase 5): plot_rate_extremes moved verbatim to time_rate_extremes.R
## Keeping removal here to reduce duplication.

#' Find highest and lowest rate periods
#'
#' Analyzes exceedance rates over time using sliding windows to identify
#' periods with highest and lowest pandemic risk. Includes advanced uncertainty
#' quantification through draw-based simulation.
#'
#' @param data The pandemic dataset
#' @param variable Variable name to analyze (deaths, deaths_scaled, etc.)
#' @param threshold Threshold value to define exceedance
#' @param start_year Starting year for the analysis (default: minimum in dataset)
#' @param end_year Ending year for the analysis (default: maximum in dataset)
#' @param window_size Size of sliding window in years (default: 75)
#' @param return_all If TRUE, returns all window data; if FALSE, returns only summary (default: FALSE)
#' @param generate_draws Whether to generate draws for uncertainty quantification (default: FALSE)
#' @param n_draws Number of draws to generate for uncertainty quantification (default: 1000)
#' @param loess_span Span parameter for loess smoothing when generate_draws=TRUE (default: 0.2)
#'
#' @details
#' This function performs a sliding window analysis to identify periods with highest and lowest
#' pandemic exceedance rates. When generate_draws=TRUE, it implements a robust uncertainty
#' quantification approach using the following process:
#'
#' 1. Fits a loess smoothing model to the observed rates to capture the underlying trend
#' 2. Estimates the residual variance from the model fit
#' 3. Generates n_draws simulations by sampling from a normal distribution centered at
#'    the fitted values with variance determined from the residuals
#' 4. For each simulation, identifies the window with maximum and minimum rates
#' 5. Summarizes the distribution of highest and lowest rates across all draws
## NOTE (Phase 5): find_rate_extremes moved; leftover internal code removed.

#' Analyze pandemic exceedance rates with time trend
#'
#' Complete time series analysis framework for pandemic exceedance rates with
#' multiple model types and uncertainty quantification through Monte Carlo simulation.
#'
#' @param data The pandemic dataset
#' @param variable Variable name to analyze (deaths, deaths_scaled, etc.)
#' @param lower_cutoff Threshold value to define exceedance
#' @param start_year Starting year for the analysis (default: minimum in dataset)
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
#' @details
#' This function serves as the main entry point for time series analysis of pandemic occurrences,
#' providing a complete framework for modeling time-varying pandemic rates with rigorous
#' uncertainty quantification. It integrates several sub-components:
#'
#' 1. **Window-based rate analysis**: Calculates exceedance rates over sliding time windows
#'    with automatic window size optimization
#'
#' 2. **Multiple model types**: Supports various temporal model options:
#'    - Constant rate (traditional Poisson process assumption)
#'    - Linear trend (for detecting gradual changes in pandemic frequency)
#'    - Quadratic trend (for capturing non-linear changes)
#'    - Negative binomial variants (for handling overdispersion)
#'    - Automatic model selection based on AIC/BIC criteria
#'
#' 3. **Uncertainty quantification**: When generate_draws=TRUE, implements:
#'    - Parameter uncertainty (distribution of model coefficients)
#'    - Process uncertainty (inherent randomness in event occurrence)
#'    - Full Monte Carlo simulation for uncertainty propagation
#'    - Visualization of uncertainty through prediction intervals
#'
#' 4. **Draw-based forecasting**: Generates draws representing possible future trajectories
#'    that can be used in downstream analyses for proper uncertainty propagation
#'
#' @return List containing:
#' \item{exceedance_rate_df}{Combined data frame with historical and forecast yearly rates}
#' \item{plots}{List of visualizations if with_plots=TRUE}
#' \item{forecast_results}{Results from the forecasting procedure}
#' \item{exceedance_results}{Results from the exceedance rate analysis}
#' \item{model_comparison}{Comparison of different temporal models}
#' \item{draws}{If generate_draws=TRUE, a list containing:
#'   \item{annual_rates}{Matrix of rate draws for each year}
#'   \item{years}{Vector of years corresponding to the rows in annual_rates}
#'   \item{parameter_draws}{Matrix of parameter draws from the selected model}
#'   \item{draw_summary}{Summary statistics for the draws}
#' }

#' Analyze exceedance rates over time
#'
#' @param data The pandemic dataset
#' @param variable Variable name to analyze (deaths, deaths_scaled, etc.)
#' @param threshold Threshold value to define exceedance
#' @param base_year Starting year for analysis
#' @param window_sizes Vector of window sizes in years to try
#' @param plot_diagnostics Whether to generate diagnostic plots
#' @param event_only_windows If TRUE, uses non-overlapping windows of the specified size,
#' creating exactly one observation per time interval; dramatically reduces computation
#' @param conf Configuration object containing time.current_year for forecasting reference
#'
#' @return List containing analysis results and diagnostics

#' Fit and compare multiple temporal models
#'
#' @param data Pandemic data
#' @param variable Variable to analyze
#' @param threshold Exceedance threshold
#' @param window_size Window size in years
#' @param base_year Starting year for analysis
#' @param event_only_windows If TRUE, uses non-overlapping windows of the specified size,
#' creating exactly one observation per time interval; dramatically reduces computation
#' @param conf Configuration object containing time.current_year for forecasting reference
#'
#' @return List of model comparison results

#' Prepare forecast data with model variables
#'
#' @param forecast_years_seq Sequence of years for forecasting (these are window end years)
#' @param window_counts Historical window data
#' @param window_size Size of the analysis window
#'
#' @return Data frame ready for predictions

#' Create prediction intervals using simulation
#'
#' @param forecast_data Forecast data frame
#' @param model Statistical model
#' @param window_size Window size for rate calculation
#' @param n_sims Number of simulations
#'
#' @return Enhanced forecast data with intervals

#' Process historical data for visualization
#'
#' @param model Statistical model
#' @param window_counts Historical window counts
#' @param window_size Analysis window size
#' @param forecast_data Forecast data
#'
#' @return Data for visualization

#' Create continuous prediction line across all time periods
#'
#' @param history_fit Historical fitted values
#' @param forecast_data Forecast data
#' @param window_counts Window counts data
#' @param gap_years Years to bridge between history and forecast
#' @param model Statistical model
#' @param window_size Window size
#'
#' @return Data frame with continuous line

#' Create forecast plot
#'
#' @param combined_data Combined historical and forecast data
#' @param continuous_line Continuous prediction line
#' @param forecast_data Forecast data with intervals
#' @param last_year Last historical year
#' @param current_year Current year reference
#' @param forecast_years Number of years forecasted
#' @param prediction_interval Whether intervals are included
#' @param window_size The size of the analysis window
#'
#' @return ggplot visualization

#' Calculate forecast statistics
#'
#' @param forecast_data Forecast data
#'
#' @return List with summary statistics

#' Create forecast of future exceedance rates using best model
#'
#' @param model_results Results from compare_temporal_models
#' @param forecast_years Number of years to forecast
#' @param window_size Window size used in original analysis
#' @param prediction_interval Whether to include prediction intervals
#' @param n_sims Number of simulations for intervals
#'
#' @return List with forecast results and visualization
