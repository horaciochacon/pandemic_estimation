# Configuration Parameter Reference

This document lists all available parameters in `config.yml` for customizing the pandemic estimation analysis.

## Data & Input Parameters

`file_path` (string) - Path to input pandemic dataset CSV file  
`pop_reference` (number) - Reference population size for scaling historical deaths  
`pop_min` (number) - Minimum population threshold value  
`base_year` (number) - Starting year for historical analysis  
`variables` (array) - Primary variables to analyze [deaths_scaled]  
`deaths_unit` (number) - Unit multiplier for death counts  

## Analysis Method Parameters

`method` (string) - Extreme value statistical method [mev, evir]  
`bulk_distribution` (string) - Distribution for bulk component [lognormal, weibull, gamma]  
`use_constrained_model` (boolean) - Enable constrained mixture model fitting  
`respiratory_only` (boolean) - Analyze only respiratory pandemics  

## Threshold Parameters

`thresholds.deaths_scaled` (number) - Population-scaled death threshold for tail distribution  
`thresholds.lower_cutoff_deaths_scaled` (number) - Lower cutoff for mixture model  
`thresholds.lower_cutoff` (number) - Legacy lower cutoff parameter  
`bulk_diagnostics.tail_threshold_fixed` (number) - Fixed tail threshold for diagnostics  
`bulk_diagnostics.lower_thresholds` (array) - Candidate lower cutoffs for testing  

## Bootstrap & Uncertainty Parameters

`draws` (number) - Number of bootstrap draws for uncertainty quantification  
`seed` (number) - Random seed for reproducibility  
`seeds.secondary` (number) - Secondary seed for rate/trend calculations  
`seeds.util` (number) - Utility seed for upload preparation  
`bootstrap.enabled` (boolean) - Enable bootstrap analysis  
`bootstrap.plot` (boolean) - Generate bootstrap diagnostic plots  
`bootstrap.plot_xi` (boolean) - Generate xi parameter plots  

## Model Fitting Parameters

`constrained_optimization.method` (string) - Optimization method (penalty-based constraints)  
`constrained_optimization.verbose` (boolean) - Enable optimization debugging output  
`constrained_optimization.penalty_weights.lognormal.lower` (number) - Lower bound penalty weight  
`constrained_optimization.penalty_weights.lognormal.upper` (number) - Upper bound penalty weight  
`constrained_optimization.penalty_weights.weibull.lower` (number) - Weibull lower penalty weight  
`constrained_optimization.penalty_weights.weibull.upper` (number) - Weibull upper penalty weight  
`constrained_optimization.max_iter` (number) - Maximum optimization iterations  
`constrained_optimization.tolerance` (number) - Convergence tolerance  
`constrained_optimization.adaptive_penalty` (boolean) - Adaptive penalty weight adjustment  

## Time Analysis Parameters

`time.future_years` (number) - Years to forecast into future  
`time.prediction_interval` (boolean) - Generate prediction intervals  
`time.use_time_trend` (boolean) - Include time trend in forecasts  
`time.model_type` (string) - Time trend model type [constant, linear, quadratic, negative_binomial]  
`time.current_year` (number) - Reference year for forecasts  
`time.generate_draws` (boolean) - Generate prediction draws  
`window.sizes` (array) - Window sizes to test for rate estimation  
`window.default_size` (number) - Default window size for analysis  

## Burden Calculation Parameters

`calc_burden` (boolean) - Enable burden calculation  
`calc_dalys` (boolean) - Enable DALY estimation  
`calc_scenarios` (boolean) - Enable best/worst scenario calculations  
`scenario_window` (number) - Window size for extreme rate periods  

## Plotting Parameters

`plots.descriptive` (boolean) - Generate descriptive plots  
`plots.tail` (boolean) - Generate tail distribution plots  
`plots.tail_mixture` (boolean) - Generate mixture model plots  
`plots.tail_mixture_y_title` (string) - Y-axis title for mixture plots  
`plots.bins` (number) - Histogram bin count  
`plots.log_scale` (boolean) - Use logarithmic scale  
`plots.forecast` (boolean) - Generate forecast plots  
`plots.tail_limit` (number) - Lower limit for tail plots  
`plots.time_trend` (boolean) - Plot time trends  
`plots.time_scaled_deaths` (boolean) - Generate time-scaled death plots  
`plots.time_scaled_deaths_labels` (boolean) - Include labels in time plots  
`plots.validation_plot_labels` (boolean) - Include labels in validation plots  
`plots.draw_plot_labels` (boolean) - Include labels in draw plots  
`plots.rate_extremes` (boolean) - Plot rate extreme periods  
`plots.scenario_comparison` (boolean) - Generate scenario comparison plots  
`plots.scenario_rate_draws` (boolean) - Generate scenario rate draw plots  
`plots.log_breaks_major` (array) - Major breaks for log-scaled axes  
`plots.point_size` (number) - Size of empirical data points  
`plots.severity_marker_height` (number) - Height of severity vertical markers  

## Output Management Parameters

`output.save_plots` (boolean) - Save plots to files  
`output.save_tables` (boolean) - Save tables to files  
`output.base_directory` (string) - Base output directory path  
`output.use_timestamps` (boolean) - Use timestamped subdirectories  
`output.plot_formats` (array) - Plot output formats [svg, png, pdf, eps]  

## Parallel Processing Parameters

`parallel.enabled` (boolean) - Enable parallel processing  
`parallel.strategy` (string) - Parallelization strategy [future, foreach, none]  
`parallel.workers` (number) - Number of parallel workers  

## Tables Parameters

`tables.projection_years` (number) - Default projection period for summary tables  

## Validation-Specific Parameters

`train_year_start` (number) - Training period start year  
`train_year_end` (number) - Training period end year  
`test_year_start` (number) - Testing period start year  
`test_year_end` (number) - Testing period end year  
`window_size` (number) - Window size for validation  
`threshold_deaths_scaled` (number) - Scaled death threshold for validation  
`observation_years` (array) - Years for observation comparison  

## Sensitivity Analysis Parameters

`sens.year.start_idx` (number) - Start index for year range sensitivity  
`sens.year.end_idx` (number) - End index for year range sensitivity  
`sens.year.cutoff` (number) - Year cutoff value  
`sens.year.base_year` (number) - Base year for analysis  
`sens.year.base_thresholds.deaths_scaled` (number) - Base threshold for year sensitivity  
`sens.threshold.quantile_limit` (number) - Upper quantile limit  
`sens.threshold.skip_first` (number) - Skip initial observations  
`sens.threshold.base_thresholds.deaths_scaled` (number) - Base threshold for sensitivity  
`sens.combined.enabled` (boolean) - Enable combined sensitivity analysis  
`sens.combined.max_combinations` (number) - Maximum parameter combinations  
`sens.combined.plots.heatmap` (boolean) - Generate heatmap visualization  

## Diagnostic Parameters

`diagnostics.threshold_selection.n_thresholds` (number) - Number of thresholds to test  
`diagnostics.threshold_selection.generate_plots` (boolean) - Generate threshold plots  
`diagnostics.threshold_selection.min_exceedances` (number) - Minimum exceedances for fitting  
`diagnostics.goodness_of_fit.generate_plots` (boolean) - Generate GOF plots  
`diagnostics.goodness_of_fit.statistical_tests` (boolean) - Run statistical tests  
`diagnostics.mixture_components.generate_plots` (boolean) - Generate component plots  
`diagnostics.mixture_components.weight_stability_test` (boolean) - Test weight stability  
`diagnostics.bootstrap_stability.generate_plots` (boolean) - Generate bootstrap plots  
`diagnostics.bootstrap_stability.convergence_threshold` (number) - Convergence rate threshold  
`diagnostics.cross_validation.method` (string) - Cross-validation method [leave_one_out, k_fold]  
`diagnostics.cross_validation.generate_plots` (boolean) - Generate CV plots  
`diagnostics.cross_validation.k_folds` (number) - Number of folds for k-fold CV  
`diagnostics.publication_tables.generate_summary` (boolean) - Generate summary table  
`diagnostics.publication_tables.generate_threshold_table` (boolean) - Generate threshold table  
`diagnostics.publication_tables.format_style` (string) - Table formatting style

## Configuration Inheritance

Configuration sections can inherit from other sections using the `inherit` parameter:
- `inherit: default` - Inherit all parameters from the default section and override specific values
- Available in validation, sensitivity, and diagnostics configurations