####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Bulk component diagnostics script for mixture model analysis of pandemic data.
##              This script generates detailed diagnostics for the truncated bulk distribution
##              component, testing different lower thresholds and distribution types (lognormal,
##              weibull, gamma) while keeping the tail threshold fixed at 2.8M deaths.
##
## Generates:   - Distribution comparison table (AIC, BIC, Log-Likelihood)
##              - Threshold sensitivity analysis
##              - Continuity constraint verification plots
##              - Goodness-of-fit diagnostics (Q-Q plots, K-S tests)
##              - Parameter stability across thresholds
##              - Publication-ready tables and figures
##
## Requires:    - Configuration file ("config.yml") with analysis parameters
##              - All R functions loaded from R/ directory
##              - Packages: tidyverse, mev, evir, gridExtra, patchwork
##
## Outputs:     - Diagnostic plots saved to output/
##              - Results saved to output/bulk_diagnostics_results.rds
##              - Console output with formatted tables
##
####################################################################################################

# Load project functions & packages explicitly
source("R/load_all.R")

# ===================================================================================
# Setup environment and configuration
# ===================================================================================
# Set active configuration
Sys.setenv(R_CONFIG_ACTIVE = "default")
conf <- load_config()

# Define test parameters (prefer config; fall back to historical literals to preserve behavior)
lower_thresholds <- conf$bulk_diagnostics$lower_thresholds %||% c(50000, 75000, 100000, 125000, 150000)
bulk_distributions <- c("lognormal", "weibull") # Distribution types to compare
tail_threshold <- conf$bulk_diagnostics$tail_threshold_fixed %||% 2.8e6
base_year <- conf$base_year %||% 1600

# ===================================================================================
# Data preparation - identical to run_PoT.R
# ===================================================================================
cat("===================================================================================\n")
cat("BULK COMPONENT DIAGNOSTICS FOR MIXTURE MODEL ANALYSIS\n")
cat("===================================================================================\n")
cat("Configuration:\n")
cat("  Dataset:", basename(conf$file_path), "\n")
cat("  Base year:", base_year, "\n")
cat("  Tail threshold:", scales::comma(tail_threshold), "deaths\n")
cat("  Lower thresholds to test:", paste(scales::comma(lower_thresholds), collapse = ", "), "\n")
cat("  Distributions to test:", paste(bulk_distributions, collapse = ", "), "\n")
cat("===================================================================================\n\n")

# Read and transform the pandemic data
data <- read_transform_data(conf)

# Filter data using base year
data <- data %>%
  filter(
    start_year >= base_year,
    deaths > conf$pop_min
  )

# Apply dual transformation
data$dual <- dual_transform(data, "deaths_scaled", h = conf$pop_reference, l = conf$pop_min)

# ===================================================================================
# Initialize results storage
# ===================================================================================
results_list <- list()
comparison_df <- data.frame()
continuity_results <- list()

# ===================================================================================
# Main analysis loop: test each combination of threshold and distribution
# ===================================================================================
cat("1. FITTING BULK COMPONENTS WITH DIFFERENT CONFIGURATIONS\n")
cat("-------------------------------------------------------------------\n")

# Determine threshold
u_best <- determine_threshold(data, conf$variables, NULL, conf)

# Create initial summary of data characteristics and thresholds
specs <- tibble(
  Variable = conf$variables,
  Year = conf$base_year,
  Threshold = inv_dual_transform(u_best, h = conf$pop_reference, l = conf$pop_min),
  Cutoff = conf$thresholds$lower_cutoff,
  n = nrow(data),
  PoT = sum(data$dual > u_best),
  PcoT = mean(data$dual > u_best),
  PoB = sum(data$dual > conf$thresholds$lower_cutoff)
)



# First fit the tail component to get exceedance probabilities
tail_results <- fit_gpd_model(
  x = data$dual,
  u_best = tail_threshold,
  h = conf$pop_reference,
  summary = specs,
  method = conf$method
)


for (lower_cutoff in lower_thresholds) {
  # Create initial summary of data characteristics and thresholds
  specs <- tibble(
    Variable = conf$variables,
    Year = conf$base_year,
    Threshold = inv_dual_transform(u_best, h = conf$pop_reference, l = conf$pop_min),
    Cutoff = lower_cutoff,
    n = nrow(data),
    PoT = sum(data$dual > u_best),
    PcoT = mean(data$dual > u_best),
    PoB = sum(data$dual > conf$thresholds$lower_cutoff)
  )

  # First fit the tail component to get exceedance probabilities
  tail_results <- fit_gpd_model(
    x = data$dual,
    u_best = tail_threshold,
    h = conf$pop_reference,
    summary = specs,
    method = conf$method
  )

  for (dist_type in bulk_distributions) {
    cat(sprintf(
      "\nFitting %s distribution with lower cutoff = %s deaths...\n",
      dist_type, scales::comma(lower_cutoff)
    ))

    # Fit the bulk component
    bulk_fit <- tryCatch(
      {
        fit_bulk_component(
          data = data,
          variable = "deaths_scaled",
          lower_cutoff = lower_cutoff,
          upper_cutoff = tail_threshold,
          dist_type = dist_type,
          constrain_survival = TRUE,
          tail_fit = tail_results,
          conf = conf
        )
      },
      error = function(e) {
        cat("  Error fitting:", e$message, "\n")
        return(NULL)
      }
    )

    if (!is.null(bulk_fit)) {
      # Store results
      result_key <- paste0(dist_type, "_", lower_cutoff)
      results_list[[result_key]] <- bulk_fit

      # Extract fit statistics
      n_bulk <- sum(data$deaths_scaled > lower_cutoff & data$deaths_scaled <= tail_threshold)

      # Calculate AIC and BIC
      n_params <- length(bulk_fit$parameters)
      aic <- 2 * n_params - 2 * bulk_fit$fit$loglik
      bic <- n_params * log(n_bulk) - 2 * bulk_fit$fit$loglik

      # Test continuity at boundaries
      # Lower boundary continuity
      emp_prob_lower <- mean(data$deaths_scaled > lower_cutoff)
      fitted_prob_lower <- bulk_fit$survival_lower
      lower_continuity_error <- abs(fitted_prob_lower - emp_prob_lower) / emp_prob_lower

      # Upper boundary continuity (would need tail fit for exact comparison)
      upper_continuity_ratio <- bulk_fit$survival_thresh / bulk_fit$survival_lower

      # Append to comparison dataframe
      comparison_df <- rbind(comparison_df, data.frame(
        Distribution = dist_type,
        Lower_Cutoff = lower_cutoff,
        N_Events = n_bulk,
        Log_Likelihood = bulk_fit$fit$loglik,
        AIC = aic,
        BIC = bic,
        Parameters = paste(names(bulk_fit$parameters), "=",
          round(unlist(bulk_fit$parameters), 3),
          collapse = ", "
        ),
        Conditional_Mean = bulk_fit$conditional_mean,
        Lower_Continuity_Error = lower_continuity_error,
        Upper_Survival_Ratio = upper_continuity_ratio
      ))

      cat(sprintf(
        "  Success: LL=%.1f, AIC=%.1f, BIC=%.1f, Mean=%.0f deaths\n",
        bulk_fit$fit$loglik, aic, bic, bulk_fit$conditional_mean
      ))
    }
  }
}

# ===================================================================================
# Display comparison table
# ===================================================================================
cat("\n2. MODEL COMPARISON TABLE\n")
cat("-------------------------------------------------------------------\n")

# Sort by AIC
comparison_df <- comparison_df %>%
  arrange(AIC)

# Create formatted table for display
formatted_comparison <- comparison_df %>%
  mutate(
    Lower_Cutoff = scales::comma(Lower_Cutoff),
    Log_Likelihood = round(Log_Likelihood, 1),
    AIC = round(AIC, 1),
    BIC = round(BIC, 1),
    Conditional_Mean = scales::comma(round(Conditional_Mean, 0)),
    Lower_Continuity_Error = sprintf("%.1f%%", Lower_Continuity_Error * 100)
  ) %>%
  dplyr::select(
    Distribution, Lower_Cutoff, N_Events, AIC, BIC, Log_Likelihood,
    Conditional_Mean, Lower_Continuity_Error
  )

print(formatted_comparison, row.names = FALSE)

# Find best model
best_idx <- which.min(comparison_df$AIC)
best_config <- comparison_df[best_idx, ]
cat("\nBest configuration (by AIC):\n")
cat("  Distribution:", best_config$Distribution, "\n")
cat("  Lower cutoff:", scales::comma(best_config$Lower_Cutoff), "deaths\n")
cat("  AIC:", round(best_config$AIC, 1), "\n")

# ===================================================================================
# Generate diagnostic plots
# ===================================================================================
cat("\n3. GENERATING DIAGNOSTIC PLOTS\n")
cat("-------------------------------------------------------------------\n")

# Create multi-panel plot showing fit quality for different configurations
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2))

# Plot 1: AIC vs Lower Threshold for each distribution
plot_data_aic <- comparison_df %>%
  group_by(Distribution) %>%
  arrange(Lower_Cutoff)

plot(1,
  type = "n", xlim = range(lower_thresholds),
  ylim = range(comparison_df$AIC, na.rm = TRUE),
  xlab = "Lower Cutoff (deaths)", ylab = "AIC",
  main = "Model AIC vs Lower Threshold", axes = FALSE
)

# Add custom axes with human-readable labels
axis(1, at = lower_thresholds, labels = scales::comma(lower_thresholds))
axis(2)
box()

# Add lines for each distribution with distinct colors
dist_colors <- c("lognormal" = "#00468B", "weibull" = "#ED0000", "gamma" = "#42B540")
for (dist in bulk_distributions) {
  dist_data <- plot_data_aic %>% filter(Distribution == dist)
  lines(dist_data$Lower_Cutoff, dist_data$AIC,
    col = dist_colors[dist], lwd = 2, type = "b", pch = 19
  )
}

# Add legend
legend("topright",
  legend = bulk_distributions, col = dist_colors[bulk_distributions],
  lwd = 2, pch = 19, bty = "n"
)

# Add grid
grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Plot 2: Number of Events vs Lower Threshold
plot(lower_thresholds, comparison_df %>%
  filter(Distribution == "lognormal") %>%
  arrange(Lower_Cutoff) %>%
  pull(N_Events),
type = "b", pch = 19, col = "#00468B", lwd = 2,
xlab = "Lower Cutoff (deaths)", ylab = "Number of Events in Bulk",
main = "Sample Size in Bulk Region", axes = FALSE
)

axis(1, at = lower_thresholds, labels = scales::comma(lower_thresholds))
axis(2)
box()
grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Plot 3: Conditional Mean vs Lower Threshold
plot_data_mean <- comparison_df %>%
  group_by(Distribution) %>%
  arrange(Lower_Cutoff)

plot(1,
  type = "n", xlim = range(lower_thresholds),
  ylim = range(comparison_df$Conditional_Mean, na.rm = TRUE),
  xlab = "Lower Cutoff (deaths)", ylab = "Conditional Mean (deaths)",
  main = "Bulk Conditional Mean vs Lower Threshold", axes = TRUE
)

axis(1, at = lower_thresholds, labels = scales::comma(lower_thresholds))
box()

for (dist in bulk_distributions) {
  dist_data <- plot_data_mean %>% filter(Distribution == dist)
  lines(dist_data$Lower_Cutoff, dist_data$Conditional_Mean,
    col = dist_colors[dist], lwd = 2, type = "b", pch = 19
  )
}

grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Plot 4: Lower Boundary Continuity Error
plot_data_cont <- comparison_df %>%
  group_by(Distribution) %>%
  arrange(Lower_Cutoff) %>%
  mutate(Continuity_Error_Pct = Lower_Continuity_Error * 100)

plot(1,
  type = "n", xlim = range(lower_thresholds),
  ylim = range(plot_data_cont$Continuity_Error_Pct, na.rm = TRUE),
  xlab = "Lower Cutoff (deaths)", ylab = "Continuity Error (%)",
  main = "Lower Boundary Continuity Error", axes = FALSE
)

axis(1, at = lower_thresholds, labels = scales::comma(lower_thresholds))
axis(2)
box()

for (dist in bulk_distributions) {
  dist_data <- plot_data_cont %>% filter(Distribution == dist)
  lines(dist_data$Lower_Cutoff, dist_data$Continuity_Error_Pct,
    col = dist_colors[dist], lwd = 2, type = "b", pch = 19
  )
}

# Add reference line at 5% error
abline(h = 5, lty = 2, col = "gray50", lwd = 1)
text(max(lower_thresholds), 5, "5% threshold", adj = c(1, -0.5), col = "gray50", cex = 0.8)

grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Plot 5: Q-Q plot for best configuration
best_key <- paste0(best_config$Distribution, "_", best_config$Lower_Cutoff)
best_fit <- results_list[[best_key]]

# Extract bulk data
bulk_data <- data$deaths_scaled[data$deaths_scaled > best_config$Lower_Cutoff &
  data$deaths_scaled <= tail_threshold]
n_bulk <- length(bulk_data)

# Debug: Print parameter structure
cat("Debug - Best fit parameters structure:\n")
print(str(best_fit$parameters))
cat("Distribution:", best_config$Distribution, "\n")

# Generate theoretical quantiles for the fitted distribution
p_seq <- ppoints(n_bulk)
dist_func_name <- switch(best_config$Distribution,
  "lognormal" = "lnorm",
  "weibull" = "weibull"
)

q_fn <- match.fun(paste0("q", dist_func_name))

# Map parameter names to R distribution function expectations
# Handle different parameter storage formats
params <- best_fit$parameters
if (is.list(params)) {
  # Parameters stored as named list
  if (best_config$Distribution == "lognormal") {
    theoretical_q <- q_fn(p_seq,
      meanlog = as.numeric(params$meanlog),
      sdlog = as.numeric(params$sdlog)
    )
  } else if (best_config$Distribution == "weibull") {
    theoretical_q <- q_fn(p_seq,
      shape = as.numeric(params$shape),
      scale = as.numeric(params$scale)
    )
  } else if (best_config$Distribution == "gamma") {
    theoretical_q <- q_fn(p_seq,
      shape = as.numeric(params$shape),
      rate = as.numeric(params$rate)
    )
  }
} else {
  # Parameters might be stored as named vector (from fitdistr)
  if (best_config$Distribution == "lognormal") {
    theoretical_q <- q_fn(p_seq,
      meanlog = as.numeric(params["meanlog"]),
      sdlog = as.numeric(params["sdlog"])
    )
  } else if (best_config$Distribution == "weibull") {
    theoretical_q <- q_fn(p_seq,
      shape = as.numeric(params["shape"]),
      scale = as.numeric(params["scale"])
    )
  } else if (best_config$Distribution == "gamma") {
    theoretical_q <- q_fn(p_seq,
      shape = as.numeric(params["shape"]),
      rate = as.numeric(params["rate"])
    )
  }
}

# Create Q-Q plot
qqplot(theoretical_q, sort(bulk_data),
  xlab = paste("Theoretical", best_config$Distribution, "quantiles"),
  ylab = "Empirical quantiles",
  main = paste(
    "Q-Q Plot: Best Model\n(", best_config$Distribution,
    ", cutoff =", scales::comma(best_config$Lower_Cutoff), ")"
  ),
  pch = 19, col = "#00468B", cex = 0.8
)

# Add 1:1 line
abline(0, 1, col = "#ED0000", lwd = 2)
grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Plot 6: Empirical vs Fitted CDF for best model
# Create empirical CDF
ecdf_bulk <- ecdf(bulk_data)
x_range <- seq(min(bulk_data), max(bulk_data), length.out = 500)
empirical_cdf <- ecdf_bulk(x_range)

# Calculate fitted CDF
p_fn <- match.fun(paste0("p", dist_func_name))

# Map parameter names for CDF calculation using same robust approach
params <- best_fit$parameters
if (is.list(params)) {
  # Parameters stored as named list
  if (best_config$Distribution == "lognormal") {
    fitted_cdf_raw <- p_fn(x_range,
      meanlog = as.numeric(params$meanlog),
      sdlog = as.numeric(params$sdlog)
    )
    p_lower <- p_fn(best_config$Lower_Cutoff,
      meanlog = as.numeric(params$meanlog),
      sdlog = as.numeric(params$sdlog)
    )
    p_upper <- p_fn(tail_threshold,
      meanlog = as.numeric(params$meanlog),
      sdlog = as.numeric(params$sdlog)
    )
  } else if (best_config$Distribution == "weibull") {
    fitted_cdf_raw <- p_fn(x_range,
      shape = as.numeric(params$shape),
      scale = as.numeric(params$scale)
    )
    p_lower <- p_fn(best_config$Lower_Cutoff,
      shape = as.numeric(params$shape),
      scale = as.numeric(params$scale)
    )
    p_upper <- p_fn(tail_threshold,
      shape = as.numeric(params$shape),
      scale = as.numeric(params$scale)
    )
  } else if (best_config$Distribution == "gamma") {
    fitted_cdf_raw <- p_fn(x_range,
      shape = as.numeric(params$shape),
      rate = as.numeric(params$rate)
    )
    p_lower <- p_fn(best_config$Lower_Cutoff,
      shape = as.numeric(params$shape),
      rate = as.numeric(params$rate)
    )
    p_upper <- p_fn(tail_threshold,
      shape = as.numeric(params$shape),
      rate = as.numeric(params$rate)
    )
  }
} else {
  # Parameters stored as named vector (from fitdistr)
  if (best_config$Distribution == "lognormal") {
    fitted_cdf_raw <- p_fn(x_range,
      meanlog = as.numeric(params["meanlog"]),
      sdlog = as.numeric(params["sdlog"])
    )
    p_lower <- p_fn(best_config$Lower_Cutoff,
      meanlog = as.numeric(params["meanlog"]),
      sdlog = as.numeric(params["sdlog"])
    )
    p_upper <- p_fn(tail_threshold,
      meanlog = as.numeric(params["meanlog"]),
      sdlog = as.numeric(params["sdlog"])
    )
  } else if (best_config$Distribution == "weibull") {
    fitted_cdf_raw <- p_fn(x_range,
      shape = as.numeric(params["shape"]),
      scale = as.numeric(params["scale"])
    )
    p_lower <- p_fn(best_config$Lower_Cutoff,
      shape = as.numeric(params["shape"]),
      scale = as.numeric(params["scale"])
    )
    p_upper <- p_fn(tail_threshold,
      shape = as.numeric(params["shape"]),
      scale = as.numeric(params["scale"])
    )
  } else if (best_config$Distribution == "gamma") {
    fitted_cdf_raw <- p_fn(x_range,
      shape = as.numeric(params["shape"]),
      rate = as.numeric(params["rate"])
    )
    p_lower <- p_fn(best_config$Lower_Cutoff,
      shape = as.numeric(params["shape"]),
      rate = as.numeric(params["rate"])
    )
    p_upper <- p_fn(tail_threshold,
      shape = as.numeric(params["shape"]),
      rate = as.numeric(params["rate"])
    )
  }
}

# Normalize to account for truncation
fitted_cdf <- (fitted_cdf_raw - p_lower) / (p_upper - p_lower)

plot(x_range, empirical_cdf,
  type = "l", lwd = 2, col = "#00468B",
  xlab = "Population-scaled deaths", ylab = "Cumulative probability",
  main = "Empirical vs Fitted CDF (Best Model)",
  log = "x", axes = TRUE
)

box()

lines(x_range, fitted_cdf, col = "#ED0000", lwd = 2, lty = 2)

legend("bottomright",
  legend = c("Empirical", "Fitted"),
  col = c("#00468B", "#ED0000"), lwd = 2, lty = c(1, 2), bty = "n"
)

grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Reset plot parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ===================================================================================
# Goodness-of-fit tests
# ===================================================================================
cat("\n4. GOODNESS-OF-FIT TESTS FOR BEST MODEL\n")
cat("-------------------------------------------------------------------\n")

# Kolmogorov-Smirnov test
# Helper to extract numeric parameter irrespective of list or named vector storage
extract_param <- function(prms, nm) {
  if (is.list(prms)) as.numeric(prms[[nm]]) else as.numeric(prms[nm])
}

# Single factory for fitted CDF (replaces duplicated blocks above; logic unchanged)
build_fitted_cdf_fn <- function(dist_type, params, p_fn) {
  switch(dist_type,
    lognormal = function(x) {
      p_fn(x,
        meanlog = extract_param(params, "meanlog"),
        sdlog = extract_param(params, "sdlog")
      )
    },
    weibull = function(x) {
      p_fn(x,
        shape = extract_param(params, "shape"),
        scale = extract_param(params, "scale")
      )
    },
    gamma = function(x) {
      p_fn(x,
        shape = extract_param(params, "shape"),
        rate = extract_param(params, "rate")
      )
    },
    stop("Unsupported distribution type: ", dist_type)
  )
}

params <- best_fit$parameters
fitted_cdf_fn <- build_fitted_cdf_fn(best_config$Distribution, params, p_fn)

ks_test <- ks.test(bulk_data, fitted_cdf_fn)
cat("Kolmogorov-Smirnov test:\n")
cat("  D-statistic:", round(ks_test$statistic, 4), "\n")
cat("  p-value:", round(ks_test$p.value, 4), "\n")
cat("  Result:", ifelse(ks_test$p.value > 0.05, "Fail to reject null (good fit)",
  "Reject null (poor fit)"
), "\n\n")

# Anderson-Darling test (if available)
if (requireNamespace("goftest", quietly = TRUE)) {
  ad_test <- goftest::ad.test(bulk_data, fitted_cdf_fn)
  cat("Anderson-Darling test:\n")
  cat("  A-statistic:", round(ad_test$statistic, 4), "\n")
  cat("  p-value:", round(ad_test$p.value, 4), "\n")
  cat("  Result:", ifelse(ad_test$p.value > 0.05, "Fail to reject null (good fit)",
    "Reject null (poor fit)"
  ), "\n")
}

# ===================================================================================
# Parameter estimates table for best model
# ===================================================================================
cat("\n5. PARAMETER ESTIMATES FOR BEST MODEL\n")
cat("-------------------------------------------------------------------\n")

param_names <- names(best_fit$parameters)
param_values <- unlist(best_fit$parameters)

cat("Best model configuration:\n")
cat("  Distribution:", best_config$Distribution, "\n")
cat("  Lower cutoff:", scales::comma(best_config$Lower_Cutoff), "deaths\n")
cat("  upper threshold:", scales::comma(tail_threshold), "deaths\n")
cat("  Number of events in bulk:", best_config$N_Events, "\n\n")

cat("Parameter estimates:\n")
for (i in seq_along(param_names)) {
  cat("  ", param_names[i], ":", round(param_values[i], 4), "\n")
}

cat("\nDerived quantities:\n")
cat("  Conditional mean:", scales::comma(round(best_fit$conditional_mean, 0)), "deaths\n")
cat("  Survival at lower boundary:", round(best_fit$survival_lower, 4), "\n")
cat("  Survival at upper boundary:", round(best_fit$survival_thresh, 6), "\n")
cat("  Log-likelihood:", round(best_fit$fit$loglik, 2), "\n")

# ===================================================================================
# Create publication-ready markdown table
# ===================================================================================
cat("\n6. PUBLICATION-READY MARKDOWN TABLE\n")
cat("-------------------------------------------------------------------\n\n")

if (exists("enhanced_comparison") && nrow(enhanced_comparison) > 0) {
  cat("### Enhanced Model Comparison (Statistical Fit + Predictive Performance)\n\n")
  cat("**Metrics:**\n")
  cat("- **AIC**: Akaike Information Criterion (lower = better statistical fit)\n")
  cat("- **RMSE**: Root Mean Squared Error (lower = better prediction accuracy)\n")
  cat("- **MAE**: Mean Absolute Error (lower = better prediction accuracy)\n")
  cat("- **MAPE**: Mean Absolute Percentage Error (lower = better relative accuracy)\n")
  cat("- **Coverage**: % of observations within 95% uncertainty intervals (higher = better calibration)\n")
  cat("- **Composite Score**: Weighted combination (30% AIC + 25% RMSE + 25% MAE + 20% Coverage)\n\n")

  cat("| Distribution | Lower Cutoff | AIC | RMSE | MAE | MAPE | Coverage | Composite Score |\n")
  cat("|--------------|--------------|-----|------|-----|------|----------|------------------|\n")

  for (i in seq_len(min(10, nrow(enhanced_comparison)))) {
    row <- enhanced_comparison[i, ]
    cat(sprintf(
      "| %s | %s | %.1f | %s | %s | %s | %s | %.2f |\n",
      row$Distribution,
      scales::comma(row$Lower_Cutoff),
      row$AIC,
      ifelse(is.na(row$RMSE), "Failed", scales::comma(round(row$RMSE, 0))),
      ifelse(is.na(row$MAE), "Failed", scales::comma(round(row$MAE, 0))),
      ifelse(is.na(row$MAPE), "Failed", paste0(round(row$MAPE, 1), "%")),
      ifelse(is.na(row$Coverage_Prob), "Failed", paste0(round(row$Coverage_Prob, 1), "%")),
      row$Composite_Score
    ))
  }
} else {
  cat("### Statistical Fit Only (No Validation Results)\n\n")
  cat("| Distribution | Lower Cutoff | AIC | BIC | Log-Likelihood | Conditional Mean | Continuity Error |\n")
  cat("|--------------|--------------|-----|-----|----------------|------------------|------------------|\n")

  for (i in seq_len(min(10, nrow(comparison_df)))) {
    row <- comparison_df[i, ]
    cat(sprintf(
      "| %s | %s | %.1f | %.1f | %.1f | %s | %.1f%% |\n",
      row$Distribution,
      scales::comma(row$Lower_Cutoff),
      row$AIC,
      row$BIC,
      row$Log_Likelihood,
      scales::comma(round(row$Conditional_Mean, 0)),
      row$Lower_Continuity_Error * 100
    ))
  }
}

# ===================================================================================
# Predictive performance validation for different bulk configurations
# ===================================================================================
cat("\n\n7. PREDICTIVE PERFORMANCE VALIDATION\n")
cat("-------------------------------------------------------------------\n")

# Define test configurations - combinations of distribution and lower cutoff
test_configs <- expand.grid(
  bulk_distribution = bulk_distributions,
  lower_cutoff = lower_thresholds,
  stringsAsFactors = FALSE
)

cat("Testing predictive performance for", nrow(test_configs), "configurations...\n")
cat("Validation setup: Training 1600-1950, Testing 1950-2025\n\n")

# Initialize storage for validation results
validation_results <- list()
validation_performance <- data.frame()

# Progress tracking
cat("Running validation for each configuration:\n")

for (i in seq_len(nrow(test_configs))) {
  config_name <- paste0(test_configs$bulk_distribution[i], "_", test_configs$lower_cutoff[i])

  cat(sprintf(
    "  [%d/%d] Testing %s distribution with %s lower cutoff...",
    i, nrow(test_configs),
    test_configs$bulk_distribution[i],
    scales::comma(test_configs$lower_cutoff[i])
  ))

  # Load fresh validation configuration for each iteration to avoid accumulation
  Sys.setenv(R_CONFIG_ACTIVE = "validation")
  validation_conf <- load_config()

  # Create modified validation configuration
  validation_conf$bulk_distribution <- test_configs$bulk_distribution[i]
  validation_conf$thresholds$lower_cutoff <- test_configs$lower_cutoff[i]
  validation_conf$thresholds$lower_cutoff_deaths_scaled <- test_configs$lower_cutoff[i]

  # Ensure consistency in configuration
  validation_conf$mixture_dist <- test_configs$bulk_distribution[i]

  # Prepare data with validation configuration (important for consistent thresholds)
  validation_data <- read_transform_data(validation_conf) %>%
    filter(
      start_year >= validation_conf$base_year,
      deaths > validation_conf$pop_min
    )

  # Apply dual transformation with validation configuration
  validation_data$dual <- dual_transform(validation_data, "deaths_scaled",
    h = validation_conf$pop_reference,
    l = validation_conf$pop_min
  )

  # Run validation
  validation_result <- tryCatch(
    {
      validate_model(validation_data, "deaths_scaled", validation_conf)
    },
    error = function(e) {
      cat(" ERROR:", e$message, "\n")
      return(NULL)
    }
  )

  if (!is.null(validation_result)) {
    # Store full result
    validation_results[[config_name]] <- validation_result

    # Extract performance metrics using the cleaned forecast data
    observed_data <- validation_result$observed
    predicted_data <- validation_result$forecast # This has cum_deaths, cum_deaths_low, cum_deaths_up

    # Calculate performance metrics for cumulative burden at observation years
    if (nrow(observed_data) > 0 && nrow(predicted_data) > 0) {
      # Merge observed and predicted data by year
      comparison_data <- observed_data %>%
        inner_join(predicted_data, by = "year") %>%
        arrange(year)

      if (nrow(comparison_data) > 0) {
        # Calculate primary performance metrics
        rmse <- sqrt(mean((comparison_data$obs - comparison_data$cum_deaths)^2, na.rm = TRUE))
        mae <- mean(abs(comparison_data$obs - comparison_data$cum_deaths), na.rm = TRUE)
        bias <- mean(comparison_data$cum_deaths - comparison_data$obs, na.rm = TRUE)

        # Calculate percentage errors
        mape <- mean(abs((comparison_data$obs - comparison_data$cum_deaths) / comparison_data$obs) * 100, na.rm = TRUE)

        # Calculate coverage probability (what % of observations fall within 95% UI)
        in_interval <- comparison_data$obs >= comparison_data$cum_deaths_low &
          comparison_data$obs <= comparison_data$cum_deaths_up
        coverage_prob <- mean(in_interval, na.rm = TRUE) * 100

        # Calculate relative metrics
        rmse_relative <- rmse / mean(comparison_data$obs) * 100
        mae_relative <- mae / mean(comparison_data$obs) * 100

        # Store detailed performance metrics
        performance_row <- data.frame(
          Distribution = test_configs$bulk_distribution[i],
          Lower_Cutoff = test_configs$lower_cutoff[i],
          RMSE = rmse,
          MAE = mae,
          Bias = bias,
          MAPE = mape,
          RMSE_Relative = rmse_relative,
          MAE_Relative = mae_relative,
          Coverage_Prob = coverage_prob,
          N_Obs = nrow(comparison_data),
          stringsAsFactors = FALSE
        )

        validation_performance <- rbind(validation_performance, performance_row)

        cat(" SUCCESS\n")
        cat(sprintf(
          "    RMSE: %s, Coverage: %.1f%%\n",
          scales::comma(round(rmse, 0)), coverage_prob
        ))
      } else {
        cat(" Warning: No matching years between observed and predicted data\n")
      }
    } else {
      cat(" Warning: Missing observed or predicted data\n")
    }
  }
}

cat("\nValidation completed for", nrow(validation_performance), "configurations.\n")

# ===================================================================================
# Enhanced model comparison combining statistical fit and predictive performance
# ===================================================================================
cat("\n8. ENHANCED MODEL COMPARISON: STATISTICAL FIT + PREDICTIVE PERFORMANCE\n")
cat("-------------------------------------------------------------------\n")

if (nrow(validation_performance) > 0) {
  # Merge statistical and predictive performance results
  enhanced_comparison <- comparison_df %>%
    left_join(
      validation_performance,
      by = c("Distribution", "Lower_Cutoff")
    ) %>%
    # Calculate composite scores (lower is better, except coverage which higher is better)
    mutate(
      AIC_rank = rank(AIC),
      RMSE_rank = rank(RMSE, na.last = TRUE),
      MAE_rank = rank(MAE, na.last = TRUE),
      MAPE_rank = rank(MAPE, na.last = TRUE),
      # For coverage probability, higher is better, so rank in descending order
      Coverage_rank = rank(-Coverage_Prob, na.last = TRUE),
      # Composite score with weighted criteria
      # 30% statistical fit (AIC), 50% prediction accuracy (RMSE+MAE), 20% uncertainty calibration (Coverage)
      Composite_Score = (0.3 * AIC_rank + 0.25 * RMSE_rank + 0.25 * MAE_rank + 0.2 * Coverage_rank)
    ) %>%
    arrange(Composite_Score)

  # Display enhanced comparison table
  cat("\nEnhanced Model Comparison (ordered by composite score):\n")
  cat("Composite Score = 30% AIC + 25% RMSE + 25% MAE + 20% Coverage\n")
  cat("(Lower composite score = better overall performance)\n\n")

  enhanced_display <- enhanced_comparison %>%
    mutate(
      Lower_Cutoff = scales::comma(Lower_Cutoff),
      AIC = round(AIC, 1),
      BIC = round(BIC, 1),
      RMSE = ifelse(is.na(RMSE), "Failed", scales::comma(round(RMSE, 0))),
      MAE = ifelse(is.na(MAE), "Failed", scales::comma(round(MAE, 0))),
      MAPE = ifelse(is.na(MAPE), "Failed", paste0(round(MAPE, 1), "%")),
      Coverage_Prob = ifelse(is.na(Coverage_Prob), "Failed", paste0(round(Coverage_Prob, 1), "%")),
      Composite_Score = round(Composite_Score, 2)
    ) %>%
    dplyr::select(Distribution, Lower_Cutoff, AIC, RMSE, MAE, MAPE, Coverage_Prob, Composite_Score)

  print(enhanced_display)

  # Identify best configuration by composite score
  best_enhanced_idx <- which.min(enhanced_comparison$Composite_Score)
  best_enhanced_config <- enhanced_comparison[best_enhanced_idx, ]

  cat("\nBest configuration by composite score:\n")
  cat("  Distribution:", best_enhanced_config$Distribution, "\n")
  cat("  Lower cutoff:", scales::comma(best_enhanced_config$Lower_Cutoff), "deaths\n")
  cat("  AIC:", round(best_enhanced_config$AIC, 1), "\n")
  cat("  RMSE:", ifelse(is.na(best_enhanced_config$RMSE), "Failed",
    scales::comma(round(best_enhanced_config$RMSE, 0))
  ), "\n")
  cat("  MAE:", ifelse(is.na(best_enhanced_config$MAE), "Failed",
    scales::comma(round(best_enhanced_config$MAE, 0))
  ), "\n")
  cat("  MAPE:", ifelse(is.na(best_enhanced_config$MAPE), "Failed",
    paste0(round(best_enhanced_config$MAPE, 1), "%")
  ), "\n")
  cat("  Coverage probability:", ifelse(is.na(best_enhanced_config$Coverage_Prob), "Failed",
    paste0(round(best_enhanced_config$Coverage_Prob, 1), "%")
  ), "\n")
  cat("  Composite score:", round(best_enhanced_config$Composite_Score, 2), "\n")

  # Compare with AIC-only best
  cat("\nComparison with AIC-only selection:\n")
  if (best_enhanced_config$Distribution != best_config$Distribution ||
    best_enhanced_config$Lower_Cutoff != best_config$Lower_Cutoff) {
    cat("  NOTE: Composite score selects different configuration than AIC alone!\n")
    cat(
      "  AIC-only best:", best_config$Distribution, "with",
      scales::comma(best_config$Lower_Cutoff), "cutoff\n"
    )
    cat(
      "  Composite best:", best_enhanced_config$Distribution, "with",
      scales::comma(best_enhanced_config$Lower_Cutoff), "cutoff\n"
    )
  } else {
    cat("  Both methods agree on the same configuration.\n")
  }
} else {
  cat("No successful validation results to combine with statistical measures.\n")
  enhanced_comparison <- comparison_df
  best_enhanced_config <- best_config
}

# ===================================================================================
# Test mixture model with best enhanced configuration
# ===================================================================================
cat("\n\n9. FULL MIXTURE MODEL WITH BEST ENHANCED CONFIGURATION\n")
cat("-------------------------------------------------------------------\n")

# Update config with best enhanced configuration
conf_test <- conf
conf_test$thresholds$lower_cutoff <- best_enhanced_config$Lower_Cutoff
conf_test$thresholds$deaths_scaled <- tail_threshold
conf_test$mixture_dist <- best_enhanced_config$Distribution

specs <- tibble(
  Variable = conf_test$variables,
  Year = conf_test$base_year,
  Threshold = inv_dual_transform(u_best, h = conf_test$pop_reference, l = conf_test$pop_min),
  Cutoff = conf_test$thresholds$lower_cutoff,
  n = nrow(data),
  PoT = sum(data$dual > u_best),
  PcoT = mean(data$dual > u_best),
  PoB = sum(data$dual > conf_test$thresholds$lower_cutoff)
)

cat("Fitting full mixture model with:\n")
cat("  Bulk distribution:", conf_test$mixture_dist, "\n")
cat("  Lower cutoff:", scales::comma(conf_test$thresholds$lower_cutoff), "\n")
cat("  Tail threshold:", scales::comma(conf_test$thresholds$deaths_scaled), "\n")

# Fit two-stage mixture model with flexible bulk and GPD tail
mixture_fit <- fit_mixture_model(
  data, variable, conf_test, conf_test$thresholds$deaths_scaled, conf_test$mixture_dist,
  constrain_survival = conf$use_constrained_model,
  summary = specs
)


cat("\nMixture model results:\n")
cat("  Bulk conditional mean:", scales::comma(round(mixture_fit$bulk$conditional_mean, 0)), "deaths\n")
cat("  Tail shadow mean:", scales::comma(round(mixture_fit$tail$shadow, 0)), "deaths\n")
cat("  Overall expected severity:", scales::comma(round(mixture_fit$mean, 0)), "deaths\n")
cat("  Bulk weight:", round(mixture_fit$bulk$bulk_prob * 100, 1), "%\n")
cat("  Tail weight:", round(mixture_fit$tail$tail_prob * 100, 1), "%\n")

# ===================================================================================
# Final Summary and Recommendations
# ===================================================================================
cat("\n===================================================================================\n")
cat("BULK DIAGNOSTICS COMPLETED - FINAL RECOMMENDATIONS\n")
cat("===================================================================================\n")

cat("STATISTICAL FIT ONLY (AIC-based selection):\n")
cat("  Distribution:", best_config$Distribution, "\n")
cat("  Lower cutoff:", scales::comma(best_config$Lower_Cutoff), "deaths\n")
cat("  AIC:", round(best_config$AIC, 1), "\n")
cat("  Conditional mean:", scales::comma(round(best_config$Conditional_Mean, 0)), "deaths\n")
cat("  Continuity error:", sprintf("%.1f%%", best_config$Lower_Continuity_Error * 100), "\n\n")

if (exists("best_enhanced_config") && !is.null(validation_performance) && nrow(validation_performance) > 0) {
  cat("ENHANCED SELECTION (Statistical Fit + Predictive Performance):\n")
  cat("  Distribution:", best_enhanced_config$Distribution, "\n")
  cat("  Lower cutoff:", scales::comma(best_enhanced_config$Lower_Cutoff), "deaths\n")
  cat("  AIC:", round(best_enhanced_config$AIC, 1), "\n")
  cat("  RMSE:", ifelse(is.na(best_enhanced_config$RMSE), "Failed",
    scales::comma(round(best_enhanced_config$RMSE, 0))
  ), "\n")
  cat("  MAE:", ifelse(is.na(best_enhanced_config$MAE), "Failed",
    scales::comma(round(best_enhanced_config$MAE, 0))
  ), "\n")
  cat("  MAPE:", ifelse(is.na(best_enhanced_config$MAPE), "Failed",
    paste0(round(best_enhanced_config$MAPE, 1), "%")
  ), "\n")
  cat("  Coverage probability:", ifelse(is.na(best_enhanced_config$Coverage_Prob), "Failed",
    paste0(round(best_enhanced_config$Coverage_Prob, 1), "%")
  ), "\n")
  cat("  Composite score:", round(best_enhanced_config$Composite_Score, 2), "\n\n")

  cat("RECOMMENDATION:\n")
  if (best_enhanced_config$Distribution == best_config$Distribution &&
    best_enhanced_config$Lower_Cutoff == best_config$Lower_Cutoff) {
    cat("  ✓ Both statistical fit and predictive performance agree.\n")
    cat(
      "  ✓ RECOMMENDED: Use", best_enhanced_config$Distribution, "distribution with",
      scales::comma(best_enhanced_config$Lower_Cutoff), "deaths lower cutoff.\n"
    )
  } else {
    cat("  ⚠ Statistical fit and predictive performance suggest different configurations.\n")
    cat(
      "  ✓ RECOMMENDED: Use", best_enhanced_config$Distribution, "distribution with",
      scales::comma(best_enhanced_config$Lower_Cutoff), "deaths lower cutoff.\n"
    )
    cat("  📊 Rationale: Prioritizes real-world predictive accuracy over in-sample fit.\n")
  }
} else {
  cat("VALIDATION RESULTS: Failed or incomplete.\n")
  cat("RECOMMENDATION: Use AIC-based selection as fallback.\n")
  cat(
    "  ✓ RECOMMENDED: Use", best_config$Distribution, "distribution with",
    scales::comma(best_config$Lower_Cutoff), "deaths lower cutoff.\n"
  )
}

cat("\nDiagnostic analysis completed. No results saved.\n")
cat("===================================================================================\n")
