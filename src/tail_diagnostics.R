####################################################################################################
## Author:      Horacio Chacon Torrico
##
## Description: Tail diagnostics script for extreme value analysis of pandemic data.
##              This script generates publication-ready diagnostic plots for extreme value
##              threshold selection and model validation using the mev package.
##
## Generates:   - Mean excess plot (mev::automrl) showing mean excess vs threshold
##              - Threshold stability plot (mev::tstab.gpd) showing parameter stability
##
## Requires:    - Configuration file ("config.yml") with analysis parameters
##              - All R functions loaded from R/ directory
##              - mev package for extreme value diagnostics
##
## Outputs:     - Two publication-ready diagnostic plots
##              - Console output with diagnostic information
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

# ===================================================================================
# Setup output management
# ===================================================================================
# Setup script-specific timestamped output directory
run_info <- setup_script_run("tail_diagnostics", conf)
output_dir <- run_info$output_dir

# ===================================================================================
# Data preparation - identical to run_PoT.R
# ===================================================================================
# Read and transform the pandemic data
data <- read_transform_data(conf)

# Filter data using same criteria as analysis
data <- data %>%
  filter(
    start_year >= conf$base_year,
    deaths > conf$pop_min
  )

# Apply dual transformation to standardize pandemic size relative to population
data$dual <- dual_transform(data, "deaths_scaled", h = conf$pop_reference, l = conf$pop_min)

# ===================================================================================
# Threshold extraction and range selection
# ===================================================================================
# Extract threshold ranges for testing
thresholds <- list(
  deaths = threshold_extract(data, "deaths", conf$base_year, conf),
  deaths_scaled = threshold_extract(data, "deaths_scaled", conf$base_year, conf)
) %>%
  map(~ .x[30:90])

# ===================================================================================
# Diagnostic plot generation with publication-ready formatting
# ===================================================================================

# Generate mean excess plot with optimized formatting for standalone presentation
cat("Generating mean excess plot...\n")

# Generate the mean excess plot computation without plotting
mrl_result <- mev::automrl(data$dual, plot = FALSE)

# Recreate the plot with custom labels using the same logic as automrl
k <- 8
xdat <- sort(data$dual[is.finite(data$dual)], decreasing = TRUE)
n <- length(xdat)
meanX <- cumsum(xdat) / seq_along(xdat)
varX <- cumsum((xdat[-1] - meanX[-1]) * (xdat[-1] - meanX[-n])) / seq.int(1L, to = n - 1, by = 1L)
excu <- (meanX[-n] - xdat[-1])[-(1:k)]
weights <- (seq_len(n - 1) / varX)[-(1:k)]
xk <- xdat[-(1:(k + 1))]

# Save mean excess plot
save_base_plot_if_enabled({
  par(
    mfrow = c(1, 1), # Reset to single plot layout
    mar = c(4, 4, 2, 2), # Reduced margins since no title/subtitle
    cex.axis = 0.8,
    cex.lab = 1.0,
    font.lab = 1,
    font.axis = 1
  )
  
  # Create the plot with custom axis labels and scaled axes
  plot(
    x = xk / 1e5, # Convert x-axis to 100,000s scale
    y = excu / 1e6, # Convert y-axis to millions scale
    ylab = "Mean Excess (Millions)",
    xlab = "Threshold (100,000s Population-Scaled Deaths)",
    pch = 16,
    col = "#00468B",
    cex = 0.8,
    bty = "l",
    axes = FALSE
  ) # Suppress default axes to add custom ones
  
  # Add custom x-axis
  axis(1, at = pretty(xk / 1e5), labels = pretty(xk / 1e5), cex.axis = 0.8)
  # Add custom y-axis
  axis(2, at = pretty(excu / 1e6), labels = pretty(excu / 1e6), cex.axis = 0.8)
  # Add box around plot
  box()
  
  # Add the 2.8 million threshold line (reference threshold from config)
  abline(v = conf$thresholds$deaths_scaled / 1e5, lty = 2, col = "#ED0000", lwd = 2)
  
  # Calculate the fitted line coefficients (simplified from automrl source)
  fit_subset <- xk > mrl_result$thresh
  if (sum(fit_subset) > 0) {
    fit <- lm(excu ~ xk, weights = weights, subset = fit_subset)
    # Convert fitted line to scaled coordinates
    # For scaled coordinates: y_scaled = (a + b*x) / 1e6, x_scaled = x / 1e5
    # So: y_scaled = (a/1e6) + (b*1e5/1e6) * x_scaled
    abline(
      a = fit$coefficients[1] / 1e6,
      b = fit$coefficients[2] * 1e5 / 1e6,
      col = "#925E9F",
      lwd = 2
    )
  }
  
  # Add grid for better readability
  grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)
}, "mean_excess_plot", output_dir, conf)

# Also display plot interactively
par(
  mfrow = c(1, 1), # Reset to single plot layout
  mar = c(4, 4, 2, 2), # Reduced margins since no title/subtitle
  cex.axis = 0.8,
  cex.lab = 1.0,
  font.lab = 1,
  font.axis = 1
)

# Create the plot with custom axis labels and scaled axes
plot(
  x = xk / 1e5, # Convert x-axis to 100,000s scale
  y = excu / 1e6, # Convert y-axis to millions scale
  ylab = "Mean Excess (Millions)",
  xlab = "Threshold (100,000s Population-Scaled Deaths)",
  pch = 16,
  col = "#00468B",
  cex = 0.8,
  bty = "l",
  axes = FALSE
) # Suppress default axes to add custom ones

# Add custom x-axis
axis(1, at = pretty(xk / 1e5), labels = pretty(xk / 1e5), cex.axis = 0.8)
# Add custom y-axis
axis(2, at = pretty(excu / 1e6), labels = pretty(excu / 1e6), cex.axis = 0.8)
# Add box around plot
box()

# Add the 2.8 million threshold line (reference threshold from config)
abline(v = conf$thresholds$deaths_scaled / 1e5, lty = 2, col = "#ED0000", lwd = 2)

# Calculate the fitted line coefficients (simplified from automrl source)
fit_subset <- xk > mrl_result$thresh
if (sum(fit_subset) > 0) {
  fit <- lm(excu ~ xk, weights = weights, subset = fit_subset)
  # Convert fitted line to scaled coordinates
  # For scaled coordinates: y_scaled = (a + b*x) / 1e6, x_scaled = x / 1e5
  # So: y_scaled = (a/1e6) + (b*1e5/1e6) * x_scaled
  abline(
    a = fit$coefficients[1] / 1e6,
    b = fit$coefficients[2] * 1e5 / 1e6,
    col = "#925E9F",
    lwd = 2
  )
}

# Add grid for better readability
grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)

# Generate threshold stability plot with publication formatting
cat("Generating threshold stability plot...\n")

# Generate threshold stability plot computation without plotting
tstab_result <- mev::tstab.gpd(data$dual, thresholds$deaths_scaled, method = "profile", plot = FALSE)

# Calculate plot limits to include confidence intervals
scale_ylim <- range(c(tstab_result$mle[, 1], tstab_result$lower[, 1], tstab_result$upper[, 1]), na.rm = TRUE)
shape_ylim <- range(c(tstab_result$mle[, 2], tstab_result$lower[, 2], tstab_result$upper[, 2]), na.rm = TRUE)

# Save threshold stability plots
save_base_plot_if_enabled({
  # Set up plots one on top of each other with optimized spacing
  par(
    mfrow = c(2, 1),
    mar = c(3, 4, 1, 2), # Reduce top and bottom margins
    oma = c(2, 0, 0, 0)
  ) # Add outer margin only at bottom for shared x-label
  
  # First plot: Scale Parameter
  plot(tstab_result$threshold / 1e5, tstab_result$mle[, 1],
    xlab = "", # Remove x-label from top plot
    ylab = "Scale Parameter",
    type = "l",
    col = "#00468B",
    lwd = 2,
    bty = "l",
    ylim = scale_ylim
  )
  
  # Add confidence intervals for scale parameter
  lines(tstab_result$threshold / 1e5, tstab_result$lower[, 1], lty = 2, col = "#00468B")
  lines(tstab_result$threshold / 1e5, tstab_result$upper[, 1], lty = 2, col = "#00468B")
  
  # Add grid and threshold line
  grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)
  abline(v = conf$thresholds$deaths_scaled / 1e5, lty = 2, col = "#ED0000", lwd = 2)
  
  # Second plot: Shape Parameter
  plot(tstab_result$threshold / 1e5, tstab_result$mle[, 2],
    xlab = "Threshold (100,000s Population-Scaled Deaths)",
    ylab = "Shape Parameter",
    type = "l",
    col = "#925E9F",
    lwd = 2,
    bty = "l",
    ylim = shape_ylim
  )
  
  # Add confidence intervals for shape parameter
  lines(tstab_result$threshold / 1e5, tstab_result$lower[, 2], lty = 2, col = "#925E9F")
  lines(tstab_result$threshold / 1e5, tstab_result$upper[, 2], lty = 2, col = "#925E9F")
  
  # Add grid and threshold line
  grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)
  abline(v = conf$thresholds$deaths_scaled / 1e5, lty = 2, col = "#ED0000", lwd = 2)
  
  # Add shared x-axis label in outer margin
  mtext("Threshold (100,000s Population-Scaled Deaths)", side = 1, outer = TRUE, line = 0.5, cex = 1.0)
}, "threshold_stability_plot", output_dir, conf)

# Also display plots interactively
par(
  mfrow = c(2, 1),
  mar = c(3, 4, 1, 2), # Reduce top and bottom margins
  oma = c(2, 0, 0, 0)
) # Add outer margin only at bottom for shared x-label

# First plot: Scale Parameter
plot(tstab_result$threshold / 1e5, tstab_result$mle[, 1],
  xlab = "", # Remove x-label from top plot
  ylab = "Scale Parameter",
  type = "l",
  col = "#00468B",
  lwd = 2,
  bty = "l",
  ylim = scale_ylim
)

# Add confidence intervals for scale parameter
lines(tstab_result$threshold / 1e5, tstab_result$lower[, 1], lty = 2, col = "#00468B")
lines(tstab_result$threshold / 1e5, tstab_result$upper[, 1], lty = 2, col = "#00468B")

# Add grid and threshold line
grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)
abline(v = conf$thresholds$deaths_scaled / 1e5, lty = 2, col = "#ED0000", lwd = 2)

# Second plot: Shape Parameter
plot(tstab_result$threshold / 1e5, tstab_result$mle[, 2],
  xlab = "Threshold (100,000s Population-Scaled Deaths)",
  ylab = "Shape Parameter",
  type = "l",
  col = "#925E9F",
  lwd = 2,
  bty = "l",
  ylim = shape_ylim
)

# Add confidence intervals for shape parameter
lines(tstab_result$threshold / 1e5, tstab_result$lower[, 2], lty = 2, col = "#925E9F")
lines(tstab_result$threshold / 1e5, tstab_result$upper[, 2], lty = 2, col = "#925E9F")

# Add grid and threshold line
grid(col = "#E5E5E5", lty = "dotted", lwd = 0.5)
abline(v = conf$thresholds$deaths_scaled / 1e5, lty = 2, col = "#ED0000", lwd = 2)

# Add shared x-axis label in outer margin
mtext("Threshold (100,000s Population-Scaled Deaths)", side = 1, outer = TRUE, line = 0.5, cex = 1.0)

# Reset to single plot layout
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

# ===================================================================================
# GPD Model Fitting and Parameter Estimation
# ===================================================================================

# Fit GPD model using the optimal threshold from automrl
optimal_thresh <- mrl_result$thresh
cat("Fitting GPD model with optimal threshold:", round(optimal_thresh, 0), "deaths\n")

# Fit the GPD model
gpd_fit <- mev::fit.gpd(data$dual, threshold = optimal_thresh)

# Extract parameter estimates and standard errors
optimal_scale <- gpd_fit$estimate[1] # Scale parameter
optimal_shape <- gpd_fit$estimate[2] # Shape parameter
scale_se <- gpd_fit$std.err[1] # Scale standard error
shape_se <- gpd_fit$std.err[2] # Shape standard error

# Calculate 95% confidence intervals using normal approximation
z_score <- qnorm(0.975) # 97.5% quantile for 95% CI
scale_ci_lower <- optimal_scale - z_score * scale_se
scale_ci_upper <- optimal_scale + z_score * scale_se
shape_ci_lower <- optimal_shape - z_score * shape_se
shape_ci_upper <- optimal_shape + z_score * shape_se

# Calculate additional diagnostics
n_total <- nrow(data)
n_exceedances <- sum(data$dual > optimal_thresh)
prop_exceedances <- n_exceedances / n_total
threshold_original_scale <- inv_dual_transform(optimal_thresh, h = conf$pop_reference, l = conf$pop_min)

# Extract model diagnostics
log_likelihood <- gpd_fit$nllh
aic <- AIC(gpd_fit)
n_params <- length(gpd_fit$estimate)
convergence <- gpd_fit$convergence == "unsuccessful"

# Print detailed summary
cat("\n=== GPD Model Fitting Results ===\n")
cat("Model Convergence:", ifelse(convergence == 0, "✓ Successful", "⚠ Warning"), "\n")
cat("Log-likelihood:", round(log_likelihood, 2), "\n")
cat("AIC:", round(aic, 2), "\n")
cat("Number of parameters:", n_params, "\n\n")

cat("Optimal Threshold Selection:\n")
cat("- Threshold (population-scaled):", round(optimal_thresh, 0), "deaths\n")
cat("- Threshold (original scale):", round(threshold_original_scale, 0), "deaths\n")
cat("- Number of exceedances:", n_exceedances, "out of", n_total, "events\n")
cat("- Proportion exceeding threshold:", round(prop_exceedances * 100, 1), "%\n\n")

cat("GPD Parameter Estimates (Maximum Likelihood):\n")
cat(
  "- Shape parameter (ξ):", round(optimal_shape, 3), "±", round(shape_se, 3),
  "(95% UI:", round(shape_ci_lower, 3), "–", round(shape_ci_upper, 3), ")\n"
)
cat(
  "- Scale parameter (β):", round(optimal_scale, 0), "±", round(scale_se, 0),
  "(95% UI:", round(scale_ci_lower, 0), "–", round(scale_ci_upper, 0), ")\n\n"
)

# Create enhanced markdown table for publication
cat("=== Publication-Ready Markdown Table ===\n")
cat("\n")
cat("| Parameter | Estimate | Standard Error | 95% Uncertainty Interval |\n")
cat("|-----------|----------|----------------|---------------------------|\n")
cat("| Threshold (population-scaled deaths) |", round(optimal_thresh, 0), "| – | – |\n")
cat("| Threshold (original scale deaths) |", round(threshold_original_scale, 0), "| – | – |\n")
cat("| Number of exceedances |", n_exceedances, "| – | – |\n")
cat("| Proportion exceeding threshold |", round(prop_exceedances * 100, 1), "% | – | – |\n")
cat("| Shape parameter (ξ) |", round(optimal_shape, 3), "|", round(shape_se, 3), "|", round(shape_ci_lower, 3), "–", round(shape_ci_upper, 3), "|\n")
cat("| Scale parameter (β) |", round(optimal_scale, 0), "|", round(scale_se, 0), "|", round(scale_ci_lower, 0), "–", round(scale_ci_upper, 0), "|\n")
cat("| Log-likelihood |", round(log_likelihood, 2), "| – | – |\n")
cat("| AIC |", round(aic, 2), "| – | – |\n")
cat("\n")

# Additional model diagnostics
cat("=== Model Diagnostics ===\n")
if (optimal_shape < 1) {
  cat("✓ Shape parameter ξ < 1: Finite mean, infinite variance (heavy-tailed)\n")
} else if (optimal_shape < 0.5) {
  cat("✓ Shape parameter ξ < 0.5: Finite mean and variance\n")
} else {
  cat("⚠ Shape parameter ξ ≥ 1: Infinite mean (extremely heavy-tailed)\n")
}

if (optimal_shape > 0) {
  cat("✓ Shape parameter ξ > 0: Fréchet domain (heavy-tailed distribution)\n")
} else if (optimal_shape == 0) {
  cat("✓ Shape parameter ξ = 0: Gumbel domain (exponential-type tail)\n")
} else {
  cat("✓ Shape parameter ξ < 0: Weibull domain (finite upper bound)\n")
}

cat("✓ Threshold selection based on automated mean residual life procedure\n")
cat("✓", round(prop_exceedances * 100, 1), "% of data used for tail fitting\n")

cat("\nCompleted tail diagnostics analysis.\n")
