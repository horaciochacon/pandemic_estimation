#' Model Diagnostics for Pandemic Estimation Mixture Models
#'
#' This module provides detailed diagnostic tools for the two-stage mixture model
#' used in pandemic severity estimation. It includes threshold selection diagnostics,
#' goodness-of-fit testing, mixture model component evaluation, bootstrap stability
#' assessment, and cross-validation procedures adapted for small sample sizes.
#'
#' Key considerations for pandemic data:
#' - Small sample size (107 observations) affects statistical power
#' - Heavy-tailed distribution requires robust extreme value diagnostics
#' - Two-stage mixture model needs component-specific validation
#' - Threshold selection is critical but constrained by data scarcity

# Helper functions for GPD distribution (if not available from packages)
if (!exists("pgpd")) {
  #' GPD Cumulative Distribution Function
  #' @param q Quantiles
  #' @param shape Shape parameter (xi)
  #' @param scale Scale parameter (beta)
  #' @return Cumulative probabilities
  pgpd <- function(q, shape, scale) {
    if (shape != 0) {
      1 - (1 + shape * q / scale)^(-1 / shape)
    } else {
      1 - exp(-q / scale)
    }
  }
}

if (!exists("qgpd")) {
  #' GPD Quantile Function
  #' @param p Probabilities
  #' @param shape Shape parameter (xi)
  #' @param scale Scale parameter (beta)
  #' @return Quantiles
  qgpd <- function(p, shape, scale) {
    if (shape != 0) {
      scale * ((1 - p)^(-shape) - 1) / shape
    } else {
      -scale * log(1 - p)
    }
  }
}

if (!exists("dgpd")) {
  #' GPD Probability Density Function
  #' @param x Values
  #' @param shape Shape parameter (xi)
  #' @param scale Scale parameter (beta)
  #' @return Densities
  dgpd <- function(x, shape, scale) {
    if (shape != 0) {
      (1 / scale) * (1 + shape * x / scale)^(-1 / shape - 1)
    } else {
      (1 / scale) * exp(-x / scale)
    }
  }
}

#' Diagnose Threshold Selection for GPD Tail Modeling
#'
#' Detailed threshold selection diagnostics that account for data scarcity.
#' Uses multiple approaches including mean excess plots, parameter stability,
#' and automated threshold selection with confidence intervals.
#'
#' @param data Data frame with dual-transformed pandemic data
#' @param variable Variable name to analyze
#' @param conf Configuration object
#' @param n_thresholds Number of thresholds to test (limited by sample size)
#' @param generate_plots Whether to generate diagnostic plots
#' @return List containing threshold diagnostics and recommendations
#' @export
diagnose_threshold_selection <- function(data, variable, conf,
                                         n_thresholds = 15, generate_plots = TRUE) {
  # Extract dual-transformed variable
  x <- data$dual

  # Create threshold candidates from upper quantiles (accounting for small sample)
  # Start from 80th percentile to ensure sufficient tail observations
  min_tail_obs <- max(5, length(x) * 0.05) # At least 5 observations or 5% of data
  max_quantile <- 1 - min_tail_obs / length(x)
  threshold_quantiles <- seq(0.70, max_quantile, length.out = n_thresholds)
  thresholds <- quantile(x, threshold_quantiles)

  # Initialize storage for diagnostics
  threshold_results <- data.frame(
    threshold = thresholds,
    n_exceedances = NA,
    xi = NA,
    xi_se = NA,
    beta = NA,
    beta_se = NA,
    loglik = NA,
    aic = NA,
    mean_excess = NA,
    mean_excess_se = NA
  )

  # Calculate diagnostics for each threshold
  for (i in seq_along(thresholds)) {
    u <- thresholds[i]
    exceedances <- x[x > u] - u
    n_exc <- length(exceedances)

    if (n_exc >= 3) { # Minimum for GPD fitting
      # Fit GPD model
      tryCatch(
        {
          gpd_fit <- mev::fit.gpd(x, threshold = u)

          # Extract parameters
          xi <- gpd_fit$estimate["shape"]
          beta <- gpd_fit$estimate["scale"]
          xi_se <- gpd_fit$std.err["shape"]
          beta_se <- gpd_fit$std.err["scale"]

          # Store results
          threshold_results$n_exceedances[i] <- n_exc
          threshold_results$xi[i] <- xi
          threshold_results$xi_se[i] <- xi_se
          threshold_results$beta[i] <- beta
          threshold_results$beta_se[i] <- beta_se
          threshold_results$loglik[i] <- -gpd_fit$nllh
          threshold_results$aic[i] <- -2 * -gpd_fit$nllh + 2 * 2 # 2 parameters

          # Mean excess statistics
          threshold_results$mean_excess[i] <- mean(exceedances)
          threshold_results$mean_excess_se[i] <- sd(exceedances) / sqrt(n_exc)
        },
        error = function(e) {
          # Handle fitting failures
          threshold_results$n_exceedances[i] <- n_exc
        }
      )
    }
  }

  # Remove failed fits
  threshold_results <- threshold_results[!is.na(threshold_results$xi), ]

  # Assess threshold stability
  stability_metrics <- assess_threshold_stability(threshold_results)

  # Recommend threshold range based on stability
  recommended_range <- recommend_threshold_range(threshold_results, stability_metrics)

  # Generate plots if requested
  plots <- list()
  if (generate_plots) {
    plots$mean_excess <- create_mean_excess_plot(data, variable, conf)
    plots$parameter_stability <- create_parameter_stability_plot(threshold_results)
    plots$threshold_choice <- create_threshold_choice_plot(
      threshold_results,
      recommended_range
    )
  }

  # Compile results
  list(
    threshold_results = threshold_results,
    stability_metrics = stability_metrics,
    recommended_range = recommended_range,
    plots = plots,
    diagnostics_summary = list(
      n_tested = nrow(threshold_results),
      stability_score = stability_metrics$overall_stability,
      recommended_threshold = recommended_range[1]
    )
  )
}

#' Assess Threshold Stability Metrics
#'
#' Evaluates the stability of GPD parameters across different threshold choices.
#' Uses coefficient of variation and trend analysis to assess stability.
#'
#' @param threshold_results Data frame with threshold diagnostic results
#' @return List containing stability metrics
assess_threshold_stability <- function(threshold_results) {
  # Calculate coefficient of variation for key parameters
  xi_cv <- sd(threshold_results$xi, na.rm = TRUE) / abs(mean(threshold_results$xi, na.rm = TRUE))
  beta_cv <- sd(threshold_results$beta, na.rm = TRUE) / abs(mean(threshold_results$beta, na.rm = TRUE))

  # Test for trends in parameters (should be stable)
  xi_trend_p <- tryCatch(
    {
      cor.test(threshold_results$threshold, threshold_results$xi)$p.value
    },
    error = function(e) NA
  )

  beta_trend_p <- tryCatch(
    {
      cor.test(threshold_results$threshold, threshold_results$beta)$p.value
    },
    error = function(e) NA
  )

  # Overall stability score (lower is better)
  overall_stability <- mean(c(xi_cv, beta_cv), na.rm = TRUE)

  list(
    xi_cv = xi_cv,
    beta_cv = beta_cv,
    xi_trend_p = xi_trend_p,
    beta_trend_p = beta_trend_p,
    overall_stability = overall_stability,
    stability_rating = ifelse(overall_stability < 0.3, "Good",
      ifelse(overall_stability < 0.5, "Moderate", "Poor")
    )
  )
}

#' Recommend Threshold Range Based on Diagnostics
#'
#' Provides threshold recommendations based on stability metrics and sample size constraints.
#'
#' @param threshold_results Data frame with threshold diagnostic results
#' @param stability_metrics List with stability assessment
#' @return Vector with recommended threshold range
recommend_threshold_range <- function(threshold_results, stability_metrics) {
  # Filter to thresholds with sufficient exceedances (at least 10 for reliability)
  reliable_thresholds <- threshold_results[threshold_results$n_exceedances >= 10, ]

  if (nrow(reliable_thresholds) == 0) {
    # Fallback to minimum acceptable sample
    reliable_thresholds <- threshold_results[threshold_results$n_exceedances >= 5, ]
  }

  if (nrow(reliable_thresholds) > 0) {
    # Select threshold with best balance of sample size and parameter stability
    # Prefer thresholds with more observations but reasonable stability
    reliability_score <- reliable_thresholds$n_exceedances / max(reliable_thresholds$n_exceedances) -
      abs(reliable_thresholds$xi - median(reliable_thresholds$xi, na.rm = TRUE))

    best_idx <- which.max(reliability_score)
    recommended_threshold <- reliable_thresholds$threshold[best_idx]

    # Create range around recommended threshold
    threshold_range <- c(
      max(recommended_threshold * 0.9, min(reliable_thresholds$threshold)),
      min(recommended_threshold * 1.1, max(reliable_thresholds$threshold))
    )
  } else {
    # Fallback to original threshold from config
    threshold_range <- c(NA, NA)
  }

  threshold_range
}

#' Create Mean Excess Plot with Confidence Intervals
#'
#' Enhanced mean excess plot that includes confidence intervals and stability assessment.
#' Adapted for small sample sizes with appropriate confidence interval adjustments.
#'
#' @param data Data frame with pandemic data
#' @param variable Variable name to analyze
#' @param conf Configuration object
#' @return ggplot object for mean excess plot
create_mean_excess_plot <- function(data, variable, conf) {
  x <- sort(data$dual)
  n <- length(x)

  # Create threshold candidates (more conservative for small sample)
  min_exceedances <- 5
  max_threshold_idx <- n - min_exceedances
  threshold_indices <- seq(min_exceedances, max_threshold_idx, by = 2)

  # Calculate mean excess for each threshold
  mean_excess_data <- data.frame(
    threshold = numeric(0),
    mean_excess = numeric(0),
    lower_ci = numeric(0),
    upper_ci = numeric(0),
    n_exceedances = numeric(0)
  )

  for (i in threshold_indices) {
    if (i <= n) {
      u <- x[i]
      exceedances <- x[x > u] - u
      n_exc <- length(exceedances)

      if (n_exc > 0) {
        mean_exc <- mean(exceedances)
        se_exc <- sd(exceedances) / sqrt(n_exc)

        # Confidence interval (wider for small samples)
        t_crit <- qt(0.975, df = max(1, n_exc - 1))

        mean_excess_data <- rbind(mean_excess_data, data.frame(
          threshold = u,
          mean_excess = mean_exc,
          lower_ci = mean_exc - t_crit * se_exc,
          upper_ci = mean_exc + t_crit * se_exc,
          n_exceedances = n_exc
        ))
      }
    }
  }

  # Transform back to original scale for interpretability
  mean_excess_data$threshold_orig <- inv_dual_transform(
    mean_excess_data$threshold,
    h = conf$pop_reference,
    l = conf$pop_min
  )
  mean_excess_data$mean_excess_orig <- inv_dual_transform(
    mean_excess_data$mean_excess,
    h = conf$pop_reference,
    l = conf$pop_min
  )

  # Create plot
  ggplot(mean_excess_data, aes(x = threshold_orig, y = mean_excess_orig)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, fill = "blue") +
    geom_line(color = "blue", size = 1) +
    geom_point(aes(size = n_exceedances), color = "darkblue", alpha = 0.7) +
    scale_x_log10(labels = scales::comma_format()) +
    scale_y_log10(labels = scales::comma_format()) +
    labs(
      title = "Mean Excess Plot for Threshold Selection",
      subtitle = "For GPD tail modeling - points show number of exceedances",
      x = "Threshold (deaths)",
      y = "Mean Excess (deaths)",
      size = "Exceedances",
      caption = "Blue ribbon shows 95% confidence interval"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
}

#' Create Parameter Stability Plot
#'
#' Plots GPD parameter estimates across different thresholds to assess stability.
#' Essential for threshold selection in extreme value analysis.
#'
#' @param threshold_results Data frame with threshold diagnostic results
#' @return ggplot object showing parameter stability
create_parameter_stability_plot <- function(threshold_results) {
  # Prepare data for plotting
  plot_data <- threshold_results %>%
    dplyr::select(threshold, xi, xi_se, beta, beta_se, n_exceedances) %>%
    tidyr::pivot_longer(cols = c(xi, beta), names_to = "parameter", values_to = "estimate") %>%
    tidyr::pivot_longer(cols = c(xi_se, beta_se), names_to = "se_param", values_to = "se") %>%
    filter(
      (parameter == "xi" & se_param == "xi_se") |
        (parameter == "beta" & se_param == "beta_se")
    ) %>%
    mutate(
      lower_ci = estimate - 1.96 * se,
      upper_ci = estimate + 1.96 * se,
      parameter_label = ifelse(parameter == "xi", "Shape (ξ)", "Scale (β)")
    )

  # Create faceted plot
  ggplot(plot_data, aes(x = threshold, y = estimate)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, fill = "red") +
    geom_line(color = "red", size = 1) +
    geom_point(aes(size = n_exceedances), color = "darkred", alpha = 0.7) +
    facet_wrap(~parameter_label, scales = "free_y", ncol = 1) +
    labs(
      title = "GPD Parameter Stability Across Thresholds",
      subtitle = "Parameters should be stable in the valid threshold range",
      x = "Threshold (dual-transformed scale)",
      y = "Parameter Estimate",
      size = "Exceedances",
      caption = "Red ribbon shows 95% confidence interval"
    ) +
    scale_x_log10() +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
}

#' Create Threshold Choice Recommendation Plot
#'
#' Visualizes the recommended threshold range with supporting diagnostics.
#'
#' @param threshold_results Data frame with threshold diagnostic results
#' @param recommended_range Vector with recommended threshold range
#' @return ggplot object for threshold recommendation
create_threshold_choice_plot <- function(threshold_results, recommended_range) {
  # Create plot showing multiple criteria
  p1 <- ggplot(threshold_results, aes(x = threshold)) +
    geom_line(aes(y = n_exceedances), color = "blue", size = 1) +
    geom_vline(xintercept = recommended_range, color = "red", linetype = "dashed", alpha = 0.7) +
    labs(
      title = "Threshold Selection Criteria",
      x = "Threshold (dual-transformed)",
      y = "Number of Exceedances"
    ) +
    scale_x_log10() +
    theme_bw()

  # AIC plot (if available)
  if (!all(is.na(threshold_results$aic))) {
    p2 <- ggplot(threshold_results, aes(x = threshold, y = aic)) +
      geom_line(color = "green", size = 1) +
      geom_point(color = "darkgreen") +
      geom_vline(xintercept = recommended_range, color = "red", linetype = "dashed", alpha = 0.7) +
      labs(
        x = "Threshold (dual-transformed)",
        y = "AIC"
      ) +
      scale_x_log10() +
      theme_bw()

    # Combine plots
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  } else {
    p1
  }
}

#' Diagnose Goodness-of-Fit for Mixture Model Components
#'
#' Detailed goodness-of-fit testing for both bulk and tail components
#' of the mixture model. Includes visual diagnostics and formal statistical tests.
#'
#' @param data Data frame with analysis data
#' @param model_results Results from mixture model fitting
#' @param conf Configuration object
#' @param generate_plots Whether to generate diagnostic plots
#' @return List containing GOF diagnostics and test results
#' @export
diagnose_goodness_of_fit <- function(data, model_results, conf, generate_plots = TRUE) {
  # Extract model components
  threshold <- inv_dual_transform(model_results$specs$Threshold,
    h = conf$pop_reference, l = conf$pop_min
  )
  lower_cutoff <- model_results$specs$Cutoff

  # Prepare data subsets
  all_data <- data[data[[conf$variables[1]]] >= lower_cutoff, ]
  bulk_data <- all_data[all_data[[conf$variables[1]]] < threshold, ]
  tail_data <- all_data[all_data[[conf$variables[1]]] >= threshold, ]

  # Initialize results storage
  test_results <- data.frame(
    component = character(0),
    test_name = character(0),
    statistic = numeric(0),
    p_value = numeric(0),
    conclusion = character(0)
  )

  plots <- list()

  # ===================================================================================
  # TAIL COMPONENT DIAGNOSTICS (GPD)
  # ===================================================================================
  if (nrow(tail_data) >= 3) {
    # Fit GPD to tail data for diagnostics
    tail_exceedances <- tail_data[[conf$variables[1]]] - threshold

    tryCatch(
      {
        # Re-fit GPD for diagnostic purposes
        gpd_fit <- mev::fit.gpd(all_data[[conf$variables[1]]], threshold = threshold)
        xi <- gpd_fit$estimate["shape"]
        beta <- gpd_fit$estimate["scale"]

        # Generate theoretical quantiles for Q-Q plot
        n_tail <- length(tail_exceedances)
        p_points <- ppoints(n_tail)
        theoretical_quantiles <- mev::qgp(p_points, shape = xi, scale = beta)

        # Kolmogorov-Smirnov test for GPD
        ks_test <- ks.test(tail_exceedances, function(x) mev::pgp(x, shape = xi, scale = beta))

        test_results <- rbind(test_results, data.frame(
          component = "Tail (GPD)",
          test_name = "Kolmogorov-Smirnov",
          statistic = ks_test$statistic,
          p_value = ks_test$p.value,
          conclusion = ifelse(ks_test$p.value > 0.05, "Adequate fit", "Poor fit")
        ))

        # Anderson-Darling test (if package available)
        tryCatch(
          {
            # Simplified AD test using empirical distribution
            ad_statistic <- calculate_ad_statistic(tail_exceedances, xi, beta, "gpd")

            test_results <- rbind(test_results, data.frame(
              component = "Tail (GPD)",
              test_name = "Anderson-Darling",
              statistic = ad_statistic,
              p_value = NA, # Critical values would need to be computed
              conclusion = ifelse(ad_statistic < 2.5, "Adequate fit", "Poor fit")
            ))
          },
          error = function(e) {}
        )

        # Create Q-Q plot for tail
        if (generate_plots) {
          qq_tail_data <- data.frame(
            theoretical = sort(theoretical_quantiles),
            empirical = sort(tail_exceedances)
          )

          plots$qq_tail <- ggplot(qq_tail_data, aes(x = theoretical, y = empirical)) +
            geom_point(alpha = 0.7, color = "blue") +
            geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
            labs(
              title = "Q-Q Plot: GPD Tail Component",
              x = "Theoretical Quantiles (GPD)",
              y = "Empirical Quantiles",
              subtitle = paste("Shape =", round(xi, 3), ", Scale =", round(beta, 0))
            ) +
            scale_y_log10() +
            scale_x_log10() +
            theme_bw()
        }
      },
      error = function(e) {
        test_results <- rbind(test_results, data.frame(
          component = "Tail (GPD)",
          test_name = "GPD Fit",
          statistic = NA,
          p_value = NA,
          conclusion = "Fitting failed"
        ))
      }
    )
  }

  # ===================================================================================
  # BULK COMPONENT DIAGNOSTICS
  # ===================================================================================
  if (nrow(bulk_data) >= 5 && !is.null(model_results$params$bulk_distribution)) {
    bulk_dist <- model_results$params$bulk_distribution

    tryCatch(
      {
        # Extract bulk distribution parameters based on type
        if (bulk_dist == "lognormal") {
          # Get parameters from model results
          bulk_values <- bulk_data[[conf$variables[1]]]

          # Fit distribution for diagnostic comparison
          bulk_fit <- MASS::fitdistr(bulk_values, "lognormal")
          meanlog <- estimate["meanlog"]
          sdlog <- bulk_fit$estimate["sdlog"]

          # KS test for bulk component
          ks_test_bulk <- ks.test(bulk_values, function(x) plnorm(x, meanlog, sdlog))

          test_results <- rbind(test_results, data.frame(
            component = "Bulk (Lognormal)",
            test_name = "Kolmogorov-Smirnov",
            statistic = ks_test_bulk$statistic,
            p_value = ks_test_bulk$p.value,
            conclusion = ifelse(ks_test_bulk$p.value > 0.05, "Adequate fit", "Poor fit")
          ))

          # Create Q-Q plot for bulk
          if (generate_plots) {
            n_bulk <- length(bulk_values)
            p_points_bulk <- ppoints(n_bulk)
            theoretical_bulk <- qlnorm(p_points_bulk, meanlog, sdlog)

            qq_bulk_data <- data.frame(
              theoretical = sort(theoretical_bulk),
              empirical = sort(bulk_values)
            )

            plots$qq_bulk <- ggplot(qq_bulk_data, aes(x = theoretical, y = empirical)) +
              geom_point(alpha = 0.7, color = "blue") +
              geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
              scale_x_log10() +
              scale_y_log10() +
              labs(
                title = "Q-Q Plot: Lognormal Bulk Component",
                x = "Theoretical Quantiles (Lognormal)",
                y = "Empirical Quantiles",
                subtitle = paste("μ =", round(meanlog, 3), ", σ =", round(sdlog, 3))
              ) +
              theme_bw()
          }
        }
      },
      error = function(e) {
        test_results <- rbind(test_results, data.frame(
          component = "Bulk",
          test_name = "Distribution Fit",
          statistic = NA,
          p_value = NA,
          conclusion = "Fitting failed"
        ))
      }
    )
  }

  # ===================================================================================
  # COMBINED PLOTS
  # ===================================================================================
  if (generate_plots) {
    # Combined Q-Q plots
    if ("qq_tail" %in% names(plots) && "qq_bulk" %in% names(plots)) {
      plots$qq_plots <- gridExtra::grid.arrange(plots$qq_bulk, plots$qq_tail, ncol = 2)
    } else if ("qq_tail" %in% names(plots)) {
      plots$qq_plots <- plots$qq_tail
    } else if ("qq_bulk" %in% names(plots)) {
      plots$qq_plots <- plots$qq_bulk
    }

    # P-P plots (empirical vs theoretical CDF)
    plots$pp_plots <- create_pp_plots(data, model_results, conf)

    # Residual plots
    plots$residual_plots <- create_residual_plots(data, model_results, conf)
  }

  # Compile results
  list(
    test_results = test_results,
    plots = plots,
    gof_summary = list(
      n_tests = nrow(test_results),
      n_adequate = sum(test_results$conclusion == "Adequate fit", na.rm = TRUE),
      overall_adequacy = ifelse(
        sum(test_results$conclusion == "Adequate fit", na.rm = TRUE) >=
          nrow(test_results) * 0.7, "Good", "Needs attention"
      )
    )
  )
}

#' Calculate Anderson-Darling Statistic
#'
#' Simplified Anderson-Darling test statistic calculation for goodness-of-fit.
#'
#' @param data Empirical data
#' @param param1 First distribution parameter
#' @param param2 Second distribution parameter
#' @param dist_type Distribution type ("gpd", "lognormal", etc.)
#' @return Anderson-Darling test statistic
calculate_ad_statistic <- function(data, param1, param2, dist_type) {
  n <- length(data)
  data_sorted <- sort(data)

  if (dist_type == "gpd") {
    # For GPD: param1 = xi (shape), param2 = beta (scale)
    F_vals <- mev::pgp(data_sorted, shape = param1, scale = param2)
  } else if (dist_type == "lognormal") {
    # For lognormal: param1 = meanlog, param2 = sdlog
    F_vals <- plnorm(data_sorted, meanlog = param1, sdlog = param2)
  }

  # Avoid numerical issues at boundaries
  F_vals <- pmax(F_vals, 1e-10)
  F_vals <- pmin(F_vals, 1 - 1e-10)

  # Calculate Anderson-Darling statistic
  i <- 1:n
  ad_stat <- -n - sum((2 * i - 1) * (log(F_vals) + log(1 - F_vals[n + 1 - i]))) / n

  return(ad_stat)
}

#' Create P-P Plots for Model Components
#'
#' Probability-probability plots comparing empirical and theoretical CDFs.
#'
#' @param data Analysis data
#' @param model_results Model fitting results
#' @param conf Configuration object
#' @return ggplot object with P-P plots
create_pp_plots <- function(data, model_results, conf) {
  # Implementation would create P-P plots for both components
  # Similar structure to Q-Q plots but comparing CDFs

  # Placeholder for now
  ggplot() +
    labs(
      title = "P-P Plots for Model Components",
      subtitle = "Empirical vs Theoretical CDFs"
    ) +
    theme_bw()
}

#' Create Residual Analysis Plots
#'
#' Generate residual plots for model diagnostic assessment.
#'
#' @param data Analysis data
#' @param model_results Model fitting results
#' @param conf Configuration object
#' @return ggplot object with residual analysis
create_residual_plots <- function(data, model_results, conf) {
  # Implementation would create standardized residual plots
  # Including plots of residuals vs fitted values, normal Q-Q of residuals, etc.

  # Placeholder for now
  ggplot() +
    labs(
      title = "Residual Analysis",
      subtitle = "Model diagnostic residuals"
    ) +
    theme_bw()
}

#' Diagnose Mixture Model Component Separation
#'
#' Evaluates how well the mixture model separates bulk and tail components
#' and assesses the appropriateness of the component weights.
#'
#' @param data Analysis data
#' @param model_results Model fitting results
#' @param conf Configuration object
#' @param generate_plots Whether to generate diagnostic plots
#' @return List containing component diagnostics
#' @export
diagnose_mixture_components <- function(data, model_results, conf, generate_plots = TRUE) {
  # Extract component information
  threshold <- inv_dual_transform(model_results$specs$Threshold,
    h = conf$pop_reference, l = conf$pop_min
  )
  lower_cutoff <- model_results$specs$Cutoff

  # Calculate empirical component weights
  total_data <- data[data[[conf$variables[1]]] >= lower_cutoff, ]
  bulk_data <- total_data[total_data[[conf$variables[1]]] < threshold, ]
  tail_data <- total_data[total_data[[conf$variables[1]]] >= threshold, ]

  empirical_bulk_weight <- nrow(bulk_data) / nrow(total_data)
  empirical_tail_weight <- nrow(tail_data) / nrow(total_data)

  # Component adequacy assessment
  bulk_adequacy <- assess_component_adequacy(
    bulk_data, "bulk",
    model_results$params$bulk_distribution
  )
  tail_adequacy <- assess_component_adequacy(tail_data, "tail", "gpd")

  # Weight stability analysis
  weight_stability <- assess_weight_stability(data, model_results, conf)

  plots <- list()
  if (generate_plots) {
    plots$component_separation <- create_component_separation_plot(data, model_results, conf)
    plots$weight_stability <- create_weight_stability_plot(weight_stability)
  }

  list(
    empirical_weights = list(bulk = empirical_bulk_weight, tail = empirical_tail_weight),
    bulk_adequacy = bulk_adequacy,
    tail_adequacy = tail_adequacy,
    weight_stability = weight_stability,
    plots = plots,
    component_summary = list(
      bulk_n = nrow(bulk_data),
      tail_n = nrow(tail_data),
      separation_quality = ifelse(empirical_tail_weight > 0.05 && empirical_tail_weight < 0.3,
        "Good", "Suboptimal"
      )
    )
  )
}

#' Assess Component Adequacy
#'
#' Evaluates whether a component has sufficient data and appropriate characteristics.
#'
#' @param component_data Data for the component
#' @param component_type Type of component ("bulk" or "tail")
#' @param distribution Distribution type
#' @return Component adequacy assessment
assess_component_adequacy <- function(component_data, component_type, distribution) {
  n <- nrow(component_data)

  if (component_type == "bulk") {
    min_required <- 10
    adequacy <- ifelse(n >= min_required, "Adequate", "Insufficient data")
  } else if (component_type == "tail") {
    min_required <- 5
    adequacy <- ifelse(n >= min_required, "Adequate", "Insufficient data")
  }

  list(
    n_observations = n,
    adequacy = adequacy,
    distribution = distribution
  )
}

#' Assess Weight Stability
#'
#' Tests stability of mixture component weights across threshold variations.
#'
#' @param data Analysis data
#' @param model_results Model fitting results
#' @param conf Configuration object
#' @return Weight stability assessment
assess_weight_stability <- function(data, model_results, conf) {
  # Test weight sensitivity to threshold changes
  base_threshold <- inv_dual_transform(model_results$specs$Threshold,
    h = conf$pop_reference, l = conf$pop_min
  )

  # Test range of thresholds around the selected one
  threshold_multipliers <- seq(0.8, 1.2, by = 0.1)
  weight_results <- data.frame(
    threshold_mult = threshold_multipliers,
    bulk_weight = NA,
    tail_weight = NA
  )

  for (i in seq_along(threshold_multipliers)) {
    test_threshold <- base_threshold * threshold_multipliers[i]

    total_data <- data[data[[conf$variables[1]]] >= model_results$specs$Cutoff, ]
    bulk_count <- sum(total_data[[conf$variables[1]]] < test_threshold)
    tail_count <- sum(total_data[[conf$variables[1]]] >= test_threshold)
    total_count <- nrow(total_data)

    weight_results$bulk_weight[i] <- bulk_count / total_count
    weight_results$tail_weight[i] <- tail_count / total_count
  }

  # Calculate weight stability metrics
  bulk_cv <- sd(weight_results$bulk_weight) / mean(weight_results$bulk_weight)
  tail_cv <- sd(weight_results$tail_weight) / mean(weight_results$tail_weight)

  list(
    weight_results = weight_results,
    bulk_cv = bulk_cv,
    tail_cv = tail_cv,
    stability_rating = ifelse(max(bulk_cv, tail_cv) < 0.2, "Stable", "Unstable")
  )
}

#' Create Component Separation Plot
#'
#' Visualizes how well the mixture model separates data into components.
#'
#' @param data Analysis data
#' @param model_results Model fitting results
#' @param conf Configuration object
#' @return ggplot object showing component separation
create_component_separation_plot <- function(data, model_results, conf) {
  threshold <- inv_dual_transform(model_results$specs$Threshold,
    h = conf$pop_reference, l = conf$pop_min
  )
  lower_cutoff <- model_results$specs$Cutoff

  # Prepare data with component labels
  plot_data <- data[data[[conf$variables[1]]] >= lower_cutoff, ] %>%
    mutate(
      component = ifelse(.data[[conf$variables[1]]] < threshold, "Bulk", "Tail")
    )

  ggplot(plot_data, aes(x = .data[[conf$variables[1]]], fill = component)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
    scale_x_log10(labels = scales::comma_format()) +
    scale_fill_manual(values = c("Bulk" = "blue", "Tail" = "red")) +
    labs(
      title = "Mixture Model Component Separation",
      subtitle = paste("Threshold at", scales::comma(threshold), "deaths"),
      x = "Deaths (log scale)",
      y = "Frequency",
      fill = "Component"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Create Weight Stability Plot
#'
#' Shows how component weights change with threshold variation.
#'
#' @param weight_stability Weight stability assessment results
#' @return ggplot object showing weight stability
create_weight_stability_plot <- function(weight_stability) {
  plot_data <- weight_stability$weight_results %>%
    tidyr::pivot_longer(
      cols = c(bulk_weight, tail_weight),
      names_to = "component", values_to = "weight"
    ) %>%
    mutate(component = stringr::str_replace(component, "_weight", ""))

  ggplot(plot_data, aes(x = threshold_mult, y = weight, color = component)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("bulk" = "blue", "tail" = "red")) +
    labs(
      title = "Component Weight Stability",
      subtitle = "Sensitivity to threshold variation",
      x = "Threshold Multiplier",
      y = "Component Weight",
      color = "Component"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Diagnose Bootstrap Stability and Uncertainty
#'
#' Evaluates the stability of model parameters and diagnostics under bootstrap
#' resampling. Critical for small sample sizes where parameter uncertainty is high.
#'
#' @param data Analysis data
#' @param model_results Model fitting results
#' @param conf Configuration object
#' @param n_bootstrap Number of bootstrap samples
#' @param generate_plots Whether to generate diagnostic plots
#' @return List containing bootstrap diagnostic results
#' @export
diagnose_bootstrap_stability <- function(data, model_results, conf,
                                         n_bootstrap = 200, generate_plots = TRUE) {
  # Storage for bootstrap results
  bootstrap_results <- data.frame(
    iteration = 1:n_bootstrap,
    xi = NA,
    beta = NA,
    bulk_mean = NA,
    tail_mean = NA,
    mixture_mean = NA,
    threshold_ok = NA,
    convergence = NA
  )

  # Original model parameters for comparison
  original_xi <- model_results$params$xi
  original_beta <- model_results$params$beta

  # Bootstrap procedure
  n <- nrow(data)

  for (i in seq_len(n_bootstrap)) {
    tryCatch(
      {
        # Bootstrap sample
        boot_indices <- sample(seq_len(n), n, replace = TRUE)
        boot_data <- data[boot_indices, ]

        # Apply dual transformation to bootstrap sample
        boot_data$dual <- dual_transform(boot_data, conf$variables[1],
          h = conf$pop_reference, l = conf$pop_min
        )

        # Fit model to bootstrap sample
        boot_results <- run_analysis(boot_data, conf$variables[1],
          year = conf$base_year,
          return_severity_draws = FALSE
        )

        # Extract parameters
        if (!is.null(boot_results$params)) {
          bootstrap_results$xi[i] <- boot_results$params$xi
          bootstrap_results$beta[i] <- boot_results$params$beta
          bootstrap_results$mixture_mean[i] <- boot_results$results$shadow
          bootstrap_results$threshold_ok[i] <- boot_results$specs$PoT >= 3
          bootstrap_results$convergence[i] <- TRUE
        } else {
          bootstrap_results$convergence[i] <- FALSE
        }
      },
      error = function(e) {
        bootstrap_results$convergence[i] <- FALSE
      }
    )
  }

  # Remove failed iterations
  successful_boots <- bootstrap_results[bootstrap_results$convergence == TRUE &
    !is.na(bootstrap_results$xi), ]

  # Calculate stability metrics
  stability_metrics <- list(
    convergence_rate = nrow(successful_boots) / n_bootstrap,
    xi_mean = mean(successful_boots$xi, na.rm = TRUE),
    xi_sd = sd(successful_boots$xi, na.rm = TRUE),
    xi_cv = sd(successful_boots$xi, na.rm = TRUE) / abs(mean(successful_boots$xi, na.rm = TRUE)),
    beta_mean = mean(successful_boots$beta, na.rm = TRUE),
    beta_sd = sd(successful_boots$beta, na.rm = TRUE),
    beta_cv = sd(successful_boots$beta, na.rm = TRUE) / abs(mean(successful_boots$beta, na.rm = TRUE)),
    xi_bias = mean(successful_boots$xi, na.rm = TRUE) - original_xi,
    beta_bias = mean(successful_boots$beta, na.rm = TRUE) - original_beta
  )

  # Assess parameter stability
  stability_assessment <- list(
    xi_stability = ifelse(stability_metrics$xi_cv < 0.3, "Stable",
      ifelse(stability_metrics$xi_cv < 0.5, "Moderate", "Unstable")
    ),
    beta_stability = ifelse(stability_metrics$beta_cv < 0.3, "Stable",
      ifelse(stability_metrics$beta_cv < 0.5, "Moderate", "Unstable")
    ),
    overall_stability = ifelse(stability_metrics$convergence_rate > 0.8, "Good", "Poor")
  )

  plots <- list()
  if (generate_plots && nrow(successful_boots) > 10) {
    plots$parameter_stability <- create_bootstrap_parameter_plot(successful_boots, original_xi, original_beta)
    plots$diagnostic_uncertainty <- create_bootstrap_uncertainty_plot(successful_boots)
  }

  list(
    bootstrap_results = successful_boots,
    stability_metrics = stability_metrics,
    stability_assessment = stability_assessment,
    plots = plots,
    bootstrap_summary = list(
      n_successful = nrow(successful_boots),
      n_attempted = n_bootstrap,
      convergence_rate = stability_metrics$convergence_rate
    )
  )
}

#' Create Bootstrap Parameter Stability Plot
#'
#' Visualizes bootstrap distributions of key model parameters.
#'
#' @param bootstrap_results Bootstrap results data frame
#' @param original_xi Original xi parameter value
#' @param original_beta Original beta parameter value
#' @return ggplot object showing parameter distributions
create_bootstrap_parameter_plot <- function(bootstrap_results, original_xi, original_beta) {
  # Prepare data for plotting
  param_data <- bootstrap_results %>%
    dplyr::select(xi, beta) %>%
    tidyr::pivot_longer(cols = c(xi, beta), names_to = "parameter", values_to = "value") %>%
    mutate(
      original_value = ifelse(parameter == "xi", original_xi, original_beta),
      parameter_label = ifelse(parameter == "xi", "Shape (ξ)", "Scale (β)")
    )

  ggplot(param_data, aes(x = value)) +
    geom_histogram(bins = 30, alpha = 0.7, fill = "skyblue", color = "black") +
    geom_vline(aes(xintercept = original_value), color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~parameter_label, scales = "free", ncol = 1) +
    labs(
      title = "Bootstrap Parameter Stability",
      subtitle = "Red line shows original estimate",
      x = "Parameter Value",
      y = "Frequency"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "bold")
    )
}

#' Create Bootstrap Uncertainty Plot
#'
#' Shows uncertainty in key model outputs from bootstrap analysis.
#'
#' @param bootstrap_results Bootstrap results data frame
#' @return ggplot object showing bootstrap uncertainty
create_bootstrap_uncertainty_plot <- function(bootstrap_results) {
  if ("mixture_mean" %in% names(bootstrap_results)) {
    ggplot(bootstrap_results, aes(x = mixture_mean)) +
      geom_histogram(bins = 30, alpha = 0.7, fill = "lightgreen", color = "black") +
      geom_vline(
        xintercept = median(bootstrap_results$mixture_mean, na.rm = TRUE),
        color = "red", linetype = "dashed", size = 1
      ) +
      labs(
        title = "Bootstrap Uncertainty in Mixture Mean",
        subtitle = "Distribution of severity estimates across bootstrap samples",
        x = "Mixture Mean Estimate",
        y = "Frequency"
      ) +
      scale_x_log10() +
      theme_bw()
  } else {
    ggplot() +
      labs(title = "Bootstrap uncertainty data not available") +
      theme_void()
  }
}

#' Diagnose Cross-Validation Robustness
#'
#' Implements leave-one-out and k-fold cross-validation adapted for small sample sizes.
#' Assesses model prediction accuracy and parameter stability.
#'
#' @param data Analysis data
#' @param variable Variable name to analyze
#' @param conf Configuration object
#' @param cv_method Cross-validation method ("leave_one_out" or "k_fold")
#' @param k Number of folds for k-fold CV (ignored for LOO)
#' @param generate_plots Whether to generate diagnostic plots
#' @return List containing cross-validation results
#' @export
diagnose_cross_validation <- function(data, variable, conf, cv_method = "leave_one_out",
                                      k = 5, generate_plots = TRUE) {
  n <- nrow(data)

  if (cv_method == "leave_one_out") {
    cv_results <- perform_loo_cv(data, variable, conf)
  } else if (cv_method == "k_fold") {
    # Adjust k for small sample size
    k_adjusted <- min(k, floor(n / 5)) # Ensure at least 5 observations per fold
    cv_results <- perform_k_fold_cv(data, variable, conf, k_adjusted)
  }

  # Calculate CV metrics
  cv_metrics <- calculate_cv_metrics(cv_results)

  plots <- list()
  if (generate_plots) {
    plots$cv_stability <- create_cv_stability_plot(cv_results)
    plots$prediction_errors <- create_prediction_error_plot(cv_results)
  }

  list(
    cv_results = cv_results,
    cv_metrics = cv_metrics,
    plots = plots,
    cv_summary = list(
      method = cv_method,
      n_folds = ifelse(cv_method == "leave_one_out", n, k),
      cv_rmse = cv_metrics$rmse,
      cv_mae = cv_metrics$mae
    )
  )
}

#' Perform Leave-One-Out Cross-Validation
#'
#' Implements LOO-CV which is appropriate for small sample sizes.
#'
#' @param data Analysis data
#' @param variable Variable name
#' @param conf Configuration object
#' @return Cross-validation results
perform_loo_cv <- function(data, variable, conf) {
  # Use configured lower cutoff for scaled deaths; fallback retains historic literal 2e5
  loo_lower_cutoff <- conf$thresholds$lower_cutoff_deaths_scaled %||% conf$thresholds$lower_cutoff %||% 2e5
  data <- data |> filter(deaths_scaled >= loo_lower_cutoff)
  n <- nrow(data)

  loo_results <- data.frame(
    fold = 1:n,
    observed = NA,
    predicted = NA,
    xi = NA,
    beta = NA,
    converged = NA
  )

  for (i in 1:n) {
    tryCatch(
      {
        # Create training set (all except observation i)
        train_data <- data[-i, ]

        # Fit model on training data
        cv_model <- run_analysis(train_data, variable,
          conf = conf, year = conf$base_year,
          return_severity_draws = FALSE
        )

        if (!is.null(cv_model$params)) {
          loo_results$xi[i] <- cv_model$params$xi
          loo_results$beta[i] <- cv_model$params$beta
          loo_results$observed[i] <- mean(train_data$deaths[train_data$deaths > loo_lower_cutoff])
          loo_results$predicted[i] <- cv_model$results$shadow
          loo_results$converged[i] <- TRUE
        } else {
          loo_results$converged[i] <- FALSE
        }
      },
      error = function(e) {
        loo_results$converged[i] <- FALSE
      }
    )
  }

  loo_results
}

#' Perform K-Fold Cross-Validation
#'
#' Implements k-fold CV adapted for small samples.
#'
#' @param data Analysis data
#' @param variable Variable name
#' @param conf Configuration object
#' @param k Number of folds
#' @return Cross-validation results
perform_k_fold_cv <- function(data, variable, conf, k) {
  n <- nrow(data)

  # Create fold assignments
  fold_assignments <- sample(rep(1:k, length.out = n))

  cv_results <- data.frame(
    fold = fold_assignments,
    observed = data[[variable]],
    predicted = NA,
    xi = NA,
    beta = NA,
    converged = NA
  )

  for (fold in 1:k) {
    tryCatch(
      {
        # Create training and test sets
        test_indices <- which(fold_assignments == fold)
        train_data <- data[-test_indices, ]

        # Apply dual transformation
        train_data$dual <- dual_transform(train_data, variable,
          h = conf$pop_reference, l = conf$pop_min
        )

        # Fit model on training data
        cv_model <- run_analysis(train_data, variable,
          year = conf$base_year,
          return_severity_draws = FALSE
        )

        if (!is.null(cv_model$params)) {
          cv_results$xi[test_indices] <- cv_model$params$xi
          cv_results$beta[test_indices] <- cv_model$params$beta
          cv_results$predicted[test_indices] <- cv_model$results$shadow
          cv_results$converged[test_indices] <- TRUE
        } else {
          cv_results$converged[test_indices] <- FALSE
        }
      },
      error = function(e) {
        test_indices <- which(fold_assignments == fold)
        cv_results$converged[test_indices] <- FALSE
      }
    )
  }

  cv_results
}

#' Calculate Cross-Validation Metrics
#'
#' Computes prediction accuracy metrics from CV results.
#'
#' @param cv_results Cross-validation results data frame
#' @return List of CV metrics
calculate_cv_metrics <- function(cv_results) {
  # Filter to successful predictions
  valid_results <- cv_results[cv_results$converged == TRUE &
    !is.na(cv_results$predicted), ]

  if (nrow(valid_results) > 0) {
    # Calculate prediction errors
    errors <- valid_results$observed - valid_results$predicted

    # Calculate metrics
    rmse <- sqrt(mean(errors^2))
    mae <- mean(abs(errors))
    mape <- mean(abs(errors / valid_results$observed)) * 100
    r_squared <- cor(valid_results$observed, valid_results$predicted)^2

    list(
      rmse = rmse,
      mae = mae,
      mape = mape,
      r_squared = r_squared,
      n_successful = nrow(valid_results),
      success_rate = nrow(valid_results) / nrow(cv_results)
    )
  } else {
    list(
      rmse = NA,
      mae = NA,
      mape = NA,
      r_squared = NA,
      n_successful = 0,
      success_rate = 0
    )
  }
}

#' Create Cross-Validation Stability Plot
#'
#' Shows parameter stability across CV folds.
#'
#' @param cv_results Cross-validation results
#' @return ggplot object showing CV stability
create_cv_stability_plot <- function(cv_results) {
  valid_results <- cv_results[cv_results$converged == TRUE &
    !is.na(cv_results$xi), ]

  if (nrow(valid_results) > 0) {
    param_data <- valid_results %>%
      dplyr::select(fold, xi, beta) %>%
      tidyr::pivot_longer(cols = c(xi, beta), names_to = "parameter", values_to = "value") %>%
      mutate(parameter_label = ifelse(parameter == "xi", "Shape (ξ)", "Scale (β)"))

    ggplot(param_data, aes(x = fold, y = value)) +
      geom_point(alpha = 0.7, color = "blue") +
      geom_smooth(method = "loess", se = TRUE, color = "red") +
      facet_wrap(~parameter_label, scales = "free_y", ncol = 1) +
      labs(
        title = "Cross-Validation Parameter Stability",
        subtitle = "Parameter estimates across CV folds",
        x = "CV Fold",
        y = "Parameter Value"
      ) +
      theme_bw() +
      theme(strip.text = element_text(size = 12, face = "bold"))
  } else {
    ggplot() +
      labs(title = "Insufficient successful CV iterations") +
      theme_void()
  }
}

#' Create Prediction Error Plot
#'
#' Visualizes prediction accuracy from cross-validation.
#'
#' @param cv_results Cross-validation results
#' @return ggplot object showing prediction errors
create_prediction_error_plot <- function(cv_results) {
  valid_results <- cv_results[cv_results$converged == TRUE &
    !is.na(cv_results$predicted), ]

  if (nrow(valid_results) > 0) {
    ggplot(valid_results, aes(x = observed, y = predicted)) +
      geom_point(alpha = 0.7, color = "blue") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = TRUE, color = "green") +
      labs(
        title = "Cross-Validation Prediction Accuracy",
        subtitle = "Observed vs Predicted Values",
        x = "Observed Values",
        y = "Predicted Values"
      ) +
      theme_bw()
  } else {
    ggplot() +
      labs(title = "Insufficient prediction data") +
      theme_void()
  }
}

#' Create Diagnostic Summary Table
#'
#' Generates a publication-ready table summarizing all diagnostic results.
#'
#' @param threshold_diag Threshold diagnostics results
#' @param gof_results Goodness-of-fit results
#' @param mixture_diag Mixture model diagnostics
#' @param bootstrap_diag Bootstrap diagnostics
#' @param cv_results Cross-validation results
#' @return Formatted diagnostic summary table
#' @export
create_diagnostic_summary_table <- function(threshold_diag, gof_results, mixture_diag,
                                            bootstrap_diag, cv_results) {
  # Create summary table
  summary_table <- data.frame(
    Diagnostic_Category = c(
      "Threshold Selection",
      "Goodness-of-Fit",
      "Component Separation",
      "Bootstrap Stability",
      "Cross-Validation"
    ),
    Test_Metric = c(
      paste("Stability Score:", round(threshold_diag$stability_metrics$overall_stability, 3)),
      paste("Tests Passed:", gof_results$gof_summary$n_adequate, "/", gof_results$gof_summary$n_tests),
      paste("Separation Quality:", mixture_diag$component_summary$separation_quality),
      paste("Convergence Rate:", round(bootstrap_diag$bootstrap_summary$convergence_rate, 3)),
      paste("Prediction R²:", round(cv_results$cv_metrics$r_squared, 3))
    ),
    Assessment = c(
      threshold_diag$stability_metrics$stability_rating,
      gof_results$gof_summary$overall_adequacy,
      mixture_diag$component_summary$separation_quality,
      bootstrap_diag$stability_assessment$overall_stability,
      ifelse(cv_results$cv_metrics$r_squared > 0.7, "Good", "Moderate")
    ),
    Details = c(
      paste("Tested", threshold_diag$diagnostics_summary$n_tested, "thresholds"),
      paste("Statistical tests:", paste(gof_results$test_results$test_name, collapse = ", ")),
      paste(
        "Bulk:", mixture_diag$component_summary$bulk_n,
        "| Tail:", mixture_diag$component_summary$tail_n
      ),
      paste("CV(ξ):", round(bootstrap_diag$stability_metrics$xi_cv, 3)),
      paste("RMSE:", round(cv_results$cv_metrics$rmse, 0))
    )
  )

  # Format as publication table
  formatted_table <- summary_table %>%
    dplyr::rename(
      "Diagnostic Category" = Diagnostic_Category,
      "Test Metric" = Test_Metric
    )

  return(formatted_table)
}

#' Create Threshold Selection Table
#'
#' Publication-ready table for threshold selection diagnostics.
#'
#' @param threshold_diagnostics Threshold diagnostics results
#' @param data_characteristics List with data characteristics
#' @return Formatted threshold selection table
#' @export
create_threshold_selection_table <- function(threshold_diagnostics, data_characteristics) {
  # Get top threshold candidates
  thresh_results <- threshold_diagnostics$threshold_results

  # Select top 5 thresholds by balance of sample size and stability
  if (nrow(thresh_results) >= 5) {
    top_thresholds <- thresh_results[1:5, ]
  } else {
    top_thresholds <- thresh_results
  }

  # Create table
  threshold_table <- top_thresholds %>%
    dplyr::select(threshold, n_exceedances, xi, xi_se, beta, beta_se) %>%
    mutate(
      threshold = round(threshold, 0),
      xi = round(xi, 4),
      xi_se = round(xi_se, 4),
      beta = round(beta, 0),
      beta_se = round(beta, 0),
      xi_ci = paste0(round(xi - 1.96 * xi_se, 4), " to ", round(xi + 1.96 * xi_se, 4))
    ) %>%
    dplyr::select(
      "Threshold" = threshold,
      "Exceedances" = n_exceedances,
      "Shape (ξ)" = xi,
      "Shape SE" = xi_se,
      "Scale (β)" = beta,
      "Scale SE" = beta_se,
      "ξ 95% CI" = xi_ci
    )

  return(threshold_table)
}

#' Assess Overall Model Adequacy
#'
#' Provides a detailed assessment of model adequacy based on all diagnostics.
#'
#' @param threshold_diag Threshold diagnostics
#' @param gof_results Goodness-of-fit results
#' @param mixture_diag Mixture diagnostics
#' @param bootstrap_diag Bootstrap diagnostics
#' @param cv_results Cross-validation results
#' @param model_results Original model results
#' @return Overall adequacy assessment
#' @export
assess_overall_model_adequacy <- function(threshold_diag, gof_results, mixture_diag,
                                          bootstrap_diag, cv_results, model_results) {
  # Individual component ratings
  threshold_rating <- threshold_diag$stability_metrics$stability_rating
  gof_rating <- gof_results$gof_summary$overall_adequacy
  mixture_rating <- mixture_diag$component_summary$separation_quality
  bootstrap_rating <- bootstrap_diag$stability_assessment$overall_stability
  cv_rating <- ifelse(cv_results$cv_metrics$r_squared > 0.7, "Good",
    ifelse(cv_results$cv_metrics$r_squared > 0.5, "Moderate", "Poor")
  )

  # Overall rating logic
  good_ratings <- sum(c(
    threshold_rating, gof_rating, mixture_rating,
    bootstrap_rating, cv_rating
  ) %in% c("Good", "Adequate"))

  overall_rating <- ifelse(good_ratings >= 4, "Excellent",
    ifelse(good_ratings >= 3, "Good",
      ifelse(good_ratings >= 2, "Adequate", "Poor")
    )
  )

  # Generate recommendations
  recommendations <- c()

  if (threshold_rating %in% c("Poor", "Moderate")) {
    recommendations <- c(recommendations, "Consider alternative threshold selection methods")
  }

  if (gof_rating == "Needs attention") {
    recommendations <- c(recommendations, "Investigate model fit adequacy with alternative distributions")
  }

  if (mixture_rating == "Suboptimal") {
    recommendations <- c(recommendations, "Review component separation and threshold choice")
  }

  if (bootstrap_rating == "Poor") {
    recommendations <- c(recommendations, "Parameter estimates show high uncertainty - consider data augmentation")
  }

  if (cv_rating == "Poor") {
    recommendations <- c(recommendations, "Poor predictive performance - model may be overfitting")
  }

  if (model_results$specs$PoT < 10) {
    recommendations <- c(recommendations, "Very few tail observations - consider lower threshold or data pooling")
  }

  if (length(recommendations) == 0) {
    recommendations <- "Model diagnostics indicate adequate performance for pandemic risk estimation"
  }

  list(
    overall_rating = overall_rating,
    threshold_rating = threshold_rating,
    gof_rating = gof_rating,
    mixture_rating = mixture_rating,
    bootstrap_rating = bootstrap_rating,
    cv_rating = cv_rating,
    recommendations = recommendations,
    summary_statistics = list(
      n_total = model_results$specs$n,
      n_tail = model_results$specs$PoT,
      tail_percentage = model_results$specs$PcoT,
      xi_estimate = model_results$params$xi,
      xi_se = model_results$params$xi_se
    )
  )
}
