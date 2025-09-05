#' Plot Empirically-Calibrated Two-Component Mixture Model with Bootstrap Parameter Uncertainty
#'
#' Creates a visualization of a two-stage mixture model fitted to extreme value data,
#' showing both the model fit and bootstrap-derived parameter uncertainty bands.
#' This represents an empirically-calibrated mixture model where the bulk component
#' (truncated lognormal) is combined with a GPD tail component, with uncertainty
#' bands derived from bootstrap resampling.
#'
#' @details
#' ## Methodological Framework
#' 
#' This function visualizes a **calibrated mixture model fit with parameter uncertainty**:
#' 
#' ### Model Structure:
#' - **Bulk Component**: Truncated lognormal distribution for moderate events (lower_cutoff ≤ x ≤ threshold)
#' - **Tail Component**: Generalized Pareto Distribution (GPD) for extreme events (x > threshold)
#' - **Empirical Calibration**: Both components are scaled to match observed exceedance frequencies
#' 
#' ### Uncertainty Representation:
#' The uncertainty bands show **parameter uncertainty** derived from bootstrap resampling:
#' - Each bootstrap sample generates new parameter estimates (μ, σ for bulk; ξ, β for tail)
#' - Uncertainty bands represent 95% quantiles across bootstrap parameter draws
#' - **Fixed calibration**: All bootstrap draws use the same empirical scaling factor
#' 
#' ### Important Interpretation Notes:
#' - **Parameter uncertainty**: Shows variability in model parameters given the observed data
#' - **NOT full predictive uncertainty**: Does not include calibration/scaling uncertainty
#' - **Appropriate for model validation**: Suitable for assessing model fit against observed data
#' - **Continuous at threshold**: Uncertainty bands are mathematically connected at the threshold
#' 
#' ### Mathematical Details:
#' - **Scaling factor (α)**: Optimized to minimize squared error between model and empirical survival
#' - **Bulk uncertainty**: S_bulk,i(x) = α_fixed × (1 - F_lognormal(x; μ_i, σ_i))
#' - **Tail uncertainty**: S_tail,i(x) = S_bulk,i(threshold) × (1 - F_GPD(x; ξ_i, β_i, threshold))
#' - **Continuity**: Each bootstrap draw maintains S_bulk,i(threshold) = S_tail,i(threshold)
#'
#' @param data Data frame containing the empirical data points for visualization
#' @param variable String, name of the variable to plot (typically "deaths")
#' @param mixture_fit Object, fitted mixture model containing bulk and tail components
#' @param pop_ref Numeric, population reference for dual transformation
#' @param title String, optional plot title
#' @param y_title String, y-axis label. Default is "Probability(Pandemic deaths > X)"
#' @param tail_limit Numeric, lower limit for display. Default is NULL (will use lower_cutoff)
#' @param log10 Logical, whether to use logarithmic scale for y-axis. Default is FALSE
#' @param y_limit Numeric, optional limit for y-axis
#' @param lancet_style Logical, whether to apply The Lancet journal style guidelines. Default is TRUE
#' @param show_means Logical, whether to show the bulk and tail conditional means. Default is FALSE
#' @param shade_regions Logical, whether to shade the bulk and tail regions. Default is TRUE
#' @param shade_alpha Numeric, transparency level for shaded regions (0-1). Default is 0.1
#' @param expected_severity Numeric, expected severity value to highlight. Default is NULL
#' @param param_draws List, bootstrap parameter draws for uncertainty bands. Must contain:
#'   - xi, beta: vectors of GPD parameter draws
#'   - bulk_params: list with distribution-specific parameter draws (e.g., meanlog, sdlog)
#' @param show_uncertainty Logical, whether to display uncertainty bands. Default is TRUE
#' @param ci_level Numeric, confidence interval level for uncertainty bands (0-1). Default is 0.95
#' @param ci_alpha Numeric, transparency level for uncertainty bands (0-1). Default is 0.2
#' @param point_size Numeric, size of empirical data points. Default is 3
#' @param severity_marker_height Numeric, height of expected severity vertical marker. Default is 0.015
#' 
#' @return A ggplot object representing the calibrated mixture model with parameter uncertainty
#' 
#' @note This visualization is designed for model validation and parameter uncertainty assessment.
#'   For predictive uncertainty that includes calibration uncertainty, consider additional
#'   approaches that account for scaling factor variability across bootstrap samples.
#'   
#' @references 
#' Cirillo, P., & Taleb, N. N. (2020). Tail risk of contagious diseases. Nature Physics, 16(6), 606-613.
#' 
#' @export
truncated_bulk_tail_plot <- function(data, variable, mixture_fit, pop_ref,
                                     title = NULL, y_title = "Probability(Pandemic deaths > X)",
                                     tail_limit = NULL, log10 = FALSE, y_limit = NULL,
                                     lancet_style = TRUE, show_means = FALSE,
                                     shade_regions = TRUE, shade_alpha = 0.1,
                                     expected_severity = NULL,
                                     param_draws = NULL, show_uncertainty = TRUE, 
                                     ci_level = 0.95, ci_alpha = 0.2, point_size = 3,
                                     severity_marker_height = 0.015,
                                     output_dir = NULL, conf = NULL) {
  # Extract parameters from mixture fit
  lower_cutoff <- mixture_fit$lower_cutoff # Lower threshold
  u_best <- mixture_fit$tail_threshold # Upper threshold (GPD threshold)

  # Bulk parameters
  dist_type <- mixture_fit$bulk_dist_type
  params <- mixture_fit$bulk$parameters

  # Function name mapping for R's distribution functions
  dist_name_map <- c(
    "lognormal" = "lnorm",
    "weibull" = "weibull",
    "gamma" = "gamma"
  )

  # Get the correct distribution function name
  dist_func_name <- dist_name_map[dist_type]
  if (is.na(dist_func_name)) dist_func_name <- dist_type

  # Define extension factor for plotting
  extend <- 1.5

  # Maximum plot range
  plot_max <- max(data[[variable]]) * extend

  # If tail_limit is not specified, use lower_cutoff
  if (is.null(tail_limit)) {
    tail_limit <- lower_cutoff
  }

  # Empirical data: sort, compute empirical survival
  all_data <- data %>%
    arrange(.data[[variable]]) %>%
    mutate(
      x = .data[[variable]],
      y = 1 - ecdf(.data$x)(.data$x)
    )

  # Avoid zero probability at minimum
  all_data$y[which.min(all_data$y)] <- 0.006

  # For compatibility with the rest of the function
  empirical_data <- all_data

  # Extract tail parameters
  xi <- mixture_fit$tail$xi
  beta <- mixture_fit$tail$beta
  fit <- mixture_fit$tail$fit
  shadow <- mixture_fit$tail$shadow

  # Generate x-values for bulk component (only up to the threshold)
  bulk_x <- seq(lower_cutoff, u_best, length.out = 500)

  # Calculate survival function for bulk distribution up to threshold
  p_fn <- match.fun(paste0("p", dist_func_name))
  bulk_surv <- 1 - do.call(p_fn, c(list(q = bulk_x), params))

  # Initial scaling to match the empirical probability at lower_cutoff
  emp_p_lower <- mean(empirical_data$x >= lower_cutoff)
  scale_factor <- emp_p_lower / (1 - do.call(p_fn, c(list(q = lower_cutoff), params)))

  # Improve scaling by optimizing across the entire bulk range
  emp_data_bulk <- empirical_data %>%
    filter(x >= lower_cutoff & x <= u_best)

  if (nrow(emp_data_bulk) > 10) {
    # Create function to evaluate error at different scaling factors
    evaluate_scaling_error <- function(sf) {
      predicted_surv <- sapply(emp_data_bulk$x, function(x_val) {
        (1 - do.call(p_fn, c(list(q = x_val), params))) * sf
      })
      # Mean squared error between predicted and empirical survival
      mean((predicted_surv - emp_data_bulk$y)^2)
    }

    # Find optimal scaling factor that minimizes overall error
    opt_result <- optimize(
      f = evaluate_scaling_error,
      interval = c(0.5 * scale_factor, 1.5 * scale_factor)
    )

    # Use the improved scaling factor
    improved_scale <- opt_result$minimum
    scale_factor <- improved_scale
  }

  # Apply the final scaling factor
  bulk_surv <- bulk_surv * scale_factor

  # Generate x-values for tail component (starting at threshold)
  tail_x <- seq(u_best, plot_max, length.out = 500)

  # Calculate fitted curve probabilities for tail (GPD)
  prob_less_thresh <- fit$p.less.thresh
  tail_surv <- (1 - prob_less_thresh) * (1 - evir::pgpd(tail_x, xi, u_best, beta))

  # Create data frames for bulk and tail curves
  bulk_data <- data.frame(x = bulk_x, y = bulk_surv)
  tail_data <- data.frame(x = tail_x, y = tail_surv)
  
  # =============================================================================
  # BOOTSTRAP PARAMETER UNCERTAINTY CALCULATION
  # =============================================================================
  # This section calculates uncertainty bands using bootstrap parameter draws.
  # Key methodological approach: FIXED CALIBRATION with PARAMETER UNCERTAINTY
  # 
  # Rationale: We show how parameter uncertainty affects the model fit while
  # keeping the empirical calibration (scaling factor) fixed. This represents
  # "parameter uncertainty given the observed data" rather than full predictive
  # uncertainty which would include calibration uncertainty.
  # =============================================================================
  
  bulk_uncertainty <- NULL
  tail_uncertainty <- NULL
  
  if (!is.null(param_draws) && show_uncertainty && length(param_draws$xi) > 1) {
    # Calculate confidence interval bounds
    alpha <- 1 - ci_level
    lower_q <- alpha / 2
    upper_q <- 1 - alpha / 2
    
    # -------------------------------------------------------------------------
    # BULK UNCERTAINTY BANDS: Parameter uncertainty with fixed calibration
    # -------------------------------------------------------------------------
    # For each bootstrap draw, we:
    # 1. Use the bootstrap lognormal parameters (μ_i, σ_i) 
    # 2. Apply the SAME scaling factor (α) from the point estimate
    # 3. This shows "if parameters were μ_i, σ_i, what would the survival curve look like?"
    #    while maintaining calibration to the observed data
    # -------------------------------------------------------------------------
    
    if (!is.null(param_draws$bulk_params) && length(param_draws$bulk_params) > 0) {
      n_draws <- length(param_draws$xi)
      
      # Matrix to store survival probabilities: rows = x-values, columns = bootstrap draws
      bulk_surv_draws <- matrix(NA, nrow = length(bulk_x), ncol = n_draws)
      
      # Calculate survival probability for each bootstrap parameter set
      for (i in 1:n_draws) {
        # Extract bulk distribution parameters for this bootstrap draw
        if (dist_type == "lognormal" && 
            "meanlog" %in% names(param_draws$bulk_params) && 
            "sdlog" %in% names(param_draws$bulk_params)) {
          
          meanlog_i <- param_draws$bulk_params$meanlog[i]
          sdlog_i <- param_draws$bulk_params$sdlog[i]
          
          if (!is.na(meanlog_i) && !is.na(sdlog_i) && sdlog_i > 0) {
            # Calculate raw survival probabilities using bootstrap parameters
            surv_i <- 1 - plnorm(bulk_x, meanlog_i, sdlog_i)
            
            # CRITICAL: Apply the FIXED scaling factor from point estimate
            # This represents parameter uncertainty while holding calibration constant
            # Formula: S_bulk,i(x) = α_fixed × (1 - F_lognormal(x; μ_i, σ_i))
            surv_i <- surv_i * scale_factor
            
            bulk_surv_draws[, i] <- surv_i
          }
        } else if (dist_type == "weibull" && 
                   "shape" %in% names(param_draws$bulk_params) && 
                   "scale" %in% names(param_draws$bulk_params)) {
          
          shape_i <- param_draws$bulk_params$shape[i]
          scale_i <- param_draws$bulk_params$scale[i]
          
          if (!is.na(shape_i) && !is.na(scale_i) && shape_i > 0 && scale_i > 0) {
            # Calculate raw survival probabilities using bootstrap Weibull parameters
            surv_i <- 1 - pweibull(bulk_x, shape_i, scale_i)
            
            # Apply the FIXED scaling factor (same approach as lognormal)
            # Formula: S_bulk,i(x) = α_fixed × (1 - F_weibull(x; α_i, β_i))
            surv_i <- surv_i * scale_factor
            
            bulk_surv_draws[, i] <- surv_i
          }
        }
      }
      
      # Calculate uncertainty bands as quantiles across bootstrap draws
      # This gives us the range of survival curves due to parameter uncertainty
      if (sum(!is.na(bulk_surv_draws[1, ])) > 2) {  # Need at least 3 valid draws
        # For each x-value, compute the specified quantiles across all bootstrap draws
        bulk_lower <- apply(bulk_surv_draws, 1, quantile, probs = lower_q, na.rm = TRUE)
        bulk_upper <- apply(bulk_surv_draws, 1, quantile, probs = upper_q, na.rm = TRUE)
        
        bulk_uncertainty <- data.frame(
          x = bulk_x,
          y_lower = bulk_lower,
          y_upper = bulk_upper
        )
      }
    }
    
    # -------------------------------------------------------------------------
    # TAIL UNCERTAINTY BANDS: Parameter uncertainty with enforced continuity
    # -------------------------------------------------------------------------
    # For the tail component, we:
    # 1. Use bootstrap GPD parameters (ξ_i, β_i) for each draw
    # 2. Calculate where each bulk curve ends at the threshold
    # 3. Scale each tail curve to start exactly where its bulk curve ends
    # 4. This ensures mathematical continuity: S_bulk,i(threshold) = S_tail,i(threshold)
    # -------------------------------------------------------------------------
    
    if (length(param_draws$xi) > 2 && length(param_draws$beta) > 2) {
      n_draws <- length(param_draws$xi)
      tail_surv_draws <- matrix(NA, nrow = length(tail_x), ncol = n_draws)
      
      # Step 1: Calculate where each bulk curve ends at the threshold
      # This ensures continuity between bulk and tail uncertainty bands
      bulk_at_threshold_draws <- numeric(n_draws)
      
      # For each bootstrap draw, calculate S_bulk,i(threshold)
      for (i in 1:n_draws) {
        if (dist_type == "lognormal" && 
            "meanlog" %in% names(param_draws$bulk_params) && 
            "sdlog" %in% names(param_draws$bulk_params)) {
          
          meanlog_i <- param_draws$bulk_params$meanlog[i]
          sdlog_i <- param_draws$bulk_params$sdlog[i]
          
          if (!is.na(meanlog_i) && !is.na(sdlog_i) && sdlog_i > 0) {
            # Calculate survival probability at threshold using bootstrap parameters + fixed scaling
            # Formula: S_bulk,i(threshold) = α_fixed × (1 - F_lognormal(threshold; μ_i, σ_i))
            bulk_at_threshold_draws[i] <- (1 - plnorm(u_best, meanlog_i, sdlog_i)) * scale_factor
          }
        } else if (dist_type == "weibull" && 
                   "shape" %in% names(param_draws$bulk_params) && 
                   "scale" %in% names(param_draws$bulk_params)) {
          
          shape_i <- param_draws$bulk_params$shape[i]
          scale_i <- param_draws$bulk_params$scale[i]
          
          if (!is.na(shape_i) && !is.na(scale_i) && shape_i > 0 && scale_i > 0) {
            # Calculate survival probability at threshold using bootstrap Weibull parameters + fixed scaling
            # Formula: S_bulk,i(threshold) = α_fixed × (1 - F_weibull(threshold; α_i, β_i))
            bulk_at_threshold_draws[i] <- (1 - pweibull(u_best, shape_i, scale_i)) * scale_factor
          }
        }
      }
      
      # Step 2: Calculate tail survival curves with enforced continuity
      for (i in 1:n_draws) {
        xi_i <- param_draws$xi[i]
        beta_i <- param_draws$beta[i]
        
        if (!is.na(xi_i) && !is.na(beta_i) && beta_i > 0 && !is.na(bulk_at_threshold_draws[i])) {
          # Calculate raw GPD survival probabilities using bootstrap parameters
          # Note: GPD survival at threshold = 1 by definition
          gpd_surv_i <- 1 - evir::pgpd(tail_x, xi_i, u_best, beta_i)
          
          # CRITICAL: Scale to ensure continuity with bulk component
          # Formula: S_tail,i(x) = S_bulk,i(threshold) × S_GPD(x; ξ_i, β_i, threshold)
          # This ensures S_tail,i(threshold) = S_bulk,i(threshold) for each draw
          tail_surv_i <- bulk_at_threshold_draws[i] * gpd_surv_i
          tail_surv_draws[, i] <- tail_surv_i
        }
      }
      
      # Step 3: Calculate uncertainty bands as quantiles across bootstrap draws
      # Result: continuous uncertainty bands that connect smoothly at threshold
      if (sum(!is.na(tail_surv_draws[1, ])) > 2) {  # Need at least 3 valid draws
        # For each x-value in tail region, compute quantiles across all bootstrap draws
        tail_lower <- apply(tail_surv_draws, 1, quantile, probs = lower_q, na.rm = TRUE)
        tail_upper <- apply(tail_surv_draws, 1, quantile, probs = upper_q, na.rm = TRUE)
        
        tail_uncertainty <- data.frame(
          x = tail_x,
          y_lower = tail_lower,
          y_upper = tail_upper
        )
      }
    }
  }

  # Labels for thresholds and means
  lower_cutoff_label <- paste0(
    scales::comma(round(lower_cutoff)),
    " deaths"
  )

  threshold_label <- paste0(
    scales::comma(round(inv_dual_transform(u_best, h = pop_ref, l = 1e4))),
    " deaths"
  )

  shadow_label <- paste0(scales::comma(round(shadow)), " deaths")
  bulk_mean_label <- paste0(
    scales::comma(round(mixture_fit$bulk$conditional_mean)),
    " deaths"
  )

  # Create label for the expected pandemic severity (use parameter if provided, otherwise use mixture mean)
  severity_value <- ifelse(is.null(expected_severity), mixture_fit$mean, expected_severity)
  expected_severity_label <- paste0(
    "Expected Pandemic\nSeverity\n",
    scales::comma(round(severity_value)),
    " deaths"
  )

  # Define Lancet color palette
  lancet_colors <- c(
    "Primary" = "#00468B", # Primary blue
    "Highlight" = "#ED0000", # Signal red
    "Secondary" = "#42B540", # Secondary green
    "Tertiary" = "#0099B4", # Tertiary blue
    "Quaternary" = "#925E9F" # Quaternary purple
  )

  # Create bulk model label based on distribution type
  bulk_model_label <- "Log-normal model"
  gpd_model_label <- "GPD model"

  # Define line weights per Lancet guidelines (≤ 0.5 pt)
  lancet_line_width <- 1

  # Construct the plot - show all data points, don't filter by tail_limit
  p <- all_data %>%
    ggplot() +

    # Optional shaded regions for bulk and tail
    {
      if (shade_regions) {
        geom_ribbon(
          data = bulk_data,
          aes(x = .data$x, ymin = 0, ymax = .data$y, fill = bulk_model_label),
          alpha = shade_alpha
        )
      }
    } +
    {
      if (shade_regions) {
        geom_ribbon(
          data = tail_data,
          aes(x = .data$x, ymin = 0, ymax = .data$y, fill = "GPD model"),
          alpha = shade_alpha
        )
      }
    } +

    # Uncertainty bands for bulk component
    {
      if (!is.null(bulk_uncertainty) && show_uncertainty) {
        geom_ribbon(
          data = bulk_uncertainty,
          aes(x = .data$x, ymin = .data$y_lower, ymax = .data$y_upper),
          fill = lancet_colors[["Quaternary"]], # Purple for log-normal
          alpha = ci_alpha
        )
      }
    } +
    
    # Uncertainty bands for tail component  
    {
      if (!is.null(tail_uncertainty) && show_uncertainty) {
        geom_ribbon(
          data = tail_uncertainty,
          aes(x = .data$x, ymin = .data$y_lower, ymax = .data$y_upper),
          fill = lancet_colors[["Secondary"]], # Green for GPD
          alpha = ci_alpha
        )
      }
    } +

    # Empirical points
    geom_point(aes(x = .data$x, y = y), size = point_size, alpha = 0.2) +

    # Bulk distribution curve (truncated at threshold)
    geom_line(
      data = bulk_data,
      aes(x = .data$x, y = y, color = bulk_model_label),
      linewidth = if (lancet_style) lancet_line_width else 1
    ) +

    # GPD tail curve
    geom_line(
      data = tail_data,
      aes(x = .data$x, y = y, color = "GPD model"),
      linewidth = if (lancet_style) lancet_line_width else 1
    ) +

    # Vertical lines for thresholds and means - use thin lines (≤ 0.25 pt) for reference lines
    geom_vline(
      aes(xintercept = lower_cutoff, color = "Lower Cutoff"),
      linetype = "dashed",
      linewidth = if (lancet_style) 0.25 else 0.8
    ) +
    geom_vline(
      aes(xintercept = u_best, color = "Tail Threshold"),
      linetype = "dashed",
      linewidth = if (lancet_style) 0.25 else 0.8
    ) +
    # Add conditional mean lines only if show_means is TRUE
    {
      if (show_means) {
        geom_vline(
          aes(xintercept = mixture_fit$bulk$conditional_mean, color = "Bulk Mean"),
          linetype = "dashed",
          linewidth = if (lancet_style) 0.25 else 0.8
        )
      }
    } +
    {
      if (show_means) {
        geom_vline(
          aes(xintercept = shadow, color = "Tail Mean"),
          linetype = "dashed",
          linewidth = if (lancet_style) 0.25 else 0.8
        )
      }
    } +

    # Log10 scale on x-axis with human-readable tick marks (multiples of 5 or 10)
    scale_x_log10(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)
    ) +

    # Axis and legend labels with proper SI units
    labs(
      title = title,
      x = "Population-adjusted pandemic-related deaths",
      y = y_title,
      color = NULL, # Remove legend title per Lancet style
      fill = NULL # Remove fill legend title
    ) +

    # Log ticks with outward pointing ticks
    # annotation_logticks(sides = "b", outside = TRUE) +

    # Threshold text annotations (always show)
    # Use smaller 8pt sans-serif font size for annotations
    annotate(
      "text",
      x = lower_cutoff * 1.2, y = 1.1,
      label = lower_cutoff_label,
      color = if (lancet_style) lancet_colors["Quaternary"] else "darkred",
      hjust = 0,
      size = if (lancet_style) 2.8 else 3,
      family = if (lancet_style) "Arial" else ""
    ) +
    annotate(
      "text",
      x = u_best * 1.2, y = 1.1,
      label = threshold_label,
      color = if (lancet_style) lancet_colors["Secondary"] else "darkgreen",
      hjust = 0,
      size = if (lancet_style) 2.8 else 3,
      family = if (lancet_style) "Arial" else ""
    ) +
    # Mean annotations (only show if show_means is TRUE)
    {
      if (show_means) {
        annotate(
          "text",
          x = mixture_fit$bulk$conditional_mean * 1.2, y = 0.85,
          label = bulk_mean_label,
          color = if (lancet_style) lancet_colors["Secondary"] else "purple",
          hjust = 0,
          size = if (lancet_style) 2.8 else 3,
          family = if (lancet_style) "Arial" else ""
        )
      }
    } +
    {
      if (show_means) {
        annotate(
          "text",
          x = shadow * 1.2, y = 0.8,
          label = shadow_label,
          color = if (lancet_style) lancet_colors["Tertiary"] else "darkblue",
          hjust = 0,
          size = if (lancet_style) 2.8 else 3,
          family = if (lancet_style) "Arial" else ""
        )
      }
    } +

    # Repel labels for extreme points
    ggrepel::geom_text_repel(
      data = subset(all_data, all_data$x > shadow),
      aes(x = .data$x, y = y, label = .data$name),
      size = if (lancet_style) 2 else 3,
      family = if (lancet_style) "Arial" else "",
      box.padding = 0.5,
      max.overlaps = 40,
      color = if (lancet_style) "#666666" else "black",
      segment.size = if (lancet_style) 0.2 else 0.5,
      segment.color = if (lancet_style) "#999999" else "black",
      min.segment.length = 0
    ) +

    # Add lollipop marker for expected pandemic severity (only if the severity value is within the plot range)
    {
      if (severity_value >= lower_cutoff && severity_value <= plot_max) {
        geom_segment(
          aes(x = severity_value, xend = severity_value, y = 0, yend = severity_marker_height),
          linewidth = 0.5,
          color = lancet_colors["Highlight"],
          linetype = 1,
          alpha = 0.3
        )
      }
    } +
    {
      if (severity_value >= lower_cutoff && severity_value <= plot_max) {
        geom_point(
          aes(x = severity_value, y = severity_marker_height),
          size = 3,
          color = lancet_colors["Highlight"],
          alpha = 0.3
        )
      }
    } +
    {
      if (severity_value >= lower_cutoff && severity_value <= plot_max) {
        annotate(
          "text",
          x = severity_value,
          y = severity_marker_height * 1.33,
          label = expected_severity_label,
          color = lancet_colors["Highlight"],
          hjust = 0.5,
          vjust = 0,
          size = if (lancet_style) 2 else 3,
          family = if (lancet_style) "Arial" else "",
          fontface = "bold",
          alpha = 0.7
        )
      }
    }

  # Apply Lancet theme if requested
  # Create dynamic color mapping using the actual distribution type - basic colors
  lancet_color_values <- c(
    "Lower Cutoff" = lancet_colors["Primary"],
    "Tail Threshold" = lancet_colors["Highlight"],
    "GPD Tail" = lancet_colors["Highlight"]
  )

  # Add mean-related colors only if show_means is TRUE
  if (show_means) {
    lancet_color_values["Bulk Mean"] <- lancet_colors["Secondary"]
    lancet_color_values["Tail Mean"] <- lancet_colors["Tertiary"]
  }
  # Add the specific bulk distribution color dynamically
  lancet_color_values[bulk_model_label] <- lancet_colors["Quaternary"]
  lancet_color_values[gpd_model_label] <- lancet_colors["Secondary"]
  if (lancet_style) {
    p <- p +
      theme_minimal() +
      theme(
        # Typography per Lancet guidelines
        plot.title = element_text(family = "Times New Roman", size = 10, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(family = "Arial", size = 9, hjust = 0.5, margin = margin(b = 20)),
        plot.caption = element_text(family = "Arial", size = 8, color = "gray30", hjust = 0, margin = margin(t = 10)),

        # Legend styling
        legend.position = "bottom",
        legend.title = element_text(family = "Arial", size = 9, face = "bold"),
        legend.text = element_text(family = "Arial", size = 8),
        legend.key.size = unit(1.2, "lines"),
        legend.key = element_rect(fill = "white", color = NA),

        # Axis styling per Lancet guidelines
        axis.title.x = element_text(family = "Arial", size = 11, margin = margin(t = 10)),
        axis.title.y = element_text(family = "Arial", size = 11, margin = margin(r = 10)),
        axis.text = element_text(family = "Arial", size = 8),
        axis.ticks = element_line(linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),

        # Grid lines - removed per Lancet guidelines
        panel.grid = element_blank(),

        # Background
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(20, 20, 20, 20)
      ) +


      # Use The Lancet's color palette with dynamically constructed values
      scale_color_manual(values = lancet_color_values) +
      # Use dynamic fill scale to match the model labels exactly
      {
        fill_values <- c()
        fill_values[bulk_model_label] <- lancet_colors["Quaternary"]
        fill_values["GPD model"] <- lancet_colors["Secondary"]
        scale_fill_manual(values = fill_values)
      }
  } else {
    # Original theme styling
    p <- p +
      theme_bw() +
      theme(
        axis.text.x = element_text(margin = margin(t = 10), size = 12),
        axis.text.y = element_text(margin = margin(r = 10), size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom"
      ) +
      # Create dynamic color mapping for non-Lancet style - basic colors
      regular_color_values <- c(
      "Lower Cutoff" = "darkred",
      "Tail Threshold" = "darkgreen",
      "GPD Tail" = "red"
    )

    # Add mean-related colors only if show_means is TRUE
    if (show_means) {
      regular_color_values["Bulk Mean"] <- "purple"
      regular_color_values["Tail Mean"] <- "darkblue"
    }

    # Add the specific bulk distribution color dynamically
    bulk_colors <- c("lognormal" = "orange", "weibull" = "darkorange", "gamma" = "goldenrod")
    default_bulk_color <- "orange"
    regular_color_values[bulk_model_label] <- bulk_colors[tolower(dist_type)]
    if (is.na(regular_color_values[bulk_model_label])) {
      regular_color_values[bulk_model_label] <- default_bulk_color
    }

    scale_color_manual(values = regular_color_values) + {
      fill_values <- c()
      fill_values[bulk_model_label] <- bulk_colors[tolower(dist_type)]
      if (is.na(fill_values[bulk_model_label])) {
        fill_values[bulk_model_label] <- default_bulk_color
      }
      fill_values["GPD model"] <- "red"
      scale_fill_manual(values = fill_values)
    }
  }

  # Optionally use log10 scale on y-axis
  if (log10) {
    p <- p + scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  }

  # Apply y-axis limit if specified
  if (!is.null(y_limit)) {
    p <- p + coord_cartesian(lim = c(0, y_limit))
  }

  # Save plot if output directory is provided
  save_plot_if_enabled(p, "mixture_model_plot", output_dir, conf, width = 12, height = 8)
  
  print(p)
  
  invisible(p)
}
