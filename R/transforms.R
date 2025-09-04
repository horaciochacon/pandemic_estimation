#' Apply Dual Transformation to a Variable
#'
#' Transforms the original data using the dual distribution transformation
#' as described in Cirillo and Taleb (2020).
#'
#' @param data Data frame containing the variable
#' @param variable String, name of variable to transform
#' @param h Numeric, high value parameter (population reference)
#' @param l Numeric, low value parameter (population minimum)
#' @return Numeric vector of transformed values
#' @references Cirillo, P., & Taleb, N. N. (2020). Tail risk of contagious diseases.
#'             Nature Physics, 16(6), 606-613.
#' @export
dual_transform <- function(data, variable, h, l = 1e4) {
  if (!is.data.frame(data) || !variable %in% names(data)) {
    stop("Invalid data or variable name")
  }
  x <- data[[variable]]
  l - h * log((h - x) / (h - l))
}

#' Calculate Shadow Mean with Smooth Transition
#'
#' Implements the shadow mean calculation from Cirillo and Taleb's (2020) methodology
#' for analyzing tail risk in contagious diseases. The function handles different cases
#' based on the shape parameter xi, with a smooth transition zone around xi = 1 to
#' prevent bootstrap spikes from extreme values when xi approaches 1 from below.
#'
#' @param h Numeric, high value parameter
#' @param l_alt Numeric, alternative low value parameter
#' @param xi Numeric, shape parameter from GPD
#' @param sigma Numeric, scale parameter from GPD
#' @param transition_width Numeric, width of transition zone around xi = 1 (default 0.05)
#' @param max_multiplier Numeric, maximum multiplier for regularization (default 50)
#' @return Numeric, the calculated shadow mean
#' @references Cirillo, P., & Taleb, N. N. (2020). Tail risk of contagious diseases.
#'             Nature Physics, 16(6), 606-613.
#' @export
shadow_mean <- function(h, l_alt, xi, sigma, transition_width = 0.01, max_multiplier = 5) {
  if (!all(is.numeric(c(h, l_alt, xi, sigma)))) {
    stop("All parameters must be numeric")
  }

  # Define transition boundaries
  xi_lower <- 1 - transition_width
  xi_upper <- 1 + transition_width
  
  if (xi >= xi_upper) {
    # Use the original Cirillo formula for xi well above 1
    (h - l_alt) *
      exp(sigma / (h * xi)) *
      (sigma / (h * xi))^(1 / xi) *
      gamma((xi - 1) / xi) *
      pgamma(sigma / (h * xi), (xi - 1) / xi, lower.tail = FALSE) + l_alt
      
  } else if (xi <= xi_lower) {
    # For xi well below 1, use regularized standard formula to prevent extreme values
    regularized_denom <- max(1 - xi, sigma / (max_multiplier * sigma + l_alt))
    sigma / regularized_denom + l_alt
    
  } else {
    # Smooth transition zone: blend both formulas
    # Calculate both estimates
    cirillo_mean <- tryCatch({
      (h - l_alt) *
        exp(sigma / (h * xi)) *
        (sigma / (h * xi))^(1 / xi) *
        gamma((xi - 1) / xi) *
        pgamma(sigma / (h * xi), (xi - 1) / xi, lower.tail = FALSE) + l_alt
    }, error = function(e) NA)
    
    # Use regularized standard formula for stability
    regularized_denom <- max(1 - xi, sigma / (max_multiplier * sigma + l_alt))
    standard_mean <- sigma / regularized_denom + l_alt
    
    # Create smooth weight using logistic function
    # Weight ranges from 0 (use standard) to 1 (use Cirillo)
    weight <- 1 / (1 + exp(-10 * (xi - 1) / transition_width))
    
    # Handle case where Cirillo calculation failed
    if (is.na(cirillo_mean) || !is.finite(cirillo_mean)) {
      return(standard_mean)
    }
    
    # Blend the two estimates
    weight * cirillo_mean + (1 - weight) * standard_mean
  }
}

#' Calculate Shadow Mean with Error Handling
#'
#' Implements the complete shadow mean calculation process from Cirillo and Taleb (2020)
#' including data transformation, GPD fitting, and parameter estimation with robust
#' error handling.
#'
#' @param data Data frame containing the data
#' @param variable String, name of the variable to analyze
#' @param threshold Numeric, threshold value
#' @param pop_ref Numeric, population reference value
#' @param fixed_u Logical, whether to use fixed threshold. Default is TRUE
#' @return Numeric, calculated shadow mean or NA if error occurs
#' @references Cirillo, P., & Taleb, N. N. (2020). Tail risk of contagious diseases.
#'             Nature Physics, 16(6), 606-613.
#' @export
calculate_shadow_mean <- function(data, variable, threshold, conf, fixed_u = TRUE) {
  tryCatch(
    {
      # Dual transformation parameters
      z <- dual_transform(data, variable, conf$pop_reference, conf$pop_min)

      # Determine threshold
      if (!fixed_u) {
        threshold <- tea::GH(z)$threshold
      }

      # Fit GPD and extract parameters
      fit <- fit_gpd_model(z, threshold, conf$pop_reference, conf)
      xi <- fit$xi
      beta <- fit$beta

      # Validate parameters
      if (is.na(xi) || is.na(beta) || xi <= -1 || beta <= 0) {
        return(NA)
      }

      # Calculate sigma
      n_u <- sum(z > threshold)
      n <- length(z)
      sigma <- beta * (n_u / n)^xi

      # Calculate and return shadow mean
      shadow_mean(conf$pop_reference, threshold, xi, sigma)
    },
    error = function(e) {
      NA
    }
  )
}

#' Extract Threshold-Exceeding Values
#'
#' Extracts and sorts the dual-transformed values from a dataset starting from
#' a specific year.
#'
#' @param data Data frame containing the data
#' @param variable String, name of variable to analyze
#' @param year Numeric, starting year for analysis
#' @param pop_ref Numeric, population reference value
#' @return Numeric vector of sorted transformed values
#' @export
threshold_extract <- function(data, variable, year, conf) {
  if (!is.data.frame(data) || !variable %in% names(data)) {
    stop("Invalid data or variable name")
  }

  new_data <- data %>%
    dplyr::filter(.data$start_year >= year)

  l <- conf$pop_min
  h <- conf$pop_reference
  z <- dual_transform(new_data, variable, h, l)

  sort(z)
}

#' Inverse Dual Transformation
#'
#' Reverses the dual transformation to recover original values.
#'
#' @param dual_value Numeric, transformed value to inverse
#' @param h Numeric, high value parameter
#' @param l Numeric, low value parameter (default 1e4)
#' @return Numeric, original value
#' @export
inv_dual_transform <- function(dual_value, h, l = 1e4) {
  if (!is.numeric(dual_value) || !is.numeric(h) || !is.numeric(l)) {
    stop("All parameters must be numeric")
  }

  h - (h - l) * exp((l - dual_value) / h)
}

#' Calculate Truncated Distribution Mean
#'
#' Generic function to calculate conditional mean E[X | L < X ≤ U] for various distributions
#' @param dist_type Distribution type (e.g., "lognormal", "weibull", "gamma")
#' @param params List of distribution parameters
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Expected value in the range
#' @export
calc_truncated_mean <- function(dist_type, params, lower, upper) {
  if (dist_type == "lognormal") {
    return(calc_truncated_lognormal_mean(
      params[["meanlog"]], params[["sdlog"]], lower, upper
    ))
  } else if (dist_type == "weibull") {
    return(calc_truncated_weibull_mean(
      params[["shape"]], params[["scale"]], lower, upper
    ))
  } else if (dist_type == "gamma") {
    return(calc_truncated_gamma_mean(
      params[["shape"]], params[["rate"]], lower, upper
    ))
  }
}

#' Calculate Truncated Lognormal Mean
#'
#' Calculates E[X | L < X ≤ U] for lognormal distribution
#' @param meanlog Mean of log(X)
#' @param sdlog SD of log(X)
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Expected value in the range
#' @export
calc_truncated_lognormal_mean <- function(meanlog, sdlog, lower, upper) {
  # For lognormal, E[X | L < X ≤ U] involves the standard normal CDF
  z1 <- (log(lower) - meanlog) / sdlog
  z2 <- (log(upper) - meanlog) / sdlog

  # Expected value formula for truncated lognormal
  full_mean <- exp(meanlog + (sdlog^2) / 2)
  phi_z1 <- pnorm(z1)
  phi_z2 <- pnorm(z2)

  norm_const <- phi_z2 - phi_z1
  if (norm_const == 0) {
    return(NA)
  }

  full_mean * (pnorm(z2 - sdlog) - pnorm(z1 - sdlog)) / norm_const
}

#' Calculate Truncated Weibull Mean
#'
#' Calculates E[X | L < X ≤ U] for Weibull distribution
#' @param shape Shape parameter
#' @param scale Scale parameter
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Expected value in the range
#' @export
calc_truncated_weibull_mean <- function(shape, scale, lower, upper) {
  # Normalizing constant (probability mass in the range)
  p_range <- pweibull(upper, shape, scale) - pweibull(lower, shape, scale)

  if (p_range == 0) {
    return(NA)
  }

  # For Weibull, we can use the incomplete gamma function relation
  # and the formula E[X] = scale * gamma(1 + 1/shape)

  # Define the truncated expectation integral
  weibull_expected <- function(a, b, shape, scale) {
    # This uses numerical integration to compute the truncated expectation
    integrand <- function(x) x * dweibull(x, shape, scale)
    stats::integrate(integrand, a, b)$value / p_range
  }

  weibull_expected(lower, upper, shape, scale)
}

#' Calculate Truncated Gamma Mean
#'
#' Calculates E[X | L < X ≤ U] for Gamma distribution
#' @param shape Shape parameter
#' @param rate Rate parameter (1/scale)
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Expected value in the range
#' @export
calc_truncated_gamma_mean <- function(shape, rate, lower, upper) {
  # Normalizing constant (probability mass in the range)
  p_range <- pgamma(upper, shape, rate) - pgamma(lower, shape, rate)

  if (p_range == 0) {
    return(NA)
  }

  # Define the truncated expectation integral
  gamma_expected <- function(a, b, shape, rate) {
    # This uses numerical integration to compute the truncated expectation
    integrand <- function(x) x * dgamma(x, shape, rate)
    stats::integrate(integrand, a, b)$value / p_range
  }

  gamma_expected(lower, upper, shape, rate)
}
