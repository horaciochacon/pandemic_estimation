#' Fits a two-stage mixture model with flexible bulk distribution and GPD tail
#'
#' Implements an advanced statistical model that combines a flexible probability distribution
#' for non-extreme pandemic events (the bulk component) with a Generalized Pareto Distribution (GPD)
#' for modeling the extreme tail behavior of pandemic severity.
#'
#' @param data Data frame containing transformed data
#' @param variable String, name of variable to analyze
#' @param conf Configuration object
#' @param u_best Upper threshold value for GPD tail
#' @param dist_type Distribution type for bulk component (default: "lognormal",
#' options: "weibull", "gamma", etc.)
#' @param constrain_survival Boolean, whether to constrain the survival function
#' continuity between bulk and tail
#' @param summary Summary statistics for the dataset
#'
#' @details
#' This function implements a two-component mixture model for pandemic severity analysis:
#'
#' 1. **Tail Component**: Uses a Generalized Pareto Distribution (GPD) to model extreme
#'    pandemic events above threshold u_best. This follows the Peaks-over-Threshold (PoT)
#'    approach from extreme value theory.
#'
#' 2. **Bulk Component**: Models non-extreme events between lower_cutoff and u_best
#'    using a flexible parametric distribution (lognormal, weibull, or gamma).
#'
#' When constrain_survival=TRUE, the function implements a constrained mixture model
#' that ensures continuity of the survival function at the threshold point u_best.
#' This is achieved through penalized maximum likelihood estimation that enforces
#' the following constraints:
#'
#' - The survival function of the bulk component at the lower cutoff matches the
#'   empirical probability of exceeding the cutoff
#' - The survival function of the bulk component at the upper threshold matches
#'   the GPD exceedance probability
#'
#' The constrained approach ensures a smooth transition between bulk and tail
#' components, avoiding discontinuities in risk estimation.
#'
#' @return List containing:
#' \item{bulk}{Bulk component fit with parameters and conditional mean}
#' \item{tail}{Tail component fit with GPD parameters and shadow mean}
#' \item{mean}{Combined conditional mean from both components}
#' \item{lower_cutoff}{Lower threshold for bulk component}
#' \item{tail_threshold}{Upper threshold for GPD tail}
#' \item{bulk_dist_type}{Distribution type used for bulk component}
#' \item{constrained}{Whether a constrained model was used}
#' @export
fit_mixture_model <- function(data, variable, conf, u_best, dist_type = "lognormal",
                              constrain_survival = FALSE, summary) {
  # First fit the tail component to get exceedance probabilities
  tail_results <- fit_gpd_model(
    x = data$dual,
    u_best = u_best,
    h = conf$pop_reference,
    summary = summary,
    method = conf$method
  )

  # Then fit bulk component with specified distribution using survival constraints if requested
  bulk_results <- fit_bulk_component(
    data = data,
    variable = variable,
    lower_cutoff = conf$thresholds$lower_cutoff,
    upper_cutoff = u_best,
    dist_type = dist_type,
    constrain_survival = constrain_survival,
    tail_fit = tail_results,
    conf = conf
  )

  # Calculate mixture mean
  mixture_mean <- calculate_conditional_mean(
    bulk_results,
    tail_results
  )

  list(
    bulk = bulk_results,
    tail = tail_results,
    mean = mixture_mean,
    lower_cutoff = conf$thresholds$lower_cutoff,
    tail_threshold = u_best,
    bulk_dist_type = dist_type,
    constrained = constrain_survival
  )
}

#' Fit GPD Model
#'
#' Fits a Generalized Pareto Distribution model using either 'mev' or 'evir' packages.
#' Based on methodology from Cirillo & Taleb (2020).
#'
#' @param data Data frame containing transformed data
#' @param u_best Threshold value
#' @param h High value parameter
#' @param summary Summary statistics
#' @param method Fitting method ('mev' or 'evir')
#' @return List containing model fit and parameters
#' @references Cirillo, P., & Taleb, N. N. (2020). Tail risk of contagious diseases.
#' @export
fit_gpd_model <- function(x, u_best, h, summary = NULL, method = "mev") {
  tryCatch(
    {
      # Fit GPD model
      fit <- if (method == "mev") {
        fit_mev <- mev::fit.gpd(x, threshold = u_best)
        list(
          par.ests = c(
            fit_mev$estimate["shape"], # xi
            fit_mev$estimate["scale"] # beta
          ),
          par.se = c(
            fit_mev$std.err["shape"],
            fit_mev$std.err["scale"]
          ),
          p.less.thresh = length(x[x <= u_best]) / length(x)
        )
      } else if (method == "evir") {
        evir::gpd(x, threshold = u_best)
      } else {
        stop("Invalid method. Use 'mev' or 'evir'")
      }

      # Extract parameters
      xi <- as.numeric(fit$par.ests[1])
      beta <- as.numeric(fit$par.ests[2])

      if (!is.null(summary)) {
        sigma <- beta * (summary$PoT / summary$n)^xi
        shadow <- shadow_mean(h, u_best, xi, sigma)

        return(
          list(
            fit = fit,
            xi = xi,
            xi_se = as.numeric(fit$par.se[1]),
            beta = beta,
            beta_se = as.numeric(fit$par.se[2]),
            tail_prob = summary$PoT / summary$PoB,
            sigma = sigma,
            shadow = shadow
          )
        )
      } else {
        return(
          list(
            fit = fit,
            xi = xi,
            beta = beta
          )
        )
      }
    },
    error = function(e) {
      warning("Error in GPD fitting: ", e$message)
      NULL
    }
  )
}

#' Fit Bulk Component with Flexible Distribution
#'
#' Fits a parameterized probability distribution to the bulk (non-extreme) portion
#' of the pandemic severity data. Supports constrained fitting to ensure survival
#' function continuity with the GPD tail component.
#'
#' @param data Data frame containing the data
#' @param variable String, name of variable to analyze
#' @param lower_cutoff Lower threshold for bulk
#' @param upper_cutoff Upper threshold for bulk (tail threshold)
#' @param dist_type Distribution type (default: "lognormal", options: "weibull", "gamma", etc.)
#' @param scale_data Boolean, whether to scale data for numerical stability
#' @param constrain_survival Boolean, whether to constrain the survival function
#' @param tail_fit GPD tail fit object (needed if constrain_survival=TRUE)
#' @param conf Configuration object
#'
#' @details
#' This function handles the bulk component of the two-stage mixture model by fitting
#' a parametric distribution to non-extreme pandemic events. It operates in two modes:
#'
#' 1. **Standard MLE Fitting**: When constrain_survival=FALSE, uses standard maximum
#'    likelihood estimation to fit the specified probability distribution to events
#'    between lower_cutoff and upper_cutoff.
#'
#' 2. **Constrained Fitting**: When constrain_survival=TRUE, uses penalized maximum
#'    likelihood estimation to ensure continuity in the survival function at both
#'    the lower cutoff and upper threshold (interface with the GPD tail).
#'    The constraints enforce:
#'    - S_bulk(lower_cutoff) = P(X > lower_cutoff)
#'    - S_bulk(upper_cutoff) = P(X > upper_cutoff) from the GPD fit
#'
#' The function supports multiple distribution types:
#' - lognormal: Works best for most pandemic data (default)
#' - weibull: Alternative heavy-tailed distribution
#' - gamma: Alternative for data with specific variance properties
#'
#' Numerical stability is enhanced through optional data scaling and multiple
#' optimization methods to improve convergence probability.
#'
#' @return List containing:
#' \item{fit}{Fitted distribution object}
#' \item{parameters}{Named list of distribution parameters}
#' \item{dist_type}{Type of distribution used}
#' \item{dual_mean}{Conditional mean in dual-transformed space}
#' \item{conditional_mean}{Conditional mean in original scale}
#' \item{bulk_prob}{Probability mass in the bulk region}
#' \item{scale_factor}{Scaling factor used (if applicable)}
#' \item{survival_lower}{Survival function value at lower cutoff}
#' \item{survival_thresh}{Survival function value at upper threshold}
#' @export
fit_bulk_component <- function(data, variable, lower_cutoff, upper_cutoff,
                               dist_type = "lognormal", scale_data = TRUE,
                               constrain_survival = FALSE, tail_fit = NULL, conf) {
  # Filter data for bulk range fitting using the dual-transformed data
  # This ensures both bulk and tail are modeled in the same dual-transformed space
  bulk_data <- data[data$dual >= lower_cutoff & data$dual <= upper_cutoff, ]

  # Check if we have enough data points
  if (nrow(bulk_data) < 5) {
    warning(
      "Too few data points (", nrow(bulk_data),
      ") between cutoffs in dual space. Falling back to lognormal."
    )
    dist_type <- "lognormal"
  }

  # Bulk probability
  bulk_prob <- nrow(bulk_data) / nrow(data[data$dual >= lower_cutoff, ])

  # If constraining survival function and tail_fit is provided
  if (constrain_survival && !is.null(tail_fit)) {
    # Get the tail exceedance probability from GPD fit
    tail_exceedance_prob <- 1 - tail_fit$fit$p.less.thresh

    # Empirical probability of exceeding lower cutoff
    p_exceed_lower <- nrow(data[data$dual >= lower_cutoff, ]) / nrow(data)

    if (dist_type == "lognormal") {
      # Use constrained optimization to enforce the survival function constraints
      objective_function <- function(params) {
        meanlog <- params[1]
        sdlog <- params[2]

        # Calculate the survival function at lower cutoff and threshold
        s_lower <- 1 - plnorm(lower_cutoff, meanlog, sdlog)
        s_thresh <- 1 - plnorm(upper_cutoff, meanlog, sdlog)

        # Calculate the negative log-likelihood for the bulk data in dual space
        neg_ll <- -sum(dlnorm(bulk_data$dual, meanlog, sdlog, log = TRUE))

        # Add penalties for constraint violations
        penalty_lower <- 1e3 * (s_lower / p_exceed_lower - 1)^2
        penalty_thresh <- 1e10 * (s_thresh / tail_exceedance_prob - 1)^2

        # Return penalized negative log-likelihood
        neg_ll + penalty_lower + penalty_thresh
      }

      # Initial parameter estimates (using MLE without constraints) in dual space
      initial_fit <- MASS::fitdistr(bulk_data$dual, "lognormal")
      initial_params <- initial_fit$estimate

      # Run constrained optimization
      opt_result <- optim(
        par = c(initial_params[["meanlog"]], initial_params[["sdlog"]]),
        fn = objective_function,
        method = "L-BFGS-B",
        hessian = TRUE,
        lower = c(-Inf, 0.001) # sdlog must be positive
      )

      vcov <- solve(opt_result$hessian)
      se <- sqrt(diag(vcov))

      # Create a "bulk_fit" object with the optimized parameters
      bulk_fit <- list(
        estimate = c(meanlog = opt_result$par[1], sdlog = opt_result$par[2]),
        sd = se,
        vcov = vcov,
        n = nrow(bulk_data),
        loglik = -opt_result$value
      )
      class(bulk_fit) <- "fitdistr"
    } else if (dist_type == "weibull") {
      # Use constrained optimization for Weibull distribution
      objective_function <- function(params) {
        shape <- params[1]
        scale <- params[2]

        # Calculate the survival function at lower cutoff and threshold
        s_lower <- 1 - pweibull(lower_cutoff, shape, scale)
        s_thresh <- 1 - pweibull(upper_cutoff, shape, scale)

        # Calculate the negative log-likelihood for the bulk data in dual space
        neg_ll <- -sum(dweibull(bulk_data$dual, shape, scale, log = TRUE))

        # Add penalties for constraint violations
        penalty_lower <- 3400 * (s_lower / p_exceed_lower - 1)^2
        penalty_thresh <- 1000 * (s_thresh / tail_exceedance_prob - 1)^2

        # Return penalized negative log-likelihood
        neg_ll + penalty_lower + penalty_thresh
      }

      # Initial parameter estimates (using MLE without constraints) in dual space
      initial_fit <- MASS::fitdistr(bulk_data$dual, "weibull")
      initial_params <- initial_fit$estimate

      # Run constrained optimization
      opt_result <- optim(
        par = c(initial_params[["shape"]], initial_params[["scale"]]),
        fn = objective_function,
        method = "L-BFGS-B",
        lower = c(0.001, 0.001) # both shape and scale must be positive
      )

      # Create a "bulk_fit" object with the optimized parameters
      bulk_fit <- list(
        estimate = c(shape = opt_result$par[1], scale = opt_result$par[2]),
        sd = initial_fit$sd, # Using original SDs as approximation
        vcov = initial_fit$vcov, # Using original variance-covariance matrix
        n = nrow(bulk_data),
        loglik = -opt_result$value
      )
      class(bulk_fit) <- "fitdistr"
    } else {
      warning(
        "Survival function constraint currently only implemented\n        for lognormal and weibull distributions."
      )
      dist_type <- "lognormal"

      # Default to lognormal implementation
      objective_function <- function(params) {
        meanlog <- params[1]
        sdlog <- params[2]

        # Calculate the survival function at lower cutoff and threshold
        s_lower <- 1 - plnorm(lower_cutoff, meanlog, sdlog)
        s_thresh <- 1 - plnorm(upper_cutoff, meanlog, sdlog)

        # Calculate the negative log-likelihood for the bulk data in dual space
        neg_ll <- -sum(dlnorm(bulk_data$dual, meanlog, sdlog, log = TRUE))

        # Add penalties for constraint violations
        penalty_lower <- 100 * (s_lower / p_exceed_lower - 1)^2
        penalty_thresh <- 100 * (s_thresh / tail_exceedance_prob - 1)^2

        # Return penalized negative log-likelihood
        neg_ll + penalty_lower + penalty_thresh
      }

      # Initial parameter estimates (using MLE without constraints) in dual space
      initial_fit <- MASS::fitdistr(bulk_data$dual, "lognormal")
      initial_params <- initial_fit$estimate

      # Run constrained optimization
      opt_result <- optim(
        par = c(initial_params[["meanlog"]], initial_params[["sdlog"]]),
        fn = objective_function,
        method = "L-BFGS-B",
        lower = c(-Inf, 0.001) # sdlog must be positive
      )

      # Create a "bulk_fit" object with the optimized parameters
      bulk_fit <- list(
        estimate = c(meanlog = opt_result$par[1], sdlog = opt_result$par[2]),
        sd = initial_fit$sd, # Using original SDs as approximation
        vcov = initial_fit$vcov, # Using original variance-covariance matrix
        n = nrow(bulk_data),
        loglik = -opt_result$value
      )
      class(bulk_fit) <- "fitdistr"
    }
  } else {
    # Original unconstrained fitting code
    # Scale data to improve numerical stability if requested
    original_data <- bulk_data$dual
    scale_factor <- 1

    if (scale_data && dist_type != "lognormal") {
      # For weibull and gamma, scaling can help convergence
      scale_factor <- median(original_data)
      fitted_data <- original_data / scale_factor
    } else {
      fitted_data <- original_data
    }

    # Try different optimization methods if the default fails
    optimization_methods <- c("Nelder-Mead", "BFGS", "L-BFGS-B", "CG")
    max_attempts <- length(optimization_methods)

    for (i in seq_len(max_attempts)) {
      current_method <- optimization_methods[i]

      bulk_fit <- tryCatch(
        {
          MASS::fitdistr(fitted_data, dist_type, method = current_method)
        },
        error = function(e) {
          if (i == max_attempts) {
            warning(
              "All optimization methods failed for ", dist_type, ": ", e$message,
              "\nFalling back to lognormal distribution"
            )
            MASS::fitdistr(original_data, "lognormal")
          } else {
            # Return NULL to try the next method
            NULL
          }
        }
      )

      # If we got a successful fit, break the loop
      if (!is.null(bulk_fit)) break
    }
  }

  # Extract distribution parameters
  params <- as.list(bulk_fit$estimate)

  # Rescale parameters if needed (for scale parameters)
  if (!constrain_survival && scale_data && dist_type != "lognormal") {
    if (dist_type == "weibull") {
      params$scale <- params$scale * scale_factor
    } else if (dist_type == "gamma") {
      params$rate <- params$rate / scale_factor
    }
  }

  # Saving distribution type for reference
  bulk_fit$dist_type <- dist_type

  # Calculate conditional mean E[X | L < X ≤ U] in dual space
  dual_cond_mean <- calc_truncated_mean(
    dist_type = dist_type,
    params = params,
    lower = lower_cutoff,
    upper = upper_cutoff
  )

  # Calculate original scale conditional mean by inverse transforming the dual mean
  # This ensures that when combined with the shadow mean (already in original scale),
  # we're working in consistent units
  original_cond_mean <- inv_dual_transform(
    dual_cond_mean,
    h = conf$pop_reference,
    l = conf$pop_min
  )

  # Calculate survival probabilities (for verification)
  # Function name mapping for R's distribution functions
  dist_name_map <- c(
    "lognormal" = "lnorm",
    "weibull" = "weibull",
    "gamma" = "gamma"
  )

  # Get the correct distribution function name
  dist_func_name <- dist_name_map[dist_type]
  if (is.na(dist_func_name)) dist_func_name <- dist_type

  # Calculate survival probabilities in dual space
  survival_lower <- 1 - do.call(paste0("p", dist_func_name), c(list(q = lower_cutoff), params))
  survival_thresh <- 1 - do.call(paste0("p", dist_func_name), c(list(q = upper_cutoff), params))

  list(
    fit = bulk_fit,
    parameters = params,
    dist_type = dist_type,
    dual_mean = dual_cond_mean,
    conditional_mean = original_cond_mean,
    bulk_prob = bulk_prob,
    scale_factor = ifelse(exists("scale_factor"), scale_factor, 1),
    survival_lower = survival_lower,
    survival_thresh = survival_thresh
  )
}

#' Calculate Combined Conditional Mean for Mixture Model
#'
#' Calculates the expected value E[X | X > lower_cutoff] by weighted combination
#' of the bulk and tail components of the mixture model. This is a key function
#' for estimating overall pandemic severity.
#'
#' @param bulk_fit Bulk component fit with conditional_mean in original scale
#' @param tail_fit Tail component fit with shadow mean in original scale
#'
#' @details
#' This function combines the conditional means of the bulk and tail components
#' to produce an overall expected value for pandemic severity. It uses a weighted
#' average approach:
#'
#' E[X | X > lower_cutoff] = P(bulk) * E[X | lower_cutoff < X ≤ threshold] + P(tail) * E[X | X > threshold]
#'
#' Where:
#' - P(bulk) is the probability mass in the bulk component (between lower_cutoff and threshold)
#' - P(tail) is the probability mass in the tail component (above threshold)
#' - E[X | lower_cutoff < X ≤ threshold] is the conditional mean of the bulk component
#' - E[X | X > threshold] is the conditional mean of the tail component (shadow mean)
#'
#' This approach properly weights the contributions of both moderate and extreme
#' pandemic events to the overall expected severity, providing a more complete
#' risk assessment than traditional approaches that focus only on the tail.
#'
#' Both component means are in the original data scale (not dual-transformed) to
#' ensure consistent units in the final calculation.
#'
#' @return Combined conditional mean in original scale
#' @export
calculate_conditional_mean <- function(bulk_fit, tail_fit) {
  # Get bulk conditional mean E[X | lower < X ≤ threshold]
  # Already inverse-transformed to original scale in fit_bulk_component
  bulk_mean <- bulk_fit$conditional_mean

  # Get tail conditional mean E[X | X > threshold] from shadow mean
  # Shadow mean is already in original scale
  tail_mean <- tail_fit$shadow

  # Total conditional mean E[X | X > lower] is the weighted sum
  # Both means are in the same scale (original data scale, not dual-transformed)
  bulk_fit$bulk_prob * bulk_mean + tail_fit$tail_prob * tail_mean
}

#' Determine Analysis Threshold
#'
#' Determines the threshold for extreme value analysis based on provided configuration
#' or manual override.
#'
#' @param data Data frame containing the data
#' @param variable Name of variable to analyze
#' @param threshold Optional manual threshold override
#' @param conf Configuration object
#' @return Numeric threshold value
#' @export
determine_threshold <- function(data, variable, threshold = NULL, conf) {
  if (!is.null(threshold)) {
    return(dual_transform(
      tibble(x = threshold),
      "x",
      h = conf$pop_reference,
      l = conf$pop_min
    ))
  }

  if (!is.null(conf$thresholds[[variable]])) {
    return(dual_transform(
      tibble(x = conf$thresholds[[variable]]),
      "x",
      h = conf$pop_reference,
      l = conf$pop_min
    ))
  }

  stop("Threshold value must be specified either in configuration or as parameter")
}
