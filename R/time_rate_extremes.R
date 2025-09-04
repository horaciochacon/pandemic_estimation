#' Time-based rate extremes utilities (extracted from time_functions.R)
#'
#' Modularized functions for identifying and visualizing periods of highest and lowest
#' pandemic exceedance rates. Logic unchanged; moved here during Phase 5 modularization.
#'
#' @keywords internal
# (Removed redundant library() calls; packages loaded via R/load_all.R)

#' Plot pandemic rate extremes with enhanced visualization
#'
#' Creates a professional visualization of pandemic exceedance rates over time highlighting
#' the periods with highest and lowest rates. Includes improved styling and visual elements
#' for better presentation quality.
#'
#' @param extreme_periods Result object from find_rate_extremes function
#' @param threshold_value Threshold used to define exceedances
#'
#' @return A ggplot object with the enhanced rate extremes visualization
#' @export
plot_rate_extremes <- function(extreme_periods, threshold_value) {
  if (is.null(extreme_periods$all_windows)) {
    stop("The extreme_periods object must include all_windows data. Set return_all=TRUE when calling find_rate_extremes.")
  }
  window_size <- extreme_periods$window_size
  rate_colors <- c(
    "highest" = "#d73027",
    "lowest" = "#4575b4",
    "line" = "#666666"
  )
  mean_rate <- mean(extreme_periods$all_windows$rate, na.rm = TRUE)
  extreme_periods$summary <- extreme_periods$summary %>%
    dplyr::mutate(
      period_label = sprintf("%d – %d", window_start, window_end),
      rate_label = sprintf("%.3f", rate)
    )
  year_min <- min(extreme_periods$all_windows$window_end, na.rm = TRUE)
  year_max <- max(extreme_periods$all_windows$window_end, na.rm = TRUE)
  year_range <- year_max - year_min
  year_padding <- ceiling(year_range * 0.05)
  windows_plot <- ggplot(extreme_periods$all_windows, aes(x = window_end, y = rate)) +
    geom_hline(yintercept = mean_rate, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
    annotate("text", x = year_min + year_padding, y = mean_rate, label = paste("Mean rate:", round(mean_rate, 3)), hjust = 0, vjust = -0.5, size = 3, color = "gray30", fontface = "italic") +
    geom_line(color = rate_colors["line"], linewidth = 0.9, alpha = 0.9) +
    geom_area(fill = "gray90", alpha = 0.3) +
    geom_point(data = extreme_periods$summary, aes(color = period_type), size = 5, shape = 21, stroke = 1.5) +
    ggrepel::geom_text_repel(data = extreme_periods$summary, aes(label = period_label, color = period_type), fontface = "bold", box.padding = 0.8, point.padding = 0.5, nudge_x = 10, nudge_y = 0.08, size = 4, segment.color = "gray50", min.segment.length = 1e100, show.legend = FALSE) +
    ggrepel::geom_text_repel(data = extreme_periods$summary, aes(label = rate_label, color = period_type), box.padding = 0.6, nudge_x = 10, nudge_y = 0.05, point.padding = 0.5, size = 3.5, min.segment.length = 0, fontface = "italic", segment.color = "gray50", show.legend = FALSE) +
    scale_color_manual(values = rate_colors, name = "Period Type", labels = c("highest" = "Highest Risk Period", "lowest" = "Lowest Risk Period")) +
    labs(title = "Historical Variation in Pandemic Risk", subtitle = paste0(window_size, "-year sliding windows • Deaths threshold: 200k population-scaled deaths"), x = "Year", y = "Annual Exceedance Rate", caption = "Higher values indicate more frequent pandemic events exceeding the threshold.") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)), plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15)), plot.caption = element_text(size = 9, color = "gray30", hjust = 1, margin = margin(t = 10)), legend.position = "bottom", legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 9), axis.title.x = element_text(margin = margin(t = 10), size = 11), axis.title.y = element_text(margin = margin(r = 10), size = 11), axis.text = element_text(size = 10), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "gray90", linewidth = 0.3), panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(20, 20, 20, 20))
  windows_plot
}

#' Find highest and lowest rate periods
#'
#' Sliding window exceedance rate analysis with optional draw-based uncertainty quantification.
#' Logic unchanged from original consolidated file.
#' @inheritParams plot_rate_extremes
#' @param variable Variable name to analyze
#' @param threshold Threshold value to define exceedance
#' @param start_year,end_year Bounds of analysis period
#' @param window_size Sliding window length in years
#' @param return_all Return all window rows
#' @param generate_draws Whether to generate simulation draws
#' @param n_draws Number of draws for uncertainty quantification
#' @param loess_span Span parameter for loess smoother
#' @export
find_rate_extremes <- function(
    data,
    variable = "deaths_scaled",
    threshold,
    start_year = NULL,
    end_year = NULL,
    window_size = 75,
    return_all = FALSE,
    generate_draws = FALSE,
    n_draws = conf$draws %||% 1000,
    loess_span = 0.2) {
  if (is.null(start_year)) start_year <- min(data$start_year)
  if (is.null(end_year)) end_year <- max(data$start_year)
  filtered_data <- data %>%
    dplyr::filter(start_year >= start_year, start_year <= end_year) %>%
    dplyr::arrange(start_year)
  window_ends <- seq(start_year + window_size - 1, end_year, by = 1)
  window_starts <- window_ends - window_size + 1
  window_counts <- tibble::tibble(window_start = window_starts, window_end = window_ends, window_size = window_size, exceedances = 0, events = 0, rate = 0)
  for (i in seq_len(nrow(window_counts))) {
    w_start <- window_counts$window_start[i]
    w_end <- window_counts$window_end[i]
    window_events <- filtered_data %>% dplyr::filter(start_year >= w_start, start_year <= w_end)
    window_counts$events[i] <- nrow(window_events)
    window_counts$exceedances[i] <- sum(window_events[[variable]] > threshold)
    window_counts$rate[i] <- window_counts$exceedances[i] / window_size
  }
  result <- list(window_size = window_size, threshold = threshold, variable = variable)
  if (return_all) result$all_windows <- window_counts
  if (generate_draws) {
    cat("Generating draws for rate extremes uncertainty quantification...\n")
    loess_model <- stats::loess(rate ~ window_end, data = window_counts, span = loess_span)
    window_counts$loess_fit <- stats::predict(loess_model, newdata = window_counts)
    window_counts$residuals <- window_counts$rate - window_counts$loess_fit
    residual_se <- stats::sd(window_counts$residuals)
    rate_draws <- matrix(0, nrow = nrow(window_counts), ncol = n_draws)
    set.seed(conf$seeds$secondary %||% 42)
    for (i in seq_len(n_draws)) {
      noise <- stats::rnorm(nrow(window_counts), mean = 0, sd = residual_se)
      rate_draws[, i] <- pmax(0, window_counts$loess_fit + noise)
    }
    draws_summary <- data.frame(window_end = window_counts$window_end, window_start = window_counts$window_start, observed_rate = window_counts$rate, loess_fit = window_counts$loess_fit, mean_draw = rowMeans(rate_draws), median = apply(rate_draws, 1, median), lower_95 = apply(rate_draws, 1, quantile, 0.025), upper_95 = apply(rate_draws, 1, quantile, 0.975), lower_80 = apply(rate_draws, 1, quantile, 0.1), upper_80 = apply(rate_draws, 1, quantile, 0.9))
    hi_idx <- which.max(draws_summary$mean_draw)
    lo_idx <- which.min(draws_summary$mean_draw)
    highest_draw_year <- draws_summary$window_end[hi_idx]
    lowest_draw_year <- draws_summary$window_end[lo_idx]
    highest_year_draws <- rate_draws[hi_idx, ]
    lowest_year_draws <- rate_draws[lo_idx, ]
    highest_year_quantiles <- quantile(highest_year_draws, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))
    lowest_year_quantiles <- quantile(lowest_year_draws, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))
    highest_period <- draws_summary[hi_idx, ] %>% dplyr::mutate(period_type = "highest", exceedances = window_counts$exceedances[hi_idx], events = window_counts$events[hi_idx], rate = mean_draw, rate_lower_95 = lower_95, rate_upper_95 = upper_95)
    lowest_period <- draws_summary[lo_idx, ] %>% dplyr::mutate(period_type = "lowest", exceedances = window_counts$exceedances[lo_idx], events = window_counts$events[lo_idx], rate = mean_draw, rate_lower_95 = lower_95, rate_upper_95 = upper_95)
    summary <- dplyr::bind_rows(highest_period, lowest_period) %>%
      dplyr::select(window_start, window_end, exceedances, events, rate, rate_lower_95, rate_upper_95, period_type) %>%
      dplyr::mutate(window_size = window_size)
    result$summary <- summary
    result$highest_period <- highest_period
    result$lowest_period <- lowest_period
    result$draws <- list(rate_draws = rate_draws, window_ends = window_counts$window_end, window_starts = window_counts$window_start, summary = draws_summary, highest_draw_year = highest_draw_year, lowest_draw_year = lowest_draw_year, highest_year_draws = highest_year_draws, lowest_year_draws = lowest_year_draws, highest_year_quantiles = highest_year_quantiles, lowest_year_quantiles = lowest_year_quantiles, loess_span = loess_span, residual_se = residual_se)
  } else {
    highest_idx <- which.max(window_counts$rate)
    lowest_idx <- which.min(window_counts$rate)
    highest_period <- window_counts[highest_idx, ] %>% dplyr::mutate(period_type = "highest")
    lowest_period <- window_counts[lowest_idx, ] %>% dplyr::mutate(period_type = "lowest")
    summary <- dplyr::bind_rows(highest_period, lowest_period) %>%
      dplyr::select(window_start, window_end, exceedances, events, rate, period_type) %>%
      dplyr::mutate(window_size = window_size)
    result$summary <- summary
    result$highest_period <- highest_period
    result$lowest_period <- lowest_period
  }
  result
}
