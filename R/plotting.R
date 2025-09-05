#' @importFrom magrittr %>%
#' @importFrom dplyr filter
utils::globalVariables(c(".data", "%>%", "dual_transform", "determine_threshold"))

#' Plot Histogram of Variable
#'
#' This function creates a histogram of the specified variable from the data.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @param bin Numeric, number of bins in the histogram. Default is 40.
#' @param log Logical, whether to use logarithmic scale for x-axis. Default is FALSE.
#' @return A ggplot object representing the histogram.
plot_histogram <- function(data, variable, bin = 40, log = FALSE) {
  p <- ggplot2::ggplot(data, aes(x = .data[[variable]])) +
    geom_histogram(
      bins = bin,
      fill = "skyblue",
      color = "black"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = paste(variable, "(in thousands)"), y = "Number of Events") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(axis.text.x = element_text(margin = margin(t = 10)))

  if (log == TRUE) {
    p <- p +
      scale_x_log10(labels = trans_format("log10", math_format())) +
      annotation_logticks(sides = "b", outside = TRUE)
  }

  return(p)
}
## (Removed unused theme_lancet helper: previously a thin wrapper over theme_minimal
## with minor grid tweak. No references in codebase; delete to reduce surface.)

#' Create Zipf Plot
#'
#' This function creates a Zipf plot (log-log survival plot) of the specified variable.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @return A ggplot object representing the Zipf plot.
zipf_plot <- function(data, variable) {
  x <- sort(as.numeric(data[[variable]]))
  ypoints <- 1 - ppoints(data[[variable]])

  temp <- tibble(x, ypoints)

  temp %>%
    ggplot(aes(x = x, y = ypoints)) +
    geom_point(
      size = 1,
      alpha = 0.8,
      col = "purple"
    ) +
    scale_x_log10(labels = trans_format("log10", math_format())) +
    scale_y_log10() +
    labs(x = paste(variable, "(in thousands)"), y = "Empirical survival (log scale)") +
    annotation_logticks(sides = "bl", outside = TRUE) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 10))
    )
}


#' Create Mean Excess Plot
#'
#' This function creates a mean excess plot of the specified variable.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @return A ggplot object representing the mean excess plot.
mean_excess_plot <- function(data, variable) {
  x <- as.numeric(data[[variable]])
  myrank <- function(x, na_last = TRUE) {
    ranks <- sort.list(sort.list(x, na.last = na_last))
    if (is.na(na_last)) {
      x <- x[!is.na(x)]
    }
    for (i in unique(x[duplicated(x)])) {
      which <- x == i & !is.na(x)
      ranks[which] <- max(ranks[which])
    }
    ranks
  }
  x <- sort(x)
  n_excess <- unique(floor(length(x) - myrank(x)))
  points <- unique(x)
  nl <- length(points)
  n_excess <- n_excess[-nl]
  points <- points[-nl]
  excess <- cumsum(rev(x))[n_excess] - n_excess * points
  y <- excess / n_excess
  omit <- 3
  xx <- points[1:(nl - omit)]
  yy <- y[1:(nl - omit)]

  ggplot(dplyr::tibble(xx, yy), aes(x = xx, y = yy)) +
    geom_point(
      size = 1,
      alpha = 0.8,
      col = "purple"
    ) +
    # scale_x_log10() +
    # scale_y_log10() +
    labs(x = "Threshold", y = "Mean Excess") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 10))
    )
}

#' Create Combined Descriptive Plots
#'
#' Creates a combined view of descriptive plots including histogram (normal and log scale),
#' Zipf plot, and mean excess plot.
#'
#' @param data Data frame containing the analysis data
#' @param variable String, name of variable to analyze
#' @param bins Integer, number of bins for histograms (default: 40)
#' @return A ggplot object containing the combined plots
#' @export
create_descriptive_plots <- function(data, variable, bins = 40, output_dir = NULL, conf = NULL) {
  # Generate individual plots
  p1 <- plot_histogram(data, variable, bins, log = FALSE) +
    labs(title = "Regular Histogram")

  p2 <- plot_histogram(data, variable, bins, log = TRUE) +
    labs(title = "Log-scale Histogram")

  p3 <- zipf_plot(data, variable) +
    labs(title = "Zipf Plot")

  p4 <- mean_excess_plot(data, variable) +
    labs(title = "Mean Excess Plot")

  # Combine plots using patchwork
  combined_plots <- (p1 + p2) / (p3 + p4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste("Descriptive Analysis:", variable),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
    )

  # Save plot if output directory is provided
  save_plot_if_enabled(combined_plots, "descriptive_plots", output_dir, conf, width = 12, height = 10)
  
  return(combined_plots)
}

#' Create Tail Plot
#'
#' This function creates a tail plot of the data with GPD fit and annotations.
#'
#' @param data Data frame containing the data.
#' @param variable String, name of the variable to plot.
#' @param fit Object, GPD fit object from evir::gpd.
#' @param xi Numeric, shape parameter from GPD.
#' @param u_best Numeric, threshold used in GPD fitting.
#' @param beta Numeric, scale parameter from GPD.
#' @param shadow Numeric, calculated shadow mean.
#' @param log10 Logical, whether to use logarithmic scale for y-axis. Default is FALSE.
#' @return None. The function prints the plot.
tail_plot <- function(
    data, variable, fit, xi, u_best, beta, pop_ref, shadow, tail_limit, log10 = FALSE,
    output_dir = NULL, conf = NULL) {
  # Define extension factor for plotting
  extend <- 1.5

  # Maximum plot range
  plotmax <- max(data[[variable]]) * extend

  # Generate probabilities and corresponding quantiles
  xx <- seq(0, 1, length.out = 1000)
  z <- evir::qgpd(xx, xi, u_best, beta)
  z <- pmax(pmin(z, plotmax), u_best)

  # Calculate fitted curve probabilities
  prob_less_thresh <- fit$p.less.thresh
  y <- (1 - prob_less_thresh) * (1 - evir::pgpd(z, xi, u_best, beta))
  curve_data <- data.frame(x = z, y = y)

  # Empirical data: sort, compute empirical survival, and ensure no zero values
  all_data <- data %>%
    arrange(.data[[variable]]) %>%
    mutate(
      x = .data[[variable]],
      y = 1 - ecdf(.data$x)(.data$x)
    )

  # Avoid zero probability at minimum
  all_data$y[which.min(all_data$y)] <- 0.006

  # Labels for threshold and shadow mean
  threshold_label <- paste0(
    scales::comma(
      round(inv_dual_transform(u_best, h = pop_ref, l = 1e4))
    ),
    " deaths"
  )

  shadow_label <- paste0(scales::comma(round(shadow)), " deaths")

  # Construct the plot
  p <- all_data %>%
    filter(.data$x > tail_limit) %>%
    ggplot() +
    # Empirical points
    geom_point(aes(x = .data$x, y = y), size = 1, alpha = 0.6) +
    # GPD fitted curve
    geom_line(
      data = curve_data, aes(x = .data$x, y = y, color = "GPD Fit"),
      linewidth = 1
    ) +
    # Vertical lines for threshold and shadow mean
    geom_vline(
      aes(xintercept = u_best, color = "Threshold"),
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(
      aes(xintercept = shadow, color = "Mean (Expected Value)"),
      linetype = "dashed",
      linewidth = 1
    ) +
    # Log10 scale on x-axis
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format())) +
    # Axis and legend labels
    labs(
      x = paste("Number of", stringr::str_to_title(stringr::str_replace_all(variable, "_", " "))),
      y = "Probability Density Function (PDF)",
      color = "Legend"
    ) +
    # Log ticks and themes
    annotation_logticks(sides = "bl", outside = TRUE) +
    # Threshold and shadow mean text
    annotate(
      "text",
      x = u_best * 1.2, y = 0.92, label = threshold_label, color = "darkgreen", hjust = 0
    ) +
    annotate(
      "text",
      x = shadow * 1.2, y = 0.92, label = shadow_label, color = "darkblue", hjust = 0
    ) +
    # Repel labels for points above shadow mean
    ggrepel::geom_text_repel(
      data = subset(all_data, all_data$x > shadow),
      aes(x = .data$x, y = y, label = .data$name),
      size = 3, box.padding = 0.5
    ) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      axis.text.x = element_text(margin = margin(t = 10), size = 16),
      axis.text.y = element_text(margin = margin(r = 10), size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    ) +
    scale_color_manual(values = c(
      "Threshold" = "darkgreen",
      "Mean (Expected Value)" = "darkblue",
      "GPD Fit" = "red"
    ))

  # Optionally use log10 scale on y-axis
  if (log10) {
    p <- p + scale_y_log10()
  }

  # Save plot if output directory is provided
  save_plot_if_enabled(p, paste0("tail_plot_", variable), output_dir, conf)
  
  # Also print for interactive use
  print(p)
  
  invisible(p)
}

#' Plot Time-Scaled Deaths
#'
#' Creates a professional visualization of pandemic and epidemic deaths over time, with deaths scaled
#' to the reference population year. Shows events through time with error bars and labels
#' for the largest events. Uses enhanced visual styling for improved readability and presentation.
#'
#' @param data Data frame containing the pandemic/epidemic data
#' @param lower_cutoff Numeric threshold for highlighting events with text labels
#' @param label_threshold Numeric threshold for applying text labels to events (defaults to 1e7)
#' @param min_year Numeric, earliest year to include in the visualization (defaults to 1600)
#' @param reference_year Numeric, population reference year for scaling (defaults to 2025)
#' @param include_labels Logical, whether to include title, subtitle, and caption (defaults to TRUE)
#' @return A ggplot object representing the enhanced time-scaled deaths plot
#' @export
plot_time_scaled_deaths <- function(data, lower_cutoff = 2e5, label_threshold = 1e7,
                                    min_year = 1600, reference_year = 2025,
                                    include_labels = TRUE,
                                    output_dir = NULL, conf = NULL) {
  # Filter data to include only events from min_year onwards
  data_time <- data |>
    dplyr::select(name, start_year, end_year, deaths, deaths_scaled, type) |>
    dplyr::filter(start_year >= min_year)

  # Add period categorization for visual encoding
  data_time <- data_time |>
    dplyr::mutate(
      # Create period groupings
      period = case_when(
        start_year < 1800 ~ "Pre-1800",
        start_year < 1900 ~ "19th Century",
        start_year < 1950 ~ "1900-1950",
        start_year < 2000 ~ "1950-2000",
        TRUE ~ "21st Century"
      ),
      # Order periods chronologically
      period = factor(period, levels = c("Pre-1800", "19th Century", "1900-1950", "1950-2000", "21st Century")),
      # Calculate deaths in millions for labels
      deaths_millions = deaths_scaled / 1e6,
      # Create formatted label with deaths
      event_label = if_else(
        deaths_scaled > label_threshold,
        paste0(name, "\n(", round(deaths_millions, 1), "M)"),
        name
      )
    )

  # Define color palette based on periods using The Lancet palette with gradient
  period_colors <- c(
    "Pre-1800" = "#925E9F", # Quaternary (purple) for oldest events
    "19th Century" = "#0099B4", # Tertiary (teal)
    "1900-1950" = "#42B540", # Secondary (green)
    "1950-2000" = "#00468B", # Primary (blue)
    "21st Century" = "#ED0000" # Highlight (red) for most recent events
  )

  # Define alpha and size based on recency and magnitude
  point_size_range <- c(2, 8)
  max_deaths <- max(data_time$deaths_scaled, na.rm = TRUE)

  # Create the enhanced plot
  p <- ggplot(data_time, aes(x = start_year, y = deaths_scaled)) +
    # Add minimal reference lines for major time periods (hair-line thin per guidelines)
    geom_vline(
      xintercept = c(1800, 1900, 1950, 2000), linetype = "dotted",
      color = "gray50", linewidth = 0.25
    ) +

    # Add threshold reference line (hair-line thin per guidelines)
    geom_hline(aes(yintercept = lower_cutoff),
      linetype = "dashed", color = "grey50", linewidth = 0.25
    ) +
    annotate("text",
      x = min(data_time$start_year), y = lower_cutoff * 0.8,
      label = "Analysis threshold", color = "grey50", hjust = 0, fontface = "italic", size = 3
    ) +

    # Removed horizontal error bars for pandemic duration

    # Add points with size proportional to deaths and color by period
    geom_point(aes(
      color = period, size = deaths_scaled, fill = period
    ), shape = 21, stroke = 1, alpha = 0.7) +

    # Add text labels for significant events
    ggrepel::geom_text_repel(
      data = data_time %>% filter(deaths_scaled > label_threshold),
      aes(label = event_label, color = period),
      box.padding = 0.8,
      point.padding = 0.5,
      min.segment.length = 0.2,
      max.overlaps = 20,
      segment.color = "gray40",
      segment.alpha = 0.7,
      fontface = "bold",
      size = 3,
      show.legend = FALSE
    ) +

    # Set scales with Lancet-compatible specifications
    scale_color_manual(values = period_colors, name = "Time Period") +
    scale_fill_manual(values = period_colors, guide = "none") +
    scale_size_continuous(range = point_size_range, guide = "none") +
    scale_y_log10(
      name = "Population-adjusted deaths (log scale)",
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(1e5, 1e6, 1e7, 1e8, 1e9)
    ) +
    scale_x_continuous(
      name = "Year",
      breaks = seq(min_year, 2025, by = 50), # Use evenly spaced, human-readable tick marks
      minor_breaks = NULL, # Remove minor breaks per Lancet guidelines
      expand = expansion(mult = c(0.02, 0.02))
    ) +

    # Set labels conditionally
    labs(
      title = if (include_labels) "Historical Pandemic Mortality Impact (1600-Present)" else NULL,
      subtitle = if (include_labels) paste0("Death tolls scaled to ", reference_year, " world population") else NULL,
      x = "Year",
      caption = if (include_labels) {
        paste0(
          "Note: Events with < ", scales::label_number(scale_cut = scales::cut_short_scale())(label_threshold),
          " deaths may not be labeled. Size of points corresponds to magnitude of impact."
        )
      } else {
        NULL
      }
    ) +

    # Set theme according to Lancet guidelines
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
    # Override the legend to use filled points rather than default
    guides(color = guide_legend(override.aes = list(
      shape = 21, # Shape with fill and border
      fill = period_colors, # Use the same colors defined for period
      size = 4, # Larger marker size for legend
      stroke = 1.2 # Border thickness
    )))

  # Save plot if output directory is provided
  save_plot_if_enabled(p, "time_scaled_deaths", output_dir, conf)
  
  return(p)
}

# Sensitivity plots

#' Plot shadow prices sensitivity analysis
#' @param results Data frame with sensitivity analysis results
#' @param var Variable analyzed ("deaths" or "deaths_scaled")
plot_shadow_sensitivity <- function(results, var, analysis_type, output_dir = NULL, conf = NULL) {
  # Define the x-axis column based on analysis type
  x_col <- case_when(
    analysis_type == "year" ~ "year",
    analysis_type == "threshold" ~ "threshold",
    analysis_type == "cutoff" ~ "cutoff",
    TRUE ~ "value"
  )

  # Only proceed if we have the required columns
  req_cols <- c(x_col, "boot", "boot_low", "boot_upp")
  if (!all(req_cols %in% colnames(results))) {
    warning(
      "Missing required columns: ",
      paste(setdiff(req_cols, colnames(results)), collapse = ", ")
    )
    return(NULL)
  }

  p <- ggplot(results, aes_string(x = x_col)) +
    geom_ribbon(aes(ymin = boot_low, ymax = boot_upp), alpha = 0.2) +
    geom_line(aes(y = boot)) +
    labs(
      y = if (var == "deaths") "Expected Deaths (millions)" else "Expected Scaled Deaths (millions)",
      x = analysis_type,
      title = paste0("Sensitivity analysis by ", analysis_type)
    ) +
    scale_y_log10() +
    theme_minimal()

  if (analysis_type %in% c("threshold", "cutoff")) {
    p <- p + scale_x_log10()
  }

  # Save plot if output directory is provided
  save_plot_if_enabled(p, paste0("shadow_sensitivity_", var, "_", analysis_type), output_dir, conf)
  
  # Also print for interactive use
  print(p)
  
  invisible(p)
}

#' Plot yearly deaths sensitivity analysis
#' @param results Data frame with sensitivity analysis results
#' @param var Variable analyzed ("deaths" or "deaths_scaled")
#' @param analysis_type Type of analysis ("year" or "threshold")
plot_yearly_deaths <- function(results, var, analysis_type, output_dir = NULL, conf = NULL) {
  # Define the x-axis column based on analysis type
  x_col <- case_when(
    analysis_type == "year" ~ "year",
    analysis_type == "threshold" ~ "threshold",
    analysis_type == "cutoff" ~ "cutoff",
    TRUE ~ "value"
  )

  # Only proceed if we have the required columns
  req_cols <- c(x_col, "yearly_deaths", "yearly_deaths_low", "yearly_deaths_up")
  if (!all(req_cols %in% colnames(results))) {
    warning(
      "Missing required columns: ",
      paste(setdiff(req_cols, colnames(results)), collapse = ", ")
    )
    return(NULL)
  }

  p <- ggplot(results, aes_string(x = x_col)) +
    geom_ribbon(
      aes(
        ymin = yearly_deaths_low,
        ymax = yearly_deaths_up
      ),
      alpha = 0.2
    ) +
    geom_line(aes(y = yearly_deaths)) +
    labs(
      y = if (var == "deaths") "Expected Annual Deaths (millions)" else "Expected Annual Scaled Deaths (millions)",
      x = analysis_type,
      title = paste0(
        "Annual",
        ifelse(var == "deaths_scaled", " scaled ", ""),
        " deaths sensitivity by ", analysis_type
      )
    ) +
    theme_minimal()

  if (analysis_type %in% c("threshold", "cutoff")) {
    p <- p + scale_x_log10()
  }
  
  # Save plot if output directory is provided
  save_plot_if_enabled(p, paste0("yearly_deaths_", var, "_", analysis_type), output_dir, conf)
  
  # Also print for interactive use
  print(p)
  
  invisible(p)
}

#' Format number as label with appropriate suffix (M for millions, k for thousands)
#' @param value Numeric value to format
#' @return Character string with formatted label
format_number_label <- function(value) {
  if (value >= 1e6) {
    # For millions, show one decimal place if not a whole number
    millions <- value / 1e6
    if (millions == round(millions)) {
      paste0(round(millions), "M")
    } else {
      paste0(round(millions, 1), "M")
    }
  } else if (value >= 1e3) {
    # For thousands, show whole numbers
    paste0(round(value / 1e3), "k")
  } else {
    as.character(round(value))
  }
}

#' Plot cumulative deaths sensitivity analysis with enhanced design
#' @param results Data frame with sensitivity analysis results
#' @param var Variable analyzed ("deaths" or "deaths_scaled")
#' @param analysis_type Type of analysis ("year", "threshold", or "cutoff")
#' @param config Configuration object containing analysis parameters
plot_cumulative_deaths <- function(results, var, analysis_type, config = NULL, output_dir = NULL, conf = NULL) {
  # Define the x-axis column based on analysis type
  x_col <- case_when(
    analysis_type == "year" ~ "year",
    analysis_type == "threshold" ~ "threshold",
    analysis_type == "cutoff" ~ "cutoff",
    TRUE ~ "value"
  )

  # Only proceed if we have the required columns
  req_cols <- c(x_col, "cum_deaths", "cum_deaths_low", "cum_deaths_up")
  if (!all(req_cols %in% colnames(results))) {
    warning(
      "Missing required columns: ",
      paste(setdiff(req_cols, colnames(results)), collapse = ", ")
    )
    return(NULL)
  }

  # Define professional color palette
  main_color <- "#377EB8" # Blue for main line
  ribbon_color <- "#C6DBEF" # Light blue for uncertainty ribbon
  accent_color <- "gray50" # Gray for vertical line indicators
  grid_color <- "gray95" # Light gray for grid lines

  # Display unit based on values
  max_value <- max(results$cum_deaths_up, na.rm = TRUE)
  if (max_value > 5000) {
    scale_factor <- "billions"
    results <- results %>%
      mutate(across(c(cum_deaths, cum_deaths_low, cum_deaths_up), ~ .x / 1000))
  } else {
    scale_factor <- "millions"
  }

  # Generate x-axis label
  x_label <- case_when(
    analysis_type == "year" ~ "Year",
    analysis_type == "threshold" ~ "Threshold (in 100,000 deaths)",
    analysis_type == "cutoff" ~ "Lower Cutoff (in 100,000 deaths)",
    TRUE ~ stringr::str_to_title(analysis_type)
  )

  # Generate title with variable formatting
  title <- case_when(
    var == "deaths" ~ "Pandemic Cumulative Death Burden Sensitivity",
    var == "deaths_scaled" ~ "Population-adjusted Pandemic Death Burden Sensitivity",
    TRUE ~ paste0("Sensitivity Analysis for ", stringr::str_to_title(var))
  )

  # Generate subtitle with analysis type information
  subtitle <- case_when(
    analysis_type == "year" ~ "Impact of varying the historical time window on cumulative death estimates",
    analysis_type == "threshold" ~ "Impact of threshold selection on cumulative death estimates",
    analysis_type == "cutoff" ~ "Impact of lower cutoff selection on cumulative death estimates",
    TRUE ~ paste0("Sensitivity by ", stringr::str_to_title(analysis_type))
  )

  # Generate formatted caption
  caption <- paste0(
    "Note: Shaded area represents 95% confidence interval. ",
    if (analysis_type %in% c("threshold", "cutoff")) "Both axes use logarithmic scale. " else "Y-axis uses logarithmic scale. ",
    "Values are shown in ", scale_factor, "."
  )

  # Build plot with professional styling
  p <- ggplot(results, aes(x = .data[[x_col]])) +
    # Add ribbon for uncertainty bounds with semi-transparent fill
    geom_ribbon(
      aes(
        ymin = cum_deaths_low,
        ymax = cum_deaths_up
      ),
      fill = ribbon_color,
      alpha = 0.6
    ) +
    # Add line with purposeful styling
    geom_line(
      aes(y = cum_deaths),
      color = main_color,
      linewidth = 1
    ) +
    # Add points for emphasis
    geom_point(
      aes(y = cum_deaths),
      color = main_color,
      size = 2,
      alpha = 0.7
    ) +
    # Set appropriate scales
    scale_y_continuous(
      name = paste0("Cumulative Deaths (", scale_factor, ", log scale)"),
      labels = scales::label_number(accuracy = 1),
      breaks = if (!is.null(config) && !is.null(config$plots$log_breaks_major)) config$plots$log_breaks_major else c(1, 10, 100, 1000),
      trans = "log1p",
      limits = c(0, max(results$cum_deaths_up, na.rm = TRUE) * 1.5)
    ) +
    # Apply conditional x-axis scaling
    {
      if (analysis_type %in% c("threshold", "cutoff")) {
        scale_x_log10(
          name = x_label,
          labels = scales::label_number(accuracy = 0.1),
          breaks = scales::breaks_log(n = 5)
        )
      } else {
        scale_x_continuous(
          name = x_label,
          breaks = scales::pretty_breaks(n = 8)
        )
      }
    } +
    # Apply professional styling
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    # Enhanced theme with attention to detail
    theme_minimal() +
    theme(
      # Typography
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20), color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray30", hjust = 0, margin = margin(t = 15)),

      # Axis styling
      axis.title.x = element_text(margin = margin(t = 10), size = 12),
      axis.title.y = element_text(margin = margin(r = 10), size = 12),
      axis.text = element_text(size = 10),

      # Panel and border styling for distinct plotting area
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),

      # Grid styling for cleaner appearance
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = grid_color, linewidth = 0.3),

      # Overall spacing
      plot.margin = margin(20, 20, 20, 20)
    )

  # Add log ticks if using log scale
  # if (analysis_type %in% c("threshold", "cutoff")) {
  #   p <- p + annotation_logticks(sides = "bl", color = "gray50", size = 0.3, short = unit(0.07, "cm"), mid = unit(0.15, "cm"), long = unit(0.2, "cm"))
  # } else {
  #   p <- p + annotation_logticks(sides = "l", color = "gray50", size = 0.3, short = unit(0.07, "cm"), mid = unit(0.15, "cm"), long = unit(0.2, "cm"))
  # }

  # Add reference line for significant value if available in config
  if (analysis_type == "year" && !is.null(results$year) && !is.null(config)) {
    baseline_year <- if (!is.null(config$sens$year$cutoff)) config$sens$year$cutoff else config$sens$year$base_year
    if (!is.null(baseline_year)) {
      p <- p +
        geom_vline(xintercept = baseline_year, linetype = "dashed", color = accent_color, linewidth = 0.5) +
        annotate("text",
          x = baseline_year, y = min(results$cum_deaths_low, na.rm = TRUE) * 1.5,
          label = baseline_year, color = accent_color, hjust = -0.2, vjust = 0, size = 3, fontface = "italic"
        )
    }
  }
  if (analysis_type == "threshold" && !is.null(results$threshold) && !is.null(config)) {
    baseline_threshold <- config$sens$threshold$base_thresholds$deaths_scaled
    if (!is.null(baseline_threshold)) {
      threshold_label <- format_number_label(baseline_threshold)
      p <- p +
        geom_vline(xintercept = baseline_threshold, linetype = "dashed", color = accent_color, linewidth = 0.5) +
        annotate("text",
          x = baseline_threshold, y = min(results$cum_deaths_low, na.rm = TRUE) * 0.5,
          label = threshold_label, color = accent_color, hjust = -0.2, vjust = 0, size = 3, fontface = "italic"
        )
    }
  }
  if (analysis_type == "cutoff" && !is.null(results$cutoff) && !is.null(config)) {
    baseline_cutoff <- config$thresholds$lower_cutoff_deaths_scaled
    if (!is.null(baseline_cutoff)) {
      cutoff_label <- format_number_label(baseline_cutoff)
      p <- p +
        geom_vline(xintercept = baseline_cutoff, linetype = "dashed", color = accent_color, linewidth = 0.5) +
        annotate("text",
          x = baseline_cutoff, y = min(results$cum_deaths_low, na.rm = TRUE) * 1.5,
          label = cutoff_label, color = accent_color, hjust = -0.2, vjust = 0, size = 3, fontface = "italic"
        )
    }
  }

  # Save plot if output directory is provided
  save_plot_if_enabled(p, paste0("cumulative_deaths_", var, "_", analysis_type), output_dir, conf)
  
  # Also print for interactive use
  print(p)
  
  invisible(p)
}


#' Plot parameter sensitivity analysis
#' @param results Data frame with sensitivity analysis results
#' @param var Variable analyzed ("deaths" or "deaths_scaled")
#' @param analysis_type Type of analysis ("year" or "threshold")
plot_par_sensitivity <- function(results, var, analysis_type, output_dir = NULL, conf = NULL) {
  # Define the x-axis column based on analysis type
  x_col <- case_when(
    analysis_type == "year" ~ "year",
    analysis_type == "threshold" ~ "threshold",
    analysis_type == "cutoff" ~ "cutoff",
    TRUE ~ "value"
  )

  # Only proceed if we have the required columns
  req_cols <- c(x_col, "Xi", "Beta")
  if (!all(req_cols %in% colnames(results))) {
    warning(
      "Missing required columns: ",
      paste(setdiff(req_cols, colnames(results)), collapse = ", ")
    )
    return(NULL)
  }

  p <- results %>%
    dplyr::select(all_of(c(x_col, "Xi", "Beta"))) %>%
    pivot_longer(
      cols = c(Xi, Beta),
      names_to = "parameter",
      values_to = "estimate"
    ) %>%
    ggplot(aes_string(x = x_col, y = "estimate", color = "parameter")) +
    geom_line() +
    facet_wrap(parameter ~ ., scales = "free_y", ncol = 1) +
    labs(
      y = "Parameter Value",
      x = analysis_type,
      title = paste0(var, " - Parameter sensitivity by ", analysis_type)
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  if (analysis_type %in% c("threshold", "cutoff")) {
    p <- p + scale_x_log10()
  }

  # Save plot if output directory is provided
  save_plot_if_enabled(p, paste0("parameter_sensitivity_", var, "_", analysis_type), output_dir, conf)
  
  # Also print for interactive use
  print(p)
  
  invisible(p)
}
