####################################################################################################
## File:        output_manager.R
## Author:      Horacio Chacon Torrico
##
## Description: Output management utilities for organizing plots, tables, and results
##              into script-specific timestamped directories.
##
## Functions:   - create_script_output_dir(): Create script-specific timestamped output directory
##              - get_output_path(): Generate file paths within script output directories
##              - setup_script_run(): Initialize script run with metadata
##              - save_plot_if_enabled(): Conditionally save plots based on configuration
##              - save_table_if_enabled(): Conditionally save tables based on configuration
##              - setup_run_all_session(): Initialize run_all session with clean directory structure
##              - run_all_scripts(): Execute all analysis scripts with unified output management
##              - generate_run_all_summary(): Generate comprehensive execution summary
####################################################################################################

#' Create Script-Specific Timestamped Output Directory
#'
#' Creates a timestamped output directory specific to the calling script.
#' Directory format: output/{script_name}_{timestamp}/
#'
#' @param script_name Character string identifying the script (e.g., "run_pot", "run_validation")
#' @param base_dir Base output directory. Default: "output"
#' @return Character string with full path to created directory
#' @examples
#' dir <- create_script_output_dir("run_pot")
#' # Returns something like: "output/run_pot_2024-09-04_14-30-25/"
create_script_output_dir <- function(script_name, base_dir = "output") {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  script_dir <- file.path(base_dir, paste0(script_name, "_", timestamp))
  
  # Create main directory and subdirectories
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(script_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(script_dir, "tables"), showWarnings = FALSE)
  
  return(script_dir)
}

#' Get Output File Path
#'
#' Generate appropriate file path within script's output directory structure.
#'
#' @param output_dir Base output directory for the script run
#' @param filename Desired filename
#' @param subdir Subdirectory within output_dir ("plots" or "tables")
#' @return Full file path
get_output_path <- function(output_dir, filename, subdir = "plots") {
  if (is.null(output_dir)) {
    return(NULL)
  }
  return(file.path(output_dir, subdir, filename))
}

#' Setup Script Run Environment
#'
#' Initialize script-specific output directory and save run metadata.
#'
#' @param script_name Name of the calling script
#' @param conf Configuration object
#' @param base_dir Base output directory
#' @return List containing output_dir path and metadata
setup_script_run <- function(script_name, conf = NULL, base_dir = "output") {
  # Check if output saving is enabled
  if (!is.null(conf) && !is.null(conf$output) && !isTRUE(conf$output$save_plots) && !isTRUE(conf$output$save_tables)) {
    return(list(output_dir = NULL, metadata = NULL))
  }
  
  # Create output directory
  output_dir <- create_script_output_dir(script_name, base_dir)
  
  # Create metadata
  metadata <- list(
    script_name = script_name,
    timestamp = Sys.time(),
    r_version = R.version.string,
    working_directory = getwd(),
    config_active = Sys.getenv("R_CONFIG_ACTIVE", "default")
  )
  
  # Save metadata
  if (!is.null(conf)) {
    metadata$configuration <- conf
  }
  
  # Write metadata to YAML file
  tryCatch({
    yaml_content <- yaml::as.yaml(metadata)
    writeLines(yaml_content, file.path(output_dir, "run_metadata.yml"))
  }, error = function(e) {
    message("Warning: Could not save metadata - ", e$message)
  })
  
  message("Output directory created: ", output_dir)
  
  return(list(
    output_dir = output_dir,
    metadata = metadata
  ))
}

#' Save Plot If Enabled
#'
#' Conditionally save a ggplot object based on configuration settings.
#'
#' @param plot ggplot object to save
#' @param filename Base filename (without extension)
#' @param output_dir Output directory path
#' @param conf Configuration object
#' @param width Plot width in inches. Default: 10
#' @param height Plot height in inches. Default: 8
#' @param dpi Plot resolution. Default: 300
#' @return Invisibly returns the plot object
save_plot_if_enabled <- function(plot, filename, output_dir = NULL, conf = NULL, 
                                 width = 10, height = 8, dpi = 300) {
  
  # Return early if no output directory or saving disabled
  if (is.null(output_dir)) {
    return(invisible(plot))
  }
  
  if (!is.null(conf) && !is.null(conf$output) && !isTRUE(conf$output$save_plots)) {
    return(invisible(plot))
  }
  
  # Determine formats to save
  formats <- c("svg")  # Default format: SVG only
  if (!is.null(conf) && !is.null(conf$output) && !is.null(conf$output$plot_formats)) {
    formats <- conf$output$plot_formats
  }
  
  # Save in each format
  for (format in formats) {
    full_filename <- paste0(filename, ".", format)
    filepath <- get_output_path(output_dir, full_filename, "plots")
    
    tryCatch({
      ggplot2::ggsave(
        filename = filepath,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        device = format
      )
      message("Plot saved: ", filepath)
    }, error = function(e) {
      message("Warning: Could not save plot ", filepath, " - ", e$message)
    })
  }
  
  return(invisible(plot))
}

#' Save Base R Plot If Enabled
#'
#' Conditionally save base R graphics plots based on configuration settings.
#' This function captures base R plotting code and saves it using appropriate graphics devices.
#'
#' @param plot_expr Expression or function containing base R plotting code
#' @param filename Base filename (without extension)
#' @param output_dir Output directory path
#' @param conf Configuration object
#' @param width Plot width in inches. Default: 10
#' @param height Plot height in inches. Default: 8
#' @param res Plot resolution (for raster formats). Default: 300
#' @return Invisibly returns TRUE if plots were saved, FALSE otherwise
save_base_plot_if_enabled <- function(plot_expr, filename, output_dir = NULL, conf = NULL,
                                     width = 10, height = 8, res = 300) {
  
  # Return early if no output directory or saving disabled
  if (is.null(output_dir)) {
    return(invisible(FALSE))
  }
  
  if (!is.null(conf) && !is.null(conf$output) && !isTRUE(conf$output$save_plots)) {
    return(invisible(FALSE))
  }
  
  # Determine formats to save
  formats <- c("svg")  # Default format: SVG only
  if (!is.null(conf) && !is.null(conf$output) && !is.null(conf$output$plot_formats)) {
    formats <- conf$output$plot_formats
  }
  
  # Save in each format
  for (format in formats) {
    full_filename <- paste0(filename, ".", format)
    filepath <- get_output_path(output_dir, full_filename, "plots")
    
    tryCatch({
      # Open appropriate graphics device
      switch(format,
        "svg" = svg(filepath, width = width, height = height),
        "png" = png(filepath, width = width * res, height = height * res, res = res),
        "pdf" = pdf(filepath, width = width, height = height),
        "eps" = postscript(filepath, width = width, height = height, paper = "special"),
        stop("Unsupported format: ", format)
      )
      
      # Execute the plotting code
      if (is.function(plot_expr)) {
        plot_expr()
      } else {
        eval(plot_expr)
      }
      
      # Close the device
      dev.off()
      message("Base R plot saved: ", filepath)
    }, error = function(e) {
      # Ensure device is closed even if there's an error
      if (dev.cur() != 1) dev.off()
      message("Warning: Could not save base R plot ", filepath, " - ", e$message)
    })
  }
  
  return(invisible(TRUE))
}

#' Save Table If Enabled
#'
#' Conditionally save table content based on configuration settings.
#'
#' @param content Character vector or data.frame to save
#' @param filename Filename with extension
#' @param output_dir Output directory path
#' @param conf Configuration object
#' @return Invisible NULL
save_table_if_enabled <- function(content, filename, output_dir = NULL, conf = NULL) {
  
  # Return early if no output directory or saving disabled
  if (is.null(output_dir)) {
    return(invisible(NULL))
  }
  
  if (!is.null(conf) && !is.null(conf$output) && !isTRUE(conf$output$save_tables)) {
    return(invisible(NULL))
  }
  
  filepath <- get_output_path(output_dir, filename, "tables")
  
  tryCatch({
    if (is.character(content)) {
      # Save character content (like markdown)
      writeLines(content, filepath)
    } else {
      # Save data.frame as CSV
      write.csv(content, filepath, row.names = FALSE)
    }
    message("Table saved: ", filepath)
  }, error = function(e) {
    message("Warning: Could not save table ", filepath, " - ", e$message)
  })
  
  return(invisible(NULL))
}

#' Save Markdown Summary If Enabled
#'
#' Conditionally save captured console output as formatted markdown based on configuration settings.
#'
#' @param output_lines Character vector of captured console output lines
#' @param filename Base filename (without .md extension)
#' @param output_dir Output directory path
#' @param conf Configuration object
#' @param title Optional title for the markdown document
#' @return Invisibly returns TRUE if saved, FALSE otherwise
save_markdown_summary <- function(output_lines, filename, output_dir = NULL, conf = NULL, title = NULL) {
  
  # Return early if no output directory or saving disabled
  if (is.null(output_dir)) {
    return(invisible(FALSE))
  }
  
  if (!is.null(conf) && !is.null(conf$output) && !isTRUE(conf$output$save_tables)) {
    return(invisible(FALSE))
  }
  
  # Ensure filename has .md extension
  if (!grepl("\\.md$", filename)) {
    filename <- paste0(filename, ".md")
  }
  
  filepath <- get_output_path(output_dir, filename, "tables")
  
  tryCatch({
    # Start building markdown content
    markdown_content <- character()
    
    # Add title if provided
    if (!is.null(title)) {
      markdown_content <- c(markdown_content, paste("#", title), "", "")
    }
    
    # Add timestamp
    markdown_content <- c(
      markdown_content,
      paste("*Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "*"),
      "",
      "---",
      ""
    )
    
    # Process output lines to format as markdown
    for (line in output_lines) {
      # Convert section headers (lines with ===)
      if (grepl("^=+", line)) {
        next # Skip separator lines
      } else if (grepl("BULK DIAGNOSTICS|FINAL RECOMMENDATIONS|ENHANCED SELECTION|STATISTICAL FIT", line)) {
        # Convert main headers
        markdown_content <- c(markdown_content, paste("##", line), "")
      } else if (grepl("^  ✓|^  ⚠|^  📊", line)) {
        # Convert recommendation lines with emojis
        markdown_content <- c(markdown_content, paste("-", line))
      } else if (grepl("^  [A-Z].*:", line)) {
        # Convert field lines to bullet points
        markdown_content <- c(markdown_content, paste("-", line))
      } else if (grepl("RECOMMENDATION:", line)) {
        # Convert recommendation header
        markdown_content <- c(markdown_content, paste("###", line), "")
      } else if (line != "" && !grepl("^\\s*$", line)) {
        # Regular content lines
        markdown_content <- c(markdown_content, line)
      } else {
        # Preserve empty lines
        markdown_content <- c(markdown_content, "")
      }
    }
    
    # Write markdown content
    writeLines(markdown_content, filepath)
    message("Markdown summary saved: ", filepath)
  }, error = function(e) {
    message("Warning: Could not save markdown summary ", filepath, " - ", e$message)
    return(invisible(FALSE))
  })
  
  return(invisible(TRUE))
}

#' Get Script Name from Call Stack
#'
#' Attempt to determine the calling script name for automatic output directory naming.
#'
#' @return Character string with script name or "unknown_script"
get_calling_script <- function() {
  # Try to get script name from command line arguments
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grepl("--file=", args)]
  
  if (length(script_arg) > 0) {
    script_path <- sub("--file=", "", script_arg[1])
    script_name <- tools::file_path_sans_ext(basename(script_path))
    return(script_name)
  }
  
  # Fallback to call stack inspection
  calls <- sys.calls()
  for (call in calls) {
    if (is.call(call) && as.character(call[[1]]) == "source") {
      if (length(call) > 1) {
        source_file <- as.character(call[[2]])
        script_name <- tools::file_path_sans_ext(basename(source_file))
        return(script_name)
      }
    }
  }
  
  return("unknown_script")
}

#' Setup Run All Session with Clean Directory Structure
#'
#' Creates a clean timestamped output directory specifically for run_all orchestration.
#' Unlike setup_script_run, this does NOT create plots/tables subdirectories or metadata files.
#'
#' @param conf Configuration object
#' @param base_dir Base output directory. Default: "output"
#' @return List containing output_dir path and metadata
setup_run_all_session <- function(conf = NULL, base_dir = "output") {
  # Check if output saving is enabled
  if (!is.null(conf) && !is.null(conf$output) && !isTRUE(conf$output$save_plots) && !isTRUE(conf$output$save_tables)) {
    return(list(output_dir = NULL, metadata = NULL))
  }
  
  # Create clean run_all directory (no subfolders)
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  output_dir <- file.path(base_dir, paste0("run_all_", timestamp))
  
  # Create ONLY the main directory (no plots/tables subfolders)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create basic metadata
  metadata <- list(
    session_type = "run_all",
    timestamp = Sys.time(),
    r_version = R.version.string,
    working_directory = getwd(),
    config_active = Sys.getenv("R_CONFIG_ACTIVE", "default")
  )
  
  message("Run all output directory created: ", output_dir)
  
  return(list(
    output_dir = output_dir,
    metadata = metadata
  ))
}

#' Execute All Analysis Scripts with Unified Output Management
#'
#' Orchestrates the execution of all pandemic estimation analysis scripts
#' with organized output directory structure and comprehensive logging.
#'
#' @param conf Configuration object
#' @param output_dir Base output directory for the run_all session
#' @return List containing execution results and summary statistics
run_all_scripts <- function(conf, output_dir = NULL) {
  
  # Define scripts to run with their configurations
  scripts_config <- list(
    list(
      name = "run_pot",
      file = "src/run_PoT.R",
      description = "Main peaks-over-threshold mixture model analysis",
      config_env = "default"
    ),
    list(
      name = "run_validation",
      file = "src/run_validation.R", 
      description = "Model validation with out-of-sample testing",
      config_env = "validation"
    ),
    list(
      name = "run_sens",
      file = "src/run_sens.R",
      description = "Comprehensive sensitivity analysis", 
      config_env = "sensitivity"
    ),
    list(
      name = "bulk_diagnostics",
      file = "src/bulk_diagnostics.R",
      description = "Bulk distribution component diagnostics",
      config_env = "diagnostics"
    ),
    list(
      name = "tail_diagnostics", 
      file = "src/tail_diagnostics.R",
      description = "Tail distribution component diagnostics",
      config_env = "diagnostics"
    )
  )
  
  # Initialize execution tracking
  execution_log <- data.frame(
    script_name = character(),
    description = character(), 
    config_env = character(),
    start_time = character(),
    end_time = character(),
    duration_minutes = numeric(),
    status = character(),
    output_dir = character(),
    stringsAsFactors = FALSE
  )
  
  # Execute each script with error handling
  overall_start_time <- Sys.time()
  successful_runs <- 0
  failed_runs <- 0
  
  message("Executing ", length(scripts_config), " analysis scripts...")
  
  for (i in seq_along(scripts_config)) {
    script_info <- scripts_config[[i]]
    
    message(sprintf("[%d/%d] Running %s...", i, length(scripts_config), script_info$name))
    
    # Create script-specific subdirectory
    script_output_dir <- file.path(output_dir, script_info$name)
    dir.create(script_output_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(script_output_dir, "plots"), showWarnings = FALSE)
    dir.create(file.path(script_output_dir, "tables"), showWarnings = FALSE)
    
    # Record start time
    script_start_time <- Sys.time()
    
    # Execute script with error handling
    tryCatch({
      # Set configuration environment for this script
      Sys.setenv(R_CONFIG_ACTIVE = script_info$config_env)
      
      # Create a modified environment that will intercept setup_script_run
      script_env <- new.env(parent = globalenv())
      
      # Load all functions into the script environment
      sys.source("R/load_all.R", envir = script_env)
      
      # Override setup_script_run to use our custom output directory
      # This prevents individual scripts from creating their own timestamped directories
      script_env$setup_script_run <- function(script_name, conf = NULL, base_dir = "output") {
        # Return pre-created script directory structure (no additional files/folders)
        list(
          output_dir = script_output_dir,
          metadata = list(
            script_name = script_name,
            timestamp = Sys.time(),
            r_version = R.version.string,
            working_directory = getwd(),
            config_active = Sys.getenv("R_CONFIG_ACTIVE", "default"),
            run_all_session = TRUE,
            master_output_dir = output_dir
          )
        )
      }
      
      # Source the script in our controlled environment
      sys.source(script_info$file, envir = script_env)
      
      # Record success
      script_end_time <- Sys.time()
      duration <- as.numeric(difftime(script_end_time, script_start_time, units = "mins"))
      
      execution_log <- rbind(execution_log, data.frame(
        script_name = script_info$name,
        description = script_info$description,
        config_env = script_info$config_env,
        start_time = format(script_start_time, "%Y-%m-%d %H:%M:%S"),
        end_time = format(script_end_time, "%Y-%m-%d %H:%M:%S"),
        duration_minutes = round(duration, 2),
        status = "SUCCESS",
        output_dir = script_output_dir,
        stringsAsFactors = FALSE
      ))
      
      successful_runs <- successful_runs + 1
      message("✓ ", script_info$name, " completed (", round(duration, 1), " min)")
      
    }, error = function(e) {
      # Record failure
      script_end_time <- Sys.time()
      duration <- as.numeric(difftime(script_end_time, script_start_time, units = "mins"))
      
      execution_log <<- rbind(execution_log, data.frame(
        script_name = script_info$name,
        description = script_info$description,
        config_env = script_info$config_env,
        start_time = format(script_start_time, "%Y-%m-%d %H:%M:%S"),
        end_time = format(script_end_time, "%Y-%m-%d %H:%M:%S"),
        duration_minutes = round(duration, 2),
        status = paste("FAILED:", e$message),
        output_dir = script_output_dir,
        stringsAsFactors = FALSE
      ))
      
      failed_runs <<- failed_runs + 1
      message("✗ ", script_info$name, " failed (", round(duration, 1), " min): ", e$message)
    })
  }
  
  # Calculate overall statistics
  overall_end_time <- Sys.time()
  total_duration <- as.numeric(difftime(overall_end_time, overall_start_time, units = "mins"))
  
  # Return execution results
  list(
    execution_log = execution_log,
    summary = list(
      total_scripts = length(scripts_config),
      successful_runs = successful_runs,
      failed_runs = failed_runs,
      success_rate = round(100 * successful_runs / length(scripts_config), 1),
      total_duration_minutes = round(total_duration, 1),
      start_time = overall_start_time,
      end_time = overall_end_time
    )
  )
}

#' Generate Comprehensive Execution Summary
#'
#' Generate summary tables and files documenting the run_all execution results.
#'
#' @param results Results list from run_all_scripts()
#' @param output_dir Output directory for summary files
#' @param conf Configuration object
#' @return Invisible NULL
generate_run_all_summary <- function(results, output_dir = NULL, conf = NULL) {
  
  if (is.null(output_dir)) {
    return(invisible(NULL))
  }
  
  # Save execution log directly to root directory
  log_file <- file.path(output_dir, "execution_log.csv")
  tryCatch({
    write.csv(results$execution_log, log_file, row.names = FALSE)
    message("Execution log saved: ", log_file)
  }, error = function(e) {
    message("Warning: Could not save execution log - ", e$message)
  })
  
  # Create and save session info directly to root directory
  session_file <- file.path(output_dir, "session_info.txt")
  tryCatch({
    capture.output(sessionInfo(), file = session_file)
    message("Session info saved: ", session_file)
  }, error = function(e) {
    message("Warning: Could not save session info - ", e$message)
  })
  
  # Create comprehensive summary
  summary_content <- c(
    "# Pandemic Estimation - Comprehensive Analysis Pipeline",
    "",
    paste("**Execution completed:** ", format(results$summary$end_time, "%Y-%m-%d %H:%M:%S")),
    paste("**Total duration:** ", results$summary$total_duration_minutes, " minutes"),
    paste("**Scripts executed:** ", results$summary$total_scripts),
    paste("**Successful:** ", results$summary$successful_runs),
    paste("**Failed:** ", results$summary$failed_runs),
    paste("**Success rate:** ", results$summary$success_rate, "%"),
    "",
    "## Output Structure",
    ""
  )
  
  # Add directory structure with status
  for (i in seq_len(nrow(results$execution_log))) {
    log_entry <- results$execution_log[i, ]
    status_symbol <- if (log_entry$status == "SUCCESS") "✓" else "✗"
    summary_content <- c(summary_content, 
      paste("- **", log_entry$script_name, "/**: ", status_symbol, " ", log_entry$status, 
            " (", log_entry$duration_minutes, " min)"))
  }
  
  summary_content <- c(summary_content, 
    "",
    "## Files Generated",
    "- `execution_log.csv`: Detailed execution log with timing and status",
    "- `session_info.txt`: R session and package information", 
    "- `run_summary.md`: This summary document",
    ""
  )
  
  # Save summary as markdown directly to root directory
  summary_file <- file.path(output_dir, "run_summary.md")
  tryCatch({
    writeLines(summary_content, summary_file)
    message("Run summary saved: ", summary_file)
  }, error = function(e) {
    message("Warning: Could not save run summary - ", e$message)
  })
  
  # Print summary to console
  message("\n", paste(rep("=", 60), collapse = ""))
  message("EXECUTION SUMMARY")
  message(paste(rep("=", 60), collapse = ""))
  message("Total duration: ", results$summary$total_duration_minutes, " minutes")
  message("Success rate: ", results$summary$success_rate, "%")
  message("Output directory: ", output_dir)
  message(paste(rep("=", 60), collapse = ""))
  
  return(invisible(NULL))
}