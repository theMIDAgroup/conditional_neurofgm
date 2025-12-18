#!/usr/bin/env Rscript

# Load required libraries
library(stringr)
library(dplyr)

# Function to extract numbers from folder names
extract_parameters <- function(folder_name) {
  # Extract p value (number after 'p')
  p_match <- str_extract(folder_name, "p(\\d+)", group = 1)
  p <- as.numeric(p_match)
  
  # Extract first n value (number after first 'n')
  n_match <- str_extract(folder_name, "n(\\d+)", group = 1)
  n <- as.numeric(n_match)
  
  return(list(p = p, n = n))
}

# Function to extract seed number
extract_seed <- function(seed_folder) {
  seed_match <- str_extract(seed_folder, "seed_(\\d+)", group = 1)
  return(as.numeric(seed_match))
}

# Function to extract computational time from log file
extract_comp_time <- function(log_file_path) {
  if (!file.exists(log_file_path)) {
    return(NA)
  }
  
  lines <- readLines(log_file_path, warn = FALSE)
  
  # Find lines with "Computational time of:"
  comp_time_lines <- grep("Computational time of:", lines, value = TRUE)
  
  if (length(comp_time_lines) == 0) {
    return(NA)
  }
  
  # Extract numbers after "Computational time of:"
  comp_times <- str_extract(comp_time_lines, "Computational time of:\\s*(\\d+\\.?\\d*)", group = 1)
  comp_times <- as.numeric(comp_times)
  comp_times <- comp_times[!is.na(comp_times)]
  
  if (length(comp_times) == 0) {
    return(NA)
  }
  
  return(max(comp_times))
}

# Main processing function
process_step1_folder <- function(step1_path) {
  # Create empty dataframe
  results_df <- data.frame(
    p = numeric(0),
    n = numeric(0),
    seed = numeric(0),
    comp_time = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Check if Step1 folder exists
  if (!dir.exists(step1_path)) {
    stop(paste("Directory", step1_path, "does not exist"))
  }
  
  # Get all folders starting with "complete_p"
  main_folders <- list.dirs(step1_path, recursive = FALSE, full.names = TRUE)
  main_folders <- main_folders[grepl("p\\d+_n\\d+_n\\d+", basename(main_folders))]
  
  cat("Found", length(main_folders), "main folders to process\n")
  
  for (main_folder in main_folders) {
    folder_name <- basename(main_folder)
    cat("Processing folder:", folder_name, "\n")
    
    # Extract p and n parameters
    params <- extract_parameters(folder_name)
    p_val <- params$p
    n_val <- params$n
    
    if (is.na(p_val) || is.na(n_val)) {
      cat("  Warning: Could not extract parameters from", folder_name, "\n")
      next
    }
    
    cat("  Parameters: p =", p_val, ", n =", n_val, "\n")
    
    # Get all seed folders
    seed_folders <- list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
    seed_folders <- seed_folders[grepl("seed_\\d+", basename(seed_folders))]
    
    cat("  Found", length(seed_folders), "seed folders\n")
    
    for (seed_folder in seed_folders) {
      seed_name <- basename(seed_folder)
      seed_val <- extract_seed(seed_name)
      
      if (is.na(seed_val)) {
        cat("    Warning: Could not extract seed from", seed_name, "\n")
        next
      }
      
      cat("    Processing seed:", seed_val, "\n")
      
      # Check for results/log folder
      log_folder <- file.path(seed_folder, "results_lit_comparison", "logs")
      
      if (!dir.exists(log_folder)) {
        cat("      Warning: Log folder not found at", log_folder, "\n")
        next
      }
      
      # Get all .log files
      log_files <- list.files(log_folder, pattern = "Comp_time_comparison.log$", full.names = TRUE)
      
      if (length(log_files) == 0) {
        cat("      Warning: No .log files found in", log_folder, "\n")
        next
      }
      
      cat("      Found", length(log_files), ".log files\n")
      
      # Extract computational times from all log files
      max_comp_time <- -Inf
      
      for (log_file in log_files) {
        comp_time <- extract_comp_time(log_file)
        
        if (!is.na(comp_time)) {
          max_comp_time <- max(max_comp_time, comp_time)
          cat("        ", basename(log_file), ": comp_time =", comp_time, "\n")
        } else {
          cat("        ", basename(log_file), ": no computational time found\n")
        }
      }
      
      # Only add row if we found at least one valid computational time
      if (max_comp_time != -Inf) {
        new_row <- data.frame(
          p = p_val,
          n = n_val,
          seed = seed_val,
          comp_time = max_comp_time
        )
        results_df <- rbind(results_df, new_row)
        cat("      Max computational time:", max_comp_time, "\n")
      } else {
        cat("      Warning: No valid computational times found\n")
      }
    }
  }
  
  return(results_df)
}

# Main execution
main <- function() {
  # Set the path to your Step1 folder
  step1_path <- "/group/diangelantonio/users/alessia_mapelli/Brain_simulations/systematic_simulations/Step1"  # Change this to your actual path
  
  # You can also accept command line argument
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    step1_path <- args[1]
  }
  
  cat("Starting processing of", step1_path, "\n")
  
  # Process the folders
  results <- process_step1_folder(step1_path)
  
  # Display results
  cat("\nProcessing complete!\n")
  cat("Total rows collected:", nrow(results), "\n")
  print(results)
  
  # Save the dataframe
  output_file <- paste0(step1_path, "/computational_times_litt_comp.csv")
  write.csv(results, output_file, row.names = FALSE)
  cat("\nResults saved to:", output_file, "\n")

  return(results)
}

# Run the script if called directly
if (!interactive()) {
  main()
}
