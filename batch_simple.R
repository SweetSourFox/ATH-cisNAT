#!/usr/bin/env Rscript
# =============================================================================
# Simple Test Script for Batch Analysis
# =============================================================================
# 
# This script provides a simple way to test the batch analysis functionality
# =============================================================================

# Load required libraries and functions without executing main
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(progress)
  library(yaml)
  library(parallel)
  library(data.table)
  library(corrplot)
  library(viridis)
  library(optparse)
  library(purrr)
  library(tidyr)
})

# Source the main script functions without executing main
source_lines <- readLines("sense_nat_correlation.R")
# Remove the execution part at the end
source_lines <- source_lines[1:(length(source_lines)-3)]
# Create a temporary file
temp_script <- tempfile(fileext = ".R")
writeLines(source_lines, temp_script)
source(temp_script)
unlink(temp_script)

# Function to find available cellranger files
find_available_files <- function() {
  cat("=== Searching for Cellranger Result Files ===\n")
  
  # Check if Cellranger_Results directory exists
  if (!dir.exists("~/Documents/test1/Cellranger_Results")) {
    cat("✗ Cellranger_Results directory not found\n")
    cat("Please ensure you have cellranger result files in the Cellranger_Results directory\n")
    return(NULL)
  }
  
  # Find all h5 files
  h5_files <- list.files(
    path = "~/Documents/test1/Cellranger_Results",
    pattern = "filtered_feature_bc_matrix.h5",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(h5_files) == 0) {
    cat("✗ No filtered_feature_bc_matrix.h5 files found\n")
    cat("Please ensure your cellranger results contain the expected h5 files\n")
    return(NULL)
  }
  
  cat(sprintf("✓ Found %d cellranger result files:\n", length(h5_files)))
  for (i in seq_along(h5_files)) {
    cat(sprintf("  %d. %s\n", i, h5_files[i]))
  }
  
  return(h5_files)
}


# Function to run batch analysis test
test_batch_analysis <- function() {
  cat("\n=== Testing Batch Analysis ===\n")
  
  h5_files <- find_available_files()
  if (is.null(h5_files)) return(FALSE)
  
  # Use all found files, and extract SRR/ERR/CRR-like sample names from their paths
  test_files <- h5_files
  # Extract sample names like "SRR8485805", "ERR1234567", or "CRR7654321" from the file path
  test_names <- sapply(test_files, function(f) {
    m <- regmatches(f, regexpr("(SRR|ERR|CRR)[0-9]+", f))
    if (length(m) > 0 && m != "") {
      return(m)
    } else {
      # fallback: use basename of parent dir
      return(basename(dirname(dirname(f))))
    }
  })
  cat(sprintf("Testing with %d files:\n", length(test_files)))
  for (i in seq_along(test_files)) {
    cat(sprintf("  %d. %s (%s)\n", i, test_files[i], test_names[i]))
  }
  
  # Create output directory
  output_dir <- "test_results/batch_analysis"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Load configuration
    config <- load_config()
    
    # Setup logging
    log_file <- file.path(output_dir, "test.log")
    setup_logging(log_file)
    
    # Run batch analysis
    result <- analyze_multiple_samples(
      input_files = test_files,
      sample_names = test_names,
      output_dir = output_dir,
      config = config,
      n_cores = 18
    )
    
    if (!is.null(result)) {
      cat("✓ Batch analysis completed successfully\n")
      cat(sprintf("  Total gene pairs: %d\n", nrow(result)))
      cat(sprintf("  Samples processed: %s\n", paste(unique(result$sample), collapse = ", ")))
      
      # Check output files
      expected_files <- c(
        "combined_sense_nat_correlation.csv",
        "batch_analysis_summary.txt",
        "sample_summary_statistics.csv"
      )
      
      missing_files <- c()
      for (file in expected_files) {
        if (file.exists(file.path(output_dir, file))) {
          cat(sprintf("  ✓ %s\n", file))
        } else {
          cat(sprintf("  ✗ %s\n", file))
          missing_files <- c(missing_files, file)
        }
      }
      
      if (length(missing_files) == 0) {
        cat("✓ All expected output files created\n")
      }
      
      return(TRUE)
    } else {
      cat("✗ Batch analysis failed\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat(sprintf("✗ Error during batch analysis: %s\n", e$message))
    return(FALSE)
  })
}


# Main test function
main <- function() {
  cat("Simple Batch Analysis Test\n")
  cat("=========================\n\n")
  
  # Check if main script exists
  if (!file.exists("sense_nat_correlation.R")) {
    cat("✗ sense_nat_correlation.R not found in current directory\n")
    cat("Please run this script from the directory containing sense_nat_correlation.R\n")
    return()
  }
  
  # Run tests
  tests <- list(
    "Batch Analysis" = test_batch_analysis,
  )
  
  results <- list()
  
  for (test_name in names(tests)) {
    cat(sprintf("\n--- %s Test ---\n", test_name))
    tryCatch({
      results[[test_name]] <- tests[[test_name]]()
    }, error = function(e) {
      cat(sprintf("✗ Test failed with error: %s\n", e$message))
      results[[test_name]] <- FALSE
    })
  }
  
  # Summary
  cat("\n=== Test Summary ===\n")
  passed <- sum(unlist(results))
  total <- length(results)
  
  for (test_name in names(results)) {
    status <- if (results[[test_name]]) "PASS" else "FAIL"
    cat(sprintf("%s: %s\n", test_name, status))
  }
  
  cat(sprintf("\nOverall: %d/%d tests passed\n", passed, total))
  
  if (passed == total) {
    cat("✓ All tests passed! Batch analysis functionality is working correctly.\n")
  } else {
    cat("✗ Some tests failed. Please check the output above for details.\n")
  }
  
  cat("\nTest results saved in 'test_results' directory\n")
}

# Execute if run as script
if (!interactive()) {
  main()
} else {
  cat("Run this script from command line: Rscript test_batch_simple.R\n")
}
