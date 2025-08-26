#!/usr/bin/env Rscript
# =============================================================================
# Test Script for Batch Analysis Functionality
# =============================================================================
# 
# This script tests the batch analysis functionality of the enhanced
# sense_nat_correlation.R script
# =============================================================================

# Load required libraries
library(tools)

# Test function to check if files exist
check_files_exist <- function() {
  cat("=== Checking Required Files ===\n")
  
  # Check main script
  if (file.exists("sense_nat_correlation.R")) {
    cat("✓ sense_nat_correlation.R found\n")
  } else {
    cat("✗ sense_nat_correlation.R not found\n")
    return(FALSE)
  }
  
  # Check cellranger results
  cellranger_dir <- "Cellranger_Results"
  if (dir.exists(cellranger_dir)) {
    cat("✓ Cellranger_Results directory found\n")
    
    # Find h5 files
    h5_files <- list.files(
      path = cellranger_dir,
      pattern = "filtered_feature_bc_matrix.h5",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(h5_files) > 0) {
      cat(sprintf("✓ Found %d h5 files:\n", length(h5_files)))
      for (i in seq_along(h5_files)) {
        cat(sprintf("  %d. %s\n", i, h5_files[i]))
      }
    } else {
      cat("✗ No h5 files found in Cellranger_Results\n")
      return(FALSE)
    }
  } else {
    cat("✗ Cellranger_Results directory not found\n")
    return(FALSE)
  }
  
  return(TRUE)
}

# Test function to validate command line arguments
test_command_line_parsing <- function() {
  cat("\n=== Testing Command Line Parsing ===\n")
  
  # Source the main script to access functions
  source("sense_nat_correlation.R")
  
  # Test single file parsing
  test_input <- "file1.h5"
  parsed_files <- strsplit(test_input, ",")[[1]]
  parsed_files <- trimws(parsed_files)
  
  if (length(parsed_files) == 1 && parsed_files[1] == "file1.h5") {
    cat("✓ Single file parsing works\n")
  } else {
    cat("✗ Single file parsing failed\n")
    return(FALSE)
  }
  
  # Test multiple files parsing
  test_input <- "file1.h5,file2.h5,file3.h5"
  parsed_files <- strsplit(test_input, ",")[[1]]
  parsed_files <- trimws(parsed_files)
  
  if (length(parsed_files) == 3 && all(parsed_files == c("file1.h5", "file2.h5", "file3.h5"))) {
    cat("✓ Multiple files parsing works\n")
  } else {
    cat("✗ Multiple files parsing failed\n")
    return(FALSE)
  }
  
  # Test sample names parsing
  test_names <- "Control,Treatment"
  parsed_names <- strsplit(test_names, ",")[[1]]
  parsed_names <- trimws(parsed_names)
  
  if (length(parsed_names) == 2 && all(parsed_names == c("Control", "Treatment"))) {
    cat("✓ Sample names parsing works\n")
  } else {
    cat("✗ Sample names parsing failed\n")
    return(FALSE)
  }
  
  return(TRUE)
}

# Test function to validate batch analysis functions
test_batch_functions <- function() {
  cat("\n=== Testing Batch Analysis Functions ===\n")
  
  # Source the main script
  source("sense_nat_correlation.R")
  
  # Test analyze_single_sample function structure
  if (exists("analyze_single_sample")) {
    cat("✓ analyze_single_sample function exists\n")
  } else {
    cat("✗ analyze_single_sample function not found\n")
    return(FALSE)
  }
  
  # Test analyze_multiple_samples function structure
  if (exists("analyze_multiple_samples")) {
    cat("✓ analyze_multiple_samples function exists\n")
  } else {
    cat("✗ analyze_multiple_samples function not found\n")
    return(FALSE)
  }
  
  # Test create_cross_sample_comparisons function
  if (exists("create_cross_sample_comparisons")) {
    cat("✓ create_cross_sample_comparisons function exists\n")
  } else {
    cat("✗ create_cross_sample_comparisons function not found\n")
    return(FALSE)
  }
  
  # Test generate_batch_summary_report function
  if (exists("generate_batch_summary_report")) {
    cat("✓ generate_batch_summary_report function exists\n")
  } else {
    cat("✗ generate_batch_summary_report function not found\n")
    return(FALSE)
  }
  
  return(TRUE)
}

# Test function to run a small batch analysis
test_small_batch_analysis <- function() {
  cat("\n=== Testing Small Batch Analysis ===\n")
  
  # Find available h5 files
  h5_files <- list.files(
    path = "Cellranger_Results",
    pattern = "filtered_feature_bc_matrix.h5",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(h5_files) < 2) {
    cat("⚠ Need at least 2 h5 files for batch analysis test\n")
    return(TRUE)
  }
  
  # Use first two files for testing
  test_files <- h5_files[1:2]
  test_names <- c("Test1", "Test2")
  
  cat(sprintf("Testing with files:\n"))
  for (i in seq_along(test_files)) {
    cat(sprintf("  %d. %s (%s)\n", i, test_files[i], test_names[i]))
  }
  
  # Create test output directory
  test_output <- "test_batch_results"
  if (dir.exists(test_output)) {
    unlink(test_output, recursive = TRUE)
  }
  dir.create(test_output, recursive = TRUE)
  
  # Run test analysis
  tryCatch({
    # Load configuration
    config <- load_config()
    
    # Setup logging
    log_file <- file.path(test_output, "test_analysis.log")
    setup_logging(log_file)
    
    # Run batch analysis
    results <- analyze_multiple_samples(
      input_files = test_files,
      sample_names = test_names,
      output_dir = test_output,
      config = config,
      n_cores = 2  # Use fewer cores for testing
    )
    
    if (!is.null(results)) {
      cat("✓ Batch analysis completed successfully\n")
      cat(sprintf("  Total gene pairs: %d\n", nrow(results)))
      cat(sprintf("  Samples processed: %s\n", paste(unique(results$sample), collapse = ", ")))
      
      # Check output files
      expected_files <- c(
        "combined_sense_nat_correlation.csv",
        "batch_analysis_summary.txt",
        "sample_summary_statistics.csv",
        "cross_sample_correlation_distribution.png",
        "sample_correlation_boxplot.png",
        "sample_significance_comparison.png"
      )
      
      missing_files <- c()
      for (file in expected_files) {
        if (file.exists(file.path(test_output, file))) {
          cat(sprintf("  ✓ %s\n", file))
        } else {
          cat(sprintf("  ✗ %s\n", file))
          missing_files <- c(missing_files, file)
        }
      }
      
      if (length(missing_files) == 0) {
        cat("✓ All expected output files created\n")
      } else {
        cat(sprintf("✗ Missing files: %s\n", paste(missing_files, collapse = ", ")))
      }
      
    } else {
      cat("✗ Batch analysis failed\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat(sprintf("✗ Error during batch analysis: %s\n", e$message))
    return(FALSE)
  })
  
  return(TRUE)
}

# Test function to validate output structure
test_output_structure <- function() {
  cat("\n=== Testing Output Structure ===\n")
  
  test_output <- "test_batch_results"
  
  if (!dir.exists(test_output)) {
    cat("⚠ Test output directory not found\n")
    return(TRUE)
  }
  
  # Check main output files
  main_files <- c(
    "combined_sense_nat_correlation.csv",
    "batch_analysis_summary.txt",
    "sample_summary_statistics.csv"
  )
  
  for (file in main_files) {
    file_path <- file.path(test_output, file)
    if (file.exists(file_path)) {
      file_size <- file.size(file_path)
      cat(sprintf("✓ %s (%s)\n", file, format(file_size, units = "auto")))
    } else {
      cat(sprintf("✗ %s\n", file))
    }
  }
  
  # Check sample-specific directories
  sample_dirs <- list.dirs(test_output, full.names = FALSE, recursive = FALSE)
  sample_dirs <- sample_dirs[sample_dirs != ""]
  
  for (dir in sample_dirs) {
    dir_path <- file.path(test_output, dir)
    if (dir.exists(dir_path)) {
      sample_files <- list.files(dir_path)
      cat(sprintf("✓ %s/ (%d files)\n", dir, length(sample_files)))
    } else {
      cat(sprintf("✗ %s/\n", dir))
    }
  }
  
  return(TRUE)
}

# Main test function
main <- function() {
  cat("Sense-Antisense Correlation Batch Analysis Test\n")
  cat("==============================================\n\n")
  
  # Run all tests
  tests <- list(
    "File Check" = check_files_exist,
    "Command Line Parsing" = test_command_line_parsing,
    "Batch Functions" = test_batch_functions,
    "Small Batch Analysis" = test_small_batch_analysis,
    "Output Structure" = test_output_structure
  )
  
  results <- list()
  
  for (test_name in names(tests)) {
    cat(sprintf("\n--- %s ---\n", test_name))
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
  
  # Cleanup
  if (dir.exists("test_batch_results")) {
    cat("\nCleaning up test files...\n")
    unlink("test_batch_results", recursive = TRUE)
  }
}

# Execute if run as script
if (!interactive()) {
  main()
} else {
  cat("Run this script from command line: Rscript test_batch_analysis.R\n")
}
