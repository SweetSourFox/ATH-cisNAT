#!/usr/bin/env Rscript
# =============================================================================
# Enhanced R Package Installer for AThNATCount
# =============================================================================
# 
# Description:
#   Automated installation of all required R packages with version checking,
#   parallel installation, and comprehensive error handling.
#
# Features:
#   - Automatic CRAN/Bioconductor detection
#   - Parallel package installation
#   - Version compatibility checking
#   - Retry mechanism for failed installations
#   - Progress tracking
# =============================================================================

# Setup
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(warn = 1)

# Enhanced logging
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Check R version
check_r_version <- function(min_version = "4.0.0") {
  current_version <- paste(R.version$major, R.version$minor, sep = ".")
  if (compareVersion(current_version, min_version) < 0) {
    stop(sprintf("R version %s or higher required. Current version: %s", 
                 min_version, current_version))
  }
  log_message(sprintf("R version %s detected", current_version))
}

# Define all required packages with sources and minimum versions
get_package_list <- function() {
  packages <- list(
    # CRAN packages
    CRAN = list(
      # Seurat and dependencies
      Seurat = list(version = "4.0.0", description = "Single-cell analysis"),
      Matrix = list(version = "1.3.0", description = "Sparse matrix support"),
      
      # Data manipulation
      dplyr = list(version = "1.0.0", description = "Data manipulation"),
      tidyr = list(version = "1.1.0", description = "Data tidying"),
      stringr = list(version = "1.4.0", description = "String manipulation"),
      data.table = list(version = "1.14.0", description = "Fast data manipulation"),
      
      # Visualization
      ggplot2 = list(version = "3.3.0", description = "Data visualization"),
      patchwork = list(version = "1.1.0", description = "Plot composition"),
      RColorBrewer = list(version = "1.1.0", description = "Color palettes"),
      gridExtra = list(version = "2.3", description = "Grid graphics"),
      scales = list(version = "1.1.0", description = "Scale functions"),
      viridis = list(version = "0.5.0", description = "Viridis color scales"),
      ggrepel = list(version = "0.9.0", description = "Text label positioning"),
      corrplot = list(version = "0.84", description = "Correlation plots"),
      pheatmap = list(version = "1.0.0", description = "Heatmaps"),
      VennDiagram = list(version = "1.6.0", description = "Venn diagrams"),
      ggpubr = list(version = "0.4.0", description = "Publication-ready plots"),
      cowplot = list(version = "1.1.0", description = "Plot composition"),
      ComplexHeatmap = list(version = "2.0.0", description = "Complex heatmaps"),
      circlize = list(version = "0.4.0", description = "Circular visualization"),
      
      # File I/O
      readr = list(version = "2.0.0", description = "Fast file reading"),
      yaml = list(version = "2.2.0", description = "YAML file support"),
      
      # Progress and parallel
      progress = list(version = "1.2.0", description = "Progress bars"),
      pbapply = list(version = "1.5.0", description = "Progress bar apply"),
      parallel = list(version = NULL, description = "Parallel processing"),
      
      # Additional utilities
      optparse = list(version = "1.7.0", description = "Command line parsing"),
      R.utils = list(version = "2.11.0", description = "R utilities")
    ),
    
    # Bioconductor packages
    Bioconductor = list(
      rtracklayer = list(version = "1.50.0", description = "GTF/GFF file support"),
      GenomicRanges = list(version = "1.42.0", description = "Genomic intervals"),
      BiocParallel = list(version = "1.24.0", description = "Parallel processing"),
      S4Vectors = list(version = "0.28.0", description = "S4 implementation"),
      IRanges = list(version = "2.24.0", description = "Integer ranges")
    )
  )
  
  return(packages)
}

# Check if package is installed with version
is_package_installed <- function(pkg, min_version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(FALSE)
  }
  
  if (!is.null(min_version)) {
    installed_version <- as.character(packageVersion(pkg))
    if (compareVersion(installed_version, min_version) < 0) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# Install package with retry mechanism
install_package_with_retry <- function(pkg, source = "CRAN", retries = 3) {
  for (attempt in 1:retries) {
    tryCatch({
      if (source == "CRAN") {
        install.packages(pkg, quiet = TRUE, dependencies = TRUE)
      } else if (source == "Bioconductor") {
        BiocManager::install(pkg, quiet = TRUE, update = FALSE, ask = FALSE)
      }
      
      # Verify installation
      if (requireNamespace(pkg, quietly = TRUE)) {
        return(TRUE)
      }
    }, error = function(e) {
      if (attempt < retries) {
        log_message(sprintf("Attempt %d failed for %s: %s", attempt, pkg, e$message), "WARNING")
        Sys.sleep(2)
      }
    })
  }
  
  return(FALSE)
}

# Install packages in parallel
install_packages_parallel <- function(packages, ncores = 2) {
  if (.Platform$OS.type == "windows") {
    ncores <- 1  # Parallel installation problematic on Windows
  }
  
  if (ncores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    
    parallel::clusterExport(cl, c("install_package_with_retry", "log_message"))
    
    results <- parallel::parLapply(cl, names(packages), function(pkg) {
      install_package_with_retry(pkg, "CRAN")
    })
    
    return(unlist(results))
  } else {
    # Sequential installation
    results <- lapply(names(packages), function(pkg) {
      install_package_with_retry(pkg, "CRAN")
    })
    return(unlist(results))
  }
}

# Main installation function
install_all_packages <- function() {
  log_message("=== AThNATCount R Package Installation ===", "INFO")
  
  # Check R version
  check_r_version()
  
  # Get package list
  all_packages <- get_package_list()
  
  # Count packages
  total_cran <- length(all_packages$CRAN)
  total_bioc <- length(all_packages$Bioconductor)
  total_packages <- total_cran + total_bioc
  
  log_message(sprintf("Total packages to check: %d (CRAN: %d, Bioconductor: %d)", 
                     total_packages, total_cran, total_bioc))
  
  # Check existing installations
  log_message("Checking installed packages...")
  
  packages_to_install <- list(CRAN = list(), Bioconductor = list())
  
  for (source in names(all_packages)) {
    for (pkg in names(all_packages[[source]])) {
      pkg_info <- all_packages[[source]][[pkg]]
      if (!is_package_installed(pkg, pkg_info$version)) {
        packages_to_install[[source]][[pkg]] <- pkg_info
      }
    }
  }
  
  # Count packages to install
  n_cran_install <- length(packages_to_install$CRAN)
  n_bioc_install <- length(packages_to_install$Bioconductor)
  n_total_install <- n_cran_install + n_bioc_install
  
  if (n_total_install == 0) {
    log_message("All packages are already installed!", "INFO")
    return(invisible(TRUE))
  }
  
  log_message(sprintf("Packages to install: %d (CRAN: %d, Bioconductor: %d)", 
                     n_total_install, n_cran_install, n_bioc_install))
  
  # Install CRAN packages
  if (n_cran_install > 0) {
    log_message("Installing CRAN packages...", "INFO")
    
    # Create progress bar
    pb <- txtProgressBar(min = 0, max = n_cran_install, style = 3)
    installed <- 0
    
    for (pkg in names(packages_to_install$CRAN)) {
      pkg_info <- packages_to_install$CRAN[[pkg]]
      log_message(sprintf("Installing %s: %s", pkg, pkg_info$description))
      
      success <- install_package_with_retry(pkg, "CRAN")
      
      if (success) {
        log_message(sprintf("Successfully installed %s", pkg), "INFO")
      } else {
        log_message(sprintf("Failed to install %s", pkg), "ERROR")
      }
      
      installed <- installed + 1
      setTxtProgressBar(pb, installed)
    }
    close(pb)
  }
  
  # Install Bioconductor packages
  if (n_bioc_install > 0) {
    log_message("Installing Bioconductor packages...", "INFO")
    
    # Ensure BiocManager is installed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      log_message("Installing BiocManager...", "INFO")
      install.packages("BiocManager", quiet = TRUE)
    }
    
    # Create progress bar
    pb <- txtProgressBar(min = 0, max = n_bioc_install, style = 3)
    installed <- 0
    
    for (pkg in names(packages_to_install$Bioconductor)) {
      pkg_info <- packages_to_install$Bioconductor[[pkg]]
      log_message(sprintf("Installing %s: %s", pkg, pkg_info$description))
      
      success <- install_package_with_retry(pkg, "Bioconductor")
      
      if (success) {
        log_message(sprintf("Successfully installed %s", pkg), "INFO")
      } else {
        log_message(sprintf("Failed to install %s", pkg), "ERROR")
      }
      
      installed <- installed + 1
      setTxtProgressBar(pb, installed)
    }
    close(pb)
  }
  
  # Verify all installations
  log_message("Verifying installations...", "INFO")
  
  failed_packages <- character()
  
  for (source in names(all_packages)) {
    for (pkg in names(all_packages[[source]])) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        failed_packages <- c(failed_packages, pkg)
      }
    }
  }
  
  if (length(failed_packages) > 0) {
    log_message(sprintf("Failed to install: %s", paste(failed_packages, collapse = ", ")), "ERROR")
    
    # Write failed packages to file
    writeLines(failed_packages, "failed_packages.txt")
    log_message("Failed packages written to failed_packages.txt", "INFO")
    
    return(invisible(FALSE))
  }
  
  log_message("All packages successfully installed!", "INFO")
  
  # Generate installation report
  generate_installation_report(all_packages)
  
  return(invisible(TRUE))
}

# Generate installation report
generate_installation_report <- function(all_packages) {
  log_message("Generating installation report...", "INFO")
  
  report_lines <- c(
    "AThNATCount R Package Installation Report",
    "=========================================",
    sprintf("Date: %s", Sys.Date()),
    sprintf("R Version: %s", R.version.string),
    sprintf("Platform: %s", R.version$platform),
    "",
    "Installed Packages:",
    "-------------------"
  )
  
  for (source in names(all_packages)) {
    report_lines <- c(report_lines, "", sprintf("%s Packages:", source))
    
    for (pkg in sort(names(all_packages[[source]]))) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        version <- as.character(packageVersion(pkg))
        desc <- all_packages[[source]][[pkg]]$description
        report_lines <- c(report_lines, 
                         sprintf("  ✓ %s (%s): %s", pkg, version, desc))
      }
    }
  }
  
  # System information
  report_lines <- c(report_lines, "", "System Information:", "-------------------")
  
  info_items <- list(
    "Total Memory" = sprintf("%.1f GB", as.numeric(system("grep MemTotal /proc/meminfo | awk '{print $2}'", intern = TRUE)) / 1024 / 1024),
    "CPU Cores" = parallel::detectCores(),
    "Library Path" = .libPaths()[1]
  )
  
  for (item in names(info_items)) {
    report_lines <- c(report_lines, sprintf("  %s: %s", item, info_items[[item]]))
  }
  
  # Write report
  report_file <- "R_packages_installation_report.txt"
  writeLines(report_lines, report_file)
  log_message(sprintf("Installation report saved to: %s", report_file), "INFO")
}

# Run installation if executed as script
if (!interactive()) {
  success <- install_all_packages()
  
  if (!success) {
    quit(status = 1)
  }
}