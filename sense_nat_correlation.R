#!/usr/bin/env Rscript
# =============================================================================
# Enhanced Sense-Antisense Correlation Analysis Pipeline
# =============================================================================
# 
# Description:
#   Computes correlation coefficients between sense and natural antisense
#   transcript (NAT) gene pairs from 10X Genomics single-cell RNA-seq data
#   Supports both single sample and batch analysis of multiple samples
#
# Features:
#   - Single sample and batch analysis modes
#   - Parallel processing for large datasets
#   - Multiple correlation methods (Pearson, Spearman, Kendall)
#   - Advanced filtering and quality control
#   - Interactive and batch visualization
#   - Comprehensive statistical analysis
#   - Cross-sample comparison and integration
# =============================================================================

# Load required libraries with error handling
load_libraries <- function() {
  required_packages <- c("Seurat", "Matrix", "dplyr", "ggplot2", "patchwork", 
                        "readr", "progress", "yaml", "parallel", "data.table",
                        "corrplot", "viridis", "optparse", "purrr", "tidyr")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Required package '%s' not installed. Run install_packages.R first.", pkg))
    }
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

suppressPackageStartupMessages(load_libraries())

# Configuration and logging setup
setup_logging <- function(log_file = NULL) {
  if (!is.null(log_file)) {
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
    
    # Create custom logging function
    log_func <- function(msg, level = "INFO") {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      formatted_msg <- sprintf("[%s] %s: %s", timestamp, level, msg)
      cat(formatted_msg, "\n")
      cat(formatted_msg, "\n", file = log_file, append = TRUE)
    }
  } else {
    log_func <- function(msg, level = "INFO") {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
    }
  }
  
  assign("log_message", log_func, envir = .GlobalEnv)
}

# Load configuration
load_config <- function(config_file = "config.yaml") {
  if (file.exists(config_file)) {
    config <- yaml::read_yaml(config_file)
  } else {
    # Default configuration
    config <- list(
      analysis = list(
        correlation = list(
          min_shared_cells = 3,
          min_expression = 0,
          correlation_method = "pearson"
        )
      ),
      performance = list(
        parallel_cores = 0,
        chunk_size = 1000
      ),
      output = list(
        figures = list(
          dpi = 300,
          format = "png",
          width = 10,
          height = 8
        )
      )
    )
  }
  return(config)
}

# Enhanced 10X data reading with validation
read_10x_data <- function(h5_file, validate = TRUE) {
  log_message(sprintf("Reading 10X data from: %s", h5_file))
  
  if (!file.exists(h5_file)) {
    stop(sprintf("H5 file not found: %s", h5_file))
  }
  
  # Read data
  raw_data <- Read10X_h5(h5_file)
  
  # Handle multiple assays
  if (is.list(raw_data)) {
    if ("Gene Expression" %in% names(raw_data)) {
      expr_matrix <- raw_data[["Gene Expression"]]
    } else {
      expr_matrix <- raw_data[[1]]
      log_message("Using first assay from multi-assay data", "WARNING")
    }
  } else {
    expr_matrix <- raw_data
  }
  
  # Validation
  if (validate) {
    n_genes <- nrow(expr_matrix)
    n_cells <- ncol(expr_matrix)
    sparsity <- 1 - (length(expr_matrix@x) / (n_genes * n_cells))
    
    log_message(sprintf("Matrix dimensions: %d genes × %d cells", n_genes, n_cells))
    log_message(sprintf("Sparsity: %.2f%%", sparsity * 100))
    
    if (n_genes < 100 || n_cells < 100) {
      warning("Very small dataset detected. Results may be unreliable.")
    }
  }
  
  return(expr_matrix)
}

# Identify sense/NAT pairs with enhanced matching
identify_sense_nat_pairs <- function(gene_names, nat_suffix = "_NAT") {
  log_message("Identifying sense/NAT gene pairs...")
  
  # Find NAT genes
  nat_mask <- grepl(paste0(nat_suffix, "$"), gene_names)
  nat_genes <- gene_names[nat_mask]
  
  # Match to sense genes
  sense_genes <- sub(paste0(nat_suffix, "$"), "", nat_genes)
  
  # Validate pairs
  valid_pairs <- sense_genes %in% gene_names
  nat_genes <- nat_genes[valid_pairs]
  sense_genes <- sense_genes[valid_pairs]
  
  # Create pair dataframe
  pairs_df <- data.frame(
    sense_gene = sense_genes,
    nat_gene = nat_genes,
    sense_idx = match(sense_genes, gene_names),
    nat_idx = match(nat_genes, gene_names),
    stringsAsFactors = FALSE
  )
  
  log_message(sprintf("Found %d valid sense/NAT pairs", nrow(pairs_df)))
  
  return(pairs_df)
}

# Process function for each chunk (defined outside to ensure proper export)
process_chunk_correlation <- function(indices, expr_matrix, pairs_df, min_cells, min_expr, method) {
  results <- list()
  
  for (i in indices) {
    sense_idx <- pairs_df$sense_idx[i]
    nat_idx <- pairs_df$nat_idx[i]
    
    # Extract expression vectors
    sense_expr <- expr_matrix[sense_idx, ]
    nat_expr <- expr_matrix[nat_idx, ]
    
    # Find cells expressing both genes
    shared_mask <- (sense_expr > min_expr) & (nat_expr > min_expr)
    n_shared <- sum(shared_mask)
    
    if (n_shared >= min_cells) {
      # Extract shared expression values
      sense_vals <- as.numeric(sense_expr[shared_mask])
      nat_vals <- as.numeric(nat_expr[shared_mask])
      
      # Compute correlation
      if (method == "pearson") {
        cor_result <- cor.test(sense_vals, nat_vals, method = "pearson")
      } else if (method == "spearman") {
        cor_result <- cor.test(sense_vals, nat_vals, method = "spearman")
      } else if (method == "kendall") {
        cor_result <- cor.test(sense_vals, nat_vals, method = "kendall")
      }
      
      # Additional statistics
      mean_sense <- mean(sense_vals)
      mean_nat <- mean(nat_vals)
      sd_sense <- sd(sense_vals)
      sd_nat <- sd(nat_vals)
      
      results[[length(results) + 1]] <- list(
        pair_idx = i,
        correlation = cor_result$estimate,
        p_value = cor_result$p.value,
        n_cells = n_shared,
        mean_sense = mean_sense,
        mean_nat = mean_nat,
        sd_sense = sd_sense,
        sd_nat = sd_nat,
        method = method
      )
    }
  }
  
  return(results)
}

# Parallel correlation computation with multiple methods
compute_correlations_parallel <- function(expr_matrix, pairs_df, config, n_cores = NULL) {
  log_message("Computing correlations with parallel processing...")
  
  # Get parameters
  min_cells <- config$analysis$correlation$min_shared_cells
  min_expr <- config$analysis$correlation$min_expression
  method <- config$analysis$correlation$correlation_method
  
  # Determine number of cores
  if (is.null(n_cores)) {
    n_cores <- config$performance$parallel_cores
    if (n_cores == 0) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
  }
  
  log_message(sprintf("Using %d cores for parallel processing", n_cores))
  
  # Prepare for parallel processing
  n_pairs <- nrow(pairs_df)
  chunk_size <- config$performance$chunk_size
  n_chunks <- ceiling(n_pairs / chunk_size)
  
  # Create chunks
  chunks <- split(1:n_pairs, cut(1:n_pairs, n_chunks, labels = FALSE))
  
  # Initialize cluster
  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    
    # Load required libraries on cluster
    clusterEvalQ(cl, {
      library(Matrix)
      library(stats)
    })
    
    # Export the function to the cluster
    clusterExport(cl, "process_chunk_correlation")
    
    # Run parallel processing with explicit variable passing
    chunk_results <- parLapply(cl, chunks, function(chunk, expr_matrix, pairs_df, min_cells, min_expr, method) {
      process_chunk_correlation(chunk, expr_matrix, pairs_df, min_cells, min_expr, method)
    }, expr_matrix, pairs_df, min_cells, min_expr, method)
  
    stopCluster(cl)
  } else {
    # Sequential processing with progress bar
    pb <- progress_bar$new(
      total = n_chunks,
      format = "Processing [:bar] :current/:total (:percent) ETA: :eta"
    )
    
    chunk_results <- lapply(chunks, function(chunk) {
      result <- process_chunk_correlation(chunk, expr_matrix, pairs_df, min_cells, min_expr, method)
      pb$tick()
      return(result)
    })
  }
  
  # Combine results
  all_results <- unlist(chunk_results, recursive = FALSE)
  
  # Convert to dataframe
  if (length(all_results) > 0) {
    results_df <- do.call(rbind, lapply(all_results, function(x) {
      data.frame(
        sense_gene = pairs_df$sense_gene[x$pair_idx],
        nat_gene = pairs_df$nat_gene[x$pair_idx],
        correlation = x$correlation,
        p_value = x$p_value,
        n_cells = x$n_cells,
        mean_sense_expr = x$mean_sense,
        mean_nat_expr = x$mean_nat,
        sd_sense_expr = x$sd_sense,
        sd_nat_expr = x$sd_nat,
        method = x$method,
        stringsAsFactors = FALSE
      )
    }))
    
    # Add multiple testing correction
    results_df$p_adjusted <- p.adjust(results_df$p_value, method = "fdr")
    
    # Add significance categories
    results_df$significance <- cut(
      results_df$p_adjusted,
      breaks = c(0, 0.001, 0.01, 0.05, 1),
      labels = c("***", "**", "*", "ns"),
      include.lowest = TRUE
    )
    
  } else {
    results_df <- data.frame()
  }
  
  log_message(sprintf("Computed correlations for %d pairs", nrow(results_df)))
  
  return(results_df)
}

# Enhanced visualization functions
create_visualizations <- function(results_df, output_dir, config, sample_name = NULL) {
  log_message("Creating visualizations...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get figure parameters
  fig_config <- config$output$figures
  
  # Add sample name to titles if provided
  title_suffix <- if (!is.null(sample_name)) sprintf(" - %s", sample_name) else ""
  
  # 1. Correlation distribution histogram with density
  p1 <- ggplot(results_df, aes(x = correlation)) +
    geom_histogram(aes(y = ..density..), bins = 50, 
                   fill = "steelblue", color = "white", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = median(results_df$correlation, na.rm = TRUE), 
               linetype = "dashed", color = "green", size = 1) +
    theme_minimal() +
    labs(
      title = paste0("Distribution of Sense/NAT Expression Correlations", title_suffix),
      subtitle = sprintf("n = %d gene pairs | Method: %s", 
                        nrow(results_df), unique(results_df$method)[1]),
      x = "Correlation Coefficient",
      y = "Density"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  
  filename <- if (!is.null(sample_name)) {
    sprintf("correlation_distribution_%s.png", sample_name)
  } else {
    "correlation_distribution.png"
  }
  
  ggsave(
    file.path(output_dir, filename),
    p1, width = fig_config$width, height = fig_config$height, 
    dpi = fig_config$dpi
  )
  
  # 2. Enhanced volcano plot
  results_df$log_p <- -log10(results_df$p_value)
  results_df$log_p[is.infinite(results_df$log_p)] <- max(results_df$log_p[!is.infinite(results_df$log_p)]) + 1
  
  p2 <- ggplot(results_df, aes(x = correlation, y = log_p)) +
    geom_point(aes(color = significance, size = n_cells), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c("***" = "darkred", "**" = "red", 
                                  "*" = "orange", "ns" = "gray50")) +
    scale_size_continuous(range = c(1, 4)) +
    theme_minimal() +
    labs(
      title = paste0("Volcano Plot of Sense/NAT Correlations", title_suffix),
      x = "Correlation Coefficient",
      y = "-log10(p-value)",
      color = "Significance",
      size = "# Cells"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  filename <- if (!is.null(sample_name)) {
    sprintf("correlation_volcano_%s.png", sample_name)
  } else {
    "correlation_volcano.png"
  }
  
  ggsave(
    file.path(output_dir, filename),
    p2, width = fig_config$width * 1.2, height = fig_config$height, 
    dpi = fig_config$dpi
  )
  
  # 3. Correlation vs expression level
  p3 <- ggplot(results_df, aes(x = log10(mean_sense_expr + 1), 
                               y = log10(mean_nat_expr + 1))) +
    geom_point(aes(color = correlation), alpha = 0.7) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste0("Sense vs NAT Expression Levels", title_suffix),
      x = "log10(Mean Sense Expression + 1)",
      y = "log10(Mean NAT Expression + 1)",
      color = "Correlation"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  filename <- if (!is.null(sample_name)) {
    sprintf("expression_scatter_%s.png", sample_name)
  } else {
    "expression_scatter.png"
  }
  
  ggsave(
    file.path(output_dir, filename),
    p3, width = fig_config$width, height = fig_config$height, 
    dpi = fig_config$dpi
  )
  
  # 4. Top correlations heatmap
  if (nrow(results_df) >= 20) {
    # Get top positive and negative correlations
    top_pos <- results_df %>% 
      filter(p_adjusted < 0.05) %>%
      arrange(desc(correlation)) %>%
      head(10)
    
    top_neg <- results_df %>% 
      filter(p_adjusted < 0.05) %>%
      arrange(correlation) %>%
      head(10)
    
    top_genes <- rbind(top_pos, top_neg)
    
    if (nrow(top_genes) > 0) {
      # Create correlation matrix for visualization
      cor_matrix <- matrix(0, nrow = nrow(top_genes), ncol = 2)
      colnames(cor_matrix) <- c("Correlation", "Significance")
      rownames(cor_matrix) <- paste(top_genes$sense_gene, "vs", 
                                   sub("_NAT$", "", top_genes$nat_gene))
      
      cor_matrix[, 1] <- top_genes$correlation
      cor_matrix[, 2] <- -log10(top_genes$p_adjusted)
      
      filename <- if (!is.null(sample_name)) {
        sprintf("top_correlations_heatmap_%s.png", sample_name)
      } else {
        "top_correlations_heatmap.png"
      }
      
      png(file.path(output_dir, filename),
          width = fig_config$width * 100, 
          height = fig_config$height * 100,
          res = fig_config$dpi)
      
      corrplot(cor_matrix, method = "color", type = "full",
               col = colorRampPalette(c("blue", "white", "red"))(100),
               tl.col = "black", tl.srt = 45, tl.cex = 0.8,
               cl.cex = 0.8, is.corr = FALSE,
               title = paste0("Top Significant Sense-NAT Correlations", title_suffix))
      
      dev.off()
    }
  }
  
  log_message(sprintf("Visualizations saved to: %s", output_dir))
}

# Export expression data for downstream analysis
export_expression_pairs <- function(expr_matrix, results_df, output_file, 
                                   sample_size = NULL) {
  log_message("Exporting expression pair data...")
  
  if (!is.null(sample_size) && sample_size < nrow(results_df)) {
    # Sample gene pairs
    sample_idx <- sample(1:nrow(results_df), sample_size)
    results_sample <- results_df[sample_idx, ]
  } else {
    results_sample <- results_df
  }
  
  # Initialize data list
  expr_data_list <- list()
  
  pb <- progress_bar$new(
    total = nrow(results_sample),
    format = "Exporting [:bar] :current/:total (:percent) ETA: :eta"
  )
  
  for (i in 1:nrow(results_sample)) {
    sense_gene <- results_sample$sense_gene[i]
    nat_gene <- results_sample$nat_gene[i]
    
    # Get expression data
    sense_idx <- which(rownames(expr_matrix) == sense_gene)
    nat_idx <- which(rownames(expr_matrix) == nat_gene)
    
    sense_expr <- expr_matrix[sense_idx, ]
    nat_expr <- expr_matrix[nat_idx, ]
    
    # Find expressing cells
    expressing_cells <- which((sense_expr > 0) | (nat_expr > 0))
    
    if (length(expressing_cells) > 0) {
      expr_df <- data.frame(
        cell_id = colnames(expr_matrix)[expressing_cells],
        sense_expr = as.numeric(sense_expr[expressing_cells]),
        nat_expr = as.numeric(nat_expr[expressing_cells]),
        gene_pair = paste0(sense_gene, "_pair"),
        stringsAsFactors = FALSE
      )
      
      expr_data_list[[i]] <- expr_df
    }
    
    pb$tick()
  }
  
  # Combine all expression data
  all_expr_data <- do.call(rbind, expr_data_list)
  
  # Save as compressed CSV
  data.table::fwrite(all_expr_data, output_file, compress = "gzip")
  
  log_message(sprintf("Expression data exported to: %s", output_file))
}

# Process LRR-RLK subset if gene list provided
process_lrr_rlk_subset <- function(results_df, lrr_rlk_file, output_file) {
  if (file.exists(lrr_rlk_file)) {
    log_message("Processing LRR-RLK gene subset...")
    
    # Read LRR-RLK gene list
    lrr_rlk_genes <- readLines(lrr_rlk_file)
    lrr_rlk_genes <- trimws(lrr_rlk_genes)
    
    # Filter results
    lrr_rlk_results <- results_df %>%
      filter(sense_gene %in% lrr_rlk_genes)
    
    if (nrow(lrr_rlk_results) > 0) {
      write_csv(lrr_rlk_results, output_file)
      log_message(sprintf("LRR-RLK results saved to: %s", output_file))
      log_message(sprintf("Found %d LRR-RLK gene pairs", nrow(lrr_rlk_results)))
    } else {
      log_message("No LRR-RLK genes found in results", "WARNING")
    }
  }
}

# Generate comprehensive summary report
generate_summary_report <- function(results_df, output_file, sample_name = NULL) {
  log_message("Generating summary report...")
  
  title_suffix <- if (!is.null(sample_name)) sprintf(" - %s", sample_name) else ""
  
  report_lines <- c(
    paste0("Sense-Antisense Correlation Analysis Summary", title_suffix),
    "===========================================",
    sprintf("Date: %s", Sys.Date()),
    sprintf("Total gene pairs analyzed: %d", nrow(results_df)),
    "",
    "Correlation Statistics:",
    sprintf("  Mean correlation: %.3f", mean(results_df$correlation, na.rm = TRUE)),
    sprintf("  Median correlation: %.3f", median(results_df$correlation, na.rm = TRUE)),
    sprintf("  SD correlation: %.3f", sd(results_df$correlation, na.rm = TRUE)),
    sprintf("  Range: [%.3f, %.3f]", 
            min(results_df$correlation, na.rm = TRUE),
            max(results_df$correlation, na.rm = TRUE)),
    "",
    "Significance Statistics:",
    sprintf("  Significant pairs (p < 0.05): %d (%.1f%%)", 
            sum(results_df$p_value < 0.05, na.rm = TRUE),
            sum(results_df$p_value < 0.05, na.rm = TRUE) / nrow(results_df) * 100),
    sprintf("  Significant pairs (FDR < 0.05): %d (%.1f%%)", 
            sum(results_df$p_adjusted < 0.05, na.rm = TRUE),
            sum(results_df$p_adjusted < 0.05, na.rm = TRUE) / nrow(results_df) * 100),
    "",
    "Correlation Direction:",
    sprintf("  Positive correlations: %d (%.1f%%)", 
            sum(results_df$correlation > 0, na.rm = TRUE),
            sum(results_df$correlation > 0, na.rm = TRUE) / nrow(results_df) * 100),
    sprintf("  Negative correlations: %d (%.1f%%)", 
            sum(results_df$correlation < 0, na.rm = TRUE),
            sum(results_df$correlation < 0, na.rm = TRUE) / nrow(results_df) * 100),
    "",
    "Top 5 Positive Correlations:"
  )
  
  # Add top correlations
  top_pos <- results_df %>% 
    arrange(desc(correlation)) %>% 
    head(5)
  
  for (i in 1:nrow(top_pos)) {
    report_lines <- c(report_lines,
      sprintf("  %d. %s: r=%.3f, p=%.2e, n=%d cells",
              i, top_pos$sense_gene[i], top_pos$correlation[i],
              top_pos$p_value[i], top_pos$n_cells[i]))
  }
  
  report_lines <- c(report_lines, "", "Top 5 Negative Correlations:")
  
  top_neg <- results_df %>% 
    arrange(correlation) %>% 
    head(5)
  
  for (i in 1:nrow(top_neg)) {
    report_lines <- c(report_lines,
      sprintf("  %d. %s: r=%.3f, p=%.2e, n=%d cells",
              i, top_neg$sense_gene[i], top_neg$correlation[i],
              top_neg$p_value[i], top_neg$n_cells[i]))
  }
  
  writeLines(report_lines, output_file)
  log_message(sprintf("Summary report saved to: %s", output_file))
}

# Batch analysis function
analyze_single_sample <- function(h5_file, output_dir, config, sample_name = NULL, 
                                 lrr_rlk_file = NULL, export_expression = FALSE, 
                                 n_cores = NULL) {
  log_message(sprintf("=== Analyzing sample: %s ===", if (!is.null(sample_name)) sample_name else basename(h5_file)))
  
  # Create sample-specific output directory
  if (!is.null(sample_name)) {
    sample_output_dir <- file.path(output_dir, sample_name)
  } else {
    sample_output_dir <- file.path(output_dir, tools::file_path_sans_ext(basename(h5_file)))
  }
  
  if (!dir.exists(sample_output_dir)) {
    dir.create(sample_output_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Read 10X data
    expr_matrix <- read_10x_data(h5_file)
    
    # Identify sense/NAT pairs
    gene_names <- rownames(expr_matrix)
    pairs_df <- identify_sense_nat_pairs(gene_names)
    
    if (nrow(pairs_df) == 0) {
      warning(sprintf("No sense/NAT pairs found in sample: %s", if (!is.null(sample_name)) sample_name else basename(h5_file)))
      return(NULL)
    }
    
    # Compute correlations
    results_df <- compute_correlations_parallel(
      expr_matrix, pairs_df, config, 
      n_cores = n_cores
    )
    
    # Add sample information
    results_df$sample <- if (!is.null(sample_name)) sample_name else basename(h5_file)
    
    # Save results
    output_file <- file.path(sample_output_dir, "sense_nat_correlation.csv")
    write_csv(results_df, output_file)
    log_message(sprintf("Results saved to: %s", output_file))
    
    # Create visualizations
    create_visualizations(results_df, sample_output_dir, config, sample_name)
    
    # Export expression data if requested
    if (export_expression) {
      expr_output <- file.path(sample_output_dir, "sense_nat_expression_pairs.csv.gz")
      export_expression_pairs(expr_matrix, results_df, expr_output)
    }
    
    # Process LRR-RLK subset if provided
    if (!is.null(lrr_rlk_file)) {
      lrr_output <- file.path(sample_output_dir, "sense_nat_correlation_LRR_RLK.csv")
      process_lrr_rlk_subset(results_df, lrr_rlk_file, lrr_output)
    }
    
    # Generate summary report
    report_file <- file.path(sample_output_dir, "analysis_summary.txt")
    generate_summary_report(results_df, report_file, sample_name)
    
    log_message(sprintf("=== Sample analysis complete: %s ===", if (!is.null(sample_name)) sample_name else basename(h5_file)))
    
    return(results_df)
    
  }, error = function(e) {
    log_message(sprintf("ERROR in sample %s: %s", 
                       if (!is.null(sample_name)) sample_name else basename(h5_file), 
                       e$message), "ERROR")
    return(NULL)
  })
}

# Batch analysis with cross-sample comparison
analyze_multiple_samples <- function(input_files, sample_names = NULL, output_dir, 
                                   config, lrr_rlk_file = NULL, export_expression = FALSE,
                                   n_cores = NULL) {
  log_message("=== Starting Batch Analysis ===")
  log_message(sprintf("Number of samples to analyze: %d", length(input_files)))
  
  # Validate input
  if (length(input_files) == 0) {
    stop("No input files provided")
  }
  
  # Generate sample names if not provided
  if (is.null(sample_names)) {
    sample_names <- tools::file_path_sans_ext(basename(input_files))
  }
  
  if (length(sample_names) != length(input_files)) {
    stop("Number of sample names must match number of input files")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Analyze each sample
  all_results <- list()
  
  for (i in seq_along(input_files)) {
    log_message(sprintf("Processing sample %d/%d: %s", i, length(input_files), sample_names[i]))
    
    result <- analyze_single_sample(
      h5_file = input_files[i],
      output_dir = output_dir,
      config = config,
      sample_name = sample_names[i],
      lrr_rlk_file = lrr_rlk_file,
      export_expression = export_expression,
      n_cores = n_cores
    )
    
    if (!is.null(result)) {
      all_results[[sample_names[i]]] <- result
    }
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    
    # Save combined results
    combined_file <- file.path(output_dir, "combined_sense_nat_correlation.csv")
    write_csv(combined_results, combined_file)
    log_message(sprintf("Combined results saved to: %s", combined_file))
    
    # Generate cross-sample comparison visualizations
    create_cross_sample_comparisons(combined_results, output_dir, config)
    
    # Generate batch summary report
    generate_batch_summary_report(combined_results, output_dir, sample_names)
    
    log_message("=== Batch Analysis Complete ===")
    
    return(combined_results)
  } else {
    log_message("No successful analyses completed", "ERROR")
    return(NULL)
  }
}

# Create cross-sample comparison visualizations
create_cross_sample_comparisons <- function(combined_results, output_dir, config) {
  log_message("Creating cross-sample comparison visualizations...")
  
  fig_config <- config$output$figures
  
  # 1. Correlation distribution comparison across samples
  p1 <- ggplot(combined_results, aes(x = correlation, fill = sample)) +
    geom_density(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    theme_minimal() +
    labs(
      title = "Correlation Distribution Comparison Across Samples",
      x = "Correlation Coefficient",
      y = "Density",
      fill = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  ggsave(
    file.path(output_dir, "cross_sample_correlation_distribution.png"),
    p1, width = fig_config$width * 1.2, height = fig_config$height, 
    dpi = fig_config$dpi
  )
  
  # 2. Sample comparison boxplot
  p2 <- ggplot(combined_results, aes(x = sample, y = correlation, fill = sample)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    theme_minimal() +
    labs(
      title = "Correlation Distribution by Sample",
      x = "Sample",
      y = "Correlation Coefficient",
      fill = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    file.path(output_dir, "sample_correlation_boxplot.png"),
    p2, width = fig_config$width * 1.5, height = fig_config$height, 
    dpi = fig_config$dpi
  )
  
  # 3. Significance comparison across samples
  sig_summary <- combined_results %>%
    group_by(sample) %>%
    summarise(
      total_pairs = n(),
      significant_pairs = sum(p_adjusted < 0.05, na.rm = TRUE),
      significant_percent = significant_pairs / total_pairs * 100,
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE)
    )
  
  p3 <- ggplot(sig_summary, aes(x = sample, y = significant_percent, fill = sample)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", significant_percent)), 
              vjust = -0.5, size = 3) +
    theme_minimal() +
    labs(
      title = "Percentage of Significant Correlations by Sample",
      x = "Sample",
      y = "Percentage of Significant Pairs (%)",
      fill = "Sample"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    file.path(output_dir, "sample_significance_comparison.png"),
    p3, width = fig_config$width * 1.5, height = fig_config$height, 
    dpi = fig_config$dpi
  )
  
  # Save summary statistics
  write_csv(sig_summary, file.path(output_dir, "sample_summary_statistics.csv"))
  
  log_message("Cross-sample comparison visualizations created")
}

# Generate batch summary report
generate_batch_summary_report <- function(combined_results, output_dir, sample_names) {
  log_message("Generating batch summary report...")
  
  report_lines <- c(
    "Batch Sense-Antisense Correlation Analysis Summary",
    "================================================",
    sprintf("Date: %s", Sys.Date()),
    sprintf("Number of samples analyzed: %d", length(sample_names)),
    sprintf("Samples: %s", paste(sample_names, collapse = ", ")),
    sprintf("Total gene pairs across all samples: %d", nrow(combined_results)),
    "",
    "Overall Statistics:",
    sprintf("  Mean correlation across all samples: %.3f", 
            mean(combined_results$correlation, na.rm = TRUE)),
    sprintf("  Median correlation across all samples: %.3f", 
            median(combined_results$correlation, na.rm = TRUE)),
    sprintf("  SD correlation across all samples: %.3f", 
            sd(combined_results$correlation, na.rm = TRUE)),
    "",
    "Sample-wise Statistics:"
  )
  
  # Add sample-wise statistics
  sample_stats <- combined_results %>%
    group_by(sample) %>%
    summarise(
      n_pairs = n(),
      mean_cor = mean(correlation, na.rm = TRUE),
      median_cor = median(correlation, na.rm = TRUE),
      sd_cor = sd(correlation, na.rm = TRUE),
      n_sig = sum(p_adjusted < 0.05, na.rm = TRUE),
      pct_sig = n_sig / n_pairs * 100
    )
  
  for (i in 1:nrow(sample_stats)) {
    report_lines <- c(report_lines,
      sprintf("  %s:", sample_stats$sample[i]),
      sprintf("    Pairs: %d, Mean cor: %.3f, Significant: %d (%.1f%%)",
              sample_stats$n_pairs[i], sample_stats$mean_cor[i],
              sample_stats$n_sig[i], sample_stats$pct_sig[i])
    )
  }
  
  # Add top correlations across all samples
  report_lines <- c(report_lines, "", "Top 5 Correlations Across All Samples:")
  
  top_all <- combined_results %>%
    arrange(desc(correlation)) %>%
    head(5)
  
  for (i in 1:nrow(top_all)) {
    report_lines <- c(report_lines,
      sprintf("  %d. %s (%s): r=%.3f, p=%.2e, n=%d cells",
              i, top_all$sense_gene[i], top_all$sample[i], 
              top_all$correlation[i], top_all$p_value[i], top_all$n_cells[i]))
  }
  
  writeLines(report_lines, file.path(output_dir, "batch_analysis_summary.txt"))
  log_message("Batch summary report generated")
}

# Main analysis pipeline
main <- function() {
  # Parse command line arguments
  parser <- OptionParser(description = "Sense-Antisense Correlation Analysis (Single or Batch)")
  
  parser <- add_option(parser, c("-i", "--input"), 
                      type = "character",
                      help = "Input 10X h5 file path(s) - single file or comma-separated list")
  
  parser <- add_option(parser, c("-n", "--names"), 
                      type = "character",
                      help = "Sample names - comma-separated list (optional, uses filenames if not provided)")
  
  parser <- add_option(parser, c("-o", "--output"), 
                      type = "character", 
                      default = "sa_correlation_results",
                      help = "Output directory [default: %default]")
  
  parser <- add_option(parser, c("-c", "--config"), 
                      type = "character", 
                      default = "config.yaml",
                      help = "Configuration file [default: %default]")
  
  parser <- add_option(parser, c("-l", "--lrr-rlk"), 
                      type = "character",
                      help = "LRR-RLK gene list file (optional)")
  
  parser <- add_option(parser, c("-e", "--export-expression"), 
                      action = "store_true", 
                      default = FALSE,
                      help = "Export expression pair data")
  
  parser <- add_option(parser, c("-m", "--method"), 
                      type = "character", 
                      default = NULL,
                      help = "Correlation method (pearson/spearman/kendall)")
  
  parser <- add_option(parser, c("-p", "--cores"), 
                      type = "integer", 
                      default = NULL,
                      help = "Number of CPU cores to use")
  
  # Parse arguments
  args <- parse_args(parser)
  
  # Validate required arguments
  if (is.null(args$input)) {
    stop("Input h5 file(s) required. Use -h for help.")
  }
  
  # Parse input files
  input_files <- strsplit(args$input, ",")[[1]]
  input_files <- trimws(input_files)
  
  # Validate input files exist
  for (file in input_files) {
    if (!file.exists(file)) {
      stop(sprintf("Input file not found: %s", file))
    }
  }
  
  # Parse sample names if provided
  sample_names <- NULL
  if (!is.null(args$names)) {
    sample_names <- strsplit(args$names, ",")[[1]]
    sample_names <- trimws(sample_names)
    
    if (length(sample_names) != length(input_files)) {
      stop("Number of sample names must match number of input files")
    }
  }
  
  # Load configuration
  config <- load_config(args$config)
  
  # Override config with command line arguments
  if (!is.null(args$method)) {
    config$analysis$correlation$correlation_method <- args$method
  }
  
  # Setup logging
  log_file <- file.path(args$output, "analysis.log")
  setup_logging(log_file)
  
  log_message("=== Sense-Antisense Correlation Analysis ===")
  log_message(sprintf("Number of input files: %d", length(input_files)))
  log_message(sprintf("Input files: %s", paste(input_files, collapse = ", ")))
  log_message(sprintf("Output directory: %s", args$output))
  
  # Create output directory
  if (!dir.exists(args$output)) {
    dir.create(args$output, recursive = TRUE)
  }
  
  tryCatch({
    if (length(input_files) == 1) {
      # Single sample analysis
      log_message("Running single sample analysis...")
      result <- analyze_single_sample(
        h5_file = input_files[1],
        output_dir = args$output,
        config = config,
        sample_name = sample_names[1],
        lrr_rlk_file = args$lrr_rlk,
        export_expression = args$export_expression,
        n_cores = args$cores
      )
    } else {
      # Batch analysis
      log_message("Running batch analysis...")
      result <- analyze_multiple_samples(
        input_files = input_files,
        sample_names = sample_names,
        output_dir = args$output,
        config = config,
        lrr_rlk_file = args$lrr_rlk,
        export_expression = args$export_expression,
        n_cores = args$cores
      )
    }
    
    if (!is.null(result)) {
      log_message("=== Analysis Complete ===")
    } else {
      log_message("Analysis completed with errors", "WARNING")
    }
    
  }, error = function(e) {
    log_message(sprintf("ERROR: %s", e$message), "ERROR")
    quit(status = 1)
  })
}

# Execute if run as script
if (!interactive()) {
  main()
}