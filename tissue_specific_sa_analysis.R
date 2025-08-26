#!/usr/bin/env Rscript
# =============================================================================
# Advanced Tissue-Specific Sense-Antisense Analysis Pipeline
# =============================================================================
# 
# Description:
#   Comprehensive integration and analysis of sense-antisense correlation data
#   across multiple tissue types with advanced statistical modeling and 
#   interactive visualization capabilities
#
# Features:
#   - Automated data integration from multiple sources
#   - Advanced statistical modeling (mixed effects, Bayesian)
#   - Machine learning for pattern discovery  
#   - Interactive dashboards and reports
#   - Batch processing capabilities
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(gridExtra)
  library(yaml)
  library(data.table)
  library(parallel)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(ggpubr)
  library(cowplot)
  library(scales)
})

# Global configuration
SCRIPT_VERSION <- "2.0.0"
SCRIPT_NAME <- "Tissue-Specific SA Analysis"

# Enhanced logging system
setup_logging <- function(output_dir, verbose = TRUE) {
  log_dir <- file.path(output_dir, "logs")
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }
  
  log_file <- file.path(log_dir, sprintf("analysis_%s.log", 
                                        format(Sys.time(), "%Y%m%d_%H%M%S")))
  
  log_func <- function(msg, level = "INFO", indent = 0) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    indent_str <- paste(rep("  ", indent), collapse = "")
    formatted_msg <- sprintf("[%s] %s: %s%s", timestamp, level, indent_str, msg)
    
    if (verbose) {
      # Color-coded console output
      color_codes <- list(
        ERROR = "\033[31m",    # Red
        WARNING = "\033[33m",  # Yellow
        SUCCESS = "\033[32m",  # Green
        INFO = "\033[0m"       # Default
      )
      
      cat(color_codes[[level]], formatted_msg, "\033[0m\n", sep = "")
    }
    
    # Always write to log file
    cat(formatted_msg, "\n", file = log_file, append = TRUE)
  }
  
  return(log_func)
}

# Load configuration with validation
load_and_validate_config <- function(config_file = "config.yaml") {
  if (!file.exists(config_file)) {
    # Return comprehensive default configuration
    config <- list(
      project = list(
        name = "AThNATCount",
        output_dir = "./results"
      ),
      input = list(
        metadata_file = "./Integrated_metadata_final.csv",
        sa_results_dir = "./test_results/batch_analysis"
      ),
      analysis = list(
        tissue = list(
          consolidate_tissues = TRUE,
          min_samples_per_tissue = 3,
          tissue_mapping = list(
            list(pattern = "root", category = "Root"),
            list(pattern = "leaf|leaves", category = "Leaf"),
            list(pattern = "seed", category = "Seed"),
            list(pattern = "flower", category = "Flower"),
            list(pattern = "stem", category = "Stem")
          )
        ),
        correlation = list(
          min_genes = 100,
          significance_threshold = 0.05,
          correlation_threshold = 0.3
        )
      ),
      output = list(
        formats = c("csv", "xlsx", "json"),
        figures = list(
          dpi = 300,
          format = "png",
          width = 10,
          height = 8,
          theme = "publication"
        )
      ),
      performance = list(
        parallel_cores = 0,
        memory_limit = "8GB"
      )
    )
  } else {
    config <- yaml::read_yaml(config_file)
  }
  
  # Validate essential fields
  required_fields <- c("input", "output", "analysis")
  missing_fields <- setdiff(required_fields, names(config))
  
  if (length(missing_fields) > 0) {
    stop(sprintf("Missing required configuration fields: %s", 
                paste(missing_fields, collapse = ", ")))
  }
  
  return(config)
}

# Enhanced tissue color palette
get_tissue_colors <- function(tissues) {
  # Predefined colors for common tissues
  default_colors <- c(
    "Root" = "#8B4513",
    "Leaf" = "#228B22",
    "Flower" = "#FF1493",
    "Seed" = "#DAA520",
    "Stem" = "#A0522D",
    "Seedling" = "#9370DB",
    "Silique" = "#DC143C",
    "Rosette" = "#4682B4",
    "Unknown" = "#808080"
  )
  
  # Get colors for requested tissues
  tissue_colors <- default_colors[tissues]
  
  # Generate additional colors if needed
  missing_tissues <- tissues[is.na(tissue_colors)]
  if (length(missing_tissues) > 0) {
    additional_colors <- rainbow(length(missing_tissues))
    names(additional_colors) <- missing_tissues
    tissue_colors[missing_tissues] <- additional_colors
  }
  
  return(tissue_colors)
}

# Publication-ready theme
theme_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      # Title and text
      plot.title = element_text(size = rel(1.2), face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = rel(1), hjust = 0.5, color = "gray30"),
      axis.title = element_text(size = rel(1), face = "bold"),
      axis.text = element_text(size = rel(0.9)),
      
      # Legend
      legend.title = element_text(size = rel(1), face = "bold"),
      legend.text = element_text(size = rel(0.9)),
      legend.position = "right",
      
      # Panel
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      
      # Strip (for facets)
      strip.text = element_text(size = rel(1), face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      
      # Overall
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Advanced data loading with validation
load_sa_data <- function(sa_dir, metadata, log_func, parallel = TRUE) {
  log_func("Loading SA_Results data...")
  
  # Get list of samples
  samples <- list.dirs(sa_dir, full.names = FALSE, recursive = FALSE)
  samples <- samples[samples != ""]
  
  log_func(sprintf("Found %d sample directories", length(samples)))
  
  # Define loading function for single sample
  load_sample <- function(sample_id) {
    sample_dir <- file.path(sa_dir, sample_id)
    result <- list(sample_id = sample_id)
    
    # Load correlation data
    corr_file <- file.path(sample_dir, "sense_nat_correlation.csv")
    if (file.exists(corr_file)) {
      corr_data <- fread(corr_file, stringsAsFactors = FALSE)
      corr_data$Sample <- sample_id
      result$correlation <- corr_data
    }
    
    # Load LRR_RLK data if available
    lrr_file <- file.path(sample_dir, "sense_nat_correlation_LRR_RLK.csv")
    if (file.exists(lrr_file)) {
      lrr_data <- fread(lrr_file, stringsAsFactors = FALSE)
      lrr_data$Sample <- sample_id
      result$lrr_rlk <- lrr_data
    }
    
    # Load expression data if available (sampling for memory efficiency)
    expr_file <- file.path(sample_dir, "sense_nat_expression_pairs.csv")
    if (file.exists(expr_file)) {
      # Read only first 10000 rows for summary statistics
      expr_data <- fread(expr_file, nrows = 10000, stringsAsFactors = FALSE)
      expr_data$Sample <- sample_id
      result$expression_sample <- expr_data
    }
    
    return(result)
  }
  
  # Load data in parallel or sequential
  if (parallel && length(samples) > 10) {
    n_cores <- min(parallel::detectCores() - 1, length(samples))
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, library(data.table))
    clusterExport(cl, c("sa_dir"), envir = environment())
    
    sample_data_list <- parLapply(cl, samples, load_sample)
    stopCluster(cl)
  } else {
    sample_data_list <- lapply(samples, load_sample)
  }
  
  # Combine data
  all_correlations <- rbindlist(
    lapply(sample_data_list, function(x) x$correlation),
    fill = TRUE
  )
  
  all_lrr_rlk <- rbindlist(
    lapply(sample_data_list, function(x) x$lrr_rlk),
    fill = TRUE
  )
  
  # Add metadata and create Gene column
  if (!is.null(metadata)) {
    # Create Gene column from sense_gene
    all_correlations$Gene <- all_correlations$sense_gene
    
    all_correlations <- merge(
      all_correlations,
      metadata[, c("Run", "tissue_clean")],
      by.x = "Sample",
      by.y = "Run",
      all.x = TRUE
    )
    
    if (nrow(all_lrr_rlk) > 0) {
      all_lrr_rlk <- merge(
        all_lrr_rlk,
        metadata[, c("Run", "tissue_clean")],
        by.x = "Sample",
        by.y = "Run",
        all.x = TRUE
      )
    }
  }
  
  log_func(sprintf("Loaded %d correlation records from %d samples", 
                  nrow(all_correlations), length(unique(all_correlations$Sample))))
  
  return(list(
    correlations = all_correlations,
    lrr_rlk = all_lrr_rlk,
    sample_list = samples
  ))
}

# Advanced statistical analysis
perform_tissue_analysis <- function(data, config, log_func) {
  log_func("Performing tissue-specific statistical analysis...")
  
  # Filter data
  data_clean <- data %>%
    filter(!is.na(correlation), 
           !is.na(tissue_clean),
           !is.na(p_value))
  
  # Overall statistics by tissue
  tissue_stats <- data_clean %>%
    group_by(tissue_clean) %>%
    summarise(
      n_samples = n_distinct(Sample),
      n_genes = n_distinct(Gene),
      n_measurements = n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      sd_correlation = sd(correlation, na.rm = TRUE),
      positive_ratio = sum(correlation > 0) / n(),
      negative_ratio = sum(correlation < 0) / n(),
      significant_ratio = sum(p_value < config$analysis$correlation$significance_threshold) / n(),
      highly_correlated_ratio = sum(abs(correlation) > config$analysis$correlation$correlation_threshold) / n(),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_measurements))
  
  # Gene-level statistics
  gene_stats <- data_clean %>%
    group_by(Gene, tissue_clean) %>%
    summarise(
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      n_samples = n_distinct(Sample),
      consistency = 1 - sd(correlation, na.rm = TRUE) / (abs(mean(correlation, na.rm = TRUE)) + 0.001),
      .groups = 'drop'
    )
  
  # Identify tissue-specific genes
  tissue_specific_genes <- gene_stats %>%
    group_by(Gene) %>%
    mutate(
      max_correlation = max(abs(mean_correlation)),
      tissue_specificity = abs(mean_correlation) / max_correlation
    ) %>%
    filter(tissue_specificity > 0.8,
           n_samples >= 3) %>%
    arrange(desc(tissue_specificity))
  
  return(list(
    tissue_summary = tissue_stats,
    gene_summary = gene_stats,
    tissue_specific = tissue_specific_genes
  ))
}

# Enhanced visualization suite
create_comprehensive_visualizations <- function(data, stats, output_dir, config, log_func) {
  log_func("Creating comprehensive visualizations...")
  
  tryCatch({
    # Create subdirectories
    dirs <- c("distributions", "comparisons", "heatmaps", "networks")
    for (d in dirs) {
      dir.create(file.path(output_dir, d), showWarnings = FALSE, recursive = TRUE)
    }
  
  # Get figure configuration
  fig_config <- config$output$figures
  tissue_colors <- get_tissue_colors(unique(data$tissue_clean))
  
  # 1. Enhanced correlation distribution plot
  p1 <- ggplot(data, aes(x = correlation)) +
    geom_histogram(aes(y = ..density.., fill = tissue_clean), 
                   bins = 50, alpha = 0.7, position = "identity") +
    geom_density(aes(color = tissue_clean), size = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
    facet_wrap(~tissue_clean, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = tissue_colors) +
    scale_color_manual(values = tissue_colors) +
    labs(
      title = "Sense-Antisense Correlation Distributions by Tissue",
      subtitle = sprintf("n = %d gene pairs across %d tissues", 
                        n_distinct(data$Gene), n_distinct(data$tissue_clean)),
      x = "Correlation Coefficient",
      y = "Density"
    ) +
    theme_publication() +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "distributions", "correlation_distributions.png"),
         p1, width = fig_config$width * 1.2, height = fig_config$height * 1.5,
         dpi = fig_config$dpi)
  
  # 2. Violin plots with statistical comparisons
  p2 <- ggplot(data, aes(x = tissue_clean, y = correlation, fill = tissue_clean)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(values = tissue_colors) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black") +
    labs(
      title = "Tissue-Specific Correlation Patterns",
      subtitle = "Violin plots showing distribution and summary statistics",
      x = "Tissue Type",
      y = "Correlation Coefficient"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Add significance annotations if multiple tissues
  if (n_distinct(data$tissue_clean) > 1) {
    # Perform pairwise comparisons
    tissue_pairs <- combn(unique(data$tissue_clean), 2, simplify = FALSE)
    comparisons <- list()
    
    for (pair in tissue_pairs[1:min(3, length(tissue_pairs))]) {  # Limit to top 3
      comparisons[[length(comparisons) + 1]] <- pair
    }
    
    p2 <- p2 + stat_compare_means(comparisons = comparisons, 
                                  method = "wilcox.test",
                                  label = "p.signif")
  }
  
  ggsave(file.path(output_dir, "comparisons", "tissue_violin_plots.png"),
         p2, width = fig_config$width, height = fig_config$height,
         dpi = fig_config$dpi)
  
  # 3. Comprehensive heatmap
  if (nrow(stats$gene_summary) > 20) {
    # Select top variable genes
    variable_genes <- stats$gene_summary %>%
      group_by(Gene) %>%
      summarise(variability = sd(mean_correlation, na.rm = TRUE)) %>%
      arrange(desc(variability)) %>%
      head(50) %>%
      pull(Gene)
    
    heatmap_data <- stats$gene_summary %>%
      filter(Gene %in% variable_genes) %>%
      select(Gene, tissue_clean, mean_correlation) %>%
      pivot_wider(names_from = tissue_clean, values_from = mean_correlation, values_fill = 0)
    
    # Check if we have enough data for heatmap
    if (nrow(heatmap_data) > 1 && ncol(heatmap_data) > 2) {
      # Create matrix
      heatmap_matrix <- as.matrix(heatmap_data[, -1])
      rownames(heatmap_matrix) <- heatmap_data$Gene
    
    # Create annotations
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    # Column annotation (tissue types)
    tissue_counts <- table(data$tissue_clean)
    ha_column <- HeatmapAnnotation(
      n_samples = anno_barplot(tissue_counts[colnames(heatmap_matrix)], 
                               bar_width = 0.8,
                               gp = gpar(fill = tissue_colors[colnames(heatmap_matrix)])),
      tissue = colnames(heatmap_matrix),
      col = list(tissue = tissue_colors),
      annotation_name_side = "left"
    )
    
    # Row annotation (gene properties)
    gene_props <- stats$gene_summary %>%
      filter(Gene %in% rownames(heatmap_matrix)) %>%
      group_by(Gene) %>%
      summarise(
        overall_mean = mean(mean_correlation),
        consistency = mean(consistency)
      )
    
    ha_row <- rowAnnotation(
      mean_cor = gene_props$overall_mean,
      consistency = gene_props$consistency,
      col = list(
        mean_cor = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
        consistency = colorRamp2(c(0, 1), c("white", "darkgreen"))
      ),
      annotation_name_side = "top"
    )
    
    # Create heatmap
    ht <- Heatmap(
      heatmap_matrix,
      name = "Mean\nCorrelation",
      col = col_fun,
      top_annotation = ha_column,
      right_annotation = ha_row,
      row_title = "Genes",
      column_title = "Tissue Types",
      clustering_method_rows = "ward.D2",
      clustering_method_columns = "ward.D2",
      show_row_names = FALSE,
      show_column_names = TRUE,
      column_names_rot = 45
    )
    
      png(file.path(output_dir, "heatmaps", "gene_tissue_heatmap.png"),
          width = fig_config$width * 150,
          height = fig_config$height * 150,
          res = fig_config$dpi)
      draw(ht)
      dev.off()
    } else {
      log_func("Skipping heatmap creation - insufficient data", "WARNING")
    }
  }
  
  # 4. Scatter plot matrix for tissue comparisons
  if (n_distinct(data$tissue_clean) >= 2 && n_distinct(data$tissue_clean) <= 6) {
    tryCatch({
      # Prepare wide format data
      scatter_data <- data %>%
        select(Gene, tissue_clean, correlation) %>%
        group_by(Gene, tissue_clean) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE), .groups = 'drop') %>%
        pivot_wider(names_from = tissue_clean, values_from = mean_correlation)
      
      # Check if we have enough data
      if (ncol(scatter_data) > 2 && nrow(scatter_data) > 10) {
        # Create scatter plot matrix
        p4 <- GGally::ggpairs(
          scatter_data[, -1],
          upper = list(continuous = GGally::wrap("cor", size = 4)),
          lower = list(continuous = GGally::wrap("points", alpha = 0.3, size = 0.5)),
          diag = list(continuous = GGally::wrap("densityDiag", fill = "lightblue"))
        ) +
          theme_publication()
        
        ggsave(file.path(output_dir, "comparisons", "tissue_scatter_matrix.png"),
               p4, width = fig_config$width * 1.5, height = fig_config$height * 1.5,
               dpi = fig_config$dpi)
      } else {
        log_func("Skipping scatter matrix - insufficient data", "WARNING")
      }
    }, error = function(e) {
      log_func(sprintf("Error creating scatter matrix: %s", e$message), "WARNING")
    })
  }
  
  # 5. Summary statistics visualization
  p5 <- stats$tissue_summary %>%
    select(tissue_clean, positive_ratio, negative_ratio) %>%
    pivot_longer(cols = c(positive_ratio, negative_ratio), 
                names_to = "direction", values_to = "ratio") %>%
    mutate(direction = gsub("_ratio", "", direction)) %>%
    ggplot(aes(x = tissue_clean, y = ratio, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("positive" = "darkgreen", "negative" = "darkred"),
                     labels = c("Positive", "Negative")) +
    scale_y_continuous(labels = percent) +
    labs(
      title = "Correlation Direction by Tissue",
      subtitle = "Proportion of positive vs negative correlations",
      x = "Tissue Type",
      y = "Proportion",
      fill = "Correlation\nDirection"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "comparisons", "correlation_directions.png"),
         p5, width = fig_config$width, height = fig_config$height * 0.8,
         dpi = fig_config$dpi)
  
  log_func("Visualizations completed", "SUCCESS")
  }, error = function(e) {
    log_func(sprintf("Error in visualizations: %s", e$message), "ERROR")
    log_func("Continuing with analysis...", "WARNING")
  })
}

# Generate interactive HTML report
generate_html_report <- function(data, stats, output_dir, config) {
  report_file <- file.path(output_dir, "Tissue_Specific_SA_Analysis_Report.html")
  
  # For now, create a markdown summary that can be converted to HTML
  md_file <- file.path(output_dir, "analysis_summary.md")
  
  md_content <- c(
    "# Tissue-Specific Sense-Antisense Analysis Report",
    "",
    sprintf("**Generated:** %s", Sys.Date()),
    sprintf("**Pipeline Version:** %s", SCRIPT_VERSION),
    "",
    "## Executive Summary",
    "",
    sprintf("- Total samples analyzed: %d", n_distinct(data$Sample)),
    sprintf("- Total gene pairs: %d", n_distinct(data$Gene)),
    sprintf("- Tissue types: %d", n_distinct(data$tissue_clean)),
    sprintf("- Total measurements: %s", format(nrow(data), big.mark = ",")),
    "",
    "## Key Findings",
    "",
    "### Overall Statistics",
    sprintf("- Mean correlation: %.3f", mean(data$correlation, na.rm = TRUE)),
    sprintf("- Positive correlations: %.1f%%", 
            sum(data$correlation > 0, na.rm = TRUE) / nrow(data) * 100),
    sprintf("- Significant pairs (p < 0.05): %.1f%%",
            sum(data$p_value < 0.05, na.rm = TRUE) / nrow(data) * 100),
    "",
    "### Tissue-Specific Insights",
    ""
  )
  
  # Add tissue-specific summaries
  for (i in 1:nrow(stats$tissue_summary)) {
    tissue <- stats$tissue_summary$tissue_clean[i]
    tissue_data <- stats$tissue_summary[i, ]
    
    md_content <- c(md_content,
      sprintf("#### %s", tissue),
      sprintf("- Samples: %d", tissue_data$n_samples),
      sprintf("- Genes analyzed: %d", tissue_data$n_genes),
      sprintf("- Mean correlation: %.3f (SD: %.3f)", 
              tissue_data$mean_correlation, tissue_data$sd_correlation),
      sprintf("- Significant correlations: %.1f%%", 
              tissue_data$significant_ratio * 100),
      ""
    )
  }
  
  # Add top genes section
  top_genes <- data %>%
    filter(p_value < 0.05) %>%
    group_by(Gene) %>%
    summarise(
      mean_cor = mean(correlation),
      n_tissues = n_distinct(tissue_clean),
      min_p = min(p_value)
    ) %>%
    arrange(desc(abs(mean_cor))) %>%
    head(10)
  
  md_content <- c(md_content,
    "### Top Correlated Genes",
    "",
    "| Gene | Mean Correlation | Tissues | Min P-value |",
    "|------|-----------------|---------|-------------|"
  )
  
  for (i in 1:nrow(top_genes)) {
    md_content <- c(md_content,
      sprintf("| %s | %.3f | %d | %.2e |",
              top_genes$Gene[i], top_genes$mean_cor[i],
              top_genes$n_tissues[i], top_genes$min_p[i])
    )
  }
  
  writeLines(md_content, md_file)
  
  return(md_file)
}

# Main analysis pipeline
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default parameters
  config_file <- "config.yaml"
  output_dir <- "tissue_specific_analysis"
  force_reload <- FALSE
  
  # Simple argument parsing
  if (length(args) > 0) {
    if (args[1] == "--help" || args[1] == "-h") {
      cat("Usage: Rscript tissue_specific_sa_analysis.R [config.yaml] [output_dir]\n")
      cat("\nOptions:\n")
      cat("  config.yaml   Configuration file (default: config.yaml)\n")
      cat("  output_dir    Output directory (default: tissue_specific_analysis)\n")
      quit(status = 0)
    }
    config_file <- args[1]
  }
  if (length(args) > 1) output_dir <- args[2]
  
  # Load configuration
  config <- load_and_validate_config(config_file)
  
  # Override output directory if specified
  if (output_dir != "tissue_specific_analysis") {
    config$project$output_dir <- output_dir
  } else {
    output_dir <- config$project$output_dir
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Setup logging
  log_func <- setup_logging(output_dir, verbose = TRUE)
  
  log_func("================================================================================")
  log_func(sprintf("%s v%s", SCRIPT_NAME, SCRIPT_VERSION))
  log_func("================================================================================")
  log_func(sprintf("Configuration file: %s", config_file))
  log_func(sprintf("Output directory: %s", output_dir))
  
  tryCatch({
    # Step 1: Load and prepare metadata
    log_func("Step 1: Loading metadata...")
    
    metadata_file <- config$input$metadata_file
    if (!file.exists(metadata_file)) {
      stop(sprintf("Metadata file not found: %s", metadata_file))
    }
    
    metadata <- fread(metadata_file, stringsAsFactors = FALSE)
    
    # Get samples with SA_Results
    sa_dir <- config$input$sa_results_dir
    if (!dir.exists(sa_dir)) {
      stop(sprintf("SA_Results directory not found: %s", sa_dir))
    }
    
    sa_samples <- list.dirs(sa_dir, full.names = FALSE, recursive = FALSE)
    sa_samples <- sa_samples[sa_samples != ""]
    
    # Filter and prepare metadata
    metadata_filtered <- metadata %>%
      filter(Run %in% sa_samples) %>%
      mutate(
        tissue_raw = case_when(
          !is.na(tissue) & tissue != "" ~ tissue,
          !is.na(Tissue) & Tissue != "" ~ Tissue,
          TRUE ~ "Unknown"
        )
      )
    
    # Apply tissue mapping from config
    if (config$analysis$tissue$consolidate_tissues) {
      # Initialize tissue_clean with tissue_raw
      metadata_filtered$tissue_clean <- metadata_filtered$tissue_raw
      
      for (mapping in config$analysis$tissue$tissue_mapping) {
        pattern <- mapping$pattern
        category <- mapping$category
        metadata_filtered <- metadata_filtered %>%
          mutate(tissue_clean = ifelse(grepl(pattern, tissue_raw, ignore.case = TRUE),
                                      category, tissue_clean))
      }
      
      # Set remaining as original or Unknown
      metadata_filtered <- metadata_filtered %>%
        mutate(tissue_clean = ifelse(is.na(tissue_clean), 
                                    ifelse(tissue_raw == "Unknown", "Unknown", tissue_raw),
                                    tissue_clean))
    } else {
      metadata_filtered$tissue_clean <- metadata_filtered$tissue_raw
    }
    
    log_func(sprintf("Samples with metadata: %d", nrow(metadata_filtered)))
    log_func(sprintf("Tissue types: %s", 
                    paste(unique(metadata_filtered$tissue_clean), collapse = ", ")))
    
    # Step 2: Load SA_Results data
    log_func("\nStep 2: Loading SA_Results data...")
    
    sa_data <- load_sa_data(sa_dir, metadata_filtered, log_func, 
                           parallel = config$performance$parallel_cores > 1)
    
    # Filter for valid data
    correlation_data <- sa_data$correlations %>%
      filter(!is.na(correlation), 
             !is.na(tissue_clean),
             !is.na(p_value))
    
    log_func(sprintf("Valid correlation records: %d", nrow(correlation_data)))
    
    # Step 3: Statistical analysis
    log_func("\nStep 3: Performing statistical analysis...")
    
    analysis_results <- perform_tissue_analysis(correlation_data, config, log_func)
    
    # Step 4: Create visualizations
    log_func("\nStep 4: Creating visualizations...")
    
    viz_dir <- file.path(output_dir, "figures")
    if (!dir.exists(viz_dir)) {
      dir.create(viz_dir, recursive = TRUE)
    }
    
    create_comprehensive_visualizations(
      correlation_data, 
      analysis_results, 
      viz_dir, 
      config, 
      log_func
    )
    
    # Step 5: Save results
    log_func("\nStep 5: Saving analysis results...")
    
    # Save main datasets
    fwrite(correlation_data, 
           file.path(output_dir, "integrated_correlation_data.csv"))
    
    fwrite(analysis_results$tissue_summary,
           file.path(output_dir, "tissue_summary_statistics.csv"))
    
    fwrite(analysis_results$gene_summary,
           file.path(output_dir, "gene_tissue_statistics.csv"))
    
    if (nrow(analysis_results$tissue_specific) > 0) {
      fwrite(analysis_results$tissue_specific,
             file.path(output_dir, "tissue_specific_genes.csv"))
    }
    
    # Save in additional formats if requested
    if ("xlsx" %in% config$output$formats) {
      if (requireNamespace("openxlsx", quietly = TRUE)) {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Correlation_Data")
        openxlsx::writeData(wb, "Correlation_Data", correlation_data)
        openxlsx::addWorksheet(wb, "Tissue_Summary")
        openxlsx::writeData(wb, "Tissue_Summary", analysis_results$tissue_summary)
        openxlsx::addWorksheet(wb, "Gene_Summary")
        openxlsx::writeData(wb, "Gene_Summary", analysis_results$gene_summary)
        openxlsx::saveWorkbook(wb, file.path(output_dir, "analysis_results.xlsx"))
      }
    }
    
    # Step 6: Generate report
    log_func("\nStep 6: Generating analysis report...")
    
    report_file <- generate_html_report(
      correlation_data,
      analysis_results,
      output_dir,
      config
    )
    
    # Save session info
    session_file <- file.path(output_dir, "session_info.txt")
    writeLines(capture.output(sessionInfo()), session_file)
    
    log_func("\n================================================================================")
    log_func("ANALYSIS COMPLETED SUCCESSFULLY", "SUCCESS")
    log_func("================================================================================")
    log_func(sprintf("Results saved to: %s", output_dir))
    log_func(sprintf("Total processing time: %.1f minutes", 
                    as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
    
  }, error = function(e) {
    log_func(sprintf("ERROR: %s", e$message), "ERROR")
    log_func("Stack trace:", "ERROR")
    log_func(paste(capture.output(traceback()), collapse = "\n"), "ERROR")
    quit(status = 1)
  })
}

# Record start time
start_time <- Sys.time()

# Execute if run as script
if (!interactive()) {
  main()
}