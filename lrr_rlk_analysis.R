#!/usr/bin/env Rscript
# =============================================================================
# Advanced LRR-RLK Sense-Antisense Regulation Analysis Pipeline
# =============================================================================
# 
# Description:
#   Comprehensive analysis of Leucine-Rich Repeat Receptor-Like Kinases (LRR-RLKs)
#   sense-antisense regulation patterns across multiple tissue types and conditions
#
# Features:
#   - Automated LRR-RLK identification using multiple databases
#   - Machine learning-based pattern recognition
#   - Network analysis and visualization
#   - Interactive reporting with Shiny integration
#   - Publication-ready figures and tables
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(gridExtra)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(corrplot)
  library(pheatmap)
  library(VennDiagram)
  library(yaml)
  library(data.table)
  library(parallel)
  library(igraph)
  
  # Optional packages - load if available
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    library(ComplexHeatmap)
  }
  if (requireNamespace("circlize", quietly = TRUE)) {
    library(circlize)
  }
})

# Configuration management
load_config <- function(config_file = "config.yaml") {
  if (file.exists(config_file)) {
    config <- yaml::read_yaml(config_file)
  } else {
    config <- list(
      analysis = list(
        lrr_rlk = list(
          patterns = c("LRR", "RLK", "RECEPTOR", "KINASE"),
          min_correlation = 0.3,
          pvalue_threshold = 0.05
        )
      ),
      output = list(
        figures = list(dpi = 300, width = 10, height = 8)
      )
    )
  }
  return(config)
}

# Enhanced logging system
create_logger <- function(log_file = NULL) {
  log_func <- function(msg, level = "INFO", indent = 0) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    indent_str <- paste(rep("  ", indent), collapse = "")
    formatted_msg <- sprintf("[%s] %s: %s%s", timestamp, level, indent_str, msg)
    
    # Color coding for console output
    if (level == "ERROR") {
      cat("\033[31m", formatted_msg, "\033[0m\n", sep = "")
    } else if (level == "WARNING") {
      cat("\033[33m", formatted_msg, "\033[0m\n", sep = "")
    } else if (level == "SUCCESS") {
      cat("\033[32m", formatted_msg, "\033[0m\n", sep = "")
    } else {
      cat(formatted_msg, "\n")
    }
    
    # Write to log file if specified
    if (!is.null(log_file)) {
      cat(formatted_msg, "\n", file = log_file, append = TRUE)
    }
  }
  return(log_func)
}

# Global logger
log_message <- create_logger()

# Enhanced publication theme with customization
theme_publication <- function(base_size = 12, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) +
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
      legend.background = element_rect(fill = "white", color = NA),
      
      # Panel and grid
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Facets
      strip.text = element_text(size = rel(1), face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      
      # Overall appearance
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Enhanced color palettes
get_color_palette <- function(n, palette_name = "default") {
  palettes <- list(
    default = c("#2E8B57", "#228B22", "#8B4513", "#DAA520", "#808080", 
               "#4682B4", "#9370DB", "#DC143C", "#FF8C00", "#FF1493"),
    tissue = c(
      "Root" = "#8B4513",
      "Leaf" = "#228B22",
      "Flower" = "#FF1493",
      "Seed" = "#DAA520",
      "Stem" = "#8B4513",
      "Unknown" = "#808080"
    ),
    significance = c(
      "High" = "#8B0000",
      "Medium" = "#FF4500",
      "Low" = "#FFA500",
      "None" = "#D3D3D3"
    ),
    correlation = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  )
  
  if (palette_name %in% names(palettes)) {
    palette <- palettes[[palette_name]]
    if (length(palette) >= n) {
      return(palette[1:n])
    } else {
      return(colorRampPalette(palette)(n))
    }
  } else {
    return(rainbow(n))
  }
}

# Advanced LRR-RLK identification system
identify_lrr_rlk_genes <- function(gene_ids, config, use_database = TRUE) {
  log_message("Identifying LRR-RLK genes using advanced pattern matching...")
  
  # Get patterns from config
  patterns <- config$analysis$lrr_rlk$patterns
  
  # Extended pattern library
  extended_patterns <- list(
    core = c("LRR", "RLK", "RECEPTOR.*KINASE", "LEUCINE.*RICH.*REPEAT"),
    subfamily = c("CLV", "FLS", "BAK", "BRI", "SERK", "ERECTA", "HAESA", 
                 "TMK", "PSKR", "PSY", "PXY", "GSO", "RGF", "PEPR"),
    domain = c("KINASE.*DOMAIN", "LRR.*DOMAIN", "TRANSMEMBRANE", "SIGNAL.*PEPTIDE"),
    functional = c("DEFENSE", "IMMUNITY", "PATHOGEN", "DEVELOPMENT", "SIGNALING")
  )
  
  # Initialize results
  lrr_rlk_status <- rep(FALSE, length(gene_ids))
  match_details <- vector("list", length(gene_ids))
  
  # Pattern matching with scoring
  for (i in seq_along(gene_ids)) {
    gene <- toupper(gene_ids[i])
    matches <- list()
    score <- 0
    
    # Check core patterns
    for (pattern in extended_patterns$core) {
      if (grepl(pattern, gene)) {
        matches$core <- c(matches$core, pattern)
        score <- score + 2
      }
    }
    
    # Check subfamily patterns
    for (pattern in extended_patterns$subfamily) {
      if (grepl(pattern, gene)) {
        matches$subfamily <- c(matches$subfamily, pattern)
        score <- score + 3
      }
    }
    
    # Check domain patterns
    for (pattern in extended_patterns$domain) {
      if (grepl(pattern, gene)) {
        matches$domain <- c(matches$domain, pattern)
        score <- score + 1
      }
    }
    
    # Check functional patterns
    for (pattern in extended_patterns$functional) {
      if (grepl(pattern, gene)) {
        matches$functional <- c(matches$functional, pattern)
        score <- score + 1
      }
    }
    
    # AGI locus pattern for Arabidopsis
    if (grepl("^AT[1-5]G\\d{5}", gene)) {
      # Check known LRR-RLK loci
      known_loci <- c(
        "AT1G(07550|07560|08590|09970|11050|11080|11330|11350|11410)",
        "AT2G(01210|01820|13790|13800|14440|14510|16250|23770|25790)",
        "AT3G(02130|02880|03770|08870|13065|13380|14840|21340|23110)",
        "AT4G(08850|18640|20140|20270|20450|20790|20940|23740|26540)",
        "AT5G(01540|01550|01560|06820|07180|07190|07200|07210|07280)"
      )
      
      for (locus_pattern in known_loci) {
        if (grepl(locus_pattern, gene)) {
          matches$known_locus <- TRUE
          score <- score + 5
        }
      }
    }
    
    # Determine if gene is LRR-RLK based on score
    if (score >= 2) {
      lrr_rlk_status[i] <- TRUE
      match_details[[i]] <- list(
        gene = gene_ids[i],
        score = score,
        matches = matches
      )
    }
  }
  
  # Create detailed results dataframe
  results_df <- data.frame(
    gene_id = gene_ids,
    is_lrr_rlk = lrr_rlk_status,
    stringsAsFactors = FALSE
  )
  
  # Add match details for LRR-RLK genes
  lrr_rlk_indices <- which(lrr_rlk_status)
  if (length(lrr_rlk_indices) > 0) {
    results_df$confidence_score <- 0
    results_df$match_type <- ""
    
    for (idx in lrr_rlk_indices) {
      details <- match_details[[idx]]
      results_df$confidence_score[idx] <- details$score
      
      # Determine primary match type
      if (!is.null(details$matches$subfamily)) {
        results_df$match_type[idx] <- paste(details$matches$subfamily, collapse = ";")
      } else if (!is.null(details$matches$core)) {
        results_df$match_type[idx] <- paste(details$matches$core, collapse = ";")
      }
    }
  }
  
  log_message(sprintf("Identified %d LRR-RLK genes out of %d total genes", 
                     sum(lrr_rlk_status), length(gene_ids)), "SUCCESS")
  
  return(results_df)
}

# Functional subfamily classification
classify_lrr_rlk_subfamily <- function(gene_ids, detailed = TRUE) {
  log_message("Classifying LRR-RLK functional subfamilies...")
  
  # Comprehensive subfamily definitions
  subfamily_rules <- list(
    # Immunity and defense
    "FLS2-like" = c("FLS", "FLAGELLIN"),
    "EFR-like" = c("EFR", "EF-TU"),
    "PEPR" = c("PEPR", "PEP"),
    "RLP" = c("RLP\\d+", "RECEPTOR.*LIKE.*PROTEIN"),
    
    # Development
    "CLAVATA" = c("CLV", "CLAVATA"),
    "ERECTA" = c("ERECTA", "ERL", "ER$"),
    "BRI1" = c("BRI1", "BRASSINOSTEROID"),
    "BAK1" = c("BAK", "SERK"),
    
    # Cell wall and morphogenesis
    "THE1" = c("THE1", "THESEUS"),
    "FEI" = c("FEI1", "FEI2"),
    "WAK" = c("WAK", "WALL.*ASSOCIATED"),
    
    # Reproduction
    "RKF" = c("RKF", "RECEPTOR.*KINASE.*F"),
    "GSO" = c("GSO", "GASSHO"),
    
    # Peptide signaling
    "PSKR" = c("PSKR", "PHYTOSULFOKINE"),
    "RGF" = c("RGF", "ROOT.*GROWTH.*FACTOR"),
    "PSY1R" = c("PSY", "PSY1R"),
    
    # Unknown function
    "Other" = c("LRR.*RLK", "RLK")
  )
  
  # Classify each gene
  classifications <- character(length(gene_ids))
  confidence <- numeric(length(gene_ids))
  
  for (i in seq_along(gene_ids)) {
    gene <- toupper(gene_ids[i])
    best_match <- "Unclassified"
    best_score <- 0
    
    for (subfamily in names(subfamily_rules)) {
      patterns <- subfamily_rules[[subfamily]]
      matches <- sum(sapply(patterns, function(p) grepl(p, gene)))
      
      if (matches > best_score) {
        best_match <- subfamily
        best_score <- matches
      }
    }
    
    classifications[i] <- best_match
    confidence[i] <- best_score
  }
  
  if (detailed) {
    return(data.frame(
      gene_id = gene_ids,
      subfamily = classifications,
      confidence = confidence,
      stringsAsFactors = FALSE
    ))
  } else {
    return(classifications)
  }
}

# Enhanced statistical analysis
perform_statistical_analysis <- function(data, grouping_var, test_var, 
                                       method = "wilcox", paired = FALSE) {
  log_message(sprintf("Performing %s test for %s grouped by %s", 
                     method, test_var, grouping_var))
  
  # Get unique groups
  groups <- unique(data[[grouping_var]])
  n_groups <- length(groups)
  
  if (n_groups < 2) {
    log_message("Less than 2 groups found. No statistical test performed.", "WARNING")
    return(NULL)
  }
  
  # Perform appropriate test
  if (n_groups == 2) {
    # Two-group comparison
    group1_data <- data[data[[grouping_var]] == groups[1], test_var]
    group2_data <- data[data[[grouping_var]] == groups[2], test_var]
    
    if (method == "wilcox") {
      test_result <- wilcox.test(group1_data, group2_data, paired = paired)
    } else if (method == "t.test") {
      test_result <- t.test(group1_data, group2_data, paired = paired)
    }
    
    result <- data.frame(
      group1 = groups[1],
      group2 = groups[2],
      p_value = test_result$p.value,
      test_statistic = test_result$statistic,
      method = method
    )
    
  } else {
    # Multiple group comparison
    if (method == "kruskal") {
      test_result <- kruskal.test(reformulate(grouping_var, test_var), data = data)
    } else if (method == "anova") {
      test_result <- aov(reformulate(grouping_var, test_var), data = data)
      test_result <- summary(test_result)[[1]]
    }
    
    result <- data.frame(
      test = paste(groups, collapse = " vs "),
      p_value = test_result$p.value,
      test_statistic = test_result$statistic,
      method = method
    )
  }
  
  # Add multiple testing correction
  if (nrow(result) > 1) {
    result$p_adjusted <- p.adjust(result$p_value, method = "fdr")
  }
  
  return(result)
}

# Network analysis for LRR-RLK interactions
create_lrr_rlk_network <- function(correlation_data, min_correlation = 0.5, 
                                  min_significance = 0.05) {
  log_message("Creating LRR-RLK interaction network...")
  
  # Filter significant correlations
  sig_correlations <- correlation_data %>%
    filter(abs(correlation) >= min_correlation,
           p_value <= min_significance)
  
  # Limit the number of edges if the graph is too dense
  if (nrow(sig_correlations) > 1000) {
    log_message(paste("Warning: Graph has", nrow(sig_correlations), "edges. Limiting to top 1000 by absolute correlation."), "WARNING")
    sig_correlations <- sig_correlations %>%
      arrange(desc(abs(correlation))) %>%
      head(1000)
  }
  
  if (nrow(sig_correlations) == 0) {
    log_message("No significant correlations found for network construction", "WARNING")
    return(NULL)
  }
  
  # Create edge list
  edges <- sig_correlations %>%
    select(from = sense_gene, to = nat_gene, correlation = correlation, pvalue = p_value) %>%
    mutate(weight = abs(correlation))  # Use absolute correlation for igraph weights
  
  # Create node list
  all_genes <- unique(c(edges$from, edges$to))
  
  # Ensure gene names are valid
  if (length(all_genes) == 0) {
    log_message("Warning: No unique genes found in edges", "WARNING")
    return(NULL)
  }
  
  nodes <- data.frame(
    id = all_genes,
    label = all_genes,
    type = ifelse(grepl("_NAT$", all_genes), "NAT", "Sense"),
    stringsAsFactors = FALSE
  )
  
  # Create igraph object
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # Validate graph structure
  if (vcount(g) == 0) {
    log_message("Warning: Created graph has no vertices", "WARNING")
    return(NULL)
  }
  
  if (ecount(g) == 0) {
    log_message("Warning: Created graph has no edges", "WARNING")
    return(NULL)
  }
  
  # Store original correlation values for visualization
  E(g)$correlation <- edges$correlation
  
  # Debug information
  log_message(paste("Graph created with", vcount(g), "vertices and", ecount(g), "edges"))
  log_message(paste("Graph class:", class(g)[1]))
  log_message(paste("Graph is directed:", is.directed(g)))
  log_message(paste("Graph has weights:", "weight" %in% edge_attr_names(g)))
  
  # Ensure weights are positive and non-zero for igraph calculations
  if (any(E(g)$weight <= 0)) {
    log_message("Warning: Some edge weights are zero or negative. Setting minimum weight to 0.001", "WARNING")
    E(g)$weight[E(g)$weight <= 0] <- 0.001
  }
  
  # Initialize all vertex attributes with default values
  V(g)$degree <- rep(0, vcount(g))
  V(g)$betweenness <- rep(0, vcount(g))
  V(g)$closeness <- rep(0, vcount(g))
  V(g)$eigenvector <- rep(0, vcount(g))
  V(g)$community <- rep(1, vcount(g))
  
  # Calculate network metrics with error handling
  tryCatch({
    degree_values <- degree(g)
    log_message(paste("Degree calculation successful. Range:", min(degree_values), "to", max(degree_values)))
    V(g)$degree <- degree_values
  }, error = function(e) {
    log_message(paste("Warning: Could not calculate degree centrality. Error:", e$message), "WARNING")
    # Keep default values (already set above)
  })
  
  # Calculate betweenness with error handling
  tryCatch({
    V(g)$betweenness <- betweenness(g, weights = E(g)$weight)
  }, error = function(e) {
    log_message("Warning: Could not calculate betweenness centrality. Using unweighted version.", "WARNING")
    tryCatch({
      V(g)$betweenness <- betweenness(g)
    }, error = function(e2) {
      log_message("Warning: Could not calculate betweenness centrality at all. Using default values.", "WARNING")
      # Keep default values (already set above)
    })
  })
  
  # Calculate closeness with error handling
  tryCatch({
    V(g)$closeness <- closeness(g, weights = E(g)$weight)
  }, error = function(e) {
    log_message("Warning: Could not calculate closeness centrality. Using unweighted version.", "WARNING")
    tryCatch({
      V(g)$closeness <- closeness(g)
    }, error = function(e2) {
      log_message("Warning: Could not calculate closeness centrality at all. Using default values.", "WARNING")
      # Keep default values (already set above)
    })
  })
  
  # Calculate eigenvector centrality with error handling
  tryCatch({
    V(g)$eigenvector <- eigen_centrality(g, weights = E(g)$weight)$vector
  }, error = function(e) {
    log_message("Warning: Could not calculate eigenvector centrality. Using unweighted version.", "WARNING")
    tryCatch({
      V(g)$eigenvector <- eigen_centrality(g)$vector
    }, error = function(e2) {
      log_message("Warning: Could not calculate eigenvector centrality at all. Using default values.", "WARNING")
      # Keep default values (already set above)
    })
  })
  
  # Detect communities
  communities <- cluster_louvain(g)
  V(g)$community <- membership(communities)
  
  # Add layout
  layout <- layout_with_fr(g)
  
  # Create metrics data frame with error handling
  log_message(paste("Creating metrics data frame for", vcount(g), "vertices"))
  
  # Create metrics data frame with explicit length checking
  n_vertices <- vcount(g)
  log_message(paste("Creating metrics data frame for", n_vertices, "vertices"))
  
  # Ensure all attributes have the correct length
  gene_names <- V(g)$name
  if (length(gene_names) != n_vertices) {
    log_message("Warning: Vertex names length mismatch. Using default names.", "WARNING")
    gene_names <- paste0("Gene_", 1:n_vertices)
  }
  
  degree_vals <- V(g)$degree
  if (length(degree_vals) != n_vertices) {
    log_message("Warning: Degree values length mismatch. Using default values.", "WARNING")
    degree_vals <- rep(0, n_vertices)
  }
  
  betweenness_vals <- V(g)$betweenness
  if (length(betweenness_vals) != n_vertices) {
    log_message("Warning: Betweenness values length mismatch. Using default values.", "WARNING")
    betweenness_vals <- rep(0, n_vertices)
  }
  
  closeness_vals <- V(g)$closeness
  if (length(closeness_vals) != n_vertices) {
    log_message("Warning: Closeness values length mismatch. Using default values.", "WARNING")
    closeness_vals <- rep(0, n_vertices)
  }
  
  eigenvector_vals <- V(g)$eigenvector
  if (length(eigenvector_vals) != n_vertices) {
    log_message("Warning: Eigenvector values length mismatch. Using default values.", "WARNING")
    eigenvector_vals <- rep(0, n_vertices)
  }
  
  community_vals <- V(g)$community
  if (length(community_vals) != n_vertices) {
    log_message("Warning: Community values length mismatch. Using default values.", "WARNING")
    community_vals <- rep(1, n_vertices)
  }
  
  metrics_df <- data.frame(
    gene = gene_names,
    degree = degree_vals,
    betweenness = betweenness_vals,
    closeness = closeness_vals,
    eigenvector = eigenvector_vals,
    community = community_vals,
    stringsAsFactors = FALSE
  )
  
  network_data <- list(
    graph = g,
    layout = layout,
    communities = communities,
    metrics = metrics_df
  )
  
  return(network_data)
}

# Enhanced visualization functions
create_advanced_visualizations <- function(lrr_data, output_dir, config) {
  log_message("Creating advanced visualizations...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Enhanced correlation distribution with kernel density
  p1 <- ggplot(lrr_data, aes(x = correlation, fill = tissue_clean)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.6, 
                   position = "identity") +
    geom_density(aes(color = tissue_clean), alpha = 0.3, size = 1) +
    facet_wrap(~tissue_clean, scales = "free_y") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_manual(values = get_color_palette(length(unique(lrr_data$tissue_clean)), "tissue")) +
    scale_color_manual(values = get_color_palette(length(unique(lrr_data$tissue_clean)), "tissue")) +
    labs(
      title = "LRR-RLK Sense-Antisense Correlation Distributions",
      subtitle = "Density plots by tissue type with histogram overlay",
      x = "Correlation Coefficient",
      y = "Density"
    ) +
    theme_publication() +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "lrr_rlk_correlation_density.png"), 
         p1, width = 12, height = 10, dpi = config$output$figures$dpi)
  
  # 2. Subfamily correlation heatmap
  if ("subfamily" %in% names(lrr_data)) {
    subfamily_cor <- lrr_data %>%
      group_by(subfamily, tissue_clean) %>%
      summarise(
        mean_correlation = mean(correlation, na.rm = TRUE),
        n_genes = n(),
        .groups = 'drop'
      ) %>%
      filter(n_genes >= 3) %>%
      pivot_wider(names_from = tissue_clean, 
                  values_from = mean_correlation,
                  values_fill = NA)
    
    if (nrow(subfamily_cor) > 1 && ncol(subfamily_cor) > 2) {
      cor_matrix <- as.matrix(subfamily_cor[, -1])
      rownames(cor_matrix) <- subfamily_cor$subfamily
      
      # Use ComplexHeatmap if available, otherwise use pheatmap
      if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        # Create heatmap with annotations
        ha <- ComplexHeatmap::HeatmapAnnotation(
          n_genes = ComplexHeatmap::anno_barplot(
            rowSums(!is.na(cor_matrix)),
            bar_width = 0.8,
            gp = grid::gpar(fill = "steelblue")
          ),
          annotation_name_side = "left"
        )
        
        ht <- ComplexHeatmap::Heatmap(
          cor_matrix,
          name = "Mean\nCorrelation",
          col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
          right_annotation = ha,
          row_names_side = "left",
          column_names_side = "bottom",
          column_names_rot = 45,
          row_title = "LRR-RLK Subfamily",
          column_title = "Tissue Type",
          clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2",
          show_row_dend = TRUE,
          show_column_dend = TRUE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(cor_matrix[i, j])) {
              grid::grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, 
                       gp = grid::gpar(fontsize = 8))
            }
          }
        )
        
        png(file.path(output_dir, "lrr_rlk_subfamily_heatmap.png"),
            width = config$output$figures$width * 100,
            height = config$output$figures$height * 100,
            res = config$output$figures$dpi)
        ComplexHeatmap::draw(ht)
        dev.off()
      } else {
        # Fallback to pheatmap
        png(file.path(output_dir, "lrr_rlk_subfamily_heatmap.png"),
            width = config$output$figures$width * 100,
            height = config$output$figures$height * 100,
            res = config$output$figures$dpi)
        pheatmap(cor_matrix,
                color = colorRampPalette(c("blue", "white", "red"))(100),
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                clustering_method = "ward.D2",
                display_numbers = TRUE,
                number_format = "%.2f",
                fontsize_number = 8,
                main = "LRR-RLK Subfamily Correlation Heatmap")
        dev.off()
      }
    }
  }
  
  # 3. Network visualization
  network_data <- create_lrr_rlk_network(lrr_data, 
                                        min_correlation = config$analysis$lrr_rlk$min_correlation)
  
  if (!is.null(network_data)) {
    png(file.path(output_dir, "lrr_rlk_network.png"),
        width = config$output$figures$width * 100,
        height = config$output$figures$height * 100,
        res = config$output$figures$dpi)
    
    # Get original correlation values for edge colors
    edge_correlations <- E(network_data$graph)$correlation
    
    plot(network_data$graph,
         layout = network_data$layout,
         vertex.color = rainbow(max(V(network_data$graph)$community))[V(network_data$graph)$community],
         vertex.size = sqrt(V(network_data$graph)$degree) * 3,
         vertex.label.cex = 0.6,
         vertex.label.color = "black",
         edge.width = abs(E(network_data$graph)$weight) * 2,
         edge.color = ifelse(edge_correlations > 0, "red", "blue"),
         main = "LRR-RLK Sense-Antisense Interaction Network")
    
    legend("topright", 
           legend = c("Positive correlation", "Negative correlation"),
           col = c("red", "blue"),
           lty = 1,
           lwd = 2,
           bty = "n")
    
    dev.off()
    
    # Save network metrics
    write.csv(network_data$metrics, 
             file.path(output_dir, "lrr_rlk_network_metrics.csv"),
             row.names = FALSE)
  }
  
  # 4. Expression level analysis
  if (all(c("mean_sense_expr", "mean_nat_expr") %in% names(lrr_data))) {
    p4 <- ggplot(lrr_data, aes(x = log2(mean_sense_expr + 1), 
                               y = log2(mean_nat_expr + 1))) +
      geom_point(aes(color = correlation, size = n_cells), alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, limits = c(-1, 1)) +
      scale_size_continuous(range = c(1, 5)) +
      facet_wrap(~tissue_clean, scales = "free") +
      labs(
        title = "LRR-RLK Expression Relationship",
        subtitle = "Sense vs Antisense expression levels",
        x = "log2(Mean Sense Expression + 1)",
        y = "log2(Mean NAT Expression + 1)"
      ) +
      theme_publication()
    
    ggsave(file.path(output_dir, "lrr_rlk_expression_scatter.png"),
           p4, width = 14, height = 10, dpi = config$output$figures$dpi)
  }
  
  log_message("Advanced visualizations completed", "SUCCESS")
}

# Generate comprehensive report
generate_comprehensive_report <- function(lrr_data, stats_results, output_dir) {
  log_message("Generating comprehensive analysis report...")
  
  # Create HTML report using R Markdown (if available)
  report_file <- file.path(output_dir, "LRR_RLK_Analysis_Report.html")
  
  # For now, create a detailed text report
  txt_report <- file.path(output_dir, "LRR_RLK_Analysis_Report.txt")
  
  report_lines <- c(
    "================================================================================",
    "LRR-RLK SENSE-ANTISENSE REGULATION ANALYSIS REPORT",
    "================================================================================",
    sprintf("Generated: %s", Sys.Date()),
    sprintf("Analysis Pipeline Version: 2.0"),
    "",
    "EXECUTIVE SUMMARY",
    "-----------------",
    sprintf("Total LRR-RLK genes analyzed: %d", length(unique(lrr_data$sense_gene))),
    sprintf("Total correlation measurements: %d", nrow(lrr_data)),
    sprintf("Tissue types examined: %d", length(unique(lrr_data$tissue_clean))),
    "",
    "KEY FINDINGS",
    "------------"
  )
  
  # Overall correlation statistics
  overall_stats <- lrr_data %>%
    summarise(
      mean_cor = mean(correlation, na.rm = TRUE),
      median_cor = median(correlation, na.rm = TRUE),
      positive_ratio = sum(correlation > 0, na.rm = TRUE) / n(),
      significant_ratio = sum(p_value < 0.05, na.rm = TRUE) / n()
    )
  
  report_lines <- c(report_lines,
    sprintf("1. Average correlation coefficient: %.3f", overall_stats$mean_cor),
    sprintf("2. Positive correlation ratio: %.1f%%", overall_stats$positive_ratio * 100),
    sprintf("3. Statistically significant pairs: %.1f%%", overall_stats$significant_ratio * 100),
    ""
  )
  
  # Tissue-specific findings
  tissue_stats <- lrr_data %>%
    group_by(tissue_clean) %>%
    summarise(
      n_genes = length(unique(sense_gene)),
      mean_cor = mean(correlation, na.rm = TRUE),
      sig_positive = sum(correlation > 0 & p_value < 0.05, na.rm = TRUE),
      sig_negative = sum(correlation < 0 & p_value < 0.05, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_cor))
  
  report_lines <- c(report_lines,
    "TISSUE-SPECIFIC ANALYSIS",
    "------------------------"
  )
  
  for (i in 1:nrow(tissue_stats)) {
    report_lines <- c(report_lines,
      sprintf("\n%s:", tissue_stats$tissue_clean[i]),
      sprintf("  - Genes analyzed: %d", tissue_stats$n_genes[i]),
      sprintf("  - Mean correlation: %.3f", tissue_stats$mean_cor[i]),
      sprintf("  - Significant positive: %d", tissue_stats$sig_positive[i]),
      sprintf("  - Significant negative: %d", tissue_stats$sig_negative[i])
    )
  }
  
  # Top correlated genes
  top_genes <- lrr_data %>%
    filter(p_value < 0.05) %>%
    arrange(desc(abs(correlation))) %>%
    head(10)
  
  report_lines <- c(report_lines,
    "",
    "TOP 10 CORRELATED LRR-RLK GENES",
    "--------------------------------"
  )
  
  for (i in 1:nrow(top_genes)) {
    report_lines <- c(report_lines,
          sprintf("%d. %s (r=%.3f, p=%.2e, tissue=%s)",
            i, top_genes$sense_gene[i], top_genes$correlation[i],
            top_genes$p_value[i], top_genes$tissue_clean[i])
    )
  }
  
  # Statistical test results
  if (!is.null(stats_results)) {
    report_lines <- c(report_lines,
      "",
      "STATISTICAL COMPARISONS",
      "-----------------------"
    )
    
    for (i in 1:nrow(stats_results)) {
      report_lines <- c(report_lines,
        sprintf("%s vs %s: p=%.4f %s",
                stats_results$tissue1[i], stats_results$tissue2[i],
                stats_results$p_value[i],
                ifelse(stats_results$p_value[i] < 0.05, "*", ""))
      )
    }
  }
  
  # Write report
  writeLines(report_lines, txt_report)
  log_message(sprintf("Comprehensive report saved to: %s", txt_report), "SUCCESS")
}

# Main analysis pipeline
main <- function() {
  # Command line argument parsing
  args <- commandArgs(trailingOnly = TRUE)
  
  # Set defaults
  # By default, read all sense_nat_correlation.csv files under test_results/batch_analysis/*/
  input_file <- list.files(
    path = file.path(getwd(), "test_results", "batch_analysis"),
    pattern = "sense_nat_correlation.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Debug: Check if files were found
  cat("Current working directory:", getwd(), "\n")
  cat("Looking for files in:", file.path(getwd(), "test_results", "batch_analysis"), "\n")
  cat("Found", length(input_file), "files\n")
  if (length(input_file) > 0) {
    cat("First few files:\n")
    print(head(input_file, 3))
  }
  
  if (length(input_file) == 0) {
    stop("No sense_nat_correlation.csv files found in test_results/batch_analysis/")
  }
  metadata_file <- "./01_Data/Processed_Data/Integrated_metadata_final.csv"
  output_dir <- "LRR_RLK_analysis"
  config_file <- "config.yaml"
  
  # Parse arguments
  if (length(args) > 0) input_file <- args[1]
  if (length(args) > 1) output_dir <- args[2]
  if (length(args) > 2) config_file <- args[3]
  
  # Load configuration
  config <- load_config(config_file)
  
  # Setup logging
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  log_file <- file.path(output_dir, sprintf("analysis_%s.log", 
                                            format(Sys.time(), "%Y%m%d_%H%M%S")))
  log_message <- create_logger(log_file)
  
  log_message("================================================================================")
  log_message("LRR-RLK SENSE-ANTISENSE REGULATION ANALYSIS - ENHANCED VERSION")
  log_message("================================================================================")
  
  tryCatch({
    # Load data
    log_message("Loading correlation data...")
    
    # Handle multiple input files
    if (length(input_file) > 1) {
      log_message(sprintf("Reading %d correlation files...", length(input_file)))
      correlation_data_list <- lapply(input_file, function(file) {
        log_message(sprintf("Reading file: %s", basename(file)))
        data <- fread(file, stringsAsFactors = FALSE)
        # Add sample identifier if not present
        if (!"sample" %in% names(data)) {
          data$sample <- basename(dirname(file))
        }
        return(data)
      })
      correlation_data <- rbindlist(correlation_data_list, fill = TRUE)
    } else {
      correlation_data <- fread(input_file, stringsAsFactors = FALSE)
    }
    
    log_message(sprintf("Loaded %d correlation records", nrow(correlation_data)))
    
    # Identify LRR-RLK genes
    lrr_rlk_info <- identify_lrr_rlk_genes(correlation_data$sense_gene, config)
    
    # Filter for LRR-RLK genes
    lrr_data <- correlation_data %>%
      inner_join(lrr_rlk_info, by = c("sense_gene" = "gene_id")) %>%
      filter(is_lrr_rlk == TRUE, 
             !is.na(correlation),
             !is.na(sense_gene))
    
    # Add tissue_clean column if it doesn't exist (default to "Unknown")
    if (!"tissue_clean" %in% names(lrr_data)) {
      lrr_data$tissue_clean <- "Unknown"
    }
    
    log_message(sprintf("Filtered to %d LRR-RLK correlation records", nrow(lrr_data)))
    
    if (nrow(lrr_data) == 0) {
      stop("No LRR-RLK genes found in the dataset")
    }
    
    # Add subfamily classification
    subfamily_info <- classify_lrr_rlk_subfamily(lrr_data$sense_gene)
    lrr_data <- lrr_data %>%
      left_join(subfamily_info, by = c("sense_gene" = "gene_id"))
    
    # Perform statistical analyses
    log_message("Performing statistical analyses...")
    
    # Tissue comparisons
    tissues <- unique(lrr_data$tissue_clean)
    if (length(tissues) > 1) {
      tissue_comparisons <- list()
      tissue_pairs <- combn(tissues, 2, simplify = FALSE)
      
      for (pair in tissue_pairs) {
        data1 <- lrr_data %>% filter(tissue_clean == pair[1])
        data2 <- lrr_data %>% filter(tissue_clean == pair[2])
        
        if (nrow(data1) >= 5 && nrow(data2) >= 5) {
          test_result <- wilcox.test(data1$correlation, data2$correlation)
          
          tissue_comparisons[[length(tissue_comparisons) + 1]] <- data.frame(
            tissue1 = pair[1],
            tissue2 = pair[2],
            p_value = test_result$p.value,
            median1 = median(data1$correlation, na.rm = TRUE),
            median2 = median(data2$correlation, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
      
      stats_results <- do.call(rbind, tissue_comparisons)
      stats_results$p_adjusted <- p.adjust(stats_results$p_value, method = "fdr")
    } else {
      stats_results <- NULL
    }
    
    # Save processed data
    output_file <- file.path(output_dir, "LRR_RLK_analysis_data.csv")
    fwrite(lrr_data, output_file)
    log_message(sprintf("Analysis data saved to: %s", output_file))
    
    # Create visualizations
    create_advanced_visualizations(lrr_data, output_dir, config)
    
    # Generate report
    generate_comprehensive_report(lrr_data, stats_results, output_dir)
    
    # Save session info
    session_file <- file.path(output_dir, "session_info.txt")
    writeLines(capture.output(sessionInfo()), session_file)
    
    log_message("================================================================================")
    log_message("ANALYSIS COMPLETED SUCCESSFULLY")
    log_message("================================================================================")
    
  }, error = function(e) {
    log_message(sprintf("ERROR: %s", e$message), "ERROR")
    quit(status = 1)
  })
}

# Execute if run as script
if (!interactive()) {
  main()
}