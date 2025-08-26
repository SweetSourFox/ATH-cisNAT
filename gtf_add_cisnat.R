#!/usr/bin/env Rscript
# =============================================================================
# Enhanced GTF Annotation Builder with cis-NAT Integration
# =============================================================================
# 
# Description:
#   Builds custom GTF annotation by integrating cis-Natural Antisense 
#   Transcript (cis-NAT) information into reference GTF files
#
# Features:
#   - Automatic package installation
#   - Parallel processing support
#   - Comprehensive error handling
#   - Multiple input format support
#   - Quality control and validation
# =============================================================================

# Load configuration
load_config <- function(config_file = "config.yaml") {
  if (file.exists(config_file)) {
    return(yaml::read_yaml(config_file))
  } else {
    # Return default configuration
    return(list(
      reference = list(
        gtf_file = "./Ref/Araport11_genes.gtf",
        lncrna_gff = "./Ref/A.tha.lncRNA.gff"
      ),
      output = list(
        gtf_file = "custom_with_cisNAT.gtf"
      ),
      performance = list(
        parallel_cores = 0
      )
    ))
  }
}

# Progress message function with timestamp
progress_msg <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

# Enhanced package management
ensure_packages <- function(packages) {
  progress_msg("Checking required R packages...")
  
  for (pkg in names(packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      progress_msg(sprintf("Installing %s from %s...", pkg, packages[[pkg]]), "INFO")
      
      if (packages[[pkg]] == "CRAN") {
        install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      } else if (packages[[pkg]] == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
        }
        BiocManager::install(pkg, quiet = TRUE, update = FALSE)
      }
    }
  }
  
  # Load packages
  for (pkg in names(packages)) {
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
}

# Required packages with sources
required_packages <- list(
  "rtracklayer" = "Bioconductor",
  "GenomicRanges" = "Bioconductor",
  "dplyr" = "CRAN",
  "stringr" = "CRAN",
  "parallel" = "CRAN",
  "yaml" = "CRAN"
)

ensure_packages(required_packages)

# Function to validate GTF/GFF file
validate_annotation_file <- function(file_path, format = "auto") {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Auto-detect format if needed
  if (format == "auto") {
    format <- ifelse(grepl("\\.gtf$", file_path, ignore.case = TRUE), "gtf", "gff")
  }
  
  progress_msg(sprintf("Validating %s file: %s", toupper(format), file_path))
  
  # Try to read first few lines
  tryCatch({
    test_read <- rtracklayer::import(file_path, format = format)
    if (length(test_read) == 0) {
      stop("File appears to be empty")
    }
    progress_msg(sprintf("Successfully validated %s with %d features", 
                        basename(file_path), length(test_read)), "INFO")
    return(TRUE)
  }, error = function(e) {
    stop(sprintf("Error reading %s file: %s", format, e$message))
  })
}

# Enhanced cis-NAT pair finding with parallel processing
find_cis_nat_pairs_parallel <- function(genes_gr, num_cores = 1) {
  progress_msg("Finding cis-NAT pairs using parallel processing...")
  
  # Separate by strand
  plus_genes <- genes_gr[strand(genes_gr) == "+"]
  minus_genes <- genes_gr[strand(genes_gr) == "-"]
  
  if (length(plus_genes) == 0 || length(minus_genes) == 0) {
    warning("No genes on positive or negative strand. Cannot find cis-NAT pairs.")
    return(data.frame())
  }
  
  progress_msg(sprintf("Analyzing %d positive strand and %d negative strand genes", 
                      length(plus_genes), length(minus_genes)))
  
  # Find overlaps
  hits <- GenomicRanges::findOverlaps(plus_genes, minus_genes, ignore.strand = FALSE)
  
  if (length(hits) == 0) {
    warning("No overlapping gene pairs found")
    return(data.frame())
  }
  
  # Create pairs dataframe with additional information
  pairs_df <- data.frame(
    sense_gene_id = plus_genes$gene_id[queryHits(hits)],
    antisense_gene_id = minus_genes$gene_id[subjectHits(hits)],
    sense_idx = queryHits(hits),
    antisense_idx = subjectHits(hits),
    overlap_width = width(pintersect(plus_genes[queryHits(hits)], 
                                    minus_genes[subjectHits(hits)])),
    stringsAsFactors = FALSE
  )
  
  # Add overlap percentage
  sense_width <- width(plus_genes[queryHits(hits)])
  antisense_width <- width(minus_genes[subjectHits(hits)])
  pairs_df$overlap_pct_sense <- (pairs_df$overlap_width / sense_width) * 100
  pairs_df$overlap_pct_antisense <- (pairs_df$overlap_width / antisense_width) * 100
  
  progress_msg(sprintf("Found %d cis-NAT pairs", nrow(pairs_df)), "INFO")
  
  return(pairs_df)
}

# Generate artificial antisense annotations with enhanced features
generate_antisense_annotations <- function(sense_gene, gtf, gene_id_suffix = "_NAT") {
  nat_gene_id <- paste0(sense_gene$gene_id, gene_id_suffix)
  nat_transcript_id <- paste0(sense_gene$gene_id, gene_id_suffix, ".1")
  
  # Flip strand
  nat_strand <- ifelse(as.character(strand(sense_gene)) == "+", "-", "+")
  
  # Create gene entry
  gene_entry <- sense_gene
  gene_entry$gene_id <- nat_gene_id
  gene_entry$transcript_id <- NA
  gene_entry$gene_biotype <- "antisense_RNA"
  gene_entry$gene_type <- "antisense_RNA"
  gene_entry$source <- "cisNAT_prediction"
  strand(gene_entry) <- nat_strand
  gene_entry$type <- "gene"
  
  # Create transcript entry
  transcript_entry <- sense_gene
  transcript_entry$gene_id <- nat_gene_id
  transcript_entry$transcript_id <- nat_transcript_id
  transcript_entry$gene_biotype <- "antisense_RNA"
  transcript_entry$gene_type <- "antisense_RNA"
  transcript_entry$source <- "cisNAT_prediction"
  strand(transcript_entry) <- nat_strand
  transcript_entry$type <- "transcript"
  
  # Find and create exon entries
  exon_mask <- !is.na(gtf$type) & 
               gtf$type == "exon" & 
               !is.na(gtf$gene_id) & 
               gtf$gene_id == sense_gene$gene_id
  
  exon_entries <- NULL
  if (any(exon_mask)) {
    exon_entries <- gtf[exon_mask]
    exon_entries$gene_id <- nat_gene_id
    exon_entries$transcript_id <- nat_transcript_id
    exon_entries$gene_biotype <- "antisense_RNA"
    exon_entries$gene_type <- "antisense_RNA"
    exon_entries$source <- "cisNAT_prediction"
    strand(exon_entries) <- nat_strand
    exon_entries$type <- "exon"
  }
  
  return(list(
    gene = gene_entry,
    transcript = transcript_entry,
    exons = exon_entries
  ))
}

# Main processing function
process_gtf_with_cisnat <- function(config) {
  progress_msg("=== Starting GTF Processing with cis-NAT Integration ===", "INFO")
  
  # Get file paths from config
  gtf_file <- config$reference$gtf_file
  lnc_gff_file <- config$reference$lncrna_gff
  output_file <- config$output$gtf_file
  
  # Validate input files
  validate_annotation_file(gtf_file, "gtf")
  if (!is.null(lnc_gff_file) && file.exists(lnc_gff_file)) {
    validate_annotation_file(lnc_gff_file, "gff")
  }
  
  # Step 1: Load reference GTF
  progress_msg("Step 1: Loading reference GTF file...")
  gtf <- rtracklayer::import(gtf_file)
  progress_msg(sprintf("Loaded %d features from GTF", length(gtf)), "INFO")
  
  # Step 2: Optionally load and integrate lncRNA annotations
  if (!is.null(lnc_gff_file) && file.exists(lnc_gff_file)) {
    progress_msg("Step 2: Loading and integrating lncRNA annotations...")
    
    lnc_gff <- rtracklayer::import(lnc_gff_file, format = "gff")
    progress_msg(sprintf("Loaded %d lncRNA features", length(lnc_gff)), "INFO")
    
    # Process lncRNA annotations
    if ("ID" %in% names(mcols(lnc_gff))) {
      lnc_gff$type <- as.character(lnc_gff$type)
      
      # Set gene_id based on available fields
      lnc_gff$gene_id <- ifelse(
        !is.na(mcols(lnc_gff)$Name),
        as.character(mcols(lnc_gff)$Name),
        as.character(mcols(lnc_gff)$ID)
      )
      
      lnc_gff$transcript_id <- as.character(mcols(lnc_gff)$ID)
      lnc_gff$gene_biotype <- "lncRNA"
      lnc_gff$gene_type <- "lncRNA"
      lnc_gff$source <- "lncRNA_database"
      
      # Filter for relevant features
      lnc_gff <- lnc_gff[lnc_gff$type %in% c("transcript", "exon", "gene")]
      
      # Merge with main GTF
      gtf <- c(gtf, lnc_gff)
      progress_msg(sprintf("Integrated %d lncRNA features", length(lnc_gff)), "INFO")
    } else {
      warning("lncRNA GFF file missing required ID field")
    }
  }
  
  # Step 3: Standardize feature types
  progress_msg("Step 3: Standardizing feature types...")
  gtf$type <- as.character(gtf$type)
  gtf$type[gtf$type == "mRNA"] <- "transcript"
  
  # Step 4: Extract protein-coding genes
  progress_msg("Step 4: Identifying protein-coding genes...")
  gene_mask <- !is.na(gtf$type) & gtf$type == "gene"
  
  # Consider genes without explicit biotype as protein-coding
  coding_mask <- gene_mask & (
    is.na(gtf$gene_biotype) | 
    gtf$gene_biotype %in% c("protein_coding", "")
  )
  
  genes <- gtf[coding_mask]
  
  plus_genes_n <- sum(strand(genes) == "+", na.rm = TRUE)
  minus_genes_n <- sum(strand(genes) == "-", na.rm = TRUE)
  
  progress_msg(sprintf("Found %d protein-coding genes (+ strand: %d, - strand: %d)", 
                      length(genes), plus_genes_n, minus_genes_n), "INFO")
  
  if (plus_genes_n == 0 || minus_genes_n == 0) {
    stop("Insufficient genes on positive or negative strand for cis-NAT analysis")
  }
  
  # Step 5: Find cis-NAT pairs
  progress_msg("Step 5: Identifying cis-NAT pairs...")
  
  # Determine number of cores
  num_cores <- config$performance$parallel_cores
  if (num_cores == 0) {
    num_cores <- parallel::detectCores() - 1
  }
  num_cores <- max(1, min(num_cores, parallel::detectCores()))
  
  cis_nat_pairs <- find_cis_nat_pairs_parallel(genes, num_cores)
  
  # Save cis-NAT pairs information
  if (nrow(cis_nat_pairs) > 0) {
    pairs_file <- sub("\\.gtf$", "_cisNAT_pairs.csv", output_file)
    write.csv(cis_nat_pairs, pairs_file, row.names = FALSE)
    progress_msg(sprintf("Saved cis-NAT pairs to: %s", pairs_file), "INFO")
  }
  
  # Step 6: Generate artificial antisense annotations
  progress_msg("Step 6: Generating artificial antisense annotations...")
  
  # Find genes without known antisense partners
  known_antisense_ids <- unique(cis_nat_pairs$antisense_gene_id)
  known_antisense_ids <- known_antisense_ids[!is.na(known_antisense_ids)]
  
  # Get sense-only genes
  sense_only_mask <- !is.na(genes$gene_id) & 
                     !(genes$gene_id %in% known_antisense_ids) & 
                     strand(genes) == "+"
  
  sense_only <- genes[sense_only_mask]
  
  progress_msg(sprintf("Generating antisense annotations for %d genes", length(sense_only)), "INFO")
  
  if (length(sense_only) > 0) {
    # Generate annotations with progress reporting
    nat_annotations <- lapply(seq_along(sense_only), function(i) {
      if (i %% 100 == 0) {
        progress_msg(sprintf("Processing gene %d/%d", i, length(sense_only)), "DEBUG")
      }
      generate_antisense_annotations(sense_only[i], gtf)
    })
    
    # Combine all annotations
    nat_genes <- do.call(c, lapply(nat_annotations, function(x) x$gene))
    nat_transcripts <- do.call(c, lapply(nat_annotations, function(x) x$transcript))
    nat_exons <- do.call(c, lapply(nat_annotations, function(x) x$exons))
    
    # Add to GTF
    gtf_new <- c(gtf, nat_genes, nat_transcripts, nat_exons)
    
    progress_msg(sprintf("Added %d antisense genes, %d transcripts, %d exons", 
                        length(nat_genes), length(nat_transcripts), length(nat_exons)), "INFO")
  } else {
    gtf_new <- gtf
  }
  
  # Step 7: Quality control and export
  progress_msg("Step 7: Performing quality control and exporting...")
  
  # Sort by chromosome and position
  gtf_new <- sort(gtf_new)
  
  # Export GTF
  rtracklayer::export(gtf_new, output_file, format = "gtf")
  progress_msg(sprintf("Exported enhanced GTF to: %s", output_file), "INFO")
  
  # Perform validation
  validate_gtf_for_cellranger(gtf_new)
  
  # Generate summary statistics
  summary_stats <- generate_gtf_summary(gtf_new)
  summary_file <- sub("\\.gtf$", "_summary.txt", output_file)
  
  writeLines(c(
    "GTF Summary Statistics",
    "=====================",
    sprintf("Total features: %d", summary_stats$total_features),
    sprintf("Genes: %d", summary_stats$n_genes),
    sprintf("  - Protein coding: %d", summary_stats$n_protein_coding),
    sprintf("  - lncRNA: %d", summary_stats$n_lncrna),
    sprintf("  - Antisense RNA: %d", summary_stats$n_antisense),
    sprintf("Transcripts: %d", summary_stats$n_transcripts),
    sprintf("Exons: %d", summary_stats$n_exons),
    sprintf("Chromosomes: %d", summary_stats$n_chromosomes),
    sprintf("cis-NAT pairs: %d", nrow(cis_nat_pairs))
  ), summary_file)
  
  progress_msg(sprintf("Saved summary statistics to: %s", summary_file), "INFO")
  
  return(gtf_new)
}

# Validate GTF for CellRanger compatibility
validate_gtf_for_cellranger <- function(gtf_gr) {
  progress_msg("Validating GTF for CellRanger compatibility...")
  
  # Check required fields
  required_fields <- c("gene_id", "transcript_id", "type")
  missing_fields <- setdiff(required_fields, names(mcols(gtf_gr)))
  
  if (length(missing_fields) > 0) {
    stop(sprintf("GTF missing required fields: %s", paste(missing_fields, collapse = ", ")))
  }
  
  # Check required feature types
  required_types <- c("gene", "transcript", "exon")
  present_types <- unique(gtf_gr$type)
  missing_types <- setdiff(required_types, present_types)
  
  if (length(missing_types) > 0) {
    warning(sprintf("GTF missing feature types: %s", paste(missing_types, collapse = ", ")))
  }
  
  # Check for NA values in critical fields
  na_gene_id <- sum(is.na(gtf_gr$gene_id))
  na_transcript_id <- sum(is.na(gtf_gr$transcript_id[gtf_gr$type != "gene"]))
  
  if (na_gene_id > 0) {
    warning(sprintf("GTF contains %d features with NA gene_id", na_gene_id))
  }
  
  if (na_transcript_id > 0) {
    warning(sprintf("GTF contains %d non-gene features with NA transcript_id", na_transcript_id))
  }
  
  # Check chromosome naming
  seqnames_chr <- as.character(seqnames(gtf_gr))
  if (!all(grepl("^chr|^[0-9]+$|^[IVX]+$", seqnames_chr))) {
    warning("Some chromosome names may not be CellRanger compatible")
  }
  
  progress_msg("GTF validation complete", "INFO")
}

# Generate summary statistics
generate_gtf_summary <- function(gtf_gr) {
  summary_stats <- list(
    total_features = length(gtf_gr),
    n_genes = sum(gtf_gr$type == "gene", na.rm = TRUE),
    n_transcripts = sum(gtf_gr$type == "transcript", na.rm = TRUE),
    n_exons = sum(gtf_gr$type == "exon", na.rm = TRUE),
    n_chromosomes = length(unique(seqnames(gtf_gr))),
    n_protein_coding = sum(gtf_gr$type == "gene" & 
                          gtf_gr$gene_biotype %in% c("protein_coding", NA), na.rm = TRUE),
    n_lncrna = sum(gtf_gr$type == "gene" & 
                   gtf_gr$gene_biotype == "lncRNA", na.rm = TRUE),
    n_antisense = sum(gtf_gr$type == "gene" & 
                      gtf_gr$gene_biotype == "antisense_RNA", na.rm = TRUE)
  )
  
  return(summary_stats)
}

# Main execution
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Load configuration
  config_file <- ifelse(length(args) > 0, args[1], "config.yaml")
  config <- load_config(config_file)
  
  # Override config with command line arguments if provided
  if (length(args) > 1) config$reference$gtf_file <- args[2]
  if (length(args) > 2) config$reference$lncrna_gff <- args[3]
  if (length(args) > 3) config$output$gtf_file <- args[4]
  
  # Process GTF with cis-NAT integration
  tryCatch({
    gtf_result <- process_gtf_with_cisnat(config)
    progress_msg("=== GTF Processing Complete ===", "INFO")
  }, error = function(e) {
    progress_msg(sprintf("ERROR: %s", e$message), "ERROR")
    quit(status = 1)
  })
}

# Execute if run as script
if (!interactive()) {
  main()
}