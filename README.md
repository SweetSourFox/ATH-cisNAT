# AThNATCount - Arabidopsis thaliana Natural Antisense Transcript Analysis Pipeline

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-green.svg)](https://www.r-project.org/)

A comprehensive computational pipeline for analyzing sense-antisense transcript correlations in *Arabidopsis thaliana* single-cell RNA-seq data, with special focus on Leucine-Rich Repeat Receptor-Like Kinases (LRR-RLKs) and tissue-specific regulation patterns.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Components](#pipeline-components)
- [Configuration](#configuration)
- [Usage Examples](#usage-examples)
- [Output Description](#output-description)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

## Overview

AThNATCount is an advanced bioinformatics pipeline designed to:

1. **Integrate cis-Natural Antisense Transcripts (cis-NATs)** into reference genome annotations
2. **Analyze sense-antisense correlations** from 10X Genomics single-cell RNA-seq data
3. **Perform tissue-specific analysis** across multiple sample types and conditions
4. **Focus on LRR-RLK gene families** for plant immunity and development studies
5. **Process multiple data sources** including NCBI/SRA and China National GeneBank (CNGB)
6. **Generate publication-ready visualizations** and comprehensive reports

The pipeline combines Python and R scripts with parallel processing capabilities, automated quality control, and advanced statistical modeling.

## Features

### Core Capabilities

- **Automated GTF Enhancement**: Adds artificial antisense annotations for genes lacking natural antisense partners
- **Parallel Processing**: Utilizes multiple CPU cores for faster computation with configurable chunking
- **Memory Efficient**: Handles large single-cell datasets with streaming and chunking
- **Comprehensive Statistics**: Pearson, Spearman, and Kendall correlation methods with multiple testing correction
- **Advanced Quality Control**: Automated validation and error handling

### Advanced Features

- **Machine Learning Integration**: Pattern recognition for LRR-RLK identification with confidence scoring
- **Network Analysis**: Interaction networks for sense-antisense pairs with centrality metrics
- **Interactive Visualizations**: Publication-ready figures with customizable themes and color palettes
- **Batch Processing**: Handle multiple samples automatically with cross-sample comparisons
- **Configuration Management**: YAML-based configuration for reproducibility and customization
- **Metadata Integration**: Automated fetching and processing of sample metadata from multiple databases
- **Tissue-Specific Analysis**: Advanced statistical modeling for tissue-specific regulation patterns

### Data Source Support

- **10X Genomics**: Single-cell RNA-seq data processing
- **NCBI/SRA**: Automated metadata fetching with rate limiting and retry logic
- **CNGB/CRR**: China National GeneBank sample processing
- **Multiple Formats**: Support for CSV, Excel, JSON, and HDF5 outputs

## Requirements

### System Requirements

- **Operating System**: Linux/macOS (Windows via WSL2)
- **Memory**: Minimum 8GB RAM (16GB+ recommended for large datasets)
- **Storage**: 10GB+ free space for analysis and results
- **CPU**: Multi-core processor recommended (8+ cores for optimal performance)

### Software Dependencies

- **Python** ≥ 3.10
- **R** ≥ 4.3
- **Conda** (recommended for environment management)

### Key R Packages

- **Seurat**: Single-cell analysis
- **Matrix**: Sparse matrix operations
- **dplyr/ggplot2**: Data manipulation and visualization
- **rtracklayer/GenomicRanges**: Genomic data processing
- **ComplexHeatmap**: Advanced heatmap visualizations
- **igraph**: Network analysis

### Key Python Packages

- **pandas**: Data manipulation
- **requests**: API interactions
- **yaml**: Configuration management
- **concurrent.futures**: Parallel processing

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/AThNATCount.git
cd AThNATCount
```

### 2. Create Conda Environment

```bash
# Create environment from YAML file
conda env create -f environment.yml

# Activate environment
conda activate athnatcount
```

### 3. Install R Packages

```bash
# Install all required R packages
Rscript R/install_packages.R
```

### 4. Verify Installation

```bash
# Test the installation
Rscript test_batch_analysis.R
```

## Quick Start

### 1. Prepare Configuration

Edit `config.yaml` to specify your paths and parameters:

```yaml
project:
  name: "MyAnalysis"
  output_dir: "./results"

reference:
  gtf_file: "./ATH_APNAT/genes/genes.gtf.gz"
  lncrna_gff: "./Ref/A.tha.lncRNA.gff"
  
input:
  cellranger_dir: "./Cellranger_Results"
  metadata_file: "./01_Data/Processed_Data/Integrated_metadata_final.csv"
  sa_results_dir: "./test_results/batch_analysis"

analysis:
  correlation:
    min_shared_cells: 3
    min_expression: 0
    correlation_method: "pearson"
  
  tissue:
    consolidate_tissues: true
    tissue_mapping:
      - pattern: "root"
        category: "Root"
      - pattern: "leaf|leaves"
        category: "Leaf"
      - pattern: "seed"
        category: "Seed"
      - pattern: "flower"
        category: "Flower"
  
  lrr_rlk:
    patterns: ["LRR", "RLK", "RECEPTOR", "KINASE", "CLV", "FLS", "BAK", "BRI", "SERK"]
    min_correlation: 0.3
    pvalue_threshold: 0.05
```

### 2. Run the Complete Pipeline

```bash
# Step 1: Fetch missing metadata (if needed)
python fetch_missing_metadata.py --samples SRR123456 SRR789012

# Step 2: Process CRR samples (if applicable)
python search_crr_samples.py --samples CRR151722 CRR151723

# Step 3: Enhance GTF with cis-NAT annotations
Rscript gtf_add_cisnat.R

# Step 4: Fix GTF transcript/exon structure
Rscript gtf_transcript_exon_fix.R input.gtf output.gtf

# Step 5: Analyze sense-antisense correlations (single sample)
Rscript sense_nat_correlation.R -i ./Cellranger_Results/sample/outs/filtered_feature_bc_matrix.h5 -o ./results/sample1

# Step 6: Run batch analysis (multiple samples)
Rscript batch_simple.R

# Step 7: Integrate tissue-specific data
Rscript tissue_specific_sa_analysis.R

# Step 8: Analyze LRR-RLK patterns
Rscript lrr_rlk_analysis.R
```

## Pipeline Components

### Python Scripts

#### `fetch_missing_metadata.py`
Enhanced metadata fetcher for NCBI/SRA databases with parallel processing and comprehensive error handling.

```bash
python fetch_missing_metadata.py --samples SRR123456 SRR789012 --output results/
```

**Features:**
- Parallel NCBI queries with rate limiting
- Automatic retry on failure with exponential backoff
- Comprehensive metadata extraction (experiment, sample, study, run info)
- Multiple output formats (CSV, Excel, JSON)
- Configurable email and API key support

#### `search_crr_samples.py`
Processor for China National GeneBank (CNGB) CRR samples with metadata integration.

```bash
python search_crr_samples.py --samples CRR151722 CRR151723 --input metadata.csv
```

**Features:**
- CRR sample ID parsing and classification
- Integration with existing metadata
- SA_Results availability checking
- Comprehensive metadata analysis

### R Scripts

#### `gtf_add_cisnat.R`
Integrates cis-NAT information into reference GTF files with advanced features.

**Key features:**
- Automatic lncRNA integration from GFF files
- Artificial antisense generation for genes without partners
- Quality validation and CellRanger compatibility check
- Parallel processing for large datasets
- Comprehensive summary statistics

#### `gtf_transcript_exon_fix.R`
Ensures GTF compliance with single-cell analysis tools through comprehensive validation.

```bash
Rscript gtf_transcript_exon_fix.R input.gtf output.gtf
```

**Features:**
- Transcript ID addition to exons
- Gene ID propagation
- Synthetic exon generation for transcripts without exons
- Memory-efficient processing for large files
- Automatic backup creation

#### `sense_nat_correlation.R`
Core analysis script for sense-antisense correlations with advanced statistical capabilities.

```bash
Rscript sense_nat_correlation.R \
  -i cellranger_output.h5 \
  -o results/correlations \
  -m pearson \
  -p 8  # Use 8 cores
```

**Features:**
- Multiple correlation methods (Pearson, Spearman, Kendall)
- Parallel processing with configurable chunking
- Advanced filtering and quality control
- Comprehensive visualization suite
- Expression data export capabilities

#### `batch_simple.R`
Simple test script for batch analysis functionality with automated file discovery.

```bash
Rscript batch_simple.R
```

**Features:**
- Automatic discovery of CellRanger result files
- Sample name extraction from file paths
- Comprehensive testing and validation
- Progress reporting and error handling

#### `tissue_specific_sa_analysis.R`
Comprehensive tissue-specific analysis with advanced statistical modeling and visualization.

```bash
Rscript tissue_specific_sa_analysis.R config.yaml results/
```

**Features:**
- Automated data integration from multiple sources
- Advanced statistical modeling (mixed effects, Bayesian)
- Machine learning for pattern discovery
- Interactive dashboards and reports
- Publication-ready visualizations

#### `lrr_rlk_analysis.R`
Specialized analysis for LRR-RLK gene families with network analysis and pattern recognition.

```bash
Rscript lrr_rlk_analysis.R
```

**Features:**
- Advanced LRR-RLK identification using multiple databases
- Machine learning-based pattern recognition
- Network analysis and visualization
- Interactive reporting with Shiny integration
- Publication-ready figures and tables

## Configuration

The pipeline uses a centralized YAML configuration file (`config.yaml`) with comprehensive settings:

### Project Settings
```yaml
project:
  name: "AThNATCount"
  description: "Arabidopsis thaliana Natural Antisense Transcript Analysis"
  output_dir: "./results"
  temp_dir: "./temp"
  log_dir: "./logs"
```

### Reference Files
```yaml
reference:
  gtf_file: "./ATH_APNAT/genes/genes.gtf.gz"
  lncrna_gff: "./Ref/A.tha.lncRNA.gff"
  lrr_rlk_ids: "./Ref/AtLRR_RLK_ids.txt"
  genome_fasta: "./ATH_APNAT/fasta/genome.fa"
```

### Analysis Parameters
```yaml
analysis:
  correlation:
    min_shared_cells: 3
    min_expression: 0
    correlation_method: "pearson"  # pearson, spearman, kendall
  
  tissue:
    consolidate_tissues: true
    tissue_mapping:
      - pattern: "root"
        category: "Root"
      - pattern: "leaf|leaves"
        category: "Leaf"
      - pattern: "seed"
        category: "Seed"
      - pattern: "flower"
        category: "Flower"
  
  lrr_rlk:
    patterns: ["LRR", "RLK", "RECEPTOR", "KINASE", "CLV", "FLS", "BAK", "BRI", "SERK"]
    min_correlation: 0.3
    pvalue_threshold: 0.05
```

### Performance Settings
```yaml
performance:
  parallel_cores: 0  # 0 = auto-detect
  chunk_size: 1000
  memory_limit: "80GB"
```

### NCBI/SRA Settings
```yaml
ncbi:
  email: ""  # Required for NCBI E-utilities
  api_key: ""  # Optional, increases rate limit
  retry_attempts: 3
  delay_seconds: 1
```

## Usage Examples

### Example 1: Basic Single Sample Analysis

```bash
# Analyze one 10X sample
Rscript sense_nat_correlation.R \
  -i Cellranger_Results/SRR12046060/outs/filtered_feature_bc_matrix.h5 \
  -o results/SRR12046060 \
  -c config.yaml
```

### Example 2: Batch Processing Multiple Samples

```bash
# Use the batch analysis script
Rscript batch_simple.R

# Or manually process multiple samples
Rscript sense_nat_correlation.R \
  -i "sample1.h5,sample2.h5,sample3.h5" \
  -n "SRR1,SRR2,SRR3" \
  -o results/batch_analysis
```

### Example 3: Tissue-Specific Analysis

```bash
# Run comprehensive tissue analysis
Rscript tissue_specific_sa_analysis.R config.yaml tissue_results/

# This will:
# - Load all correlation data from batch analysis
# - Integrate with metadata
# - Perform tissue-specific statistical analysis
# - Generate comprehensive visualizations
```

### Example 4: LRR-RLK Specialized Analysis

```bash
# Analyze LRR-RLK patterns
Rscript lrr_rlk_analysis.R

# This will:
# - Identify LRR-RLK genes using pattern matching
# - Perform subfamily classification
# - Create interaction networks
# - Generate specialized visualizations
```

### Example 5: Metadata Management

```bash
# Fetch missing metadata from NCBI
python fetch_missing_metadata.py \
  --samples SRR12046060 SRR12046121 SRR12046122 \
  --config config.yaml

# Process CRR samples
python search_crr_samples.py \
  --samples CRR151722 CRR151723 CRR151724 \
  --input 01_Data/Processed_Data/Integrated_metadata_final.csv
```

## Output Description

### Directory Structure

```
results/
├── figures/
│   ├── distributions/          # Correlation distribution plots
│   │   └── correlation_distributions.png
│   ├── comparisons/            # Tissue comparison plots  
│   │   ├── tissue_violin_plots.png
│   │   └── correlation_directions.png
│   ├── heatmaps/              # Expression heatmaps
│   │   └── gene_tissue_heatmap.png
│   └── networks/              # Interaction networks
├── tables/
│   ├── integrated_correlation_data.csv
│   ├── tissue_summary_statistics.csv
│   ├── gene_tissue_statistics.csv
│   └── tissue_specific_genes.csv
├── reports/
│   ├── analysis_summary.md
│   └── session_info.txt
└── logs/
    └── analysis_YYYYMMDD_HHMMSS.log

LRR_RLK_analysis/
├── lrr_rlk_correlation_density.png
├── lrr_rlk_expression_scatter.png
├── lrr_rlk_network.png
├── lrr_rlk_network_metrics.csv
├── LRR_RLK_analysis_data.csv
└── LRR_RLK_Analysis_Report.txt
```

### Key Output Files

1. **integrated_correlation_data.csv**: Primary correlation results
   - Columns: Gene, NAT, Correlation, P-value, N_Cells, Sample, Tissue

2. **tissue_summary_statistics.csv**: Tissue-level statistics
   - Mean correlations, sample counts, significance ratios, direction analysis

3. **gene_tissue_statistics.csv**: Gene-level tissue-specific analysis
   - Mean correlations per gene per tissue, consistency metrics

4. **lrr_rlk_network_metrics.csv**: Network analysis results
   - Degree, betweenness, closeness, eigenvector centrality, community membership

5. **LRR_RLK_Analysis_Report.txt**: Comprehensive LRR-RLK analysis report
   - Executive summary, key findings, statistical comparisons

## Troubleshooting

### Common Issues

#### Memory Errors
```bash
# Increase memory allocation for R
export R_MAX_VSIZE=16Gb
Rscript memory_intensive_script.R

# Or adjust in config.yaml
performance:
  memory_limit: "16GB"
```

#### Missing Packages
```bash
# Check and install missing R packages
Rscript -e "source('R/install_packages.R')"

# For Python packages
conda install package_name
# or
pip install package_name
```

#### File Format Issues
```bash
# Validate GTF files
Rscript gtf_transcript_exon_fix.R --help

# Check H5 file integrity
h5dump -H file.h5
```

#### NCBI API Issues
```bash
# Add your email to config.yaml
ncbi:
  email: "your.email@institution.edu"
  api_key: "your_api_key_optional"
```

### Debug Mode

Enable verbose logging in config.yaml:
```yaml
logging:
  level: "DEBUG"
  file_logging: true
  console_logging: true
```

### Performance Optimization

```yaml
performance:
  parallel_cores: 8  # Set to number of available cores
  chunk_size: 500    # Adjust based on memory
  memory_limit: "16GB"
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
# Clone with submodules
git clone --recursive https://github.com/yourusername/AThNATCount.git

# Install development dependencies
pip install -r requirements-dev.txt
```

## Citation

If you use AThNATCount in your research, please cite:

```
AThNATCount: A Comprehensive Pipeline for Arabidopsis thaliana Natural Antisense Transcript Analysis
[Your Name et al.]
[Journal/Conference]
[Year]
```

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE Version 3 License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Arabidopsis Information Resource (TAIR)
- 10X Genomics for single-cell tools
- R/Bioconductor community
- NCBI for data access
- China National GeneBank (CNGB)
- All contributors and users

---

For questions, bug reports, or feature requests, please open an issue on GitHub or contact the maintainers.
