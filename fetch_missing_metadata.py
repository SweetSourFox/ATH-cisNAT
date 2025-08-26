#!/usr/bin/env python3
"""
Enhanced NCBI/SRA Metadata Fetcher for Missing Samples

This script fetches missing metadata from NCBI/SRA databases with improved
error handling, parallel processing, and configuration support.
"""

import os
import sys
import yaml
import logging
import argparse
import requests
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import time
import json

# Setup logging
def setup_logging(config: Dict[str, Any]) -> logging.Logger:
    """Setup logging configuration"""
    log_config = config.get('logging', {})
    log_level = getattr(logging, log_config.get('level', 'INFO'))
    log_format = log_config.get('format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Create logs directory if needed
    log_dir = Path(config.get('project', {}).get('log_dir', './logs'))
    log_dir.mkdir(exist_ok=True)
    
    # Configure logging
    logger = logging.getLogger('metadata_fetcher')
    logger.setLevel(log_level)
    
    # Console handler
    if log_config.get('console_logging', True):
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(console_handler)
    
    # File handler
    if log_config.get('file_logging', True):
        file_handler = logging.FileHandler(
            log_dir / f'metadata_fetch_{datetime.now():%Y%m%d_%H%M%S}.log'
        )
        file_handler.setFormatter(logging.Formatter(log_format))
        logger.addHandler(file_handler)
    
    return logger

class NCBIFetcher:
    """Enhanced NCBI/SRA metadata fetcher with retry logic and rate limiting"""
    
    def __init__(self, config: Dict[str, Any], logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.ncbi_config = config.get('ncbi', {})
        
        # NCBI E-utilities base URLs
        self.esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        # Rate limiting
        self.delay = self.ncbi_config.get('delay_seconds', 1)
        self.retry_attempts = self.ncbi_config.get('retry_attempts', 3)
        
        # Session for connection pooling
        self.session = requests.Session()
        
    def search_sample(self, accession: str) -> Optional[Dict[str, Any]]:
        """Search for a sample in NCBI/SRA with retry logic"""
        self.logger.info(f"Searching for sample: {accession}")
        
        for attempt in range(self.retry_attempts):
            try:
                # Search for the accession
                search_params = {
                    'db': 'sra',
                    'term': accession,
                    'retmode': 'json',
                    'retmax': 1
                }
                
                # Add email and API key if provided
                if self.ncbi_config.get('email'):
                    search_params['email'] = self.ncbi_config['email']
                if self.ncbi_config.get('api_key'):
                    search_params['api_key'] = self.ncbi_config['api_key']
                
                response = self.session.get(
                    self.esearch_url, 
                    params=search_params, 
                    timeout=30
                )
                response.raise_for_status()
                
                search_data = response.json()
                
                if not search_data.get('esearchresult', {}).get('idlist'):
                    self.logger.warning(f"No results found for {accession}")
                    return None
                
                uid = search_data['esearchresult']['idlist'][0]
                
                # Fetch detailed information
                fetch_params = {
                    'db': 'sra',
                    'id': uid,
                    'rettype': 'xml',
                    'retmode': 'text'
                }
                
                if self.ncbi_config.get('email'):
                    fetch_params['email'] = self.ncbi_config['email']
                if self.ncbi_config.get('api_key'):
                    fetch_params['api_key'] = self.ncbi_config['api_key']
                
                response = self.session.get(
                    self.efetch_url, 
                    params=fetch_params, 
                    timeout=30
                )
                response.raise_for_status()
                
                # Parse XML response
                root = ET.fromstring(response.content)
                metadata = self._extract_metadata(root, accession)
                
                # Rate limiting
                time.sleep(self.delay)
                
                return metadata
                
            except requests.exceptions.RequestException as e:
                self.logger.error(f"Request error for {accession} (attempt {attempt + 1}): {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(self.delay * (attempt + 1))  # Exponential backoff
                else:
                    return None
                    
            except Exception as e:
                self.logger.error(f"Unexpected error for {accession}: {e}")
                return None
    
    def _extract_metadata(self, root: ET.Element, accession: str) -> Dict[str, Any]:
        """Extract comprehensive metadata from SRA XML"""
        metadata = {
            'Run': accession,
            'fetch_date': datetime.now().isoformat()
        }
        
        # Navigate through XML structure with error handling
        for package in root.findall('.//EXPERIMENT_PACKAGE'):
            # Extract experiment information
            experiment = package.find('.//EXPERIMENT')
            if experiment is not None:
                self._extract_experiment_info(experiment, metadata)
            
            # Extract sample information
            sample = package.find('.//SAMPLE')
            if sample is not None:
                self._extract_sample_info(sample, metadata)
            
            # Extract study information
            study = package.find('.//STUDY')
            if study is not None:
                self._extract_study_info(study, metadata)
            
            # Extract run information
            run = package.find('.//RUN')
            if run is not None:
                self._extract_run_info(run, metadata)
        
        self.logger.debug(f"Extracted metadata fields: {list(metadata.keys())}")
        return metadata
    
    def _extract_experiment_info(self, experiment: ET.Element, metadata: Dict[str, Any]):
        """Extract experiment-related metadata"""
        # Title
        title = experiment.find('.//TITLE')
        if title is not None and title.text:
            metadata['experiment_title'] = title.text
        
        # Library information
        library_desc = experiment.find('.//LIBRARY_DESCRIPTOR')
        if library_desc is not None:
            for field in ['LIBRARY_NAME', 'LIBRARY_STRATEGY', 'LIBRARY_SOURCE', 
                         'LIBRARY_SELECTION', 'LIBRARY_LAYOUT']:
                elem = library_desc.find(f'.//{field}')
                if elem is not None:
                    if field == 'LIBRARY_LAYOUT':
                        # Handle paired/single layout
                        layout_type = elem.find('.//PAIRED') if elem.find('.//PAIRED') is not None else elem.find('.//SINGLE')
                        if layout_type is not None:
                            metadata['library_layout'] = layout_type.tag.lower()
                    else:
                        metadata[field.lower()] = elem.text if elem.text else ''
        
        # Platform information
        platform = experiment.find('.//PLATFORM')
        if platform is not None:
            # Get platform type (ILLUMINA, PACBIO, etc.)
            for child in platform:
                if child.tag != 'PLATFORM':
                    metadata['platform'] = child.tag
                    # Get instrument model
                    instrument = child.find('.//INSTRUMENT_MODEL')
                    if instrument is not None and instrument.text:
                        metadata['instrument_model'] = instrument.text
                    break
    
    def _extract_sample_info(self, sample: ET.Element, metadata: Dict[str, Any]):
        """Extract sample-related metadata"""
        # Scientific name
        sci_name = sample.find('.//SCIENTIFIC_NAME')
        if sci_name is not None and sci_name.text:
            metadata['organism'] = sci_name.text
        
        # Common name
        common_name = sample.find('.//COMMON_NAME')
        if common_name is not None and common_name.text:
            metadata['common_name'] = common_name.text
        
        # Sample attributes
        for attr in sample.findall('.//SAMPLE_ATTRIBUTE'):
            tag = attr.find('TAG')
            value = attr.find('VALUE')
            if tag is not None and value is not None and tag.text and value.text:
                # Standardize common attribute names
                tag_text = tag.text.lower().replace(' ', '_')
                metadata[f'sample_{tag_text}'] = value.text
    
    def _extract_study_info(self, study: ET.Element, metadata: Dict[str, Any]):
        """Extract study-related metadata"""
        # Study title
        study_title = study.find('.//STUDY_TITLE')
        if study_title is not None and study_title.text:
            metadata['study_title'] = study_title.text
        
        # Study abstract
        study_abstract = study.find('.//STUDY_ABSTRACT')
        if study_abstract is not None and study_abstract.text:
            metadata['study_abstract'] = study_abstract.text
    
    def _extract_run_info(self, run: ET.Element, metadata: Dict[str, Any]):
        """Extract run-related metadata"""
        # Run statistics
        for stat in ['total_spots', 'total_bases', 'size', 'published']:
            elem = run.find(f'.//{stat.upper()}')
            if elem is not None and elem.text:
                metadata[stat] = elem.text

def fetch_samples_parallel(samples: List[str], config: Dict[str, Any], 
                          logger: logging.Logger) -> List[Dict[str, Any]]:
    """Fetch multiple samples in parallel"""
    fetcher = NCBIFetcher(config, logger)
    all_metadata = []
    
    # Determine number of workers
    max_workers = min(len(samples), config.get('performance', {}).get('parallel_cores', 4))
    if max_workers == 0:  # Auto-detect
        max_workers = min(len(samples), os.cpu_count() or 4)
    
    logger.info(f"Fetching {len(samples)} samples using {max_workers} parallel workers")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_sample = {
            executor.submit(fetcher.search_sample, sample): sample 
            for sample in samples
        }
        
        # Process completed tasks
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                metadata = future.result()
                if metadata:
                    all_metadata.append(metadata)
                    logger.info(f"Successfully fetched metadata for {sample}")
                else:
                    logger.warning(f"No metadata retrieved for {sample}")
            except Exception as e:
                logger.error(f"Error processing {sample}: {e}")
    
    return all_metadata

def update_metadata_file(new_metadata: List[Dict[str, Any]], 
                        config: Dict[str, Any], 
                        logger: logging.Logger) -> pd.DataFrame:
    """Update or create metadata file with new information"""
    metadata_file = Path(config.get('input', {}).get('metadata_file', './Integrated_metadata.csv'))
    output_dir = Path(config.get('project', {}).get('output_dir', './results'))
    output_dir.mkdir(exist_ok=True)
    
    # Convert to DataFrame
    new_df = pd.DataFrame(new_metadata)
    
    if metadata_file.exists():
        logger.info(f"Updating existing metadata file: {metadata_file}")
        existing_df = pd.read_csv(metadata_file)
        
        # Merge with existing data, avoiding duplicates
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
        combined_df = combined_df.drop_duplicates(subset=['Run'], keep='last')
    else:
        logger.info(f"Creating new metadata file: {metadata_file}")
        combined_df = new_df
    
    # Save in multiple formats
    output_formats = config.get('output', {}).get('formats', ['csv'])
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    for fmt in output_formats:
        output_file = output_dir / f'metadata_updated_{timestamp}.{fmt}'
        if fmt == 'csv':
            combined_df.to_csv(output_file, index=False)
        elif fmt == 'xlsx':
            combined_df.to_excel(output_file, index=False, engine='openpyxl')
        elif fmt == 'json':
            combined_df.to_json(output_file, orient='records', indent=2)
        logger.info(f"Saved metadata to: {output_file}")
    
    # Also save to the standard location
    combined_df.to_csv(metadata_file.parent / 'Integrated_metadata_updated.csv', index=False)
    
    return combined_df

def main():
    """Main execution function"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Fetch missing metadata from NCBI/SRA databases'
    )
    parser.add_argument(
        '--config', '-c',
        type=str,
        default='config.yaml',
        help='Configuration file path (default: config.yaml)'
    )
    parser.add_argument(
        '--samples', '-s',
        type=str,
        nargs='+',
        help='Sample accessions to fetch (overrides default list)'
    )
    parser.add_argument(
        '--input-file', '-i',
        type=str,
        help='File containing sample accessions (one per line)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        help='Output directory (overrides config)'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config = {}
    if Path(args.config).exists():
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    else:
        print(f"Warning: Config file {args.config} not found, using defaults")
    
    # Override config with command line arguments
    if args.output_dir:
        config.setdefault('project', {})['output_dir'] = args.output_dir
    
    # Setup logging
    logger = setup_logging(config)
    
    logger.info("="*70)
    logger.info("NCBI/SRA METADATA FETCHER - ENHANCED VERSION")
    logger.info("="*70)
    
    # Determine samples to fetch
    if args.samples:
        samples = args.samples
    elif args.input_file:
        with open(args.input_file, 'r') as f:
            samples = [line.strip() for line in f if line.strip()]
    else:
        # Default sample list
        samples = [
            'CRR151722', 'CRR151723', 'CRR151724', 'CRR151725',
            'CRR995220', 'CRR995221', 'CRR995222',
            'SRR15992287', 'SRR15992288', 'SRR21622665', 'SRR21622666',
            'SRR28356471', 'SRR32380306'
        ]
    
    logger.info(f"Samples to fetch: {len(samples)}")
    
    # Check NCBI email requirement
    if not config.get('ncbi', {}).get('email'):
        logger.warning("NCBI email not configured. Please add your email to config.yaml")
        logger.warning("This is required by NCBI E-utilities terms of use")
    
    # Fetch metadata
    start_time = time.time()
    metadata_list = fetch_samples_parallel(samples, config, logger)
    elapsed_time = time.time() - start_time
    
    if metadata_list:
        logger.info(f"Successfully fetched metadata for {len(metadata_list)} samples")
        logger.info(f"Time elapsed: {elapsed_time:.2f} seconds")
        
        # Update metadata file
        updated_df = update_metadata_file(metadata_list, config, logger)
        
        # Print summary statistics
        logger.info("\nMetadata Summary:")
        logger.info(f"Total records: {len(updated_df)}")
        
        # Tissue distribution
        if 'sample_tissue' in updated_df.columns:
            tissue_counts = updated_df['sample_tissue'].value_counts()
            logger.info("\nTissue distribution:")
            for tissue, count in tissue_counts.items():
                logger.info(f"  {tissue}: {count}")
        
        # Organism distribution
        if 'organism' in updated_df.columns:
            organism_counts = updated_df['organism'].value_counts()
            logger.info("\nOrganism distribution:")
            for organism, count in organism_counts.items():
                logger.info(f"  {organism}: {count}")
    else:
        logger.warning("No metadata could be fetched")
    
    logger.info("\n" + "="*70)
    logger.info("METADATA FETCHING COMPLETE")
    logger.info("="*70)

if __name__ == "__main__":
    main()