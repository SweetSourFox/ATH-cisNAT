#!/usr/bin/env python3
"""
Enhanced CRR Sample Processor

Processes China National GeneBank (CNGB) CRR samples and integrates them
with existing metadata. Supports multiple data sources and formats.
"""

import os
import sys
import yaml
import logging
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime
import requests
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

class CRRProcessor:
    """Process and integrate CRR sample metadata"""
    
    def __init__(self, config: Dict[str, Any], logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.output_dir = Path(config.get('project', {}).get('output_dir', './results'))
        self.output_dir.mkdir(exist_ok=True)
        
        # CNGB API endpoints (if available in future)
        self.cngb_base_url = "https://db.cngb.org/api/v1"
        
    def generate_crr_metadata(self, sample_ids: List[str]) -> List[Dict[str, Any]]:
        """Generate metadata for CRR samples with enhanced information"""
        self.logger.info(f"Generating metadata for {len(sample_ids)} CRR samples")
        
        metadata_list = []
        
        for sample_id in sample_ids:
            # Extract information from sample ID pattern
            metadata = self._parse_sample_id(sample_id)
            
            # Add default values
            metadata.update({
                'Run': sample_id,
                'organism': 'Arabidopsis thaliana',
                'library_source': 'TRANSCRIPTOMIC',
                'platform': 'ILLUMINA',
                'sample_name': f'{sample_id}_sample',
                'database': 'GSA/CNGB',
                'data_source': 'CRR',
                'fetch_date': datetime.now().isoformat()
            })
            
            # Try to infer additional information
            if sample_id.startswith('CRR15'):
                metadata['project_hint'] = 'Early submission batch'
            elif sample_id.startswith('CRR99'):
                metadata['project_hint'] = 'Recent submission batch'
            
            metadata_list.append(metadata)
        
        return metadata_list
    
    def _parse_sample_id(self, sample_id: str) -> Dict[str, Any]:
        """Parse sample ID to extract potential information"""
        metadata = {}
        
        # CRR ID pattern analysis
        if sample_id.startswith('CRR'):
            id_number = sample_id[3:]
            if id_number.isdigit():
                id_num = int(id_number)
                # Estimate submission year based on ID ranges (hypothetical)
                if id_num < 200000:
                    metadata['estimated_year'] = '2018-2019'
                elif id_num < 1000000:
                    metadata['estimated_year'] = '2020-2021'
                else:
                    metadata['estimated_year'] = '2022-present'
        
        return metadata
    
    def check_sa_results(self, sample_ids: List[str]) -> Dict[str, bool]:
        """Check which samples have SA_Results data"""
        sa_dir = Path(self.config.get('input', {}).get('sa_results_dir', './SA_Results'))
        
        results = {}
        for sample_id in sample_ids:
            sample_path = sa_dir / sample_id
            results[sample_id] = sample_path.exists() and sample_path.is_dir()
        
        return results
    
    def integrate_with_existing(self, crr_metadata: List[Dict[str, Any]], 
                               existing_file: Optional[Path] = None) -> pd.DataFrame:
        """Integrate CRR metadata with existing metadata"""
        # Convert to DataFrame
        crr_df = pd.DataFrame(crr_metadata)
        
        # Check for existing metadata
        if existing_file and existing_file.exists():
            self.logger.info(f"Integrating with existing metadata: {existing_file}")
            existing_df = pd.read_csv(existing_file)
            
            # Remove duplicates if any
            existing_df = existing_df[~existing_df['Run'].isin(crr_df['Run'])]
            
            # Combine dataframes
            combined_df = pd.concat([existing_df, crr_df], ignore_index=True, sort=False)
        else:
            self.logger.info("Creating new metadata file")
            combined_df = crr_df
        
        # Sort by Run ID
        combined_df = combined_df.sort_values('Run')
        
        return combined_df
    
    def analyze_metadata_completeness(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Analyze metadata completeness and quality"""
        analysis = {
            'total_samples': len(df),
            'completeness': {},
            'tissue_distribution': {},
            'platform_distribution': {},
            'year_distribution': {}
        }
        
        # Check completeness for key fields
        key_fields = ['organism', 'tissue', 'platform', 'library_source']
        for field in key_fields:
            if field in df.columns:
                non_empty = df[field].notna() & (df[field] != '') & (df[field] != 'Unknown')
                analysis['completeness'][field] = {
                    'count': non_empty.sum(),
                    'percentage': (non_empty.sum() / len(df)) * 100
                }
        
        # Analyze distributions
        if 'tissue' in df.columns:
            analysis['tissue_distribution'] = df['tissue'].value_counts().to_dict()
        
        if 'platform' in df.columns:
            analysis['platform_distribution'] = df['platform'].value_counts().to_dict()
        
        return analysis

def setup_logging(log_level: str = 'INFO') -> logging.Logger:
    """Setup logging configuration"""
    logger = logging.getLogger('crr_processor')
    logger.setLevel(getattr(logging, log_level))
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    )
    logger.addHandler(console_handler)
    
    return logger

def main():
    """Main execution function"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Process CRR samples and integrate metadata'
    )
    parser.add_argument(
        '--config', '-c',
        type=str,
        default='config.yaml',
        help='Configuration file path'
    )
    parser.add_argument(
        '--input-metadata', '-i',
        type=str,
        help='Input metadata file to update'
    )
    parser.add_argument(
        '--samples', '-s',
        type=str,
        nargs='+',
        help='CRR sample IDs to process'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        help='Output file path'
    )
    parser.add_argument(
        '--analyze-only', '-a',
        action='store_true',
        help='Only analyze existing metadata'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config = {}
    if Path(args.config).exists():
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    
    # Setup logging
    logger = setup_logging(config.get('logging', {}).get('level', 'INFO'))
    
    logger.info("="*70)
    logger.info("CRR SAMPLE PROCESSOR - ENHANCED VERSION")
    logger.info("="*70)
    
    # Initialize processor
    processor = CRRProcessor(config, logger)
    
    # Analyze only mode
    if args.analyze_only:
        if args.input_metadata and Path(args.input_metadata).exists():
            df = pd.read_csv(args.input_metadata)
            analysis = processor.analyze_metadata_completeness(df)
            
            logger.info("\nMetadata Analysis Results:")
            logger.info(f"Total samples: {analysis['total_samples']}")
            
            logger.info("\nField completeness:")
            for field, stats in analysis['completeness'].items():
                logger.info(f"  {field}: {stats['count']}/{analysis['total_samples']} "
                          f"({stats['percentage']:.1f}%)")
            
            logger.info("\nTissue distribution:")
            for tissue, count in analysis['tissue_distribution'].items():
                logger.info(f"  {tissue}: {count}")
            
            # Save analysis report
            report_path = processor.output_dir / f'metadata_analysis_{datetime.now():%Y%m%d_%H%M%S}.json'
            with open(report_path, 'w') as f:
                json.dump(analysis, f, indent=2)
            logger.info(f"\nAnalysis report saved to: {report_path}")
        else:
            logger.error("Input metadata file required for analysis")
        return
    
    # Process CRR samples
    if args.samples:
        crr_samples = args.samples
    else:
        # Default CRR samples
        crr_samples = [
            'CRR151722', 'CRR151723', 'CRR151724', 'CRR151725',
            'CRR995220', 'CRR995221', 'CRR995222'
        ]
    
    logger.info(f"Processing {len(crr_samples)} CRR samples")
    
    # Generate metadata
    crr_metadata = processor.generate_crr_metadata(crr_samples)
    
    # Check SA_Results availability
    sa_results = processor.check_sa_results(crr_samples)
    for metadata in crr_metadata:
        metadata['has_sa_results'] = sa_results.get(metadata['Run'], False)
    
    # Log SA_Results status
    sa_count = sum(sa_results.values())
    logger.info(f"SA_Results available for {sa_count}/{len(crr_samples)} samples")
    
    # Integrate with existing metadata
    input_file = Path(args.input_metadata) if args.input_metadata else None
    if not input_file:
        # Try default locations
        for default_path in ['./Integrated_metadata_updated.csv', './Integrated_metadata.csv']:
            if Path(default_path).exists():
                input_file = Path(default_path)
                break
    
    combined_df = processor.integrate_with_existing(crr_metadata, input_file)
    
    # Add tissue information if missing
    if 'tissue' not in combined_df.columns or combined_df['tissue'].isna().any():
        # For CRR samples without tissue info, mark as Unknown
        mask = combined_df['Run'].str.startswith('CRR') & combined_df.get('tissue', pd.Series()).isna()
        if 'tissue' not in combined_df.columns:
            combined_df['tissue'] = 'Unknown'
        else:
            combined_df.loc[mask, 'tissue'] = 'Unknown'
    
    # Save results
    output_path = Path(args.output) if args.output else processor.output_dir / 'Integrated_metadata_final.csv'
    combined_df.to_csv(output_path, index=False)
    logger.info(f"Saved integrated metadata to: {output_path}")
    
    # Save in additional formats
    output_formats = config.get('output', {}).get('formats', ['csv'])
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    for fmt in output_formats:
        if fmt != 'csv':  # Already saved as CSV
            output_file = processor.output_dir / f'metadata_final_{timestamp}.{fmt}'
            if fmt == 'xlsx':
                combined_df.to_excel(output_file, index=False, engine='openpyxl')
            elif fmt == 'json':
                combined_df.to_json(output_file, orient='records', indent=2)
            logger.info(f"Saved metadata to: {output_file}")
    
    # Print summary
    logger.info("\nIntegration Summary:")
    logger.info(f"Total samples: {len(combined_df)}")
    
    # Sample type distribution
    sample_types = combined_df['Run'].str[:3].value_counts()
    logger.info("\nSample type distribution:")
    for prefix, count in sample_types.items():
        logger.info(f"  {prefix}xxx: {count}")
    
    # Samples with SA_Results
    if 'has_sa_results' in combined_df.columns:
        sa_samples = combined_df[combined_df['has_sa_results']]['Run'].tolist()
        logger.info(f"\nSamples with SA_Results: {len(sa_samples)}")
        if len(sa_samples) <= 30:  # Only show if reasonable number
            logger.info(f"SA_Results samples: {', '.join(sa_samples)}")
    
    logger.info("\n" + "="*70)
    logger.info("CRR PROCESSING COMPLETE")
    logger.info("="*70)

if __name__ == "__main__":
    main()