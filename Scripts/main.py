#!/usr/bin/env python3
"""
Rotavirus A VP4 Host Adaptation / Spillover Prediction - Preprocessing Pipeline
Publication-quality preprocessing for machine learning model training and evaluation.

This pipeline:
1. Loads and validates training and evaluation datasets
2. Performs sequence quality control (length, ambiguity, ORF detection)
3. Detects and extracts VP8* regions through alignment
4. Removes exact duplicates and clusters similar sequences
5. Normalizes metadata
6. Detects and reports host-genotype conflicts
7. Generates comprehensive QC reports
8. Produces clean, analysis-ready outputs

Author: Automated ML Pipeline
Date: 2024
"""

import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
from datetime import datetime

# Import pipeline components
from config.settings import (
    CONFIG,
    setup_logging,
    validate_config
)
from preprocessing.data_loader import DataLoader
from preprocessing.sequence_validator import SequenceValidator
from preprocessing.vp8_extractor import VP8Extractor
from preprocessing.deduplicator import Deduplicator
from preprocessing.metadata_normalizer import MetadataNormalizer
from preprocessing.conflict_detector import ConflictDetector
from preprocessing.report_generator import ReportGenerator


class RotavirusPipeline:
    """Main preprocessing pipeline orchestrator."""
    
    def __init__(self, config: Dict):
        """
        Initialize the pipeline with configuration.
        
        Args:
            config: Configuration dictionary from settings.py
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.logger.info("="*80)
        self.logger.info("ROTAVIRUS VP4 PREPROCESSING PIPELINE")
        self.logger.info("="*80)
        self.logger.info(f"Pipeline started at {datetime.now().isoformat()}")
        
        # Initialize directories
        self._setup_directories()
        
        # Initialize components
        self.data_loader = DataLoader(config)
        self.sequence_validator = SequenceValidator(config)
        self.vp8_extractor = VP8Extractor(config)
        self.deduplicator = Deduplicator(config)
        self.metadata_normalizer = MetadataNormalizer(config)
        self.conflict_detector = ConflictDetector(config)
        self.report_generator = ReportGenerator(config)
        
        # Track statistics
        self.stats = {
            'training': {},
            'evaluation': {}
        }
    
    def _setup_directories(self) -> None:
        """Create necessary output directories."""
        dirs_to_create = [
            self.config['output_dir'],
            self.config['references_dir'],
            self.config['analysis_dir'],
            self.config['reports_dir'],
            self.config['composition_dir'],
            self.config['logs_dir'],
        ]
        
        for dir_path in dirs_to_create:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Directory ready: {dir_path}")
    
    def run(self) -> bool:
        """
        Execute the complete preprocessing pipeline.
        
        Returns:
            bool: True if pipeline completed successfully
        """
        try:
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 1: LOADING DATA")
            self.logger.info("="*80)
            training_data, evaluation_data = self._load_data()
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 2: SEQUENCE VALIDATION AND FILTERING")
            self.logger.info("="*80)
            training_data, evaluation_data = self._validate_sequences(
                training_data, evaluation_data
            )
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 3: VP8* EXTRACTION AND MAPPING")
            self.logger.info("="*80)
            training_data, evaluation_data = self._extract_vp8(
                training_data, evaluation_data
            )
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 4: EXACT DEDUPLICATION")
            self.logger.info("="*80)
            training_data = self._deduplicate(training_data, 'training')
            evaluation_data = self._deduplicate(evaluation_data, 'evaluation')
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 5: REDUNDANCY CLUSTERING")
            self.logger.info("="*80)
            training_data = self._cluster_redundancy(training_data, 'training')
            evaluation_data = self._cluster_redundancy(evaluation_data, 'evaluation')
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 6: METADATA NORMALIZATION")
            self.logger.info("="*80)
            training_data = self._normalize_metadata(training_data, 'training')
            evaluation_data = self._normalize_metadata(evaluation_data, 'evaluation')
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 7: CONFLICT DETECTION")
            self.logger.info("="*80)
            self._detect_conflicts(training_data, 'training')
            self._detect_conflicts(evaluation_data, 'evaluation')
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 8: SAVING OUTPUTS")
            self.logger.info("="*80)
            self._save_outputs(training_data, evaluation_data)
            
            self.logger.info("\n" + "="*80)
            self.logger.info("STEP 9: GENERATING REPORTS")
            self.logger.info("="*80)
            self._generate_reports(training_data, evaluation_data)
            
            self.logger.info("\n" + "="*80)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info(f"Completed at {datetime.now().isoformat()}")
            self.logger.info("="*80 + "\n")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Pipeline failed with error: {str(e)}", exc_info=True)
            return False
    
    def _load_data(self) -> Tuple[Dict, Dict]:
        """Load training and evaluation data."""
        self.logger.info("Loading training data...")
        training_data = self.data_loader.load_training_data()
        self.logger.info(f"Loaded {len(training_data['metadata'])} training sequences")
        if len(training_data['metadata']) > 0:
            hosts = training_data['metadata']['host'].unique() if 'host' in training_data['metadata'].columns else []
            self.logger.info(f"Training hosts: {hosts}")
        
        self.logger.info("\nLoading evaluation data...")
        evaluation_data = self.data_loader.load_evaluation_data()
        self.logger.info(f"Loaded {len(evaluation_data['metadata'])} evaluation sequences")
        if len(evaluation_data['metadata']) > 0:
            hosts = evaluation_data['metadata']['host'].unique() if 'host' in evaluation_data['metadata'].columns else []
            self.logger.info(f"Evaluation hosts: {hosts}")
        
        return training_data, evaluation_data
    
    def _validate_sequences(
        self, 
        training_data: Dict, 
        evaluation_data: Dict
    ) -> Tuple[Dict, Dict]:
        """Perform sequence quality control."""
        self.logger.info("\nValidating training sequences...")
        training_data = self.sequence_validator.validate(training_data, 'training')
        self.stats['training']['after_validation'] = len(training_data['metadata'])
        self.logger.info(
            f"After validation: {self.stats['training']['after_validation']} sequences"
        )
        
        self.logger.info("\nValidating evaluation sequences...")
        evaluation_data = self.sequence_validator.validate(evaluation_data, 'evaluation')
        self.stats['evaluation']['after_validation'] = len(evaluation_data['metadata'])
        self.logger.info(
            f"After validation: {self.stats['evaluation']['after_validation']} sequences"
        )
        
        return training_data, evaluation_data
    
    def _extract_vp8(
        self,
        training_data: Dict,
        evaluation_data: Dict
    ) -> Tuple[Dict, Dict]:
        """Extract VP8* regions."""
        self.logger.info("\nExtracting VP8* from training sequences...")
        training_data = self.vp8_extractor.extract(training_data, 'training')
        self.stats['training']['after_vp8_extraction'] = len(training_data['metadata'])
        self.logger.info(
            f"After VP8* extraction: {self.stats['training']['after_vp8_extraction']} sequences"
        )
        
        self.logger.info("\nExtracting VP8* from evaluation sequences...")
        evaluation_data = self.vp8_extractor.extract(evaluation_data, 'evaluation')
        self.stats['evaluation']['after_vp8_extraction'] = len(evaluation_data['metadata'])
        self.logger.info(
            f"After VP8* extraction: {self.stats['evaluation']['after_vp8_extraction']} sequences"
        )
        
        return training_data, evaluation_data
    
    def _deduplicate(self, data: Dict, dataset_type: str) -> Dict:
        """Remove exact duplicates."""
        initial_count = len(data['metadata'])
        data = self.deduplicator.deduplicate(data, dataset_type)
        final_count = len(data['metadata'])
        self.logger.info(
            f"After deduplication: {final_count} sequences "
            f"(removed {initial_count - final_count})"
        )
        self.stats[dataset_type]['after_deduplication'] = final_count
        return data
    
    def _cluster_redundancy(self, data: Dict, dataset_type: str) -> Dict:
        """Cluster similar sequences to prevent memorization."""
        initial_count = len(data['metadata'])
        data = self.deduplicator.cluster_redundancy(data, dataset_type)
        final_count = len(data['metadata'])
        self.logger.info(
            f"After redundancy clustering: {final_count} sequences "
            f"(representative sequences)"
        )
        self.stats[dataset_type]['after_clustering'] = final_count
        return data
    
    def _normalize_metadata(self, data: Dict, dataset_type: str) -> Dict:
        """Normalize metadata fields."""
        self.logger.info(f"Normalizing {dataset_type} metadata...")
        data = self.metadata_normalizer.normalize(data, dataset_type)
        self.logger.info("Metadata normalization complete")
        return data
    
    def _detect_conflicts(self, data: Dict, dataset_type: str) -> None:
        """Detect and report host-genotype conflicts."""
        conflicts = self.conflict_detector.detect(data, dataset_type)
        if len(conflicts) > 0:
            self.logger.warning(
                f"Found {len(conflicts)} host-genotype conflicts in {dataset_type} data"
            )
        else:
            self.logger.info(f"No host-genotype conflicts found in {dataset_type} data")
    
    def _save_outputs(self, training_data: Dict, evaluation_data: Dict) -> None:
        """Save cleaned outputs to disk."""
        self.logger.info("Saving training outputs...")
        self._save_dataset_outputs(training_data, 'training')
        
        self.logger.info("Saving evaluation outputs...")
        self._save_dataset_outputs(evaluation_data, 'evaluation')
        
        self.logger.info("All outputs saved successfully")
    
    def _save_dataset_outputs(self, data: Dict, dataset_type: str) -> None:
        """Save a single dataset's outputs."""
        # Metadata
        metadata_path = Path(self.config['analysis_dir']) / f"cleaned_metadata_{dataset_type}.csv"
        data['metadata'].to_csv(metadata_path, index=False)
        self.logger.info(f"Saved metadata: {metadata_path}")
        
        # Nucleotide sequences
        nuc_path = Path(self.config['analysis_dir']) / f"cleaned_nucleotide_{dataset_type}.fasta"
        self._write_fasta(nuc_path, data['nucleotide'])
        self.logger.info(f"Saved nucleotide sequences: {nuc_path}")
        
        # Protein sequences
        prot_path = Path(self.config['analysis_dir']) / f"cleaned_protein_{dataset_type}.fasta"
        self._write_fasta(prot_path, data['protein'])
        self.logger.info(f"Saved protein sequences: {prot_path}")
        
        # VP8* sequences
        vp8_path = Path(self.config['analysis_dir']) / f"cleaned_vp8_{dataset_type}.fasta"
        self._write_fasta(vp8_path, data['vp8'])
        self.logger.info(f"Saved VP8* sequences: {vp8_path}")
    
    def _write_fasta(self, path: Path, sequences: Dict[str, str]) -> None:
        """Write sequences to FASTA file."""
        with open(path, 'w', encoding='utf-8') as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n")
                # Wrap at 80 characters
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
    
    def _generate_reports(self, training_data: Dict, evaluation_data: Dict) -> None:
        """Generate comprehensive QC and composition reports."""
        self.logger.info("Generating QC reports...")
        self.report_generator.generate_all_reports(
            training_data,
            evaluation_data,
            self.stats
        )
        self.logger.info("Report generation complete")


def main():
    """Main entry point."""
    # Setup logging
    setup_logging()
    logger = logging.getLogger(__name__)
    
    try:
        # Validate configuration
        config = CONFIG
        validate_config(config)
        
        # Run pipeline
        pipeline = RotavirusPipeline(config)
        success = pipeline.run()
        
        if not success:
            sys.exit(1)
        
        logger.info("\n" + "="*80)
        logger.info("All preprocessing steps completed successfully!")
        logger.info("Analysis-ready outputs are in: " + config['analysis_dir'])
        logger.info("QC reports are in: " + config['reports_dir'])
        logger.info("="*80 + "\n")
        
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
