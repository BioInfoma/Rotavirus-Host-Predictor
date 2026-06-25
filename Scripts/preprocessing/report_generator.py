"""
Report generation module.
Generates comprehensive QC and composition reports.
"""

import logging
from typing import Dict
from pathlib import Path
import pandas as pd


class ReportGenerator:
    """Generate QC and composition reports."""
    
    def __init__(self, config: Dict):
        """
        Initialize report generator.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
    
    def generate_all_reports(
        self,
        training_data: Dict,
        evaluation_data: Dict,
        stats: Dict
    ) -> None:
        """
        Generate all QC and composition reports.
        
        Args:
            training_data: Training dataset dictionary
            evaluation_data: Evaluation dataset dictionary
            stats: Pipeline statistics
        """
        self.logger.info("Generating QC Summary Report...")
        self._generate_qc_summary(stats)
        
        self.logger.info("Generating Integrity Reports...")
        self._generate_integrity_report(training_data, 'training')
        self._generate_integrity_report(evaluation_data, 'evaluation')
        
        self.logger.info("Generating Composition Reports...")
        self._generate_host_distribution(training_data, 'training')
        self._generate_host_distribution(evaluation_data, 'evaluation')
        
        self._generate_genotype_distribution(training_data, 'training')
        self._generate_genotype_distribution(evaluation_data, 'evaluation')
        
        self._generate_label_distribution(training_data, 'training')
        self._generate_label_distribution(evaluation_data, 'evaluation')
        
        self._generate_adaptation_distribution(training_data, 'training')
        self._generate_adaptation_distribution(evaluation_data, 'evaluation')
        
        self.logger.info("Report generation complete!")
    
    def _generate_qc_summary(self, stats: Dict) -> None:
        """Generate QC summary report."""
        report_path = Path(self.config['reports_dir']) / "qc_summary.csv"
        
        # Create summary rows
        summary_rows = []
        
        for dataset_type in ['training', 'evaluation']:
            if dataset_type in stats:
                dataset_stats = stats[dataset_type]
                summary_rows.append({
                    'dataset': dataset_type,
                    'after_validation': dataset_stats.get('after_validation', 0),
                    'after_vp8_extraction': dataset_stats.get('after_vp8_extraction', 0),
                    'after_deduplication': dataset_stats.get('after_deduplication', 0),
                    'after_clustering': dataset_stats.get('after_clustering', 0),
                })
        
        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_csv(report_path, index=False)
            self.logger.info(f"Saved QC summary: {report_path}")
    
    def _generate_integrity_report(self, data: Dict, dataset_type: str) -> None:
        """
        Generate data integrity report.
        
        Args:
            data: Dataset dictionary
            dataset_type: 'training' or 'evaluation'
        """
        report_path = Path(self.config['reports_dir']) / f"integrity_summary_{dataset_type}.txt"
        
        metadata = data['metadata']
        nucleotide = data['nucleotide']
        protein = data['protein']
        vp8_seqs = data['vp8']
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(f"Data Integrity Report - {dataset_type.upper()}\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Total sequences: {len(metadata)}\n")
            f.write(f"Nucleotide sequences: {len(nucleotide)}\n")
            f.write(f"Protein sequences: {len(protein)}\n")
            f.write(f"VP8* sequences: {len(vp8_seqs)}\n\n")
            
            # Check completeness
            nucleotide_complete = len(nucleotide) == len(metadata)
            protein_complete = len(protein) == len(metadata)
            vp8_complete = len(vp8_seqs) == len(metadata)
            
            f.write("Completeness Check:\n")
            f.write(f"  Nucleotide: {'PASS' if nucleotide_complete else 'FAIL'}\n")
            f.write(f"  Protein: {'PASS' if protein_complete else 'FAIL'}\n")
            f.write(f"  VP8*: {'PASS' if vp8_complete else 'FAIL'}\n\n")
            
            # Metadata completeness
            f.write("Metadata Completeness:\n")
            required_cols = ['accession', 'full_name', 'host', 'genotype', 'label',
                           'adaptation_group', 'country', 'collection_year']
            for col in required_cols:
                missing = metadata[col].isna().sum()
                f.write(f"  {col}: {len(metadata) - missing}/{len(metadata)} complete\n")
            
            # VP8* coverage statistics
            f.write("\nVP8* Coverage Statistics:\n")
            if len(metadata) > 0:
                coverage = metadata['vp8_coverage'].dropna()
                f.write(f"  Mean: {coverage.mean():.2%}\n")
                f.write(f"  Min: {coverage.min():.2%}\n")
                f.write(f"  Max: {coverage.max():.2%}\n")
            
            f.write(f"\nReport generated at {pd.Timestamp.now()}\n")
        
        self.logger.info(f"Saved integrity report: {report_path}")
    
    def _generate_host_distribution(self, data: Dict, dataset_type: str) -> None:
        """Generate host distribution report."""
        report_path = Path(self.config['composition_dir']) / f"host_distribution_{dataset_type}.csv"
        
        metadata = data['metadata']
        host_dist = metadata['host'].value_counts().reset_index()
        host_dist.columns = ['host', 'count']
        host_dist['percentage'] = (host_dist['count'] / host_dist['count'].sum() * 100).round(2)
        
        host_dist.to_csv(report_path, index=False)
        self.logger.info(f"Saved host distribution: {report_path}")
    
    def _generate_genotype_distribution(self, data: Dict, dataset_type: str) -> None:
        """Generate genotype distribution report."""
        report_path = Path(self.config['composition_dir']) / f"genotype_distribution_{dataset_type}.csv"
        
        metadata = data['metadata']
        geno_dist = metadata['genotype'].value_counts().reset_index()
        geno_dist.columns = ['genotype', 'count']
        geno_dist['percentage'] = (geno_dist['count'] / geno_dist['count'].sum() * 100).round(2)
        
        geno_dist.to_csv(report_path, index=False)
        self.logger.info(f"Saved genotype distribution: {report_path}")
    
    def _generate_label_distribution(self, data: Dict, dataset_type: str) -> None:
        """Generate label distribution report."""
        report_path = Path(self.config['composition_dir']) / f"label_distribution_{dataset_type}.csv"
        
        metadata = data['metadata']
        if 'label' in metadata.columns:
            label_dist = metadata['label'].value_counts().reset_index()
            label_dist.columns = ['label', 'count']
            label_dist['percentage'] = (label_dist['count'] / label_dist['count'].sum() * 100).round(2)
            
            label_dist.to_csv(report_path, index=False)
            self.logger.info(f"Saved label distribution: {report_path}")
    
    def _generate_adaptation_distribution(self, data: Dict, dataset_type: str) -> None:
        """Generate adaptation group distribution report."""
        report_path = Path(self.config['composition_dir']) / f"adaptation_distribution_{dataset_type}.csv"
        
        metadata = data['metadata']
        if 'adaptation_group' in metadata.columns:
            adapt_dist = metadata['adaptation_group'].value_counts().reset_index()
            adapt_dist.columns = ['adaptation_group', 'count']
            adapt_dist['percentage'] = (adapt_dist['count'] / adapt_dist['count'].sum() * 100).round(2)
            
            adapt_dist.to_csv(report_path, index=False)
            self.logger.info(f"Saved adaptation distribution: {report_path}")
