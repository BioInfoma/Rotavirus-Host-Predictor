"""
Conflict detection module.
Detects and reports host-genotype combinations that may be biologically ambiguous.
"""

import logging
from typing import Dict, List
from pathlib import Path
import pandas as pd


class ConflictDetector:
    """Detect host-genotype conflicts."""
    
    def __init__(self, config: Dict):
        """
        Initialize conflict detector.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.ambiguous_combinations = config['ambiguous_combinations']
    
    def detect(self, data: Dict, dataset_type: str) -> List[Dict]:
        """
        Detect host-genotype conflicts.
        
        Args:
            data: Dictionary with metadata
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            List of conflict records
        """
        metadata = data['metadata']
        conflicts = []
        
        for idx, row in metadata.iterrows():
            host = row['host']
            genotype = row['genotype']
            accession = row['accession']
            
            # Check if this is an ambiguous combination
            for ambig_host, ambig_geno in self.ambiguous_combinations:
                # Fuzzy matching for genotypes
                if host == ambig_host and self._genotype_matches(genotype, ambig_geno):
                    conflicts.append({
                        'accession': accession,
                        'full_name': row['full_name'],
                        'host': host,
                        'genotype': genotype,
                        'label': row['label'],
                        'adaptation_group': row['adaptation_group'],
                        'conflict_type': f'Ambiguous_{ambig_host}_{ambig_geno}',
                        'notes': self._get_conflict_notes(ambig_host, ambig_geno)
                    })
        
        # Save conflict report
        self._save_conflict_report(conflicts, dataset_type)
        
        return conflicts
    
    def _genotype_matches(self, genotype_str: str, pattern: str) -> bool:
        """
        Check if genotype matches pattern.
        
        Args:
            genotype_str: Genotype string (e.g., 'P[8]')
            pattern: Pattern to match (e.g., 'P[9]')
            
        Returns:
            True if matches
        """
        if pd.isna(genotype_str) or pd.isna(pattern):
            return False
        
        genotype_str = str(genotype_str).strip()
        pattern = str(pattern).strip()
        
        # Exact match
        if genotype_str == pattern:
            return True
        
        # Check without brackets
        gen_clean = genotype_str.replace('[', '').replace(']', '')
        pat_clean = pattern.replace('[', '').replace(']', '')
        
        return gen_clean == pat_clean
    
    def _get_conflict_notes(self, host: str, genotype: str) -> str:
        """
        Get explanatory notes for conflict type.
        
        Args:
            host: Host species
            genotype: Genotype
            
        Returns:
            Explanatory note
        """
        notes = {
            ('Human', 'P[9]'): 'Human P[9] - often zoonotic from animal sources',
            ('Human', 'P[14]'): 'Human P[14] - zoonotic, potentially animal-derived',
            ('Human', 'P[19]'): 'Human P[19] - zoonotic origin uncertain',
            ('Human', 'P[25]'): 'Human P[25] - zoonotic origin uncertain',
            ('Porcine', 'P[6]'): 'Porcine P[6] - P[6] genotype not typical for swine',
            ('Bat', 'P[6]'): 'Bat P[6] - P[6] genotype not typical for bats',
        }
        
        return notes.get((host, genotype), 'Ambiguous host-genotype combination')
    
    def _save_conflict_report(self, conflicts: List[Dict], dataset_type: str) -> None:
        """
        Save conflict report to CSV.
        
        Args:
            conflicts: List of conflict records
            dataset_type: 'training' or 'evaluation'
        """
        report_path = Path(self.config['reports_dir']) / f"host_genotype_conflicts_{dataset_type}.csv"
        
        if conflicts:
            conflicts_df = pd.DataFrame(conflicts)
            conflicts_df.to_csv(report_path, index=False)
            self.logger.info(f"Saved conflicts report ({len(conflicts)} conflicts): {report_path}")
        else:
            # Create empty report file
            pd.DataFrame(columns=['accession', 'full_name', 'host', 'genotype',
                                 'label', 'adaptation_group', 'conflict_type', 'notes'
            ]).to_csv(report_path, index=False)
            self.logger.info(f"Saved conflicts report (0 conflicts): {report_path}")
