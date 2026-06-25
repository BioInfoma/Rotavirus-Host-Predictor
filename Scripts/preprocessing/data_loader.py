"""
Data loading module for Rotavirus VP4 preprocessing pipeline.
Loads training and evaluation datasets from metadata and FASTA files.
"""

import logging
from pathlib import Path
from typing import Dict, Tuple
import pandas as pd


class SimpleFASTAParser:
    """Simple FASTA parser without Biopython dependency."""
    
    @staticmethod
    def parse(fasta_path: str) -> Dict[str, str]:
        """
        Parse FASTA file.
        
        Args:
            fasta_path: Path to FASTA file
            
        Returns:
            Dictionary mapping sequence IDs to sequences
        """
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.rstrip('\r\n')  # Handle both Unix and Windows line endings
                
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id and current_seq:
                        seq_str = ''.join(current_seq)
                        sequences[current_id] = seq_str
                    
                    # Parse new ID (take first part before space/pipe)
                    header = line[1:]  # Remove '>'
                    current_id = header.split('|')[0].split()[0].strip()
                    current_seq = []
                elif line:  # Only add non-empty lines
                    if current_id:
                        current_seq.append(line.upper())
            
            # Save last sequence
            if current_id and current_seq:
                seq_str = ''.join(current_seq)
                sequences[current_id] = seq_str
        
        return sequences


class DataLoader:
    """Load training and evaluation datasets."""
    
    def __init__(self, config: Dict):
        """
        Initialize data loader.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.fasta_parser = SimpleFASTAParser()
    
    def load_training_data(self) -> Dict:
        """
        Load training dataset from metadata and FASTA files.
        
        Returns:
            Dictionary with 'metadata', 'nucleotide', and 'protein' keys
        """
        metadata_path = Path(self.config['training_metadata'])
        fasta_path = Path(self.config['training_fasta'])
        
        self.logger.info(f"Loading training metadata from: {metadata_path}")
        metadata = pd.read_excel(metadata_path)
        
        self.logger.info(f"Loading training sequences from: {fasta_path}")
        sequences = self.fasta_parser.parse(str(fasta_path))
        
        # Match sequences to metadata by accession
        matched_data = self._match_sequences_to_metadata(
            metadata, sequences, 'training'
        )
        
        self.logger.info(f"Successfully loaded {len(matched_data['metadata'])} training sequences")
        
        return matched_data
    
    def load_evaluation_data(self) -> Dict:
        """
        Load evaluation dataset from metadata and FASTA files.
        
        Returns:
            Dictionary with 'metadata', 'nucleotide', and 'protein' keys
        """
        metadata_path = Path(self.config['evaluation_metadata'])
        fasta_path = Path(self.config['evaluation_fasta'])
        
        self.logger.info(f"Loading evaluation metadata from: {metadata_path}")
        metadata = pd.read_excel(metadata_path)
        
        self.logger.info(f"Loading evaluation sequences from: {fasta_path}")
        sequences = self.fasta_parser.parse(str(fasta_path))
        
        # Match sequences to metadata by accession
        matched_data = self._match_sequences_to_metadata(
            metadata, sequences, 'evaluation'
        )
        
        self.logger.info(f"Successfully loaded {len(matched_data['metadata'])} evaluation sequences")
        
        return matched_data
    
    def _match_sequences_to_metadata(
        self,
        metadata: pd.DataFrame,
        sequences: Dict[str, str],
        dataset_type: str
    ) -> Dict:
        """
        Match loaded sequences to metadata records.
        
        Args:
            metadata: DataFrame with metadata
            sequences: Dictionary of sequences
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Dictionary with matched metadata and sequences
        """
        # Ensure metadata has required columns
        required_cols = ['Accession', 'Full_name', 'Host', 'Genotype', 'label', 
                        'adaptation_group', 'Length', 'Country', 'Collection_Date']
        
        missing_cols = [col for col in required_cols if col not in metadata.columns]
        if missing_cols:
            self.logger.warning(
                f"Missing columns in {dataset_type} metadata: {missing_cols}"
            )
        
        # Standardize column names to lowercase
        metadata.columns = [col.lower() for col in metadata.columns]
        
        # Filter to only records where sequences exist
        initial_count = len(metadata)
        metadata = metadata[metadata['accession'].isin(sequences.keys())].copy()
        final_count = len(metadata)
        
        if initial_count != final_count:
            self.logger.warning(
                f"Dropped {initial_count - final_count} {dataset_type} records "
                f"without matching sequences"
            )
        
        # Initialize sequence dictionaries
        nucleotide_seqs = {}
        protein_seqs = {}
        
        # Extract sequences for matched records
        for accession in metadata['accession']:
            if accession in sequences:
                nucleotide_seqs[accession] = sequences[accession]
                # Protein will be computed later
                protein_seqs[accession] = None
        
        # Add metadata columns for tracking
        if 'dataset_split' not in metadata.columns:
            metadata['dataset_split'] = dataset_type
        
        if 'sequence_md5' not in metadata.columns:
            metadata['sequence_md5'] = None
        
        if 'protein_md5' not in metadata.columns:
            metadata['protein_md5'] = None
        
        if 'selected_frame' not in metadata.columns:
            metadata['selected_frame'] = None
        
        if 'protein_length' not in metadata.columns:
            metadata['protein_length'] = None
        
        if 'vp8_reference' not in metadata.columns:
            metadata['vp8_reference'] = None
        
        if 'vp8_coverage' not in metadata.columns:
            metadata['vp8_coverage'] = None
        
        if 'vp8_status' not in metadata.columns:
            metadata['vp8_status'] = None
        
        # Reset index
        metadata = metadata.reset_index(drop=True)
        
        return {
            'metadata': metadata,
            'nucleotide': nucleotide_seqs,
            'protein': protein_seqs,
            'vp8': {},  # Will be populated during VP8 extraction
        }
