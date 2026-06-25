"""
Sequence validation module for quality control filtering.
Performs:
- Sequence length filtering
- Ambiguity (N bases) filtering
- ORF detection and translation
"""

import logging
from typing import Dict, Tuple, Optional
from pathlib import Path
import pandas as pd
import hashlib


class SequenceValidator:
    """Validate sequence quality and translate to protein."""
    
    # Genetic code (standard)
    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    def __init__(self, config: Dict):
        """
        Initialize sequence validator.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.length_min = config['sequence_length_min']
        self.length_max = config['sequence_length_max']
        self.ambiguity_max = config['ambiguity_max_percent'] / 100.0
        self.orf_min_length = config['orf_min_length']
        self.stop_codons = config['stop_codons']
    
    def validate(self, data: Dict, dataset_type: str) -> Dict:
        """
        Validate sequences and remove low-quality entries.
        
        Args:
            data: Dictionary with metadata, nucleotide, and protein sequences
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Filtered data dictionary
        """
        metadata = data['metadata'].copy()
        nucleotide = data['nucleotide'].copy()
        protein = data['protein'].copy()
        
        initial_count = len(metadata)
        self.logger.info(f"Starting with {initial_count} {dataset_type} sequences")
        
        # Track reasons for filtering
        filters_applied = {
            'total_filtered': 0,
            'length_too_short': 0,
            'length_too_long': 0,
            'high_ambiguity': 0,
            'orf_failed': 0,
            'internal_stops': 0,
        }
        
        sequences_to_keep = []
        
        for idx, row in metadata.iterrows():
            accession = row['accession']
            nuc_seq = nucleotide[accession]
            
            # Length check
            if len(nuc_seq) < self.length_min:
                filters_applied['length_too_short'] += 1
                filters_applied['total_filtered'] += 1
                continue
            
            if len(nuc_seq) > self.length_max:
                filters_applied['length_too_long'] += 1
                filters_applied['total_filtered'] += 1
                continue
            
            # Ambiguity check
            n_count = nuc_seq.count('N') + nuc_seq.count('n')
            ambiguity = n_count / len(nuc_seq)
            
            if ambiguity > self.ambiguity_max:
                filters_applied['high_ambiguity'] += 1
                filters_applied['total_filtered'] += 1
                continue
            
            # ORF detection and translation
            translation_result = self._detect_and_translate_orf(nuc_seq)
            
            if translation_result is None:
                filters_applied['orf_failed'] += 1
                filters_applied['total_filtered'] += 1
                continue
            
            prot_seq, frame_number = translation_result
            
            # Check for internal stop codons (except terminal)
            if len(prot_seq) > 1 and '*' in prot_seq[:-1]:
                filters_applied['internal_stops'] += 1
                filters_applied['total_filtered'] += 1
                continue
            
            # Sequence passes validation
            sequences_to_keep.append(accession)
            protein[accession] = prot_seq
            
            # Store frame info
            metadata.loc[idx, 'selected_frame'] = frame_number
            metadata.loc[idx, 'protein_length'] = len(prot_seq)
            metadata.loc[idx, 'sequence_md5'] = self._calculate_md5(nuc_seq)
            metadata.loc[idx, 'protein_md5'] = self._calculate_md5(prot_seq)
        
        # Filter metadata and sequences
        metadata = metadata[metadata['accession'].isin(sequences_to_keep)].copy()
        nucleotide = {k: v for k, v in nucleotide.items() if k in sequences_to_keep}
        protein = {k: v for k, v in protein.items() if k in sequences_to_keep}
        
        # Log filtering results
        self.logger.info(f"After length filtering: {len(metadata)} sequences")
        self.logger.info(f"  - Too short (<{self.length_min} bp): {filters_applied['length_too_short']}")
        self.logger.info(f"  - Too long (>{self.length_max} bp): {filters_applied['length_too_long']}")
        self.logger.info(f"  - High ambiguity (>{self.ambiguity_max*100:.1f}%): {filters_applied['high_ambiguity']}")
        self.logger.info(f"  - ORF detection failed: {filters_applied['orf_failed']}")
        self.logger.info(f"  - Internal stop codons: {filters_applied['internal_stops']}")
        self.logger.info(f"Total removed: {filters_applied['total_filtered']}")
        
        # Save QC report
        self._save_validation_report(filters_applied, dataset_type)
        
        metadata = metadata.reset_index(drop=True)
        
        return {
            'metadata': metadata,
            'nucleotide': nucleotide,
            'protein': protein,
            'vp8': data['vp8'],
        }
    
    def _detect_and_translate_orf(self, nuc_seq: str) -> Optional[Tuple[str, int]]:
        """
        Detect the longest ORF and translate it.
        
        Args:
            nuc_seq: Nucleotide sequence
            
        Returns:
            Tuple of (protein_sequence, frame_number) or None if no valid ORF
        """
        longest_prot = None
        best_frame = None
        longest_length = 0
        
        # Try all three reading frames
        for frame in range(3):
            prot = self._translate_frame(nuc_seq[frame:], frame)
            
            if prot and len(prot) > longest_length:
                longest_prot = prot
                best_frame = frame
                longest_length = len(prot)
        
        if longest_length >= self.orf_min_length:
            return (longest_prot, best_frame)
        
        return None
    
    def _translate_frame(self, nuc_seq: str, frame: int) -> Optional[str]:
        """
        Translate a nucleotide sequence in a specific frame.
        Extracts the longest sequence between stop codons.
        
        Args:
            nuc_seq: Nucleotide sequence
            frame: Reading frame (0, 1, or 2)
            
        Returns:
            Protein sequence or None if no valid ORF
        """
        # Pad sequence to multiple of 3
        remainder = len(nuc_seq) % 3
        if remainder:
            nuc_seq = nuc_seq[:-remainder]
        
        # Translate to codons
        codons = [nuc_seq[i:i+3] for i in range(0, len(nuc_seq), 3)]
        
        # Translate to protein
        prot_parts = []
        current_part = []
        
        for codon in codons:
            if len(codon) == 3:
                # Handle degenerate codons (N's)
                if 'N' in codon:
                    aa = 'X'  # Unknown amino acid
                else:
                    aa = self.CODON_TABLE.get(codon, 'X')
                
                if aa == '*':  # Stop codon
                    if current_part:
                        prot_parts.append(''.join(current_part))
                    current_part = []
                else:
                    current_part.append(aa)
        
        # Don't forget the last part
        if current_part:
            prot_parts.append(''.join(current_part))
        
        # Return longest ORF
        if prot_parts:
            return max(prot_parts, key=len)
        
        return None
    
    def _calculate_md5(self, sequence: str) -> str:
        """Calculate MD5 hash of sequence."""
        return hashlib.md5(sequence.upper().encode()).hexdigest()
    
    def _save_validation_report(
        self,
        filters_applied: Dict,
        dataset_type: str
    ) -> None:
        """Save validation report to CSV."""
        report_dir = Path(self.config['reports_dir'])
        report_path = report_dir / f"sequence_validation_{dataset_type}.csv"
        
        report_df = pd.DataFrame([{
            'filter': k,
            'count': v
        } for k, v in filters_applied.items()])
        
        report_df.to_csv(report_path, index=False)
        self.logger.info(f"Saved validation report: {report_path}")
