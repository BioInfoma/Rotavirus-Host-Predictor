"""
Deduplication and redundancy clustering module.
Removes exact duplicates and clusters similar sequences to prevent memorization.
"""

import logging
from typing import Dict, List, Set, Optional
from pathlib import Path
import pandas as pd
import hashlib
import subprocess
import tempfile


class Deduplicator:
    """Handle sequence deduplication and redundancy clustering."""
    
    def __init__(self, config: Dict):
        """
        Initialize deduplicator.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.cdhit_identity = config['cdhit_identity']
        self.cdhit_wordsize = config['cdhit_wordsize']
    
    def deduplicate(self, data: Dict, dataset_type: str) -> Dict:
        """
        Remove exact duplicate sequences by MD5 hash.
        
        Args:
            data: Dictionary with metadata and sequences
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Filtered data with duplicates removed
        """
        metadata = data['metadata'].copy()
        nucleotide = data['nucleotide'].copy()
        protein = data['protein'].copy()
        vp8_seqs = data['vp8'].copy()
        
        initial_count = len(metadata)
        
        # Calculate MD5 hashes if not already done
        for idx, row in metadata.iterrows():
            accession = row['accession']
            if pd.isna(row['sequence_md5']) or not row['sequence_md5']:
                nuc_seq = nucleotide[accession]
                metadata.loc[idx, 'sequence_md5'] = self._calculate_md5(nuc_seq)
        
        # Find duplicates
        duplicates_df = []
        seen_hashes = {}
        sequences_to_keep = []
        
        for idx, row in metadata.iterrows():
            accession = row['accession']
            seq_hash = row['sequence_md5']
            
            if seq_hash in seen_hashes:
                # This is a duplicate
                original_acc = seen_hashes[seq_hash]
                duplicates_df.append({
                    'original': original_acc,
                    'duplicate': accession,
                    'sequence_md5': seq_hash,
                })
            else:
                # First occurrence, keep it
                seen_hashes[seq_hash] = accession
                sequences_to_keep.append(accession)
        
        # Filter data
        metadata = metadata[metadata['accession'].isin(sequences_to_keep)].copy()
        nucleotide = {k: v for k, v in nucleotide.items() if k in sequences_to_keep}
        protein = {k: v for k, v in protein.items() if k in sequences_to_keep}
        vp8_seqs = {k: v for k, v in vp8_seqs.items() if k in sequences_to_keep}
        
        # Log results
        duplicates_removed = initial_count - len(metadata)
        self.logger.info(f"Exact deduplication: removed {duplicates_removed} duplicates")
        
        # Save duplicate report
        if duplicates_df:
            dup_report = pd.DataFrame(duplicates_df)
            report_path = Path(self.config['reports_dir']) / f"duplicate_report_{dataset_type}.csv"
            dup_report.to_csv(report_path, index=False)
            self.logger.info(f"Saved duplicate report: {report_path}")
        
        metadata = metadata.reset_index(drop=True)
        
        return {
            'metadata': metadata,
            'nucleotide': nucleotide,
            'protein': protein,
            'vp8': vp8_seqs,
        }
    
    def cluster_redundancy(self, data: Dict, dataset_type: str) -> Dict:
        """
        Cluster similar sequences using CD-HIT-like approach.
        Keep representative sequence from each cluster.
        
        Args:
            data: Dictionary with metadata and sequences
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Data with cluster representatives only
        """
        metadata = data['metadata'].copy()
        nucleotide = data['nucleotide'].copy()
        protein = data['protein'].copy()
        vp8_seqs = data['vp8'].copy()
        
        initial_count = len(metadata)
        self.logger.info(f"Starting clustering with {initial_count} sequences")
        
        # Try to use CD-HIT if available
        cluster_result = self._cluster_with_cdhit(nucleotide, dataset_type)
        
        if cluster_result is not None:
            representatives = cluster_result
        else:
            # Fallback to simple clustering
            self.logger.info("Using fallback clustering method")
            representatives = self._simple_clustering(nucleotide)
        
        # Filter to representatives only
        metadata = metadata[metadata['accession'].isin(representatives)].copy()
        nucleotide = {k: v for k, v in nucleotide.items() if k in representatives}
        protein = {k: v for k, v in protein.items() if k in representatives}
        vp8_seqs = {k: v for k, v in vp8_seqs.items() if k in representatives}
        
        # Log results
        clusters_removed = initial_count - len(metadata)
        self.logger.info(
            f"Redundancy clustering: {clusters_removed} clustered sequences removed, "
            f"{len(metadata)} representatives kept"
        )
        
        # Save clustering report
        self._save_clustering_report(metadata, dataset_type)
        
        metadata = metadata.reset_index(drop=True)
        
        return {
            'metadata': metadata,
            'nucleotide': nucleotide,
            'protein': protein,
            'vp8': vp8_seqs,
        }
    
    def _cluster_with_cdhit(
        self,
        nucleotide: Dict[str, str],
        dataset_type: str
    ) -> Optional[Set[str]]:
        """
        Cluster sequences using CD-HIT.
        
        Args:
            nucleotide: Dictionary of nucleotide sequences
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Set of representative sequence accessions or None if CD-HIT failed
        """
        try:
            # Create temporary FASTA file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                temp_fasta = f.name
                for accession, seq in nucleotide.items():
                    f.write(f">{accession}\n")
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\n")
            
            # Run CD-HIT
            temp_dir = tempfile.gettempdir()
            output_prefix = Path(temp_dir) / f"cdhit_{dataset_type}"
            
            cmd = [
                'cd-hit',
                '-i', temp_fasta,
                '-o', str(output_prefix),
                '-c', str(self.cdhit_identity),
                '-n', str(self.cdhit_wordsize),
                '-d', '0',  # Don't retain sequence descriptions
            ]
            
            # Run CD-HIT (may not be available on all systems)
            result = subprocess.run(cmd, capture_output=True, timeout=300)
            
            if result.returncode != 0:
                self.logger.warning("CD-HIT failed, using fallback clustering")
                return None
            
            # Parse .clstr file
            clstr_file = Path(str(output_prefix) + '.clstr')
            representatives = self._parse_cdhit_clusters(clstr_file)
            
            # Clean up temporary files
            Path(temp_fasta).unlink()
            Path(output_prefix).unlink(missing_ok=True)
            clstr_file.unlink(missing_ok=True)
            
            self.logger.info(f"CD-HIT clustering successful: {len(representatives)} clusters")
            return representatives
            
        except Exception as e:
            self.logger.warning(f"CD-HIT clustering failed: {str(e)}")
            return None
    
    def _parse_cdhit_clusters(self, clstr_file: Path) -> Set[str]:
        """
        Parse CD-HIT .clstr file and extract representative sequences.
        
        Args:
            clstr_file: Path to .clstr file
            
        Returns:
            Set of representative accessions
        """
        representatives = set()
        
        try:
            with open(clstr_file, encoding='utf-8') as f:
                for line in f:
                    if line.startswith('>Cluster'):
                        continue
                    # Line format: 0	123aa, >acc_id... *
                    if '*' in line:  # Representative sequence
                        # Extract accession
                        match = line.split('>')[1].split('...')[0]
                        representatives.add(match)
            return representatives
        except Exception as e:
            self.logger.error(f"Error parsing CD-HIT clusters: {str(e)}")
            return set()
    
    def _simple_clustering(self, nucleotide: Dict[str, str]) -> Set[str]:
        """
        Simple clustering by approximate sequence similarity.
        Keeps first sequence in order, removes similar ones.
        
        Args:
            nucleotide: Dictionary of nucleotide sequences
            
        Returns:
            Set of representative accessions
        """
        representatives = set()
        representatives_seqs = []
        
        min_similarity = self.cdhit_identity
        
        for accession, seq in nucleotide.items():
            # Check similarity to existing representatives
            is_similar = False
            for rep_seq in representatives_seqs:
                similarity = self._calculate_similarity(seq, rep_seq)
                if similarity >= min_similarity:
                    is_similar = True
                    break
            
            if not is_similar:
                # Add as new representative
                representatives.add(accession)
                representatives_seqs.append(seq)
        
        self.logger.info(f"Simple clustering: {len(representatives)} clusters from "
                        f"{len(nucleotide)} sequences")
        return representatives
    
    def _calculate_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate sequence similarity (identity).
        
        Args:
            seq1: Sequence 1
            seq2: Sequence 2
            
        Returns:
            Similarity score (0-1)
        """
        min_len = min(len(seq1), len(seq2))
        if min_len == 0:
            return 0.0
        
        matches = sum(1 for i in range(min_len)
                     if seq1[i].upper() == seq2[i].upper())
        
        return matches / max(len(seq1), len(seq2))
    
    def _calculate_md5(self, sequence: str) -> str:
        """Calculate MD5 hash of sequence."""
        return hashlib.md5(sequence.upper().encode()).hexdigest()
    
    def _save_clustering_report(self, metadata: pd.DataFrame, dataset_type: str) -> None:
        """Save clustering summary report."""
        report_path = Path(self.config['reports_dir']) / f"cluster_summary_{dataset_type}.csv"
        
        summary = pd.DataFrame({
            'accession': metadata['accession'],
            'host': metadata['host'],
            'genotype': metadata['genotype'],
            'vp8_coverage': metadata['vp8_coverage'],
        })
        
        summary.to_csv(report_path, index=False)
        self.logger.info(f"Saved clustering report: {report_path}")
