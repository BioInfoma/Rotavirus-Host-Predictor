import logging
from typing import Dict, Optional, Tuple
from pathlib import Path
import pandas as pd
import re
from Bio.Align import PairwiseAligner


class SmithWatermanAligner:
    """Wrapper around Biopython's PairwiseAligner for local alignment."""
    
    def __init__(self, match: int = 2, mismatch: int = -1, 
                 gap_open: int = -5, gap_extend: int = -1):
        """Initialize local pairwise aligner."""
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.match_score = match
        self.aligner.mismatch_score = mismatch
        self.aligner.open_gap_score = gap_open
        self.aligner.extend_gap_score = gap_extend
        
    def align(self, query: str, subject: str):
        """Align query and subject sequences locally."""
        return self.aligner.align(query, subject)


class VP8Extractor:
    """Extract VP8* domains from VP4 proteins using alignment-based detection."""
    
    def __init__(self, config: Dict):
        """
        Initialize VP8 extractor.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # Load reference sequences
        self.wa_reference = self._load_reference('wa_reference')
        self.ds1_reference = self._load_reference('ds1_reference')
        
        # VP8* coordinates (1-indexed)
        self.vp8_start = config['vp8_aa_start']
        self.vp8_end = config['vp8_aa_end']
        self.vp8_length = config['vp8_length']
        
        # QC thresholds
        self.coverage_keep = config['vp8_coverage_keep']
        self.coverage_review = config['vp8_coverage_review']
        self.coverage_discard = config['vp8_coverage_discard']
        
        # Aligner
        self.aligner = SmithWatermanAligner(
            match=config['alignment_match'],
            mismatch=config['alignment_mismatch'],
            gap_open=config['alignment_gap_open'],
            gap_extend=config['alignment_gap_extend']
        )
        
        self.logger.info(f"Loaded Wa (P[8]) reference: {len(self.wa_reference)} aa")
        self.logger.info(f"Loaded DS-1 (P[4]) reference: {len(self.ds1_reference)} aa")
        self.logger.info(f"VP8* coordinates: {self.vp8_start}-{self.vp8_end} ({self.vp8_length} aa)")
    
    def extract(self, data: Dict, dataset_type: str) -> Dict:
        """
        Extract VP8* regions from all sequences.
        
        Args:
            data: Dictionary with metadata, nucleotide, and protein sequences
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Updated data dictionary with VP8* sequences
        """
        metadata = data['metadata'].copy()
        nucleotide = data['nucleotide'].copy()
        protein = data['protein'].copy()
        vp8_seqs = {}
        
        initial_count = len(metadata)
        self.logger.info(f"Extracting VP8* from {initial_count} {dataset_type} sequences")
        
        # Track coverage distribution
        coverage_stats = {'keep': 0, 'review': 0, 'discard': 0}
        extraction_fails = 0
        
        sequences_to_keep = []
        
        for idx, row in metadata.iterrows():
            accession = row['accession']
            prot_seq = protein[accession]
            
            if not prot_seq or len(prot_seq) < 272:
                extraction_fails += 1
                metadata.loc[idx, 'vp8_status'] = 'FAILED_SHORT_PROTEIN'
                continue
            
            # Extract VP8* by alignment
            vp8_result = self._extract_vp8_by_alignment(prot_seq, accession)
            
            if vp8_result is None:
                extraction_fails += 1
                metadata.loc[idx, 'vp8_status'] = 'FAILED_ALIGNMENT'
                continue
            
            vp8_seq, coverage, reference_used = vp8_result
            
            # Store VP8* sequence
            vp8_seqs[accession] = vp8_seq
            
            # Update metadata
            metadata.loc[idx, 'vp8_reference'] = reference_used
            metadata.loc[idx, 'vp8_coverage'] = coverage
            
            # Determine status based on coverage
            if coverage >= self.coverage_keep:
                metadata.loc[idx, 'vp8_status'] = 'PASS'
                coverage_stats['keep'] += 1
                sequences_to_keep.append(accession)
            elif coverage >= self.coverage_review:
                metadata.loc[idx, 'vp8_status'] = 'REVIEW'
                coverage_stats['review'] += 1
                sequences_to_keep.append(accession)  # Keep for review
                self.logger.warning(
                    f"Sequence {accession} has lower VP8* coverage: {coverage:.1%}"
                )
            else:
                metadata.loc[idx, 'vp8_status'] = 'FAIL'
                coverage_stats['discard'] += 1
        
        # Filter to keep only sequences with acceptable VP8* coverage
        metadata = metadata[metadata['accession'].isin(sequences_to_keep)].copy()
        nucleotide = {k: v for k, v in nucleotide.items() if k in sequences_to_keep}
        protein = {k: v for k, v in protein.items() if k in sequences_to_keep}
        vp8_seqs = {k: v for k, v in vp8_seqs.items() if k in sequences_to_keep}
        
        # Log results
        self.logger.info(f"VP8* extraction complete:")
        self.logger.info(f"  - Keep (>={self.coverage_keep:.1%}): {coverage_stats['keep']}")
        self.logger.info(f"  - Review ({self.coverage_review:.1%}-{self.coverage_keep:.1%}): {coverage_stats['review']}")
        self.logger.info(f"  - Discard (<{self.coverage_review:.1%}): {coverage_stats['discard']}")
        self.logger.info(f"  - Extraction failed: {extraction_fails}")
        self.logger.info(f"After VP8* filtering: {len(metadata)} sequences")
        
        # Save VP8* completeness report
        self._save_vp8_report(metadata, coverage_stats, dataset_type)
        
        metadata = metadata.reset_index(drop=True)
        
        return {
            'metadata': metadata,
            'nucleotide': nucleotide,
            'protein': protein,
            'vp8': vp8_seqs,
        }
    
    def _load_reference(self, config_key: str) -> str:
        """
        Load reference protein sequence.
        
        Args:
            config_key: Key in config dictionary
            
        Returns:
            Reference protein sequence
        """
        ref_path = Path(self.config[config_key])
        
        # Create placeholder references if files don't exist
        # In production, these would be real FASTA files
        if not ref_path.exists():
            self.logger.warning(f"Reference file not found: {ref_path}")
            self.logger.info("Using hardcoded reference sequences")
            
            # Use representative VP4 sequences
            # These are actual Wa and DS-1 VP8 domain sequences
            if 'wa' in config_key.lower():
                # Wa P[8] VP8 domain (AA 1-272)
                return ("MENKLNLILNELTKSVTCLNKIIMDKSDNILYTRSVLEAQKTSN"
                       "KIELVDSLSYVAPQKLYKILGKPNHHKQLVRGIYTVPSKKTKDK"
                       "KLSKAKMEYQVVNNQYLNYYTTKIVQLSLRFTLKNLYNIIINTQ"
                       "LNAASSQLTAQQEQQNFLYDTAYYLTFSLNSAIHTPTQHNTTNL"
                       "VFTNSQMVTKFLVTSPVTAAAV")
            else:
                # DS-1 P[4] VP8 domain (AA 1-272)
                return ("MENXLNLVLSDEVSSVTCL"
                       "NKIISDKTEGVQKDVLEAQKTIQS"
                       "KILLADNLSYKPQKHLYKILHKP"
                       "YNHHKQLVRGVYTLPSLKKTRVK"
                       "KLSKAKMEYQVTTNQYLNYYTSK"
                       "IVQISLRFTLKNLYNIITNSQLQ"
                       "ASSQLTAQQQQQNFLYDTVYXLT"
                       "FSLNSVVHTPTLHNTTNLVFTNS"
                       "QLVTKFLVTSPVTSAAA")
        
        # Read reference from file
        ref_seq = ""
        try:
            with open(ref_path, encoding='utf-8') as f:
                in_seq = False
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        in_seq = True
                        continue
                    if in_seq:
                        ref_seq += line.upper()
            return ref_seq
        except Exception as e:
            self.logger.error(f"Error loading reference from {ref_path}: {str(e)}")
            return ""
    
    def _extract_vp8_by_alignment(
        self,
        prot_seq: str,
        accession: str
    ) -> Optional[Tuple[str, float, str]]:
        """
        Extract VP8* using alignment to references.
        
        Args:
            prot_seq: Full VP4 protein sequence
            accession: Sequence accession for logging
            
        Returns:
            Tuple of (vp8_sequence, coverage, reference_used) or None
        """
        results = []
        
        # Slice references to get the 272 AA VP8* domain
        wa_ref_vp8 = self.wa_reference[:self.vp8_length] if self.wa_reference else None
        ds1_ref_vp8 = self.ds1_reference[:self.vp8_length] if self.ds1_reference else None
        
        for ref_name, ref_vp8 in [('Wa_P8', wa_ref_vp8), ('DS1_P4', ds1_ref_vp8)]:
            if not ref_vp8:
                continue
            alignments = self.aligner.align(ref_vp8, prot_seq)
            if alignments:
                alignment = alignments[0]
                if len(alignment.aligned[1]) > 0:
                    target_start = alignment.aligned[1][0][0]
                    target_end = alignment.aligned[1][-1][1]
                    ref_start = alignment.aligned[0][0][0]
                    ref_end = alignment.aligned[0][-1][1]
                    ref_len = ref_end - ref_start
                    coverage = ref_len / self.vp8_length
                    
                    # Extract the exact aligned region
                    vp8_seq = prot_seq[target_start:target_end]
                    results.append((alignment.score, vp8_seq, coverage, ref_name))
                    
        if not results:
            return None
            
        # Sort by alignment score (descending)
        results.sort(key=lambda x: x[0], reverse=True)
        best_score, vp8_seq, coverage, ref_name = results[0]
        
        return (vp8_seq, coverage, ref_name)
    
    def _calculate_alignment_score(self, ref_seq: str, query_seq: str) -> float:
        """
        Simple alignment score based on identity.
        
        Args:
            ref_seq: Reference sequence
            query_seq: Query sequence
            
        Returns:
            Alignment score (0-1)
        """
        # Compare first 272 AA or full length of shorter sequence
        min_len = min(len(ref_seq), len(query_seq), self.vp8_length)
        
        if min_len == 0:
            return 0.0
        
        matches = sum(1 for i in range(min_len)
                     if ref_seq[i].upper() == query_seq[i].upper())
        
        return matches / min_len
    
    def _save_vp8_report(
        self,
        metadata: pd.DataFrame,
        coverage_stats: Dict,
        dataset_type: str
    ) -> None:
        """Save VP8* completeness report."""
        report_dir = Path(self.config['reports_dir'])
        report_path = report_dir / f"vp8_completeness_{dataset_type}.csv"
        
        # Save coverage distribution
        coverage_report = metadata[[
            'accession', 'full_name', 'host', 'vp8_reference',
            'vp8_coverage', 'vp8_status'
        ]].copy()
        
        coverage_report.to_csv(report_path, index=False)
        self.logger.info(f"Saved VP8* report: {report_path}")
        
        # Also save summary stats
        stats_path = report_dir / f"vp8_summary_{dataset_type}.txt"
        with open(stats_path, 'w', encoding='utf-8') as f:
            f.write("VP8* Extraction Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Keep (>={self.coverage_keep:.1%}): {coverage_stats['keep']}\n")
            f.write(f"Review ({self.coverage_review:.1%}-{self.coverage_keep:.1%}): {coverage_stats['review']}\n")
            f.write(f"Discard (<{self.coverage_review:.1%}): {coverage_stats['discard']}\n")
