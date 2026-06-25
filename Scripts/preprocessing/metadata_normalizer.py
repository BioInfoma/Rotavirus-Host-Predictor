"""
Metadata normalization module.
Standardizes host, genotype, and other metadata fields.
"""

import logging
from typing import Dict, Optional
import pandas as pd
import re


class MetadataNormalizer:
    """Normalize metadata fields."""
    
    def __init__(self, config: Dict):
        """
        Initialize normalizer.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.host_map = config['host_normalization']
        self.genotype_map = config['genotype_normalization']
    
    def normalize(self, data: Dict, dataset_type: str) -> Dict:
        """
        Normalize metadata fields.
        
        Args:
            data: Dictionary with metadata
            dataset_type: 'training' or 'evaluation'
            
        Returns:
            Data with normalized metadata
        """
        metadata = data['metadata'].copy()
        
        self.logger.info(f"Normalizing {dataset_type} metadata...")
        
        # Normalize hosts
        metadata['host'] = metadata['host'].apply(self._normalize_host)
        self.logger.info(f"Host normalization complete")
        self.logger.info(f"Unique hosts: {metadata['host'].unique()}")
        
        # Normalize genotypes
        metadata['genotype'] = metadata['genotype'].apply(self._normalize_genotype)
        self.logger.info(f"Genotype normalization complete")
        self.logger.info(f"Unique genotypes (sample): {metadata['genotype'].unique()[:5]}")
        
        # Normalize collection dates
        metadata['collection_year'] = metadata['collection_date'].apply(self._extract_year)
        self.logger.info(f"Collection year extraction complete")
        
        # Standardize column names
        metadata.columns = [col.lower() for col in metadata.columns]
        
        # Ensure required columns exist
        required_cols = [
            'accession', 'full_name', 'host', 'genotype', 'label',
            'adaptation_group', 'country', 'collection_year', 'dataset_split'
        ]
        
        for col in required_cols:
            if col not in metadata.columns:
                metadata[col] = None
        
        # Reorder columns
        metadata = metadata[[
            'accession', 'full_name', 'host', 'genotype', 'label',
            'adaptation_group', 'country', 'collection_year', 'dataset_split',
            'selected_frame', 'protein_length', 'vp8_reference', 'vp8_coverage',
            'vp8_status', 'sequence_md5', 'protein_md5'
        ]]
        
        data['metadata'] = metadata
        return data
    
    def _normalize_host(self, host: str) -> str:
        """
        Normalize host name.
        
        Args:
            host: Original host name
            
        Returns:
            Normalized host name
        """
        if pd.isna(host):
            return 'Unknown'
        
        host = str(host).strip()
        
        # Check direct mapping
        if host in self.host_map:
            return self.host_map[host]
        
        # Try case-insensitive matching
        for key, value in self.host_map.items():
            if key.lower() == host.lower():
                return value
        
        # If no match, return as-is
        self.logger.warning(f"Unknown host: {host}")
        return host
    
    def _normalize_genotype(self, genotype: str) -> str:
        """
        Normalize genotype naming.
        Extract P genotype from full genotype string.
        
        Args:
            genotype: Original genotype string
            
        Returns:
            Normalized P genotype
        """
        if pd.isna(genotype):
            return 'Unknown'
        
        genotype = str(genotype).strip()
        
        # Check direct mapping
        if genotype in self.genotype_map:
            return self.genotype_map[genotype]
        
        # Try to extract P genotype using regex
        # Match patterns like P[8], P8, P[12], P12, etc.
        p_match = re.search(r'P\[?(\d+)\]?', genotype, re.IGNORECASE)
        if p_match:
            p_num = p_match.group(1)
            normalized = f'P[{p_num}]'
            return normalized
        
        # If no clear P genotype, return simplified version
        if genotype.startswith('P'):
            return genotype
        
        self.logger.warning(f"Could not normalize genotype: {genotype}")
        return genotype
    
    def _extract_year(self, date_value) -> Optional[int]:
        """
        Extract year from collection date.
        
        Args:
            date_value: Collection date (various formats possible)
            
        Returns:
            Year as integer or None
        """
        if pd.isna(date_value):
            return None
        
        # Convert to string
        date_str = str(date_value).strip()
        
        # Try to extract year (4 digits)
        year_match = re.search(r'(\d{4})', date_str)
        if year_match:
            try:
                year = int(year_match.group(1))
                # Sanity check for reasonable years
                if 1900 <= year <= 2100:
                    return year
            except ValueError:
                pass
        
        return None
