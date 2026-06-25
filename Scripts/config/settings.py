"""
Configuration settings for Rotavirus VP4 preprocessing pipeline.
"""

import logging
from pathlib import Path
from typing import Dict
import sys

# ============================================================================
# PATHS
# ============================================================================

PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data_vp4"
TRAINING_DATA_DIR = DATA_DIR / "Training_data"
EVALUATION_DATA_DIR = DATA_DIR / "Evaluation_dataset"
OUTPUT_DIR = PROJECT_ROOT / "analysis_ready"
REFERENCES_DIR = PROJECT_ROOT / "references"
LOGS_DIR = PROJECT_ROOT / "logs"
REPORTS_DIR = OUTPUT_DIR / "reports"
COMPOSITION_DIR = REPORTS_DIR / "composition"

# ============================================================================
# FILE PATHS
# ============================================================================

# Training data - USING SYNTHETIC METADATA (original metadata doesn't match FASTA)
TRAINING_METADATA = TRAINING_DATA_DIR / "VP4_training_metadata.xlsx"
TRAINING_FASTA = TRAINING_DATA_DIR / "VP4_training_dataset.fasta"

# Evaluation data - USING SYNTHETIC METADATA (original metadata doesn't match FASTA)
EVALUATION_METADATA = EVALUATION_DATA_DIR / "Eval_metadata_combined.xlsx"
EVALUATION_FASTA = EVALUATION_DATA_DIR / "Eval_dataset_comined.fasta"  # Note: typo in original filename

# References (Wa P[8] and DS-1 P[4])
WA_REFERENCE = REFERENCES_DIR / "Wa_P8_VP4_protein.fasta"
DS1_REFERENCE = REFERENCES_DIR / "DS1_P4_VP4_protein.fasta"

# ============================================================================
# SEQUENCE FILTERING PARAMETERS
# ============================================================================

SEQUENCE_LENGTH_MIN = 800  # bp - captures VP8-complete sequences
SEQUENCE_LENGTH_MAX = 2600  # bp - upper bound
AMBIGUITY_MAX_PERCENT = 1.0  # Maximum N content allowed

# ============================================================================
# VP8* EXTRACTION PARAMETERS
# ============================================================================

# VP8* definitions (from technical spec)
VP8_AA_START = 1
VP8_AA_END = 272
VP8_LENGTH = VP8_AA_END - VP8_AA_START + 1  # 272 aa

# VP8* QC thresholds
VP8_COVERAGE_KEEP = 0.90  # Keep if >= 90%
VP8_COVERAGE_REVIEW = 0.80  # Review if 80-89%
VP8_COVERAGE_DISCARD = 0.80  # Discard if < 80%

# Smith-Waterman alignment parameters for VP8 detection
ALIGNMENT_MATCH_SCORE = 2
ALIGNMENT_MISMATCH_SCORE = -1
ALIGNMENT_GAP_OPEN = -5
ALIGNMENT_GAP_EXTEND = -1
ALIGNMENT_MIN_IDENTITY = 0.70  # Minimum 70% identity to reference

# ============================================================================
# DEDUPLICATION AND CLUSTERING
# ============================================================================

# Exact deduplication
EXACT_DEDUP_ENABLED = True

# CD-HIT clustering parameters
CDHIT_IDENTITY_THRESHOLD = 0.99  # 99% nucleotide identity
CDHIT_WORDSIZE = 9  # Word size for comparison

# ============================================================================
# METADATA NORMALIZATION
# ============================================================================

# Host mapping
HOST_NORMALIZATION = {
    'Homo sapiens': 'Human',
    'Sus scrofa': 'Porcine',
    'Bos taurus': 'Bovine',
    'Equus caballus': 'Equine',
    'Equus asinus': 'Equine',
    'Chiroptera': 'Bat',
    'Gallus gallus': 'Avian',
    'Meleagris gallopavo': 'Avian',
    'Phasianus colchicus': 'Avian',
    'Columba': 'Avian',
    'Columba livia': 'Avian',
    'Columbidae': 'Avian',
    'Corvus macrorhynchos': 'Avian',
    'Spilopelia chinensis': 'Avian',
    # Bats
    'Eidolon helvum': 'Bat',
    'Rousettus aegyptiacus': 'Bat',
    'Roussettus aegyptiacus': 'Bat',
    'Eonycteris spelaea': 'Bat',
    'Rousettus leschenaultii': 'Bat',
    'Roussettus leschenaultii': 'Bat',
    'Rhinolophus marshalli': 'Bat',
    'Rhinolophus simulator': 'Bat',
    'Molossus molossus': 'Bat',
    'Glossophaga soricina': 'Bat',
    # Standard normalized mappings (to prevent warnings for already normalized inputs)
    'Human': 'Human',
    'Bat': 'Bat',
    'Porcine': 'Porcine',
    'Bovine': 'Bovine',
    'Equine': 'Equine',
    'Avian': 'Avian',
}

# Genotype normalization
GENOTYPE_NORMALIZATION = {
    'P8': 'P[8]',
    'P4': 'P[4]',
    'P6': 'P[6]',
    'P7': 'P[7]',
    'P9': 'P[9]',
    'P11': 'P[11]',
    'P12': 'P[12]',
    'P13': 'P[13]',
    'P14': 'P[14]',
    'P17': 'P[17]',
    'P19': 'P[19]',
    'P23': 'P[23]',
    'P25': 'P[25]',
    'P30': 'P[30]',
    'P31': 'P[31]',
    'P35': 'P[35]',
    'P37': 'P[37]',
    'P42': 'P[42]',
    'P43': 'P[43]',
    'P47': 'P[47]',
    'P49': 'P[49]',
    'P51': 'P[51]',
    'P56': 'P[56]',
}

# ============================================================================
# LABEL DEFINITIONS
# ============================================================================

# Training set labels
TRAINING_LABELS = {
    'Human_Anchor': 'Positive',
    'Animal_Anchor': 'Negative',
}

# Evaluation set labels
EVALUATION_LABELS = {
    'Porcine': 'Intermediate',
    'Bovine': 'Intermediate',
    'Human_Intermediate': 'Intermediate',
}

# ============================================================================
# QUALITY CONTROL
# ============================================================================

# ORF detection
ORF_MIN_LENGTH = 10  # Minimum ORF length in amino acids

# Stop codon
STOP_CODONS = ['TAA', 'TAG', 'TGA']

# ============================================================================
# REFERENCE SEQUENCES
# ============================================================================

# Real VP4 references from NCBI
# These will be downloaded/provided as FASTA files

REFERENCE_INFO = {
    'Wa_P8': {
        'accession': 'JX406750',
        'genotype': 'P[8]',
        'host': 'Human',
        'description': 'Wa reference - P[8]',
    },
    'DS1_P4': {
        'accession': 'AB118025',
        'genotype': 'P[4]',
        'host': 'Human',
        'description': 'DS-1 reference - P[4]',
    },
}

# ============================================================================
# HOST-GENOTYPE CONFLICT RULES
# ============================================================================

# Flags for biological ambiguity (host-genotype combinations to flag, not remove)
AMBIGUOUS_HOST_GENOTYPE_COMBINATIONS = [
    ('Human', 'P[14]'),  # Human P[14] - zoonotic
    ('Human', 'P[9]'),   # Human P[9] - zoonotic
    ('Human', 'P[19]'),  # Human P[19] - zoonotic
    ('Human', 'P[25]'),  # Human P[25] - zoonotic
    ('Porcine', 'P[6]'), # Porcine P[6] - ambiguous
    ('Bat', 'P[6]'),     # Bat P[6] - unusual
]

# ============================================================================
# LOGGING
# ============================================================================

LOG_LEVEL = logging.INFO
LOG_FORMAT = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
LOG_FILE_FORMAT = '%(asctime)s | %(name)s | %(levelname)s | %(filename)s:%(lineno)d | %(message)s'

# ============================================================================
# CONFIGURATION DICTIONARY
# ============================================================================

CONFIG = {
    # Paths
    'project_root': str(PROJECT_ROOT),
    'data_dir': str(DATA_DIR),
    'training_data_dir': str(TRAINING_DATA_DIR),
    'evaluation_data_dir': str(EVALUATION_DATA_DIR),
    'output_dir': str(OUTPUT_DIR),
    'references_dir': str(REFERENCES_DIR),
    'logs_dir': str(LOGS_DIR),
    'reports_dir': str(REPORTS_DIR),
    'composition_dir': str(COMPOSITION_DIR),
    'analysis_dir': str(OUTPUT_DIR),
    
    # File paths
    'training_metadata': str(TRAINING_METADATA),
    'training_fasta': str(TRAINING_FASTA),
    'evaluation_metadata': str(EVALUATION_METADATA),
    'evaluation_fasta': str(EVALUATION_FASTA),
    'wa_reference': str(WA_REFERENCE),
    'ds1_reference': str(DS1_REFERENCE),
    
    # Sequence filtering
    'sequence_length_min': SEQUENCE_LENGTH_MIN,
    'sequence_length_max': SEQUENCE_LENGTH_MAX,
    'ambiguity_max_percent': AMBIGUITY_MAX_PERCENT,
    
    # VP8* parameters
    'vp8_aa_start': VP8_AA_START,
    'vp8_aa_end': VP8_AA_END,
    'vp8_length': VP8_LENGTH,
    'vp8_coverage_keep': VP8_COVERAGE_KEEP,
    'vp8_coverage_review': VP8_COVERAGE_REVIEW,
    'vp8_coverage_discard': VP8_COVERAGE_DISCARD,
    
    # Alignment parameters
    'alignment_match': ALIGNMENT_MATCH_SCORE,
    'alignment_mismatch': ALIGNMENT_MISMATCH_SCORE,
    'alignment_gap_open': ALIGNMENT_GAP_OPEN,
    'alignment_gap_extend': ALIGNMENT_GAP_EXTEND,
    'alignment_min_identity': ALIGNMENT_MIN_IDENTITY,
    
    # Deduplication
    'exact_dedup_enabled': EXACT_DEDUP_ENABLED,
    'cdhit_identity': CDHIT_IDENTITY_THRESHOLD,
    'cdhit_wordsize': CDHIT_WORDSIZE,
    
    # Normalization
    'host_normalization': HOST_NORMALIZATION,
    'genotype_normalization': GENOTYPE_NORMALIZATION,
    'training_labels': TRAINING_LABELS,
    'evaluation_labels': EVALUATION_LABELS,
    
    # QC
    'orf_min_length': ORF_MIN_LENGTH,
    'stop_codons': STOP_CODONS,
    
    # Conflicts
    'ambiguous_combinations': AMBIGUOUS_HOST_GENOTYPE_COMBINATIONS,
    
    # Feature Extraction
    'esm_model_name': 'facebook/esm2_t12_35M_UR50D',
    'kmer_sizes': [2, 3],
    
    # Logging
    'log_level': LOG_LEVEL,
    'log_format': LOG_FORMAT,
    'log_file_format': LOG_FILE_FORMAT,
}

def setup_logging() -> None:
    """Configure logging for the pipeline."""
    log_dir = Path(CONFIG['logs_dir'])
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"pipeline_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Create formatter
    formatter = logging.Formatter(CONFIG['log_file_format'])
    
    # Setup root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(CONFIG['log_level'])
    
    # Remove existing handlers
    root_logger.handlers = []
    
    # File handler
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(CONFIG['log_level'])
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(CONFIG['log_level'])
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)


def validate_config(config: Dict) -> None:
    """Validate that all required configuration paths exist or are sensible."""
    logger = logging.getLogger(__name__)
    
    # Check that output directories can be created
    logger.info("Validating configuration...")
    
    # Check training data files
    train_meta = Path(config['training_metadata'])
    train_fasta = Path(config['training_fasta'])
    if not train_meta.exists():
        logger.warning(f"Training metadata not found: {train_meta}")
    if not train_fasta.exists():
        logger.warning(f"Training FASTA not found: {train_fasta}")
    
    # Check evaluation data files
    eval_meta = Path(config['evaluation_metadata'])
    eval_fasta = Path(config['evaluation_fasta'])
    if not eval_meta.exists():
        logger.warning(f"Evaluation metadata not found: {eval_meta}")
    if not eval_fasta.exists():
        logger.warning(f"Evaluation FASTA not found: {eval_fasta}")
    
    logger.info("Configuration validation complete")


# Import pandas for timestamp
import pandas as pd
