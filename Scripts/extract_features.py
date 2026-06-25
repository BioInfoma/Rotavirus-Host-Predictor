#!/usr/bin/env python3
"""
Rotavirus A VP4 Host Adaptation / Spillover Prediction - Feature Extraction Pipeline

This script extracts:
1. ESM-2 Protein Embeddings: Dense evolutionary representations from a pre-trained protein language model.
2. k-mer Motif Frequencies: Standardized local sequence motif counts (2-mers and 3-mers) trained strictly on
   the training set to prevent data leakage.
3. Combined Features: A consolidated representation of both feature sets.

All features are aligned with their metadata and saved in the analysis_ready/ folder.

Author: Automated ML Pipeline
Date: 2026
"""

import os
import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO
from tqdm import tqdm

# Add project root and Scripts directory to Python path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(SCRIPT_DIR))

# Import settings and logging
from config.settings import CONFIG

def setup_extraction_logging() -> logging.Logger:
    """Set up loggers for the feature extraction process."""
    log_dir = Path(CONFIG['logs_dir'])
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "feature_extraction.log"
    
    logger = logging.getLogger("feature_extraction")
    logger.setLevel(logging.INFO)
    
    # Avoid duplicate handlers if script is re-run in same session
    if logger.handlers:
        return logger
        
    formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
    
    # File handler
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(ch)
    
    return logger

def load_sequences(fasta_path: Path) -> dict:
    """Load sequences from a FASTA file into a dictionary of {accession: sequence_string}."""
    seq_dict = {}
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
        
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Header accessions might contain extra spaces or comments, split at whitespace
        accession = record.id.strip()
        seq_dict[accession] = str(record.seq).upper().replace("-", "") # Strip any gaps if present
    return seq_dict

def extract_esm2_embeddings(seq_dict: dict, model_name: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Extract dense sequence-level protein embeddings from a pre-trained ESM-2 model.
    Performs mean pooling across the sequence length, excluding special tokens (<cls>, <eos>).
    """
    logger.info(f"Loading ESM-2 tokenizer and model: {model_name}...")
    import torch
    from transformers import AutoTokenizer, AutoModel
    
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    
    # Check device (GPU vs CPU)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    model.eval()
    logger.info(f"Model successfully loaded on device: {device}")
    
    embeddings = {}
    
    logger.info("Extracting embeddings for sequences...")
    for accession, seq in tqdm(seq_dict.items(), desc="ESM-2 Embeddings"):
        # ESM-2 requires uppercase amino acids
        seq_clean = seq.upper()
        
        # Tokenize sequence
        inputs = tokenizer(seq_clean, return_tensors="pt")
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        with torch.no_grad():
            outputs = model(**inputs)
            
        # Shape: (batch_size, sequence_length, hidden_size) -> (1, L_special, hidden_size)
        last_hidden_state = outputs.last_hidden_state[0] # Shape: (L_special, hidden_size)
        
        # Slice off <cls> (index 0) and <eos> (index L_special - 1)
        residue_states = last_hidden_state[1:-1, :] # Shape: (L, hidden_size)
        
        # Compute mean representation (mean pooling)
        mean_embedding = residue_states.mean(dim=0).cpu().numpy()
        embeddings[accession] = mean_embedding
        
    # Convert to DataFrame
    emb_dim = len(next(iter(embeddings.values())))
    cols = [f"esm_dim_{i}" for i in range(emb_dim)]
    
    df_esm = pd.DataFrame.from_dict(embeddings, orient="index", columns=cols)
    df_esm.index.name = "accession"
    df_esm = df_esm.reset_index()
    
    logger.info(f"Extracted ESM embeddings of shape {df_esm.shape}")
    return df_esm

def extract_kmer_frequencies(
    train_seqs: dict, 
    eval_seqs: dict, 
    k_sizes: list, 
    logger: logging.Logger
) -> tuple:
    """
    Extract k-mer relative frequency features.
    Learns the k-mer vocabulary STRICTLY on the training set to prevent leakage.
    Normalizes frequency counts to sum to 1.0 (L1-normalization).
    """
    logger.info(f"Fitting CountVectorizer for k-mer sizes {k_sizes} strictly on training sequences...")
    from sklearn.feature_extraction.text import CountVectorizer
    from sklearn.preprocessing import normalize
    
    # Extract training and evaluation accessions and sequences
    train_accs = list(train_seqs.keys())
    train_raw = list(train_seqs.values())
    
    eval_accs = list(eval_seqs.keys())
    eval_raw = list(eval_seqs.values())
    
    # We fit a CountVectorizer. For amino acids, we analyze characters.
    # ngram_range=(min_k, max_k) will fit all lengths between min(k_sizes) and max(k_sizes)
    min_k = min(k_sizes)
    max_k = max(k_sizes)
    
    vectorizer = CountVectorizer(
        analyzer='char',
        ngram_range=(min_k, max_k),
        lowercase=False
    )
    
    # Fit strictly on training set
    vectorizer.fit(train_raw)
    vocab_size = len(vectorizer.vocabulary_)
    logger.info(f"K-mer vocabulary learned. Total unique k-mers (k={min_k}-{max_k}): {vocab_size}")
    
    # Transform both datasets to raw counts
    train_counts = vectorizer.transform(train_raw).toarray()
    eval_counts = vectorizer.transform(eval_raw).toarray()
    
    # Normalize count vectors to sum to 1.0 (relative frequencies) to handle length variations
    logger.info("Normalizing k-mer counts to relative frequencies (L1-normalization)...")
    train_freqs = normalize(train_counts, norm='l1')
    eval_freqs = normalize(eval_counts, norm='l1')
    
    # Create DataFrames
    features = [f"kmer_{feature}" for feature in vectorizer.get_feature_names_out()]
    
    df_kmer_train = pd.DataFrame(train_freqs, columns=features, index=train_accs)
    df_kmer_train.index.name = "accession"
    df_kmer_train = df_kmer_train.reset_index()
    
    df_kmer_eval = pd.DataFrame(eval_freqs, columns=features, index=eval_accs)
    df_kmer_eval.index.name = "accession"
    df_kmer_eval = df_kmer_eval.reset_index()
    
    logger.info(f"Training k-mer features shape: {df_kmer_train.shape}")
    logger.info(f"Evaluation k-mer features shape: {df_kmer_eval.shape}")
    
    return df_kmer_train, df_kmer_eval

def main():
    logger = setup_extraction_logging()
    logger.info("=" * 80)
    logger.info("ROTAVIRUS VP4 FEATURE EXTRACTION PIPELINE")
    logger.info("=" * 80)
    
    try:
        # 1. Define paths
        analysis_dir = Path(CONFIG['analysis_dir'])
        
        train_metadata_path = analysis_dir / "cleaned_metadata_training.csv"
        eval_metadata_path = analysis_dir / "cleaned_metadata_evaluation.csv"
        
        train_fasta_path = analysis_dir / "cleaned_vp8_training.fasta"
        eval_fasta_path = analysis_dir / "cleaned_vp8_evaluation.fasta"
        
        # Validate input paths
        for path in [train_metadata_path, eval_metadata_path, train_fasta_path, eval_fasta_path]:
            if not path.exists():
                raise FileNotFoundError(f"Required preprocessing output is missing: {path}")
                
        # 2. Load sequence datasets
        logger.info("Loading preprocessed training and evaluation sequences...")
        train_seqs = load_sequences(train_fasta_path)
        eval_seqs = load_sequences(eval_fasta_path)
        logger.info(f"Loaded {len(train_seqs)} training and {len(eval_seqs)} evaluation sequences.")
        
        # 3. Load metadata tables
        logger.info("Loading preprocessed metadata tables...")
        df_meta_train = pd.read_csv(train_metadata_path)
        df_meta_eval = pd.read_csv(eval_metadata_path)
        logger.info(f"Metadata rows: training={len(df_meta_train)}, evaluation={len(df_meta_eval)}")
        
        # Verify accessions match exactly between fasta and metadata
        missing_train_seqs = set(df_meta_train['accession']) - set(train_seqs.keys())
        missing_eval_seqs = set(df_meta_eval['accession']) - set(eval_seqs.keys())
        if missing_train_seqs:
            logger.warning(f"Training metadata accessions missing from FASTA: {missing_train_seqs}")
        if missing_eval_seqs:
            logger.warning(f"Evaluation metadata accessions missing from FASTA: {missing_eval_seqs}")
            
        # 4. Extract k-mer frequencies (unsupervised, fit on train, transform train & eval)
        df_kmer_train, df_kmer_eval = extract_kmer_frequencies(
            train_seqs, 
            eval_seqs, 
            CONFIG['kmer_sizes'], 
            logger
        )
        
        # 5. Extract ESM-2 Embeddings
        df_esm_train = extract_esm2_embeddings(train_seqs, CONFIG['esm_model_name'], logger)
        df_esm_eval = extract_esm2_embeddings(eval_seqs, CONFIG['esm_model_name'], logger)
        
        # 6. Merge features with metadata
        logger.info("Merging features with metadata tables...")
        
        # Helper to merge and check integrity
        def merge_features(df_meta, df_feat, name):
            merged = pd.merge(df_meta, df_feat, on="accession", how="inner")
            if len(merged) != len(df_meta):
                logger.warning(f"Merge mismatch for {name}: metadata={len(df_meta)}, merged={len(merged)}")
            return merged
            
        # ESM-2 Features Datasets
        df_esm_train_final = merge_features(df_meta_train, df_esm_train, "ESM Training")
        df_esm_eval_final = merge_features(df_meta_eval, df_esm_eval, "ESM Evaluation")
        
        # k-mer Features Datasets
        df_kmer_train_final = merge_features(df_meta_train, df_kmer_train, "k-mer Training")
        df_kmer_eval_final = merge_features(df_meta_eval, df_kmer_eval, "k-mer Evaluation")
        
        # Combined Features Datasets (ESM + k-mer)
        df_combined_train = pd.merge(df_esm_train, df_kmer_train, on="accession", how="inner")
        df_combined_eval = pd.merge(df_esm_eval, df_kmer_eval, on="accession", how="inner")
        
        df_combined_train_final = merge_features(df_meta_train, df_combined_train, "Combined Training")
        df_combined_eval_final = merge_features(df_meta_eval, df_combined_eval, "Combined Evaluation")
        
        # 7. Save outputs
        logger.info("Saving extracted feature files to analysis_ready/ ...")
        
        # ESM outputs
        df_esm_train_final.to_csv(analysis_dir / "features_esm2_training.csv", index=False)
        df_esm_eval_final.to_csv(analysis_dir / "features_esm2_evaluation.csv", index=False)
        logger.info(f"Saved ESM-2 CSVs: {analysis_dir / 'features_esm2_[training/evaluation].csv'}")
        
        # k-mer outputs
        df_kmer_train_final.to_csv(analysis_dir / "features_kmer_training.csv", index=False)
        df_kmer_eval_final.to_csv(analysis_dir / "features_kmer_evaluation.csv", index=False)
        logger.info(f"Saved k-mer CSVs: {analysis_dir / 'features_kmer_[training/evaluation].csv'}")
        
        # Combined outputs
        df_combined_train_final.to_csv(analysis_dir / "features_combined_training.csv", index=False)
        df_combined_eval_final.to_csv(analysis_dir / "features_combined_evaluation.csv", index=False)
        logger.info(f"Saved combined CSVs: {analysis_dir / 'features_combined_[training/evaluation].csv'}")
        
        logger.info("=" * 80)
        logger.info("FEATURE EXTRACTION COMPLETED SUCCESSFULLY!")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Feature extraction failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
