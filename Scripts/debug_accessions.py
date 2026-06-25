#!/usr/bin/env python3
"""
Debug script to identify accession matching issues between metadata and FASTA files.
Run this FIRST to understand your data structure before running the main pipeline.
"""

import pandas as pd
from pathlib import Path
import sys


def parse_fasta_accessions(fasta_path):
    """Extract all accessions from FASTA file."""
    accessions = []
    with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].rstrip('\r\n')
                # Try multiple parsing strategies
                acc = header.split('|')[0].split()[0].strip()
                accessions.append(acc)
    return accessions


def find_similar_accessions(acc1, acc2_list):
    """Find similar accessions (for cases like 'PP211180' vs 'PP211180.1')."""
    similar = []
    base1 = acc1.split('.')[0]  # Remove version number if present
    
    for acc2 in acc2_list:
        base2 = acc2.split('.')[0]
        
        # Exact match
        if acc1 == acc2:
            similar.append((acc2, 'EXACT'))
        # Match without version number
        elif base1 == base2:
            similar.append((acc2, 'VERSION_DIFF'))
        # Substring match
        elif acc1 in acc2 or acc2 in acc1:
            similar.append((acc2, 'PARTIAL'))
    
    return similar


def analyze_dataset(metadata_path, fasta_path, dataset_name):
    """Analyze a single dataset for accession matching."""
    print(f"\n{'='*80}")
    print(f"ANALYZING: {dataset_name}")
    print(f"{'='*80}\n")
    
    # Load metadata
    try:
        metadata = pd.read_excel(metadata_path)
        meta_accessions = set(metadata['Accession'].astype(str).str.strip())
        print(f"✓ Metadata loaded: {len(meta_accessions)} accessions")
        print(f"  First 5: {sorted(list(meta_accessions))[:5]}")
    except Exception as e:
        print(f"✗ Error loading metadata: {e}")
        return False
    
    # Load FASTA
    try:
        fasta_accessions = set(parse_fasta_accessions(fasta_path))
        print(f"✓ FASTA loaded: {len(fasta_accessions)} accessions")
        print(f"  First 5: {sorted(list(fasta_accessions))[:5]}")
    except Exception as e:
        print(f"✗ Error loading FASTA: {e}")
        return False
    
    # Find matches
    matches = meta_accessions & fasta_accessions
    only_meta = meta_accessions - fasta_accessions
    only_fasta = fasta_accessions - meta_accessions
    
    print(f"\n{'─'*80}")
    print(f"MATCHING RESULTS:")
    print(f"{'─'*80}")
    print(f"✓ Exact matches: {len(matches)}")
    print(f"✗ In metadata only: {len(only_meta)}")
    print(f"✗ In FASTA only: {len(only_fasta)}")
    
    if len(matches) == 0:
        print("\n⚠️  WARNING: NO EXACT MATCHES FOUND!")
        print("\nTrying to find similar accessions...\n")
        
        # Try to find patterns
        fasta_list = sorted(list(fasta_accessions))
        
        found_similar = False
        for meta_acc in sorted(list(only_meta))[:3]:  # Show first 3
            similar = find_similar_accessions(meta_acc, fasta_list)
            if similar:
                found_similar = True
                print(f"Metadata: {meta_acc}")
                for fasta_acc, match_type in similar[:3]:
                    print(f"  → {fasta_acc} ({match_type})")
        
        if not found_similar:
            print("No similar patterns found. Accessions appear to be from different datasets.")
    
    if len(only_meta) > 0:
        print(f"\n📋 Metadata accessions NOT in FASTA (showing first 10 of {len(only_meta)}):")
        for acc in sorted(list(only_meta))[:10]:
            print(f"   {acc}")
    
    if len(only_fasta) > 0:
        print(f"\n📋 FASTA accessions NOT in metadata (showing first 10 of {len(only_fasta)}):")
        for acc in sorted(list(only_fasta))[:10]:
            print(f"   {acc}")
    
    # Statistics
    print(f"\n{'─'*80}")
    print(f"STATISTICS:")
    print(f"{'─'*80}")
    match_percentage = (len(matches) / len(meta_accessions) * 100) if len(meta_accessions) > 0 else 0
    print(f"Match rate: {match_percentage:.1f}% of metadata")
    
    if match_percentage >= 90:
        print("✓ Good match rate - pipeline can proceed")
        return True
    elif match_percentage > 0:
        print("⚠️  Partial match - some sequences will be skipped")
        return True
    else:
        print("✗ No matches - accessions appear to be from different datasets")
        print("\nTo fix this:")
        print("1. Check if accessions have version numbers (.1, .2) that differ")
        print("2. Edit FASTA headers to match metadata accessions")
        print("3. Or edit metadata to match FASTA accessions")
        print("4. Then re-run this debug script to verify")
        return False


def main():
    """Main debug function."""
    print("\n" + "="*80)
    print("ROTAVIRUS VP4 ACCESSION MATCHING DEBUG TOOL")
    print("="*80)
    print("\nThis tool identifies mismatches between metadata and FASTA files.")
    print("Run this BEFORE main.py to diagnose any accession issues.\n")
    
    # Paths
    training_meta = Path("data_vp4/Training_data/VP4_training_metadata.xlsx")
    training_fasta = Path("data_vp4/Training_data/VP4_training_dataset.fasta")
    eval_meta = Path("data_vp4/Evaluation_dataset/Eval_metadata_combined.xlsx")
    eval_fasta = Path("data_vp4/Evaluation_dataset/Eval_dataset_comined.fasta")
    
    # Check training data
    success_train = analyze_dataset(
        training_meta, training_fasta, "TRAINING DATASET"
    )
    
    # Check evaluation data
    success_eval = analyze_dataset(
        eval_meta, eval_fasta, "EVALUATION DATASET"
    )
    
    # Final verdict
    print(f"\n{'='*80}")
    if success_train and success_eval:
        print("✓ READY: You can proceed with main.py")
        print("="*80)
        return 0
    else:
        print("✗ ISSUES FOUND: Fix accession mismatches before running main.py")
        print("="*80)
        print("\nQuick fixes:")
        print("1. Add/remove version numbers (.1) from FASTA or metadata")
        print("2. Use a text editor to batch-replace accession prefixes")
        print("3. Re-run this script to verify changes")
        return 1


if __name__ == "__main__":
    sys.exit(main())
