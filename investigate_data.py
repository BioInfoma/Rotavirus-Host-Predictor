#!/usr/bin/env python3
"""Investigate accession mismatches between individual CSVs, combined Excel, and FASTA files."""

import pandas as pd
from Bio import SeqIO
from pathlib import Path
import sys

BASE = Path(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4")

def get_fasta_ids(path, max_show=5):
    """Read FASTA, return list of (raw_id, raw_description) and full set of IDs."""
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        records.append((rec.id, rec.description))
    return records

def show_separator(title):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

def compare_sets(name_a, set_a, name_b, set_b):
    overlap = set_a & set_b
    only_a = set_a - set_b
    only_b = set_b - set_a
    print(f"\n  --- {name_a} vs {name_b} ---")
    print(f"  {name_a} total: {len(set_a)}")
    print(f"  {name_b} total: {len(set_b)}")
    print(f"  Overlap: {len(overlap)}")
    print(f"  Only in {name_a}: {len(only_a)}")
    if only_a:
        print(f"    Examples: {sorted(only_a)[:8]}")
    print(f"  Only in {name_b}: {len(only_b)}")
    if only_b:
        print(f"    Examples: {sorted(only_b)[:8]}")
    return overlap, only_a, only_b

# ============================================================
# TRAINING DATA
# ============================================================
show_separator("TRAINING DATA INVESTIGATION")

# 1. Individual CSVs
train_csvs = {
    "human": BASE / "Training_data" / "human_metadata.csv",
    "bat": BASE / "Training_data" / "bat_metadata.csv",
    "avian": BASE / "Training_data" / "Avian_metadata.csv",
    "equine": BASE / "Training_data" / "equine_metadata.csv",
}

all_csv_accessions = set()
for name, path in train_csvs.items():
    df = pd.read_csv(path)
    print(f"\n  [{name.upper()}] {path.name}: {len(df)} rows")
    print(f"    Columns: {list(df.columns)}")
    # Find accession column
    acc_col = None
    for c in df.columns:
        if 'accession' in c.lower() or 'acc' == c.lower():
            acc_col = c
            break
    if acc_col is None:
        acc_col = df.columns[0]  # assume first column
    print(f"    Accession column: '{acc_col}'")
    accessions = set(df[acc_col].astype(str).str.strip())
    print(f"    Accession count: {len(accessions)}")
    print(f"    First 5: {sorted(accessions)[:5]}")
    all_csv_accessions.update(accessions)

print(f"\n  TOTAL individual CSV accessions (combined): {len(all_csv_accessions)}")

# 2. Combined Excel
xlsx_path = BASE / "Training_data" / "VP4_training_metadata.xlsx"
xlsx_df = pd.read_excel(xlsx_path, engine="openpyxl")
print(f"\n  [COMBINED EXCEL] {xlsx_path.name}: {len(xlsx_df)} rows")
print(f"    Columns: {list(xlsx_df.columns)}")
print(f"    First 5 rows:")
print(xlsx_df.head().to_string(index=False))
acc_col_xlsx = None
for c in xlsx_df.columns:
    if 'accession' in c.lower() or 'acc' == c.lower():
        acc_col_xlsx = c
        break
if acc_col_xlsx is None:
    acc_col_xlsx = xlsx_df.columns[0]
print(f"    Accession column: '{acc_col_xlsx}'")
xlsx_accessions = set(xlsx_df[acc_col_xlsx].dropna().astype(str).str.strip()) - {'nan', 'None', ''}
print(f"    Accession count: {len(xlsx_accessions)}")
print(f"    First 10: {sorted(xlsx_accessions)[:10]}")

# 3. Individual FASTAs
train_fastas = {
    "human": BASE / "Training_data" / "human_seq.fasta",
    "bat": BASE / "Training_data" / "bat_seq.fasta",
    "avian": BASE / "Training_data" / "Avian_seq.fasta",
    "equine": BASE / "Training_data" / "equine_seq.fasta",
}

all_fasta_ids = set()
all_fasta_ids_base = set()  # without version
for name, path in train_fastas.items():
    records = get_fasta_ids(path)
    ids = set()
    ids_base = set()
    print(f"\n  [{name.upper()} FASTA] {path.name}: {len(records)} sequences")
    print(f"    First 5 headers:")
    for rid, rdesc in records[:5]:
        print(f"      ID='{rid}' | DESC='{rdesc}'")
        ids.add(rid)
        ids_base.add(rid.split('.')[0])
    for rid, rdesc in records[5:]:
        ids.add(rid)
        ids_base.add(rid.split('.')[0])
    all_fasta_ids.update(ids)
    all_fasta_ids_base.update(ids_base)

print(f"\n  TOTAL individual FASTA IDs (combined): {len(all_fasta_ids)}")
print(f"  TOTAL individual FASTA base IDs (no version): {len(all_fasta_ids_base)}")

# 4. Combined FASTA
combined_fasta_path = BASE / "Training_data" / "VP4_training_dataset.fasta"
combined_records = get_fasta_ids(combined_fasta_path)
combined_fasta_ids = set()
combined_fasta_ids_base = set()
print(f"\n  [COMBINED FASTA] {combined_fasta_path.name}: {len(combined_records)} sequences")
print(f"    First 10 headers:")
for rid, rdesc in combined_records[:10]:
    print(f"      ID='{rid}' | DESC='{rdesc}'")
    combined_fasta_ids.add(rid)
    combined_fasta_ids_base.add(rid.split('.')[0])
for rid, rdesc in combined_records[10:]:
    combined_fasta_ids.add(rid)
    combined_fasta_ids_base.add(rid.split('.')[0])

# 5. Comparisons
show_separator("TRAINING COMPARISONS")

# CSV accessions - strip version numbers for comparison too
csv_acc_base = set(a.split('.')[0] for a in all_csv_accessions)
xlsx_acc_base = set(a.split('.')[0] for a in xlsx_accessions)

print("\n>> RAW accession comparisons:")
compare_sets("Individual_CSVs", all_csv_accessions, "Combined_Excel", xlsx_accessions)
compare_sets("Individual_FASTAs", all_fasta_ids, "Combined_FASTA", combined_fasta_ids)
compare_sets("Combined_Excel", xlsx_accessions, "Combined_FASTA_IDs", combined_fasta_ids)

print("\n>> BASE accession comparisons (version stripped):")
compare_sets("Individual_CSVs_base", csv_acc_base, "Combined_Excel_base", xlsx_acc_base)
compare_sets("Individual_FASTAs_base", all_fasta_ids_base, "Combined_FASTA_base", combined_fasta_ids_base)
compare_sets("Combined_Excel_base", xlsx_acc_base, "Combined_FASTA_base", combined_fasta_ids_base)
compare_sets("Individual_CSVs_base", csv_acc_base, "Individual_FASTAs_base", all_fasta_ids_base)
compare_sets("Individual_CSVs_base", csv_acc_base, "Combined_FASTA_base", combined_fasta_ids_base)

# ============================================================
# EVALUATION DATA
# ============================================================
show_separator("EVALUATION DATA INVESTIGATION")

eval_csvs = {
    "pig": BASE / "Evaluation_dataset" / "pig_metadata.csv",
    "bovine": BASE / "Evaluation_dataset" / "Bovine_metadata.csv",
    "human_eval": BASE / "Evaluation_dataset" / "Human_eval_metadata.csv",
}

all_eval_csv_accessions = set()
for name, path in eval_csvs.items():
    df = pd.read_csv(path)
    print(f"\n  [{name.upper()}] {path.name}: {len(df)} rows")
    print(f"    Columns: {list(df.columns)}")
    acc_col = None
    for c in df.columns:
        if 'accession' in c.lower() or 'acc' == c.lower():
            acc_col = c
            break
    if acc_col is None:
        acc_col = df.columns[0]
    print(f"    Accession column: '{acc_col}'")
    accessions = set(df[acc_col].astype(str).str.strip())
    print(f"    Accession count: {len(accessions)}")
    print(f"    First 5: {sorted(accessions)[:5]}")
    all_eval_csv_accessions.update(accessions)

print(f"\n  TOTAL individual eval CSV accessions: {len(all_eval_csv_accessions)}")

# Combined Excel
eval_xlsx_path = BASE / "Evaluation_dataset" / "Eval_metadata_combined.xlsx"
eval_xlsx_df = pd.read_excel(eval_xlsx_path, engine="openpyxl")
print(f"\n  [COMBINED EXCEL] {eval_xlsx_path.name}: {len(eval_xlsx_df)} rows")
print(f"    Columns: {list(eval_xlsx_df.columns)}")
print(f"    First 5 rows:")
print(eval_xlsx_df.head().to_string(index=False))
acc_col_eval = None
for c in eval_xlsx_df.columns:
    if 'accession' in c.lower() or 'acc' == c.lower():
        acc_col_eval = c
        break
if acc_col_eval is None:
    acc_col_eval = eval_xlsx_df.columns[0]
print(f"    Accession column: '{acc_col_eval}'")
eval_xlsx_accessions = set(eval_xlsx_df[acc_col_eval].dropna().astype(str).str.strip()) - {'nan', 'None', ''}
print(f"    Accession count: {len(eval_xlsx_accessions)}")
print(f"    First 10: {sorted(eval_xlsx_accessions)[:10]}")

# Individual FASTAs
eval_fastas = {
    "pig": BASE / "Evaluation_dataset" / "pig_seq.fasta",
    "bovine": BASE / "Evaluation_dataset" / "Bovine_seq.fasta",
    "human_eval": BASE / "Evaluation_dataset" / "Human_eval_seq.fasta",
}

all_eval_fasta_ids = set()
all_eval_fasta_ids_base = set()
for name, path in eval_fastas.items():
    records = get_fasta_ids(path)
    ids = set()
    ids_base = set()
    print(f"\n  [{name.upper()} FASTA] {path.name}: {len(records)} sequences")
    print(f"    First 5 headers:")
    for rid, rdesc in records[:5]:
        print(f"      ID='{rid}' | DESC='{rdesc}'")
        ids.add(rid)
        ids_base.add(rid.split('.')[0])
    for rid, rdesc in records[5:]:
        ids.add(rid)
        ids_base.add(rid.split('.')[0])
    all_eval_fasta_ids.update(ids)
    all_eval_fasta_ids_base.update(ids_base)

print(f"\n  TOTAL individual eval FASTA IDs: {len(all_eval_fasta_ids)}")

# Combined FASTA
eval_combined_fasta = BASE / "Evaluation_dataset" / "Eval_dataset_comined.fasta"
eval_combined_records = get_fasta_ids(eval_combined_fasta)
eval_combined_fasta_ids = set()
eval_combined_fasta_ids_base = set()
print(f"\n  [COMBINED FASTA] {eval_combined_fasta.name}: {len(eval_combined_records)} sequences")
print(f"    First 10 headers:")
for rid, rdesc in eval_combined_records[:10]:
    print(f"      ID='{rid}' | DESC='{rdesc}'")
    eval_combined_fasta_ids.add(rid)
    eval_combined_fasta_ids_base.add(rid.split('.')[0])
for rid, rdesc in eval_combined_records[10:]:
    eval_combined_fasta_ids.add(rid)
    eval_combined_fasta_ids_base.add(rid.split('.')[0])

# Comparisons
show_separator("EVALUATION COMPARISONS")

eval_csv_base = set(a.split('.')[0] for a in all_eval_csv_accessions)
eval_xlsx_base = set(a.split('.')[0] for a in eval_xlsx_accessions)

print("\n>> RAW accession comparisons:")
compare_sets("Individual_CSVs", all_eval_csv_accessions, "Combined_Excel", eval_xlsx_accessions)
compare_sets("Individual_FASTAs", all_eval_fasta_ids, "Combined_FASTA", eval_combined_fasta_ids)
compare_sets("Combined_Excel", eval_xlsx_accessions, "Combined_FASTA_IDs", eval_combined_fasta_ids)

print("\n>> BASE accession comparisons (version stripped):")
compare_sets("Individual_CSVs_base", eval_csv_base, "Combined_Excel_base", eval_xlsx_base)
compare_sets("Individual_FASTAs_base", all_eval_fasta_ids_base, "Combined_FASTA_base", eval_combined_fasta_ids_base)
compare_sets("Combined_Excel_base", eval_xlsx_base, "Combined_FASTA_base", eval_combined_fasta_ids_base)
compare_sets("Individual_CSVs_base", eval_csv_base, "Individual_FASTAs_base", all_eval_fasta_ids_base)
compare_sets("Individual_CSVs_base", eval_csv_base, "Combined_FASTA_base", eval_combined_fasta_ids_base)

print("\n\nDONE.")
