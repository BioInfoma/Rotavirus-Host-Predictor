#!/usr/bin/env python3
"""
Apply data fixes to the combined training and evaluation files:
1. Remove 22 porcine/bat sequences from the combined training FASTA
2. Harmonise accession naming (add version suffixes to Excel to match FASTA)
3. Verify evaluation files are correct
4. Back up originals before overwriting
"""

import shutil
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

BASE = Path(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4")
BACKUP_DIR = BASE / "_backups" / datetime.now().strftime("%Y%m%d_%H%M%S")

# === The 22 accessions to REMOVE from training (base IDs, no version) ===
REMOVE_FROM_TRAINING = {
    "OM471829", "OM471830", "OM471832", "OM471833", "OM471834",
    "OM471835", "OM471836", "OM471837", "OM471838",
    "MW718936", "MW718937", "MW718938", "MW718939", "MW718940",
    "MW718941", "MW718942", "MW718943", "MW718944", "MW718945",
    "MW718946", "KM820719", "GU983674",
}

# === File paths ===
TRAIN_XLSX = BASE / "Training_data" / "VP4_training_metadata.xlsx"
TRAIN_FASTA = BASE / "Training_data" / "VP4_training_dataset.fasta"
EVAL_XLSX = BASE / "Evaluation_dataset" / "Eval_metadata_combined.xlsx"
EVAL_FASTA = BASE / "Evaluation_dataset" / "Eval_dataset_comined.fasta"


def backup(path: Path):
    """Back up a file before modifying it."""
    if path.exists():
        dest = BACKUP_DIR / path.relative_to(BASE)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
        print(f"  Backed up: {path.name} -> {dest}")


def base_acc(acc: str) -> str:
    """Strip version suffix: DQ525192.1 -> DQ525192"""
    return acc.split(".")[0].strip()


def main():
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("  APPLYING DATA FIXES")
    print("=" * 70)

    # ------------------------------------------------------------------
    # STEP 1: Fix the combined training FASTA
    #         Remove the 22 porcine/bat sequences
    # ------------------------------------------------------------------
    print("\n--- STEP 1: Fix combined training FASTA ---")
    backup(TRAIN_FASTA)

    train_records = list(SeqIO.parse(TRAIN_FASTA, "fasta"))
    print(f"  Before: {len(train_records)} sequences")

    # Build a lookup: base_accession -> full FASTA ID (with version)
    train_fasta_lookup = {}
    for rec in train_records:
        train_fasta_lookup[base_acc(rec.id)] = rec.id

    kept = []
    removed = []
    for rec in train_records:
        if base_acc(rec.id) in REMOVE_FROM_TRAINING:
            removed.append(rec.id)
        else:
            kept.append(rec)

    print(f"  Removed: {len(removed)} sequences")
    for rid in sorted(removed):
        print(f"    - {rid}")
    print(f"  After: {len(kept)} sequences")

    SeqIO.write(kept, TRAIN_FASTA, "fasta")
    print(f"  Saved: {TRAIN_FASTA.name}")

    # Rebuild the lookup from kept records only
    train_fasta_lookup_clean = {}
    for rec in kept:
        train_fasta_lookup_clean[base_acc(rec.id)] = rec.id

    # ------------------------------------------------------------------
    # STEP 2: Fix accession naming in the training Excel
    #         Add version suffixes to match FASTA IDs
    #         Also remove any of the 22 sequences if present
    # ------------------------------------------------------------------
    print("\n--- STEP 2: Fix combined training Excel ---")
    backup(TRAIN_XLSX)

    train_df = pd.read_excel(TRAIN_XLSX, engine="openpyxl")
    print(f"  Before: {len(train_df)} rows")
    print(f"  Columns: {list(train_df.columns)}")

    # Remove the 22 porcine/bat if any slipped in
    before_len = len(train_df)
    train_df = train_df[~train_df["Accession"].astype(str).str.strip().isin(REMOVE_FROM_TRAINING)].reset_index(drop=True)
    removed_from_xlsx = before_len - len(train_df)
    if removed_from_xlsx > 0:
        print(f"  Removed {removed_from_xlsx} porcine/bat rows from Excel")

    # Add version suffixes to accessions to match FASTA
    matched = 0
    unmatched = []
    new_accessions = []
    for _, row in train_df.iterrows():
        acc = str(row["Accession"]).strip()
        if acc in train_fasta_lookup_clean:
            new_accessions.append(train_fasta_lookup_clean[acc])
            matched += 1
        else:
            # Keep original — these won't match a FASTA record
            new_accessions.append(acc)
            unmatched.append(acc)

    train_df["Accession"] = new_accessions
    print(f"  Accessions matched to FASTA IDs: {matched}")
    if unmatched:
        print(f"  WARNING: {len(unmatched)} accessions had no FASTA match:")
        for u in unmatched[:10]:
            print(f"    - {u}")

    print(f"  After: {len(train_df)} rows")

    # Save as xlsx
    train_df.to_excel(TRAIN_XLSX, index=False, engine="openpyxl")
    print(f"  Saved: {TRAIN_XLSX.name}")

    # ------------------------------------------------------------------
    # STEP 3: Fix accession naming in the evaluation Excel
    #         Add version suffixes to match FASTA IDs
    #         Verify the 22 porcine/bat sequences are present
    # ------------------------------------------------------------------
    print("\n--- STEP 3: Fix combined evaluation Excel ---")
    backup(EVAL_XLSX)

    eval_df = pd.read_excel(EVAL_XLSX, engine="openpyxl")
    print(f"  Before: {len(eval_df)} rows")

    # Build eval FASTA lookup
    eval_records = list(SeqIO.parse(EVAL_FASTA, "fasta"))
    eval_fasta_lookup = {}
    for rec in eval_records:
        eval_fasta_lookup[base_acc(rec.id)] = rec.id

    # Verify the 22 sequences are in eval
    eval_accs = set(eval_df["Accession"].dropna().astype(str).str.strip())
    present = REMOVE_FROM_TRAINING & eval_accs
    missing = REMOVE_FROM_TRAINING - eval_accs
    print(f"  Porcine/bat sequences present in eval Excel: {len(present)}/22")
    if missing:
        print(f"  WARNING: {len(missing)} still missing from eval Excel: {sorted(missing)}")

    # Add version suffixes to eval accessions
    matched_eval = 0
    unmatched_eval = []
    new_eval_accessions = []
    for _, row in eval_df.iterrows():
        acc = str(row["Accession"]).strip()
        if acc in eval_fasta_lookup:
            new_eval_accessions.append(eval_fasta_lookup[acc])
            matched_eval += 1
        else:
            new_eval_accessions.append(acc)
            unmatched_eval.append(acc)

    eval_df["Accession"] = new_eval_accessions
    print(f"  Accessions matched to FASTA IDs: {matched_eval}")
    if unmatched_eval:
        print(f"  WARNING: {len(unmatched_eval)} accessions had no FASTA match:")
        for u in unmatched_eval[:10]:
            print(f"    - {u}")

    # Also fix the label: the audit showed eval porcine sequences were labelled
    # "Negative" — they should be "Intermediate" for evaluation
    porcine_mask = eval_df["adaptation_group"] == "Porcine"
    bovine_mask = eval_df["adaptation_group"] == "Bovine"
    bat_eval_mask = eval_df["adaptation_group"] == "Bat"
    human_int_mask = eval_df["adaptation_group"].isin(["Human_Intermediate", "Human_Eval"])

    wrong_labels = eval_df.loc[
        (porcine_mask | bovine_mask | bat_eval_mask | human_int_mask) &
        (eval_df["label"] != "Intermediate"),
        "label"
    ]
    if len(wrong_labels) > 0:
        print(f"  Fixing {len(wrong_labels)} eval labels from '{wrong_labels.iloc[0]}' to 'Intermediate'")
        eval_df.loc[porcine_mask | bovine_mask | bat_eval_mask | human_int_mask, "label"] = "Intermediate"

    print(f"  After: {len(eval_df)} rows")
    eval_df.to_excel(EVAL_XLSX, index=False, engine="openpyxl")
    print(f"  Saved: {EVAL_XLSX.name}")

    # ------------------------------------------------------------------
    # STEP 4: Verify no training/eval overlap in final files
    # ------------------------------------------------------------------
    print("\n--- STEP 4: Final Verification ---")

    # Re-read the saved files
    final_train_df = pd.read_excel(TRAIN_XLSX, engine="openpyxl")
    final_eval_df = pd.read_excel(EVAL_XLSX, engine="openpyxl")
    final_train_fasta = list(SeqIO.parse(TRAIN_FASTA, "fasta"))
    final_eval_fasta = list(SeqIO.parse(EVAL_FASTA, "fasta"))

    train_xl_accs = set(final_train_df["Accession"].dropna().astype(str).str.strip())
    eval_xl_accs = set(final_eval_df["Accession"].dropna().astype(str).str.strip())
    train_fa_accs = set(rec.id for rec in final_train_fasta)
    eval_fa_accs = set(rec.id for rec in final_eval_fasta)

    # Check Excel-FASTA match within each dataset
    train_xl_base = set(base_acc(a) for a in train_xl_accs)
    train_fa_base = set(base_acc(a) for a in train_fa_accs)
    eval_xl_base = set(base_acc(a) for a in eval_xl_accs)
    eval_fa_base = set(base_acc(a) for a in eval_fa_accs)

    print(f"\n  Training Excel: {len(train_xl_accs)} accessions")
    print(f"  Training FASTA: {len(train_fa_accs)} sequences")
    train_xl_only = train_xl_base - train_fa_base
    train_fa_only = train_fa_base - train_xl_base
    print(f"  Excel↔FASTA match: {len(train_xl_base & train_fa_base)}")
    if train_xl_only:
        print(f"  In Excel but not FASTA: {len(train_xl_only)} -> {sorted(train_xl_only)[:5]}")
    if train_fa_only:
        print(f"  In FASTA but not Excel: {len(train_fa_only)} -> {sorted(train_fa_only)[:5]}")

    print(f"\n  Evaluation Excel: {len(eval_xl_accs)} accessions")
    print(f"  Evaluation FASTA: {len(eval_fa_accs)} sequences")
    eval_xl_only = eval_xl_base - eval_fa_base
    eval_fa_only = eval_fa_base - eval_xl_base
    print(f"  Excel↔FASTA match: {len(eval_xl_base & eval_fa_base)}")
    if eval_xl_only:
        print(f"  In Excel but not FASTA: {len(eval_xl_only)} -> {sorted(eval_xl_only)[:5]}")
    if eval_fa_only:
        print(f"  In FASTA but not Excel: {len(eval_fa_only)} -> {sorted(eval_fa_only)[:5]}")

    # Cross-contamination check
    train_eval_overlap = train_fa_base & eval_fa_base
    print(f"\n  Train↔Eval FASTA overlap: {len(train_eval_overlap)}")
    if train_eval_overlap:
        print(f"  LEAKING accessions: {sorted(train_eval_overlap)[:10]}")
    else:
        print("  ✓ No data leakage between training and evaluation!")

    train_eval_xl_overlap = train_xl_base & eval_xl_base
    print(f"  Train↔Eval Excel overlap: {len(train_eval_xl_overlap)}")
    if train_eval_xl_overlap:
        print(f"  LEAKING accessions: {sorted(train_eval_xl_overlap)[:10]}")
    else:
        print("  ✓ No data leakage between training and evaluation!")

    # Exact accession format match (with version)
    train_exact = train_xl_accs & train_fa_accs
    eval_exact = eval_xl_accs & eval_fa_accs
    print(f"\n  Training exact accession match (with version): {len(train_exact)}/{len(train_xl_accs)}")
    print(f"  Evaluation exact accession match (with version): {len(eval_exact)}/{len(eval_xl_accs)}")

    print("\n" + "=" * 70)
    print("  ALL FIXES APPLIED SUCCESSFULLY")
    print("=" * 70)
    print(f"\n  Backups saved to: {BACKUP_DIR}")


if __name__ == "__main__":
    main()
