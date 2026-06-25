import pandas as pd
from Bio import SeqIO

# Check what host/species these 21 accessions have in human_metadata.csv
df = pd.read_csv(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\human_metadata.csv")
problem_accs = ["KM820719","MW718936","MW718937","MW718938","MW718939","MW718940","MW718941","MW718942","MW718943","MW718944","MW718945","MW718946","OM471829","OM471830","OM471832","OM471833","OM471834","OM471835","OM471836","OM471837","OM471838"]
problem = df[df["Accession"].isin(problem_accs)]
print("=== These 21 accessions in human_metadata.csv ===")
for _, r in problem.iterrows():
    acc = r["Accession"]
    host = r.get("Host", "?")
    species = r.get("Species", "?")
    genotype = r.get("Genotype", "?")
    print(f"  {acc}  Host={host}  Species={species}  Genotype={genotype}")

# Check GU983674 in bat_metadata.csv
bat = pd.read_csv(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\bat_metadata.csv")
gu = bat[bat["Accession"] == "GU983674"]
print("\n=== GU983674 in bat_metadata.csv ===")
for _, r in gu.iterrows():
    print(f"  {r['Accession']}  Host={r.get('Host','?')}  Species={r.get('Species','?')}")

# Where are the 22 accessions in FASTA files?
all_22 = set(problem_accs + ["GU983674"])

# Training combined FASTA
train_fasta_ids = set()
for rec in SeqIO.parse(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\VP4_training_dataset.fasta", "fasta"):
    train_fasta_ids.add(rec.id.split(".")[0])

# Eval combined FASTA
eval_fasta_ids = set()
for rec in SeqIO.parse(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Evaluation_dataset\Eval_dataset_comined.fasta", "fasta"):
    eval_fasta_ids.add(rec.id.split(".")[0])

in_train = all_22 & train_fasta_ids
in_eval = all_22 & eval_fasta_ids

print(f"\n=== Where are the 22 missing accessions in combined FASTAs? ===")
print(f"  In Training FASTA: {len(in_train)}")
if in_train:
    print(f"    {sorted(in_train)[:10]}")
print(f"  In Eval FASTA: {len(in_eval)}")
if in_eval:
    print(f"    {sorted(in_eval)[:10]}")
print(f"  In Neither: {len(all_22 - in_train - in_eval)}")

# Check individual FASTAs for these
for name, path in [
    ("human_seq.fasta", r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\human_seq.fasta"),
    ("bat_seq.fasta", r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\bat_seq.fasta"),
    ("pig_seq.fasta", r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Evaluation_dataset\pig_seq.fasta"),
]:
    ids = set()
    for rec in SeqIO.parse(path, "fasta"):
        ids.add(rec.id.split(".")[0])
    found = all_22 & ids
    if found:
        print(f"  In {name}: {sorted(found)}")

# Cross-contamination check: are any training CSV accessions also in eval Excel?
train_csvs_all = set()
for path in [
    r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\human_metadata.csv",
    r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\bat_metadata.csv",
    r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\Avian_metadata.csv",
    r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\equine_metadata.csv",
]:
    d = pd.read_csv(path)
    train_csvs_all.update(d["Accession"].astype(str).str.strip())

eval_xlsx = pd.read_excel(r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Evaluation_dataset\Eval_metadata_combined.xlsx", engine="openpyxl")
eval_xlsx_acc = set(eval_xlsx["Accession"].dropna().astype(str).str.strip())

cross = train_csvs_all & eval_xlsx_acc
print(f"\n=== Cross-contamination: accessions in BOTH training CSVs AND eval Excel ===")
print(f"  Count: {len(cross)}")
if cross:
    for acc in sorted(cross):
        # Which training CSV?
        for name, path in [
            ("human", r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\human_metadata.csv"),
            ("bat", r"C:\Users\USER\Documents\persornal_project\Zoonotic_ML\data_vp4\Training_data\bat_metadata.csv"),
        ]:
            d = pd.read_csv(path)
            if acc in d["Accession"].values:
                # Get host from eval Excel
                eval_row = eval_xlsx[eval_xlsx["Accession"] == acc]
                eval_host = eval_row["Host"].values[0] if len(eval_row) > 0 else "?"
                eval_group = eval_row["adaptation_group"].values[0] if len(eval_row) > 0 else "?"
                train_row = d[d["Accession"] == acc]
                train_host = train_row["Host"].values[0] if len(train_row) > 0 else "?"
                print(f"  {acc}: training CSV={name} (Host={train_host}), eval Excel (Host={eval_host}, Group={eval_group})")
