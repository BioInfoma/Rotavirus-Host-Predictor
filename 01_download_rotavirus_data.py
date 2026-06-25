#!/usr/bin/env python3
"""Rotavirus A VP4 preprocessing pipeline for ML-ready, publication-grade datasets.

This pipeline uses the aligned per-group CSV/FASTA pairs already present in
`data_vp4/` and performs:
- metadata/FASTA integrity validation
- exact sequence deduplication
- optional redundancy clustering (CD-HIT-like, 99% nt identity)
- ORF detection across frames 1/2/3
- VP8* completeness scoring against VP4 references
- host/genotype conflict flagging
- dataset composition summaries
- a comprehensive QC summary
- export of cleaned metadata, nucleotide FASTA, protein FASTA, and reports
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Config:
    DATA_ROOT = Path("data_vp4")
    TRAINING_SOURCES = [
        (DATA_ROOT / "Training_data" / "human_metadata.csv", DATA_ROOT / "Training_data" / "human_seq.fasta", "training", "Human_Anchor", "Positive"),
        (DATA_ROOT / "Training_data" / "bat_metadata.csv", DATA_ROOT / "Training_data" / "bat_seq.fasta", "training", "Animal_Anchor", "Negative"),
        (DATA_ROOT / "Training_data" / "Avian_metadata.csv", DATA_ROOT / "Training_data" / "Avian_seq.fasta", "training", "Animal_Anchor", "Negative"),
        (DATA_ROOT / "Training_data" / "equine_metadata.csv", DATA_ROOT / "Training_data" / "equine_seq.fasta", "training", "Animal_Anchor", "Negative"),
    ]
    EVAL_SOURCES = [
        (DATA_ROOT / "Evaluation_dataset" / "pig_metadata.csv", DATA_ROOT / "Evaluation_dataset" / "pig_seq.fasta", "evaluation", "Porcine", "Intermediate"),
        (DATA_ROOT / "Evaluation_dataset" / "Bovine_metadata.csv", DATA_ROOT / "Evaluation_dataset" / "Bovine_seq.fasta", "evaluation", "Bovine", "Intermediate"),
        (DATA_ROOT / "Evaluation_dataset" / "Human_eval_metadata.csv", DATA_ROOT / "Evaluation_dataset" / "Human_eval_seq.fasta", "evaluation", "Human_Intermediate", "Intermediate"),
    ]

    OUTPUT_BASE = DATA_ROOT / "analysis_ready"
    TEMP_DIR = OUTPUT_BASE / "tmp"
    CLEAN_METADATA = OUTPUT_BASE / "cleaned_metadata.csv"
    CLEAN_NT_FASTA = OUTPUT_BASE / "cleaned_nucleotide.fasta"
    CLEAN_AA_FASTA = OUTPUT_BASE / "cleaned_protein.fasta"
    QC_SUMMARY = OUTPUT_BASE / "qc_summary.csv"
    INTEGRITY_REPORT = OUTPUT_BASE / "integrity_report.csv"
    VP8_REPORT = OUTPUT_BASE / "vp8_completeness_report.csv"
    DUPLICATE_REPORT = OUTPUT_BASE / "duplicate_qc_report.csv"
    CLUSTER_REPORT = OUTPUT_BASE / "cluster_summary.csv"
    COMPOSITION_DIR = OUTPUT_BASE / "composition"
    LOG_DIR = OUTPUT_BASE / "logs"

    MIN_NT_LENGTH = 800
    MAX_NT_LENGTH = 2600
    MAX_AMBIGUOUS_PCT = 1.0
    VP8_REFERENCE_LENGTH = 272
    VP8_KEEP_THRESHOLD = 0.90
    VP8_REVIEW_THRESHOLD = 0.80
    CLUSTER_IDENTITY = 0.99
    ENABLE_CLUSTERING_DEFAULT = True
    LOG_LEVEL = logging.INFO

    HUMAN_POSITIVE_GENOTYPES = {"P[8]", "P[4]", "P[6]"}
    ANIMAL_POSITIVE_HOSTS = {"Bat", "Avian", "Equine"}
    REVIEW_GENOTYPES = {"P[6]", "P[7]", "P[13]", "P[9]", "P[14]", "P[19]", "P[25]"}


@dataclass
class SequenceCall:
    accession: str
    sequence_nt: str
    selected_frame: int
    protein_sequence: str
    protein_length: int
    stop_count: int
    orf_qc_pass: bool
    orf_reason: str


class Logger:
    def __init__(self, name: str):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(Config.LOG_LEVEL)
        self.logger.handlers.clear()
        formatter = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        file_handler = logging.FileHandler(Config.LOG_DIR / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log")
        stream_handler = logging.StreamHandler()
        file_handler.setFormatter(formatter)
        stream_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        self.logger.addHandler(stream_handler)

    def info(self, msg: str, *args, **kwargs):
        self.logger.info(msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs):
        self.logger.warning(msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs):
        self.logger.error(msg, *args, **kwargs)


class PipelineError(RuntimeError):
    pass


HOST_MAP = {
    "human": "Human",
    "homo sapiens": "Human",
    "pig": "Porcine",
    "pigs": "Porcine",
    "sus scrofa": "Porcine",
    "sus scrofa domesticus": "Porcine",
    "porcine": "Porcine",
    "bovine": "Bovine",
    "bos taurus": "Bovine",
    "cow": "Bovine",
    "cattle": "Bovine",
    "equine": "Equine",
    "equus caballus": "Equine",
    "horse": "Equine",
    "bat": "Bat",
    "chiroptera": "Bat",
    "avian": "Avian",
    "chicken": "Avian",
    "gallus gallus": "Avian",
    "columba": "Avian",
    "columba livia": "Avian",
}


def ensure_output_dirs() -> None:
    for path in [Config.OUTPUT_BASE, Config.COMPOSITION_DIR, Config.LOG_DIR, Config.TEMP_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def read_csv_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def parse_year(value) -> Optional[int]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"unknown", "nan", "none"}:
        return None
    match = re.search(r"\b(19\d{2}|20\d{2})\b", text)
    if match:
        return int(match.group(1))
    try:
        numeric = float(text)
        if numeric > 30000:
            return int((datetime(1899, 12, 30) + pd.to_timedelta(numeric, unit="D")).year)
        if numeric > 1000:
            return int(numeric)
    except Exception:
        pass
    return None


def normalize_host(value: str) -> str:
    if value is None:
        return "Unknown"
    clean = str(value).strip()
    if not clean or clean.lower() in {"unknown", "nan", "none"}:
        return "Unknown"
    lowered = clean.lower()
    for key, mapped in HOST_MAP.items():
        if key in lowered:
            return mapped
    return clean.title()


def normalize_genotype(value: str, description: str = "") -> str:
    text = f"{value or ''} {description or ''}"
    match = re.search(r"P\[(\d+|x|X)\]", text)
    if match:
        return f"P[{match.group(1)}]"
    match = re.search(r"(?:^|[^\w])P(\d+)(?:[^\w]|$)", text)
    if match:
        return f"P[{match.group(1)}]"
    clean = str(value).strip()
    return clean if clean and clean.lower() not in {"nan", "none", "unknown"} else "Unknown"


def standardize_metadata(frame: pd.DataFrame, dataset_split: str, adaptation_group: str, label: str) -> pd.DataFrame:
    rename_map = {
        "Accession": "accession",
        "Full_name": "full_name",
        "Host": "host",
        "Genotype": "genotype",
        "Length": "reported_length",
        "Country": "country",
        "Geo_Location": "geo_location",
        "Collection_Date": "collection_date",
        "Release_Date": "release_date",
    }
    frame = frame.rename(columns=rename_map).copy()
    if "accession" not in frame.columns:
        raise PipelineError("metadata table is missing an Accession column")
    if "full_name" not in frame.columns:
        fallback_columns = [column for column in ["Isolate", "Organism_Name", "Species"] if column in frame.columns]
        if fallback_columns:
            frame["full_name"] = frame[fallback_columns].astype(str).bfill(axis=1).iloc[:, 0]
        else:
            frame["full_name"] = frame["accession"]
    frame["accession"] = frame["accession"].astype(str).str.strip()
    frame["host"] = frame.get("host", pd.Series(["Unknown"] * len(frame))).map(normalize_host)
    frame["genotype"] = frame.apply(lambda row: normalize_genotype(row.get("genotype", ""), row.get("full_name", "")), axis=1)
    if "country" in frame.columns and "geo_location" in frame.columns:
        country_source = frame["country"].fillna("").astype(str)
        geo_source = frame["geo_location"].fillna("").astype(str)
        frame["country"] = country_source.where(country_source.str.strip() != "", geo_source)
    elif "country" not in frame.columns and "geo_location" in frame.columns:
        frame["country"] = frame["geo_location"]
    elif "country" not in frame.columns:
        frame["country"] = "Unknown"
    frame["country"] = frame["country"].astype(str).replace({"nan": "Unknown", "None": "Unknown", "": "Unknown"})
    frame["collection_year"] = frame.get("collection_date", pd.Series([None] * len(frame))).map(parse_year)
    frame["dataset_split"] = dataset_split
    frame["adaptation_group"] = adaptation_group
    frame["label"] = label
    if "reported_length" in frame.columns:
        frame["reported_length"] = pd.to_numeric(frame["reported_length"], errors="coerce")
    frame = frame.dropna(subset=["accession"])
    frame = frame[frame["accession"] != ""]
    frame = frame.drop_duplicates(subset=["accession"], keep="first")
    frame = frame.loc[:, ~frame.columns.duplicated()].copy()
    return frame.reset_index(drop=True)


def read_fasta(path: Path) -> Dict[str, SeqRecord]:
    records: Dict[str, SeqRecord] = {}
    for record in SeqIO.parse(path, "fasta"):
        accession = record.id.split("|")[0].strip()
        records[accession] = record
        records[accession.split(".")[0]] = record
    return records


def normalize_nt_sequence(sequence: str) -> str:
    cleaned = re.sub(r"\s+", "", sequence).upper().replace("U", "T")
    cleaned = cleaned.replace("-", "").replace(".", "")
    cleaned = re.sub(r"[^ATGCN]", "N", cleaned)
    return cleaned


def translate_frame(sequence: str, frame: int) -> Tuple[str, int]:
    trimmed = sequence[frame - 1 :]
    usable = len(trimmed) - (len(trimmed) % 3)
    if usable <= 0:
        return "", 0
    protein = str(Seq(trimmed[:usable]).translate(to_stop=False))
    stop_count = protein.count("*")
    protein = protein.rstrip("*")
    return protein, stop_count


def select_best_orf(sequence: str) -> SequenceCall:
    best: Optional[SequenceCall] = None
    for frame in (1, 2, 3):
        protein, stop_count = translate_frame(sequence, frame)
        protein_length = len(protein.replace("*", ""))
        orf_qc_pass = protein_length > 0 and stop_count <= 1 and "*" not in protein[:-1]
        candidate = SequenceCall(
            accession="",
            sequence_nt=sequence,
            selected_frame=frame,
            protein_sequence=protein,
            protein_length=protein_length,
            stop_count=stop_count,
            orf_qc_pass=orf_qc_pass,
            orf_reason="pass" if orf_qc_pass else "internal_stop_or_short",
        )
        if best is None:
            best = candidate
            continue
        if candidate.protein_length > best.protein_length:
            best = candidate
            continue
        if candidate.protein_length == best.protein_length and candidate.stop_count < best.stop_count:
            best = candidate
    assert best is not None
    return best


def load_reference_vp8_proteins() -> Dict[str, str]:
    reference_candidates = [
        (Config.DATA_ROOT / "Training_data" / "human_seq.fasta", "Wa"),
        (Config.DATA_ROOT / "Evaluation_dataset" / "Bovine_seq.fasta", "DS-1"),
    ]
    references: Dict[str, str] = {}
    for fasta_path, key in reference_candidates:
        if not fasta_path.exists():
            continue
        best_record = None
        for record in SeqIO.parse(fasta_path, "fasta"):
            description = record.description.lower()
            if key.lower() in description:
                best_record = record
                break
        if best_record is None:
            continue
        nt_sequence = normalize_nt_sequence(str(best_record.seq))
        protein = select_best_orf(nt_sequence).protein_sequence
        if protein:
            references[key] = protein
    return references


def assess_vp8_completeness(protein_seq: str, reference_vp4: str) -> Dict[str, object]:
    """Determine whether a sequence contains a near-complete VP8* domain."""
    vp8_start = 1
    vp8_end = 272
    min_vp8_coverage = 0.90

    aln = pairwise2.align.globalms(
        reference_vp4,
        protein_seq,
        2,
        -1,
        -5,
        -0.5,
        one_alignment_only=True,
    )
    if not aln:
        return {"vp8_coverage": 0, "vp8_status": "Discard"}

    ref_aln, query_aln, score, start, end = aln[0]
    ref_pos = 0
    aligned_vp8 = 0
    for ref_char, query_char in zip(ref_aln, query_aln):
        if ref_char != "-":
            ref_pos += 1
        if vp8_start <= ref_pos <= vp8_end and query_char != "-":
            aligned_vp8 += 1

    coverage = aligned_vp8 / (vp8_end - vp8_start + 1)
    if coverage >= min_vp8_coverage:
        status = "Keep"
    elif coverage >= 0.80:
        status = "Review"
    else:
        status = "Discard"
    return {"vp8_coverage": round(coverage, 3), "vp8_status": status}


def load_bundle(metadata_path: Path, fasta_path: Path, dataset_split: str, adaptation_group: str, label: str) -> Tuple[pd.DataFrame, Dict[str, SeqRecord]]:
    metadata = standardize_metadata(read_csv_table(metadata_path), dataset_split, adaptation_group, label)
    fasta = read_fasta(fasta_path)
    return metadata, fasta


def reconcile_metadata_and_fasta(metadata_df: pd.DataFrame, fasta_records: Dict[str, SeqRecord]) -> Tuple[pd.DataFrame, List[SeqRecord], pd.DataFrame]:
    fasta_lookup = {record.id.split("|")[0].split(".")[0]: record for record in fasta_records.values()}
    metadata_ids = set(metadata_df["accession"].astype(str))
    fasta_ids = set(fasta_lookup.keys())
    matched = metadata_ids & fasta_ids
    metadata_only = metadata_ids - fasta_ids
    fasta_only = fasta_ids - metadata_ids

    report_rows = []
    for accession in sorted(metadata_only):
        report_rows.append({"accession": accession, "status": "metadata_only"})
    for accession in sorted(fasta_only):
        report_rows.append({"accession": accession, "status": "fasta_only"})

    reconciled_metadata = metadata_df[metadata_df["accession"].isin(matched)].copy().reset_index(drop=True)
    reconciled_fasta = [fasta_lookup[accession] for accession in sorted(matched)]
    summary = pd.DataFrame([
        {"metric": "total_metadata_records", "count": len(metadata_ids)},
        {"metric": "total_fasta_records", "count": len(fasta_ids)},
        {"metric": "matched_records", "count": len(matched)},
        {"metric": "metadata_only_records", "count": len(metadata_only)},
        {"metric": "fasta_only_records", "count": len(fasta_only)},
    ])
    return reconciled_metadata, reconciled_fasta, pd.DataFrame(report_rows), summary


def sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("utf-8")).hexdigest()


def host_genotype_conflict(host: str, genotype: str) -> Optional[str]:
    host = (host or "Unknown").strip()
    genotype = (genotype or "Unknown").strip()
    if genotype in {"P[6]", "P[7]", "P[13]", "P[9]", "P[14]", "P[19]", "P[25]"} and host in {"Human", "Porcine", "Bat", "Bovine", "Avian", "Equine"}:
        return "Host_Genotype_Ambiguous"
    return None


def vp8_reference_sequences() -> Dict[str, str]:
    references = {
        "Wa": "MKTKILLLLAVALATLSTVDA...",
        "DS-1": "MKTKILLLLAVALATLSTVDA...",
    }
    # If exact reference sequences are unavailable, use the best in-workspace complete VP4-like candidates.
    # These are only used for relative VP8 coverage estimation, not for training labels.
    candidates = {
        "Wa": "",
        "DS-1": "",
    }
    return candidates if any(candidates.values()) else {"Wa": "", "DS-1": ""}


def estimate_vp8_coverage(protein: str) -> Tuple[float, str]:
    if not protein:
        return 0.0, "discard"
    coverage = min(1.0, len(protein) / Config.VP8_REFERENCE_LENGTH)
    if coverage >= Config.VP8_KEEP_THRESHOLD:
        return coverage, "keep"
    if coverage >= Config.VP8_REVIEW_THRESHOLD:
        return coverage, "review"
    return coverage, "discard"


def run_cd_hit(input_fasta: Path, output_prefix: Path, identity: float = 0.99, threads: int = 8) -> None:
    executable = shutil.which("cd-hit")
    if executable is None:
        raise PipelineError("cd-hit is not available on PATH")
    cmd = [executable, "-i", str(input_fasta), "-o", str(output_prefix), "-c", str(identity), "-T", str(threads), "-M", "0"]
    subprocess.run(cmd, check=True)


def parse_cd_hit_clstr(clstr_path: Path) -> pd.DataFrame:
    cluster_id = None
    member_rows = []
    representative = None
    cluster_size = 0
    with clstr_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id = line.replace(">", "").replace(" ", "_")
                representative = None
                cluster_size = 0
                continue
            match = re.search(r">([^\.\s>]+(?:\.[^\.\s>]+)?)\.\.\.", line)
            if not match:
                match = re.search(r">([^\s>]+)", line)
            accession = match.group(1) if match else ""
            cluster_size += 1
            if line.endswith("*"):
                representative = accession
            member_rows.append({"accession": accession, "cluster_id": cluster_id, "cluster_size": None, "representative_accession": representative, "is_representative": line.endswith("*")})
    cluster_df = pd.DataFrame(member_rows)
    if not cluster_df.empty:
        cluster_df["cluster_size"] = cluster_df.groupby("cluster_id")["accession"].transform("count")
        reps = cluster_df.loc[cluster_df["is_representative"], ["cluster_id", "accession"]].rename(columns={"accession": "representative_accession"})
        cluster_df = cluster_df.drop(columns=["is_representative"]).merge(reps, on="cluster_id", how="left", suffixes=("", "_from_star"))
        cluster_df["representative_accession"] = cluster_df["representative_accession_from_star"].fillna(cluster_df["representative_accession"])
        cluster_df = cluster_df.drop(columns=["representative_accession_from_star"])
    return cluster_df


def greedy_cluster_fallback(records: pd.DataFrame, identity: float = Config.CLUSTER_IDENTITY) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if records.empty:
        return records.copy(), pd.DataFrame(columns=["accession", "cluster_id", "cluster_size", "representative_accession"])
    clusters: List[Dict[str, object]] = []
    used = set()
    cluster_index = 0
    sequences = records.sort_values(["sequence_nt_length", "accession"], ascending=[False, True]).reset_index(drop=True)
    for idx, row in sequences.iterrows():
        accession = row["accession"]
        if accession in used:
            continue
        cluster_index += 1
        cluster_id = f"C{cluster_index:05d}"
        representative_accession = accession
        member_indices = [idx]
        used.add(accession)
        rep_seq = row["sequence_nt"]
        for jdx in range(idx + 1, len(sequences)):
            candidate = sequences.iloc[jdx]
            candidate_acc = candidate["accession"]
            if candidate_acc in used:
                continue
            alignment = pairwise2.align.globalxx(rep_seq, candidate["sequence_nt"], one_alignment_only=True)
            if not alignment:
                continue
            aln = alignment[0]
            aligned_len = max(len(aln.seqA), len(aln.seqB))
            matches = sum(1 for a, b in zip(aln.seqA, aln.seqB) if a == b and a != "-")
            if aligned_len and (matches / aligned_len) >= identity:
                member_indices.append(jdx)
                used.add(candidate_acc)
        for member_index in member_indices:
            clusters.append({"accession": sequences.iloc[member_index]["accession"], "cluster_id": cluster_id, "cluster_size": len(member_indices), "representative_accession": representative_accession})
    cluster_df = pd.DataFrame(clusters)
    reduced = records.merge(cluster_df, on="accession", how="left")
    reduced = reduced[reduced["representative_accession"] == reduced["accession"]].copy().reset_index(drop=True)
    return reduced, cluster_df


def cluster_sequences(records: pd.DataFrame, identity: float = Config.CLUSTER_IDENTITY, threads: int = 8) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if records.empty:
        return records.copy(), pd.DataFrame(columns=["accession", "cluster_id", "cluster_size", "representative_accession"])

    input_fasta = Config.TEMP_DIR / "cdhit_input.fasta"
    output_prefix = Config.TEMP_DIR / "cdhit_output"
    SeqIO.write(
        [SeqRecord(Seq(row.sequence_nt), id=row.accession, description="") for row in records.itertuples(index=False)],
        input_fasta,
        "fasta",
    )
    try:
        run_cd_hit(input_fasta, output_prefix, identity=identity, threads=threads)
        cluster_df = parse_cd_hit_clstr(output_prefix.with_suffix(".clstr"))
        if cluster_df.empty:
            raise PipelineError("cd-hit produced no cluster assignments")
        reduced = records.merge(cluster_df[["accession", "cluster_id", "cluster_size", "representative_accession"]], on="accession", how="left")
        reduced = reduced[reduced["representative_accession"] == reduced["accession"]].copy().reset_index(drop=True)
        return reduced, cluster_df[["accession", "cluster_id", "cluster_size", "representative_accession"]]
    except Exception as exc:
        logging.getLogger("VP4Pipeline").warning(f"CD-HIT unavailable or failed ({exc}); using greedy fallback clustering.")
        return greedy_cluster_fallback(records, identity=identity)


def dataset_composition_tables(frame: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    def table(column: str, name: str) -> pd.DataFrame:
        counts = frame[column].fillna("Unknown").value_counts().reset_index()
        counts.columns = [name, "count"]
        return counts

    return {
        "host": table("host", "Host"),
        "genotype": table("genotype", "Genotype"),
        "label": table("label", "Label"),
        "adaptation_group": table("adaptation_group", "Adaptation_Group"),
    }


def find_best_wa_reference(records: pd.DataFrame) -> Optional[pd.Series]:
    wa = records[records["full_name"].astype(str).str.contains("Wa", case=False, na=False)]
    if wa.empty:
        return None
    complete = wa[wa["sequence_nt_length"] >= 2200]
    if complete.empty:
        complete = wa
    return complete.sort_values(["sequence_nt_length", "accession"], ascending=[False, True]).iloc[0]


def build_host_label(adaptation_group: str, label: str, host: str) -> Tuple[str, str]:
    return adaptation_group, label


class Pipeline:
    def __init__(self, enable_clustering: bool = Config.ENABLE_CLUSTERING_DEFAULT):
        ensure_output_dirs()
        self.logger = Logger("VP4Pipeline")
        self.enable_clustering = enable_clustering
        self.metadata_frames: List[pd.DataFrame] = []
        self.sequence_frames: List[pd.DataFrame] = []
        self.qc_reports: Dict[str, pd.DataFrame] = {}
        self.qc_counts = Counter()
        self.final_records = pd.DataFrame()
        self.reference_vp8 = load_reference_vp8_proteins()
        if not self.reference_vp8:
            self.logger.warning("No in-workspace VP8 reference proteins were found; completeness scoring will fall back to length-based estimates.")

    def calculate_vp8_coverage(self, protein: str) -> Tuple[float, str, str]:
        if not protein:
            return 0.0, "discard", "none"

        if not self.reference_vp8:
            coverage, status = estimate_vp8_coverage(protein)
            return coverage, status, "length_proxy"

        best_coverage = 0.0
        best_status = "discard"
        best_ref = "none"
        for ref_name, ref_protein in self.reference_vp8.items():
            alignment = pairwise2.align.localms(ref_protein, protein, 2, -1, -2, -0.5, one_alignment_only=True)
            if not alignment:
                continue
            aln = alignment[0]
            aligned_ref_len = sum(1 for char in aln.seqA if char != "-")
            coverage = aligned_ref_len / Config.VP8_REFERENCE_LENGTH
            if coverage >= Config.VP8_KEEP_THRESHOLD:
                status = "keep"
            elif coverage >= Config.VP8_REVIEW_THRESHOLD:
                status = "review"
            else:
                status = "discard"
            if coverage > best_coverage:
                best_coverage = coverage
                best_status = status
                best_ref = ref_name
        return best_coverage, best_status, best_ref

    def load_sources(self):
        for metadata_path, fasta_path, dataset_split, adaptation_group, label in Config.TRAINING_SOURCES + Config.EVAL_SOURCES:
            self.logger.info(f"Loading {metadata_path.name} and {fasta_path.name}")
            metadata, fasta = load_bundle(metadata_path, fasta_path, dataset_split, adaptation_group, label)
            matched_metadata, matched_fasta, integrity_report, integrity_summary = reconcile_metadata_and_fasta(metadata, fasta)
            integrity_summary["source"] = metadata_path.name
            integrity_report["source"] = metadata_path.name
            self.qc_reports[f"integrity_{metadata_path.stem}"] = integrity_summary
            self.qc_reports[f"reconcile_{metadata_path.stem}"] = integrity_report
            self.logger.info(
                f"  total metadata={int(integrity_summary.loc[integrity_summary['metric']=='total_metadata_records', 'count'].iloc[0])}, "
                f"total fasta={int(integrity_summary.loc[integrity_summary['metric']=='total_fasta_records', 'count'].iloc[0])}, "
                f"matched={int(integrity_summary.loc[integrity_summary['metric']=='matched_records', 'count'].iloc[0])}, "
                f"metadata-only={int(integrity_summary.loc[integrity_summary['metric']=='metadata_only_records', 'count'].iloc[0])}, "
                f"fasta-only={int(integrity_summary.loc[integrity_summary['metric']=='fasta_only_records', 'count'].iloc[0])}"
            )
            self.qc_counts["metadata_only_removed"] += int(integrity_summary.loc[integrity_summary['metric']=='metadata_only_records', 'count'].iloc[0])
            self.qc_counts["fasta_only_removed"] += int(integrity_summary.loc[integrity_summary['metric']=='fasta_only_records', 'count'].iloc[0])
            self.metadata_frames.append(matched_metadata)
            self.sequence_frames.append(self.build_sequence_frame(matched_metadata, fasta, metadata_path.stem))

    def build_sequence_frame(self, metadata: pd.DataFrame, fasta: Dict[str, SeqRecord], source_name: str) -> pd.DataFrame:
        rows: List[Dict[str, object]] = []
        duplicate_tracker: Dict[str, str] = {}
        duplicate_rows: List[Dict[str, object]] = []
        for row in metadata.to_dict(orient="records"):
            accession = row["accession"]
            record = fasta.get(accession) or fasta.get(accession.split(".")[0])
            if record is None:
                continue
            nt_sequence = normalize_nt_sequence(str(record.seq))
            nt_length = len(nt_sequence)
            ambiguous_pct = 100 * sum(1 for base in nt_sequence if base not in {"A", "T", "G", "C"}) / nt_length if nt_length else 0.0
            length_fail = not (Config.MIN_NT_LENGTH <= nt_length <= Config.MAX_NT_LENGTH)
            ambiguity_fail = ambiguous_pct > Config.MAX_AMBIGUOUS_PCT
            if length_fail:
                self.qc_counts["length_filter_failures"] += 1
                continue
            if ambiguity_fail:
                self.qc_counts["ambiguity_filter_failures"] += 1
                continue

            selected = select_best_orf(nt_sequence)
            selected.accession = accession
            if not selected.orf_qc_pass:
                self.qc_counts["orf_failures"] += 1
                continue

            protein_seq = selected.protein_sequence
            protein_length = selected.protein_length
            vp8_coverage, vp8_status, vp8_reference = self.calculate_vp8_coverage(protein_seq)
            if vp8_status == "discard":
                self.qc_counts["vp8_completeness_failures"] += 1
                continue
            if vp8_status == "review":
                self.qc_counts["vp8_review_records"] += 1

            md5 = sequence_md5(nt_sequence)
            if md5 in duplicate_tracker:
                duplicate_rows.append(
                    {
                        "accession": accession,
                        "duplicate_of": duplicate_tracker[md5],
                        "sequence_md5": md5,
                        "source": source_name,
                    }
                )
                self.qc_counts["exact_duplicates_removed"] += 1
                continue
            duplicate_tracker[md5] = accession

            conflict_flag = host_genotype_conflict(row.get("host", "Unknown"), row.get("genotype", "Unknown"))
            rows.append(
                {
                    **row,
                    "sequence_nt": nt_sequence,
                    "sequence_nt_length": nt_length,
                    "sequence_md5": md5,
                    "seq_hash": md5,
                    "selected_frame": selected.selected_frame,
                    "protein_sequence": protein_seq,
                    "protein_length": protein_length,
                    "vp8_coverage": round(vp8_coverage, 4),
                    "vp8_status": vp8_status,
                    "vp8_reference": vp8_reference,
                    "conflict_flag": conflict_flag or "None",
                    "source_file": source_name,
                }
            )

        self.qc_reports[f"duplicates_{source_name}"] = pd.DataFrame(duplicate_rows)
        return pd.DataFrame(rows)

    def apply_sequence_deduplication(self, frame: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if frame.empty:
            return frame, pd.DataFrame(columns=["accession", "duplicate_of", "sequence_md5"])
        dedup_rows: List[Dict[str, object]] = []
        kept_rows: List[pd.Series] = []
        seen: Dict[str, str] = {}
        for _, row in frame.iterrows():
            md5 = row["sequence_md5"]
            if md5 in seen:
                dedup_rows.append(
                    {
                        "accession": row["accession"],
                        "duplicate_of": seen[md5],
                        "sequence_md5": md5,
                    }
                )
                self.qc_counts["exact_duplicates_removed"] += 1
                continue
            seen[md5] = row["accession"]
            kept_rows.append(row)
        return pd.DataFrame(kept_rows).reset_index(drop=True), pd.DataFrame(dedup_rows)

    def run(self):
        self.logger.info("\n" + "=" * 88)
        self.logger.info("Rotavirus A VP4 preprocessing pipeline")
        self.logger.info("=" * 88)
        self.load_sources()

        metadata = pd.concat(self.metadata_frames, ignore_index=True).drop_duplicates(subset=["accession"], keep="first").reset_index(drop=True)
        sequences = pd.concat(self.sequence_frames, ignore_index=True).drop_duplicates(subset=["accession"], keep="first").reset_index(drop=True)

        matched_accessions = set(metadata["accession"]) & set(sequences["accession"])
        metadata = metadata[metadata["accession"].isin(matched_accessions)].reset_index(drop=True)
        sequences = sequences[sequences["accession"].isin(matched_accessions)].reset_index(drop=True)

        sequence_deduped, duplicate_report = self.apply_sequence_deduplication(sequences)
        self.qc_reports["duplicate_report"] = duplicate_report
        self.qc_counts["initial_sequence_count"] = len(sequences)
        self.qc_counts["duplicate_removals"] = len(duplicate_report)

        if self.enable_clustering and not sequence_deduped.empty:
            clustered, cluster_membership = cluster_sequences(sequence_deduped, identity=Config.CLUSTER_IDENTITY)
            self.qc_reports["cluster_summary"] = cluster_membership.groupby(["cluster_id", "representative_accession"]).size().reset_index(name="cluster_size")
            self.qc_counts["redundancy_cluster_removals"] = len(sequence_deduped) - len(clustered)
            sequence_deduped = clustered
        else:
            self.qc_reports["cluster_summary"] = pd.DataFrame(columns=["cluster_id", "cluster_size", "representative_accession"])
            self.qc_counts["redundancy_cluster_removals"] = 0

        merged = metadata.merge(
            sequence_deduped[["accession", "sequence_nt", "sequence_nt_length", "sequence_md5", "selected_frame", "protein_sequence", "protein_length", "vp8_coverage", "vp8_status", "vp8_reference", "conflict_flag", "source_file"]],
            on="accession",
            how="inner",
        )

        merged["collection_year"] = merged["collection_year"].fillna(pd.NA)
        merged["seq_hash"] = merged["sequence_md5"]
        merged["protein_md5"] = merged["protein_sequence"].map(lambda x: hashlib.md5(str(x).encode("utf-8")).hexdigest())
        merged["source_record_type"] = "cleaned"

        merged = merged[merged["vp8_status"].isin(["keep", "review"])].reset_index(drop=True)
        merged = merged.drop_duplicates(subset=["sequence_md5"], keep="first").reset_index(drop=True)

        # Label harmonization for downstream ML.
        merged["label"] = merged["label"].replace({"Intermediate": "Intermediate", "Positive": "Positive", "Negative": "Negative"})

        self.final_records = merged.copy()
        self.qc_counts["final_retained_sequences"] = len(self.final_records)
        self.qc_counts["vp8_completeness_failures"] += int((sequence_deduped["vp8_status"] == "discard").sum())

        self.write_outputs(metadata, sequence_deduped, duplicate_report)
        self.write_reports(metadata, sequence_deduped)
        self.write_qc_summary()
        self.log_console_summary()

    def write_outputs(self, metadata: pd.DataFrame, sequence_frame: pd.DataFrame, duplicate_report: pd.DataFrame):
        cleaned_metadata = self.final_records[
            [
                "accession",
                "full_name",
                "host",
                "genotype",
                "label",
                "adaptation_group",
                "collection_year",
                "country",
                "dataset_split",
                "selected_frame",
                "protein_length",
                "vp8_coverage",
                "vp8_status",
                "vp8_reference",
                "conflict_flag",
                "sequence_md5",
                "seq_hash",
                "protein_md5",
                "source_file",
            ]
        ].copy()
        cleaned_metadata.to_csv(Config.CLEAN_METADATA, index=False)

        nt_records = [SeqRecord(Seq(row.sequence_nt), id=row.accession, description=f"{row.host} | {row.genotype} | {row.vp8_status}") for row in self.final_records.itertuples(index=False)]
        aa_records = [SeqRecord(Seq(row.protein_sequence), id=row.accession, description=f"{row.host} | {row.genotype} | frame={row.selected_frame}") for row in self.final_records.itertuples(index=False)]
        SeqIO.write(nt_records, Config.CLEAN_NT_FASTA, "fasta")
        SeqIO.write(aa_records, Config.CLEAN_AA_FASTA, "fasta")
        duplicate_report.to_csv(Config.DUPLICATE_REPORT, index=False)

    def write_reports(self, metadata: pd.DataFrame, sequence_frame: pd.DataFrame):
        if not self.final_records.empty:
            vp8_report = self.final_records[["accession", "genotype", "host", "vp8_coverage", "vp8_status", "selected_frame", "protein_length", "conflict_flag"]].copy()
        else:
            vp8_report = pd.DataFrame(columns=["accession", "genotype", "host", "vp8_coverage", "vp8_status", "selected_frame", "protein_length", "conflict_flag"])
        vp8_report.to_csv(Config.VP8_REPORT, index=False)

        comp = dataset_composition_tables(self.final_records)
        for name, table in comp.items():
            table.to_csv(Config.COMPOSITION_DIR / f"{name}_distribution.csv", index=False)

        conflict_report = self.final_records[self.final_records["conflict_flag"] != "None"].copy()
        conflict_report.to_csv(Config.OUTPUT_BASE / "host_genotype_conflicts.csv", index=False)

        cluster_summary = self.qc_reports.get("cluster_summary", pd.DataFrame())
        cluster_summary.to_csv(Config.CLUSTER_REPORT, index=False)

        integrity_frames = []
        reconcile_frames = []
        for name, report in self.qc_reports.items():
            if name.startswith("integrity_"):
                tmp = report.copy()
                if "source" in tmp.columns:
                    tmp["source"] = tmp["source"].fillna(name.replace("integrity_", ""))
                else:
                    tmp.insert(0, "source", name.replace("integrity_", ""))
                integrity_frames.append(tmp)
            if name.startswith("reconcile_"):
                tmp = report.copy()
                if "source" not in tmp.columns:
                    tmp.insert(0, "source", name.replace("reconcile_", ""))
                reconcile_frames.append(tmp)
        if integrity_frames:
            pd.concat(integrity_frames, ignore_index=True).to_csv(Config.OUTPUT_BASE / "integrity_summary.csv", index=False)
        if reconcile_frames:
            pd.concat(reconcile_frames, ignore_index=True).to_csv(Config.INTEGRITY_REPORT, index=False)

    def write_qc_summary(self):
        summary = pd.DataFrame([
            {
                "initial_sequences": int(self.qc_counts.get("initial_sequence_count", 0)),
                "metadata_only_removed": int(self.qc_counts.get("metadata_only_removed", 0)),
                "fasta_only_removed": int(self.qc_counts.get("fasta_only_removed", 0)),
                "length_failures": int(self.qc_counts.get("length_filter_failures", 0)),
                "ambiguity_failures": int(self.qc_counts.get("ambiguity_filter_failures", 0)),
                "orf_failures": int(self.qc_counts.get("orf_failures", 0)),
                "vp8_failures": int(self.qc_counts.get("vp8_completeness_failures", 0)),
                "duplicate_removals": int(self.qc_counts.get("duplicate_removals", 0)),
                "cluster_removals": int(self.qc_counts.get("redundancy_cluster_removals", 0)),
                "final_sequences": int(self.qc_counts.get("final_retained_sequences", 0)),
            }
        ])
        summary.to_csv(Config.QC_SUMMARY, index=False)

    def log_console_summary(self):
        self.logger.info("\nQC summary")
        self.logger.info("----------")
        qc = pd.read_csv(Config.QC_SUMMARY)
        row = qc.iloc[0]
        for column in qc.columns:
            self.logger.info(f"{column}: {row[column]}")
        self.logger.info("\nDataset composition summaries written to %s", Config.COMPOSITION_DIR)
        self.logger.info("Cleaned metadata: %s", Config.CLEAN_METADATA)
        self.logger.info("Cleaned nucleotide FASTA: %s", Config.CLEAN_NT_FASTA)
        self.logger.info("Cleaned protein FASTA: %s", Config.CLEAN_AA_FASTA)
        self.logger.info("VP8 report: %s", Config.VP8_REPORT)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Rotavirus A VP4 preprocessing pipeline")
    parser.add_argument("--no-clustering", action="store_true", help="Disable 99% redundancy clustering")
    args = parser.parse_args(argv)

    pipeline = Pipeline(enable_clustering=not args.no_clustering)
    pipeline.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
