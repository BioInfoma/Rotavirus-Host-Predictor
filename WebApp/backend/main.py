import os
import sys
import logging
import shutil
import subprocess
import tempfile
import re
import json
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from fastapi import FastAPI, UploadFile, File, HTTPException, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel
import pandas as pd
import numpy as np
import xgboost as xgb
from Bio import SeqIO
from io import StringIO
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

# Add project root to path
backend_dir = Path(__file__).resolve().parent
project_root = backend_dir.parent.parent
sys.path.append(str(project_root))
from Scripts.config.settings import CONFIG

app = FastAPI(title="Rotavirus VP4 Host Adaptation API")

# Enable CORS for Vite frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

model_state = {}

@app.on_event("startup")
def load_models():
    import torch
    from transformers import AutoTokenizer, AutoModel
    
    logger.info("Starting API server, loading models...")
    # 1. Load ESM-2
    model_name = CONFIG['esm_model_name']
    model_state['tokenizer'] = AutoTokenizer.from_pretrained(model_name)
    model_state['esm2'] = AutoModel.from_pretrained(model_name)
    model_state['device'] = torch.device("cpu")
    model_state['esm2'] = model_state['esm2'].to(model_state['device'])
    model_state['esm2'].eval()
    logger.info("ESM-2 loaded.")
    
    # 2. Fit CountVectorizer
    train_fasta = Path(CONFIG['analysis_dir']) / "cleaned_vp8_training.fasta"
    train_seqs = []
    for record in SeqIO.parse(train_fasta, "fasta"):
        train_seqs.append(str(record.seq).upper().replace("-", ""))
    
    k_sizes = CONFIG['kmer_sizes']
    vectorizer = CountVectorizer(analyzer='char', ngram_range=(min(k_sizes), max(k_sizes)), lowercase=False)
    vectorizer.fit(train_seqs)
    model_state['vectorizer'] = vectorizer
    model_state['kmer_features'] = [f"kmer_{feature}" for feature in vectorizer.get_feature_names_out()]
    logger.info("K-mer vectorizer fitted.")
    
    # 3. Load XGBoost Model
    xgb_model = xgb.XGBClassifier()
    model_path = Path(CONFIG['analysis_dir']) / "models" / "xgboost_vp4_model.json"
    xgb_model.load_model(model_path)
    model_state['xgb'] = xgb_model
    logger.info("XGBoost model loaded.")

    # 4. Load VP8 reference for alignment
    wa_ref_path = Path(CONFIG['wa_reference'])
    record = next(SeqIO.parse(wa_ref_path, "fasta"))
    # The reference is full length, but settings say VP8 is aa 1-272.
    # Actually Wa reference is DNA or Protein? Let's check. 
    # extract_features.py uses the VP8 reference used during preprocessing, but actually during extraction from scratch
    # we need the VP8 sequence.
    # The preprocessing pipeline extracts VP8 using Smith-Waterman.
    # For inference, if we upload protein, we align. Wait, if it's protein we use Wa_P8_VP4_protein.fasta
    model_state['vp8_ref_seq'] = str(record.seq)[CONFIG['vp8_aa_start']-1:CONFIG['vp8_aa_end']]
    logger.info("Initialization complete.")

def get_tool_path(tool_name: str) -> str:
    # 1. Check system PATH
    which_path = shutil.which(tool_name)
    if which_path:
        return which_path
    
    # 2. Check standard alias
    if tool_name == "iqtree2":
        which_path = shutil.which("iqtree")
        if which_path:
            return which_path

    # 3. Check local bin directory inside project
    local_bin = project_root / "bin"
    if local_bin.exists():
        if tool_name == "mafft":
            for p in local_bin.glob("**/mafft.bat"):
                return str(p)
            for p in local_bin.glob("**/mafft.exe"):
                return str(p)
        elif tool_name in ["iqtree", "iqtree2"]:
            for p in local_bin.glob("**/iqtree2.exe"):
                return str(p)
            for p in local_bin.glob("**/iqtree.exe"):
                return str(p)
                
    return tool_name  # Fallback to name and let OS try to run it

def clean_header(header: str) -> str:
    return re.sub(r'[^a-zA-Z0-9_]', '_', header)

def get_pdb_sequence(pdb_path: Path) -> Tuple[str, List[int]]:
    res_names = []
    res_nums = []
    seen_residues = set()
    
    d3to1 = {
        'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
        'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
        'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
        'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
        'ASX':'B', 'GLX':'Z', 'CSO':'C', 'HIP':'H', 'MSE':'M'
    }
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                chain = line[21]
                if chain != 'A':
                    continue
                res_name = line[17:20].strip()
                res_num = int(line[22:26].strip())
                
                res_id = (chain, res_num)
                if res_id not in seen_residues:
                    seen_residues.add(res_id)
                    res_names.append(d3to1.get(res_name, 'X'))
                    res_nums.append(res_num)
                    
    return "".join(res_names), res_nums

def is_nucleotide(seq: str) -> bool:
    valid_chars = set("ACGTNU")
    seq_upper = seq.upper()
    matches = sum(1 for c in seq_upper if c in valid_chars)
    return (matches / len(seq_upper)) > 0.85

def extract_vp8(sequence: str):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    
    ref_seq = model_state['vp8_ref_seq']
    
    candidate_seqs = [sequence]
    if is_nucleotide(sequence):
        from Bio.Seq import Seq
        seq_obj = Seq(sequence)
        candidate_seqs = []
        for i in range(3):
            candidate_seqs.append(str(seq_obj[i:].translate(to_stop=False)).replace("*", "X"))
        rev_seq = seq_obj.reverse_complement()
        for i in range(3):
            candidate_seqs.append(str(rev_seq[i:].translate(to_stop=False)).replace("*", "X"))
            
    best_score = -1
    best_extraction = None
    
    for c_seq in candidate_seqs:
        alignments = aligner.align(ref_seq, c_seq)
        if alignments:
            best_al = alignments[0]
            if best_al.score > best_score:
                best_score = best_al.score
                target_start = best_al.aligned[1][0][0]
                target_end = best_al.aligned[1][-1][-1]
                extracted_seq = c_seq[target_start:target_end].replace("-", "")
                
                ref_aligned = str(best_al[0])
                query_aligned = str(best_al[1])
                
                best_extraction = (extracted_seq, target_start, target_end, ref_aligned, query_aligned)
                
    if not best_extraction or best_score < 50:
        raise ValueError("Could not align sequence to VP8 reference. Are you sure this is a Rotavirus VP4 sequence?")
        
    return best_extraction

def extract_esm2(sequence: str):
    import torch
    tokenizer = model_state['tokenizer']
    model = model_state['esm2']
    device = model_state['device']
    
    inputs = tokenizer(sequence.upper(), return_tensors="pt")
    inputs = {k: v.to(device) for k, v in inputs.items()}
    
    with torch.no_grad():
        outputs = model(**inputs)
        
    last_hidden_state = outputs.last_hidden_state[0]
    residue_states = last_hidden_state[1:-1, :]
    mean_embedding = residue_states.mean(dim=0).cpu().numpy()
    
    cols = [f"esm_dim_{i}" for i in range(len(mean_embedding))]
    return pd.DataFrame([mean_embedding], columns=cols)

def extract_kmer(sequence: str):
    vectorizer = model_state['vectorizer']
    features = model_state['kmer_features']
    
    counts = vectorizer.transform([sequence]).toarray()
    freqs = normalize(counts, norm='l1')
    return pd.DataFrame(freqs, columns=features)

class ShapFeature(BaseModel):
    feature: str
    impact: float

class PredictionResponse(BaseModel):
    accession: str
    zoonotic_potential: float
    vp8_length: int
    alignment_start: int
    alignment_end: int
    is_human_adapted: bool
    ref_aligned: str
    query_aligned: str
    top_human_features: List[ShapFeature]
    top_animal_features: List[ShapFeature]
    interpretation: str

class PhyloNodePrediction(BaseModel):
    score: float
    is_human_adapted: bool
    vp8_length: int
    alignment_start: int
    alignment_end: int
    interpretation: str

class PhyloResponse(BaseModel):
    newick: str
    predictions: Dict[str, PhyloNodePrediction]
    model_used: str
    alignment_length: int
    num_sequences: int

class StructureResidue(BaseModel):
    pdb_res_num: int
    input_res_num: Optional[int] = None
    input_aa: Optional[str] = None
    shap_val: float
    is_mutation: bool

class StructureResponse(BaseModel):
    pdb_content: str
    residues: List[StructureResidue]
    zoonotic_potential: float
    is_human_adapted: bool

def generate_interpretation(prob: float, is_human: bool, sorted_shap: list, vp8_seq: str) -> str:
    """Translate SHAP values into a plain-English biological interpretation."""
    direction = "human-adapted" if is_human else "animal-adapted"
    confidence = "high" if abs(prob - 0.5) > 0.35 else "moderate" if abs(prob - 0.5) > 0.15 else "borderline"
    
    # Separate k-mer features (biologically named) from ESM-2 features
    top_kmers_human = []
    top_kmers_animal = []
    esm_human_count = 0
    esm_animal_count = 0
    
    for feat_name, shap_val in sorted_shap:
        if feat_name.startswith("kmer_"):
            motif = feat_name.replace("kmer_", "")
            # Format as hyphenated single-letter amino acids: GSE -> G-S-E
            motif_display = "-".join(list(motif))
            if shap_val > 0 and len(top_kmers_human) < 3:
                top_kmers_human.append((motif_display, motif, shap_val))
            elif shap_val < 0 and len(top_kmers_animal) < 3:
                top_kmers_animal.append((motif_display, motif, shap_val))
        elif feat_name.startswith("esm_dim_"):
            if shap_val > 0:
                esm_human_count += 1
            elif shap_val < 0:
                esm_animal_count += 1
    
    # Check which k-mers are actually present in the query sequence
    def check_presence(motif_raw, seq):
        return motif_raw in seq
    
    sentences = []
    
    # Opening sentence
    sentences.append(
        f"The model predicts this sequence is {direction} with {confidence} confidence ({prob*100:.1f}%)."
    )
    
    if is_human:
        # Explain human-adapted drivers
        if top_kmers_human:
            present = [m for m in top_kmers_human if check_presence(m[1], vp8_seq)]
            absent = [m for m in top_kmers_human if not check_presence(m[1], vp8_seq)]
            if present:
                motif_list = ", ".join([f"'{m[0]}'" for m in present])
                sentences.append(
                    f"The VP8* domain contains the amino acid motif(s) {motif_list}, "
                    f"which are characteristic of human-tropic rotavirus strains and strongly "
                    f"pushed the prediction towards human adaptation."
                )
            if absent:
                motif_list = ", ".join([f"'{m[0]}'" for m in absent])
                sentences.append(
                    f"Additionally, the absence or low frequency of motif(s) {motif_list} "
                    f"contributed to the overall scoring."
                )
        if top_kmers_animal:
            present_animal = [m for m in top_kmers_animal if check_presence(m[1], vp8_seq)]
            if present_animal:
                motif_list = ", ".join([f"'{m[0]}'" for m in present_animal])
                sentences.append(
                    f"However, the sequence also retains animal-associated motif(s) {motif_list}, "
                    f"which partially counteracted the human signal."
                )
    else:
        # Explain animal-adapted drivers
        if top_kmers_animal:
            present = [m for m in top_kmers_animal if check_presence(m[1], vp8_seq)]
            if present:
                motif_list = ", ".join([f"'{m[0]}'" for m in present])
                sentences.append(
                    f"The VP8* domain contains the amino acid motif(s) {motif_list}, "
                    f"which are characteristic of animal-origin rotavirus strains and strongly "
                    f"pushed the prediction away from human adaptation."
                )
        if top_kmers_human:
            present_human = [m for m in top_kmers_human if check_presence(m[1], vp8_seq)]
            if present_human:
                motif_list = ", ".join([f"'{m[0]}'" for m in present_human])
                sentences.append(
                    f"Notably, the sequence does contain human-associated motif(s) {motif_list}, "
                    f"suggesting partial but incomplete adaptation to human hosts."
                )
    
    # ESM-2 structural context
    esm_total = esm_human_count + esm_animal_count
    if esm_total > 0:
        dominant = "human" if esm_human_count > esm_animal_count else "animal"
        sentences.append(
            f"Beyond sequence motifs, {esm_total} ESM-2 structural embedding dimensions contributed "
            f"to the prediction, with the majority encoding {dominant}-associated protein folding patterns."
        )
    
    return " ".join(sentences)

@app.post("/predict", response_model=PredictionResponse)
async def predict_sequence(
    file: Optional[UploadFile] = File(None),
    raw_sequence: Optional[str] = Form(None)
):
    if file:
        contents = await file.read()
        fasta_str = contents.decode("utf-8")
    elif raw_sequence:
        fasta_str = raw_sequence
    else:
        raise HTTPException(status_code=400, detail="Must provide either a file or a pasted sequence.")
        
    records = list(SeqIO.parse(StringIO(fasta_str), "fasta"))
    if not records:
        # Fallback for raw text without a FASTA header
        sequence = "".join(fasta_str.upper().split())
        record_id = "Pasted_Sequence"
    else:
        record = records[0]
        sequence = str(record.seq).upper().replace("-", "")
        record_id = record.id
    
    # 1. Extract VP8
    try:
        vp8_seq, start, end, ref_aligned, query_aligned = extract_vp8(sequence)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
        
    # 2. Extract Features
    df_esm = extract_esm2(vp8_seq)
    df_kmer = extract_kmer(vp8_seq)
    
    # 3. Combine Features
    X_input = pd.concat([df_esm, df_kmer], axis=1)
    
    # 4. Predict
    xgb_model = model_state['xgb']
    prob = float(xgb_model.predict_proba(X_input)[0, 1])
    is_human = prob >= 0.5
    
    import shap
    explainer = shap.TreeExplainer(xgb_model)
    shap_values = explainer.shap_values(X_input, approximate=True)
    
    feature_names = X_input.columns.tolist()
    shap_dict = {feature_names[i]: float(shap_values[0][i]) for i in range(len(feature_names))}
    sorted_shap = sorted(shap_dict.items(), key=lambda x: x[1], reverse=True)
    
    top_human = [{"feature": k, "impact": v} for k, v in sorted_shap[:5] if v > 0]
    top_animal = [{"feature": k, "impact": v} for k, v in sorted_shap[-5:] if v < 0]
    top_animal = sorted(top_animal, key=lambda x: x["impact"]) # Most negative first
    
    interpretation = generate_interpretation(prob, is_human, sorted_shap, vp8_seq)
    
    return PredictionResponse(
        accession=record_id,
        zoonotic_potential=prob * 100.0,
        vp8_length=len(vp8_seq),
        alignment_start=start + 1,  # 1-indexed for display
        alignment_end=end,
        is_human_adapted=is_human,
        ref_aligned=ref_aligned,
        query_aligned=query_aligned,
        top_human_features=top_human,
        top_animal_features=top_animal,
        interpretation=interpretation
    )

@app.post("/predict/phylo", response_model=PhyloResponse)
async def predict_phylo(file: UploadFile = File(...)):
    # 1. Read FASTA and check sequence limit
    contents = await file.read()
    fasta_str = contents.decode("utf-8")
    records = list(SeqIO.parse(StringIO(fasta_str), "fasta"))
    
    if len(records) < 2:
        raise HTTPException(status_code=400, detail="Must provide at least 2 sequences for phylogenetic analysis.")
    if len(records) > 100:
        raise HTTPException(status_code=400, detail="Maximum of 100 sequences is allowed for performance limits.")

    # 2. Extract VP8* domain for each sequence
    vp8_records = []
    skipped_records = []
    
    for r in records:
        sequence = str(r.seq).upper().replace("-", "")
        clean_id = clean_header(r.id)
        try:
            vp8_seq, start, end, ref_aligned, query_aligned = extract_vp8(sequence)
            vp8_records.append((clean_id, vp8_seq, start, end))
        except Exception as e:
            skipped_records.append((r.id, str(e)))
            
    if len(vp8_records) < 2:
        error_details = "; ".join([f"{rid}: {err}" for rid, err in skipped_records])
        raise HTTPException(status_code=400, detail=f"Failed to align sufficient sequences to VP8 reference. Details: {error_details}")

    # 3. Create temp workspace directory inside scratch
    scratch_dir = project_root / "scratch"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    
    predictions = {}
    
    with tempfile.TemporaryDirectory(dir=str(scratch_dir)) as tmpdir:
        tmp_path = Path(tmpdir)
        in_fasta_path = tmp_path / "input.fasta"
        aligned_fasta_path = tmp_path / "aligned.fasta"
        
        # Write extracted VP8* domains to fasta
        with open(in_fasta_path, "w") as f:
            for clean_id, vp8_seq, start, end in vp8_records:
                f.write(f">{clean_id}\n{vp8_seq}\n")
                
        # Get tool paths
        mafft_path = get_tool_path("mafft")
        iqtree_path = get_tool_path("iqtree2")
        
        # 4. Run MAFFT Alignment
        try:
            is_win = os.name == 'nt'
            mafft_cmd = [mafft_path, "--auto", str(in_fasta_path)]
            res = subprocess.run(mafft_cmd, capture_output=True, text=True, check=True, shell=is_win)
            with open(aligned_fasta_path, "w") as f:
                f.write(res.stdout)
        except Exception as e:
            logger.error(f"MAFFT alignment failed: {e}")
            raise HTTPException(status_code=500, detail=f"MAFFT alignment failed: {str(e)}")
            
        # 5. Run IQ-TREE
        try:
            prefix = tmp_path / "iqtree_run"
            iqtree_cmd = [
                iqtree_path, 
                "-s", str(aligned_fasta_path), 
                "-m", "MFP", 
                "-bb", "1000", 
                "-nt", "2", 
                "-pre", str(prefix),
                "-redo"
            ]
            subprocess.run(iqtree_cmd, capture_output=True, check=True, shell=is_win)
        except Exception as e:
            logger.error(f"IQ-TREE tree building failed: {e}")
            raise HTTPException(status_code=500, detail=f"IQ-TREE tree building failed: {str(e)}")
            
        # 6. Read output tree
        tree_file = prefix.with_suffix(".treefile")
        if not tree_file.exists():
            raise HTTPException(status_code=500, detail="IQ-TREE did not produce a treefile.")
            
        with open(tree_file, "r") as f:
            newick_str = f.read().strip()
            
        # Parse selected model from .iqtree report or log file
        model_used = "Unknown"
        
        # Check .iqtree file first (primary source)
        iqtree_report = prefix.with_suffix(".iqtree")
        if iqtree_report.exists():
            with open(iqtree_report, "r", errors="ignore") as f:
                for line in f:
                    if "Best-fit model according to" in line or "Best-fit model:" in line:
                        parts = line.split(":")
                        if len(parts) > 1:
                            model_used = parts[1].split("chosen")[0].strip()
                            break
                            
        # Fallback to .log file
        if model_used == "Unknown":
            log_file = prefix.with_suffix(".log")
            if log_file.exists():
                with open(log_file, "r", errors="ignore") as f:
                    for line in f:
                        if "Best-fit model according to" in line or "Best-fit model:" in line:
                            parts = line.split(":")
                            if len(parts) > 1:
                                model_used = parts[1].split("chosen")[0].strip()
                                break

        # 7. Run predictions for all leaves
        xgb_model = model_state['xgb']
        import shap
        explainer = shap.TreeExplainer(xgb_model)
        
        for clean_id, vp8_seq, start, end in vp8_records:
            df_esm = extract_esm2(vp8_seq)
            df_kmer = extract_kmer(vp8_seq)
            X_input = pd.concat([df_esm, df_kmer], axis=1)
            
            prob = float(xgb_model.predict_proba(X_input)[0, 1])
            is_human = prob >= 0.5
            
            shap_values = explainer.shap_values(X_input, approximate=True)
            feature_names = X_input.columns.tolist()
            shap_dict = {feature_names[i]: float(shap_values[0][i]) for i in range(len(feature_names))}
            sorted_shap = sorted(shap_dict.items(), key=lambda x: x[1], reverse=True)
            
            interpretation = generate_interpretation(prob, is_human, sorted_shap, vp8_seq)
            
            predictions[clean_id] = PhyloNodePrediction(
                score=prob * 100.0,
                is_human_adapted=is_human,
                vp8_length=len(vp8_seq),
                alignment_start=start + 1,
                alignment_end=end,
                interpretation=interpretation
            )
            
        # Parse alignment length
        from Bio import AlignIO
        alignment = AlignIO.read(aligned_fasta_path, "fasta")
        alignment_length = alignment.get_alignment_length()

        return PhyloResponse(
            newick=newick_str,
            predictions=predictions,
            model_used=model_used,
            alignment_length=alignment_length,
            num_sequences=len(vp8_records)
        )

@app.post("/predict/structure", response_model=StructureResponse)
async def predict_structure(
    file: Optional[UploadFile] = File(None),
    raw_sequence: Optional[str] = Form(None)
):
    if file:
        contents = await file.read()
        fasta_str = contents.decode("utf-8")
    elif raw_sequence:
        fasta_str = raw_sequence
    else:
        raise HTTPException(status_code=400, detail="Must provide either a file or a pasted sequence.")
        
    records = list(SeqIO.parse(StringIO(fasta_str), "fasta"))
    if not records:
        sequence = "".join(fasta_str.upper().split())
    else:
        sequence = str(records[0].seq).upper().replace("-", "")
        
    # 1. Extract VP8
    try:
        vp8_seq, start, end, ref_aligned, query_aligned = extract_vp8(sequence)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
        
    # 2. Predict & SHAP
    df_esm = extract_esm2(vp8_seq)
    df_kmer = extract_kmer(vp8_seq)
    X_input = pd.concat([df_esm, df_kmer], axis=1)
    
    xgb_model = model_state['xgb']
    prob = float(xgb_model.predict_proba(X_input)[0, 1])
    is_human = prob >= 0.5
    
    import shap
    explainer = shap.TreeExplainer(xgb_model)
    shap_values = explainer.shap_values(X_input, approximate=True)
    feature_names = X_input.columns.tolist()
    
    # 3. Map SHAP values to residues (K-mer based mapping)
    query_len = len(vp8_seq)
    residue_shap = [0.0] * query_len
    
    for feat_name, val in zip(feature_names, shap_values[0]):
        if feat_name.startswith("kmer_"):
            kmer = feat_name[5:]
            k_len = len(kmer)
            for m in re.finditer(f"(?={kmer})", vp8_seq):
                idx = m.start()
                for offset in range(k_len):
                    residue_shap[idx + offset] += val / k_len
                    
    # 4. Load static PDB
    pdb_path = Path(__file__).resolve().parent / "static" / "2dwr.pdb"
    if not pdb_path.exists():
        pdb_dir = pdb_path.parent
        pdb_dir.mkdir(parents=True, exist_ok=True)
        pdb_url = "https://files.rcsb.org/download/2DWR.pdb"
        try:
            import urllib.request
            urllib.request.urlretrieve(pdb_url, pdb_path)
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"PDB file 2dwr.pdb not found and download failed: {e}")
            
    with open(pdb_path, "r") as f:
        pdb_content = f.read()
        
    pdb_seq, pdb_res_nums = get_pdb_sequence(pdb_path)
    
    # 5. Align query to PDB sequence
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    
    alignments = aligner.align(pdb_seq, vp8_seq)
    if not alignments:
        raise HTTPException(status_code=500, detail="Could not align query sequence to reference PDB sequence.")
        
    alignment = alignments[0]
    
    pdb_aligned = str(alignment[0])
    query_aligned = str(alignment[1])
    
    mapped_residues = []
    pdb_idx = 0
    query_idx = 0
    
    for char_pdb, char_query in zip(pdb_aligned, query_aligned):
        is_pdb_gap = char_pdb == '-'
        is_query_gap = char_query == '-'
        
        if not is_pdb_gap:
            res_num_pdb = pdb_res_nums[pdb_idx]
            pdb_idx += 1
            
            if not is_query_gap:
                res_num_query = query_idx
                aa_query = char_query
                val_shap = residue_shap[query_idx]
                is_mut = char_pdb != char_query
                query_idx += 1
                
                mapped_residues.append(StructureResidue(
                    pdb_res_num=res_num_pdb,
                    input_res_num=res_num_query + 1,
                    input_aa=aa_query,
                    shap_val=float(val_shap),
                    is_mutation=is_mut
                ))
            else:
                mapped_residues.append(StructureResidue(
                    pdb_res_num=res_num_pdb,
                    input_res_num=None,
                    input_aa=None,
                    shap_val=0.0,
                    is_mutation=True
                ))
        else:
            if not is_query_gap:
                query_idx += 1
                
    return StructureResponse(
        pdb_content=pdb_content,
        residues=mapped_residues,
        zoonotic_potential=prob * 100.0,
        is_human_adapted=is_human
    )

@app.get("/health")
def health():
    return {"status": "ok"}

# Serve static files and frontend
static_dir = Path(__file__).resolve().parent / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

frontend_dist = project_root / "WebApp" / "frontend" / "dist"

if frontend_dist.exists():
    app.mount("/assets", StaticFiles(directory=str(frontend_dist / "assets")), name="assets")

    @app.get("/")
    async def serve_index():
        return FileResponse(str(frontend_dist / "index.html"))
        
    @app.get("/{full_path:path}")
    async def serve_frontend(full_path: str):
        # Prevent intercepting API/static routes
        if full_path.startswith("predict") or full_path.startswith("health") or full_path.startswith("static"):
            raise HTTPException(status_code=404, detail="Not Found")
        return FileResponse(str(frontend_dist / "index.html"))
