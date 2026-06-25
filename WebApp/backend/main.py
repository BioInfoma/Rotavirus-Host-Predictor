import os
import sys
import logging
from pathlib import Path
from typing import Optional, List, Dict
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
        top_animal_features=top_animal
    )

@app.get("/health")
def health():
    return {"status": "ok"}

# Serve frontend static files
frontend_dist = project_root / "WebApp" / "frontend" / "dist"

if frontend_dist.exists():
    app.mount("/assets", StaticFiles(directory=str(frontend_dist / "assets")), name="assets")

    @app.get("/")
    async def serve_index():
        return FileResponse(str(frontend_dist / "index.html"))
        
    @app.get("/{full_path:path}")
    async def serve_frontend(full_path: str):
        # Prevent intercepting API routes
        if full_path.startswith("predict") or full_path.startswith("health"):
            raise HTTPException(status_code=404, detail="Not Found")
        return FileResponse(str(frontend_dist / "index.html"))
