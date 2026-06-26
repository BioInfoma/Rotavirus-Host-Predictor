# Rotavirus VP4 Zoonotic Potential Predictor - Comprehensive Documentation

## 1. Project Goals & Objectives
Rotavirus A is a leading cause of severe gastroenteritis worldwide. While most strains are host-specific (adapted to humans or specific animal species), "spillover" events occur where animal strains infect humans, potentially leading to novel pandemic strains.

**The Goal:** The primary objective of this project is to develop an Interpretable Machine Learning pipeline that predicts the zoonotic potential (host adaptation) of Rotavirus A sequences directly from their raw FASTA sequences. By analyzing the VP8* domain of the VP4 spike protein (which is responsible for host cell receptor binding), the model calculates a "Zoonotic Risk Score" and highlights exactly which amino acid mutations drive that risk.

---

## 2. Dataset Collection & Preprocessing
The model relies on meticulously curated and rigorously preprocessed data to prevent data leakage and ensure biological validity.

### Training Dataset (Ground Truth Anchors)
- **Total Sequences:** 427 cleaned sequences
- **Positive Labels (Human-Adapted):** 213 sequences (primarily Human P[8], P[4] genotypes)
- **Negative Labels (Animal-Adapted):** 214 sequences (deeply divergent animal strains: Equine, Bat, Avian, etc.)
- *Purpose:* Used to train the model to confidently distinguish between pure human adaptation and pure animal adaptation.

### Evaluation Dataset (Intermediate & Zoonotic Hosts)
- **Total Sequences:** 229 cleaned sequences
- **Composition:** Predominantly Bovine (101) and Porcine (63) strains, along with atypical human zoonotic strains (P[6], P[9]).
- *Purpose:* Kept strictly isolated from training. Used to evaluate if the model correctly assigns "intermediate" (30-70%) risk scores to strains known to jump between species.

### Preprocessing Steps
1. **Length & Quality Filtering:** Retained sequences of 800-2600 bp; discarded any with >1% ambiguous nucleotides.
2. **Translation & Open Reading Frame (ORF):** 6-frame translation used to find the longest valid ORF.
3. **VP8* Domain Extraction:** Sequences aligned to reference Wa (P[8]) and DS-1 (P[4]) using Smith-Waterman to isolate the exact VP8* binding domain (amino acids 1-272), ensuring all downstream ML features align structurally.
4. **Redundancy Clustering:** CD-HIT clustering (99% identity) applied to prevent the model from artificially memorizing localized outbreak clusters.

---

## 3. Methodology & Feature Engineering
The pipeline uses a dual-feature representation strategy, combining deep learning structural context with evolutionary motifs.

### A. ESM-2 Protein Language Model Embeddings
- **Model Used:** Facebook/Meta `esm2_t12_35M_UR50D` (35 Million parameters).
- **Process:** The raw VP8* amino acid sequence is passed through the ESM-2 transformer.
- **Output:** A 480-dimensional mean-pooled vector representing the structural, biochemical, and biophysical context of the folded protein domain.

### B. K-mer Frequency Analysis
- **Process:** Extraction of continuous amino acid subsequences (3-mers and 4-mers).
- **Purpose:** Captures specific binding motifs and localized evolutionary signatures that are critical for host receptor attachment.

### C. Machine Learning Classifier (XGBoost)
- The concatenated ESM-2 and K-mer features are fed into a highly optimized **XGBoost Gradient Boosting Classifier**.
- Evaluated via Leave-One-P-Genotype-Out cross-validation to ensure the model learns generalized adaptation rules rather than just memorizing specific P-types.
- **Accuracy & Results:** The model achieves extremely high cross-validation accuracy on the training anchors (>95%). More importantly, when tested on the strictly isolated evaluation set of Porcine/Bovine sequences, the model outputs scores largely in the 30-70% range, correctly identifying them as having intermediate/zoonotic characteristics.

### D. Explainable AI (SHAP)
- **Method:** TreeSHAP (SHapley Additive exPlanations) is calculated for every prediction.
- **Purpose:** Unboxes the "black box" of the XGBoost model. It maps the mathematical risk back to specific amino acid features, quantifying exactly how much a particular mutation pushed the score toward "Human" or "Animal".

---

## 4. Web Dashboard & User Interface
The project features a modern, glassmorphic React frontend communicating with a FastAPI Python backend, designed to be accessible to basic researchers and bioinformaticians alike.

### Single Sequence Analysis Tab
- **Input:** Users can paste raw DNA or Protein FASTA sequences.
- **Backend Pipeline:** The backend automatically handles translation, VP8* extraction via pairwise alignment, ESM-2 embedding, and XGBoost inference on-the-fly.
- **Plain English Interpretation:** Instead of overwhelming biologists with raw SHAP mathematics, the AI automatically translates the SHAP arrays into a readable "Interpretation" paragraph (e.g., *"The model predicts this sequence is animal-adapted. The motif D-I-T strongly pulled the prediction toward animal adaptation..."*).
- **Advanced Math:** Researchers can toggle to view the raw SHAP feature weights if desired.

### Interactive 3D Structure Visualizer Tab
- **Engine:** Built using `3dmol.js` embedded in React.
- **Mapping:** The user's sequence is globally aligned to the crystal structure of the rotavirus VP8* domain (PDB: `2DWR`).
- **Coloring Modes:** 
  1. **SHAP Risk:** Colors the 3D protein structure based on SHAP values. Red residues indicate areas driving zoonotic risk; Blue indicates human-adapted regions.
  2. **Mutations:** Highlights exactly where the user's sequence physically differs from the reference crystal structure.
- **Interaction:** Hovering over or clicking a 3D atom opens a sidebar detailing the exact amino acid, its position, and its exact mathematical contribution to the zoonotic risk score.
