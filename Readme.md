---
title: Rotavirus Host Predictor
emoji: 🦠
colorFrom: indigo
colorTo: purple
sdk: docker
pinned: false
---

# Rotavirus VP4 Zoonotic Potential Predictor 🦠

A deep-learning-powered web application for predicting the zoonotic potential (host adaptation) of Rotavirus A sequences directly from raw FASTA files.

## Overview
Rotaviruses are a major cause of severe gastroenteritis in infants and young animals worldwide. Understanding whether an animal rotavirus strain has the potential to jump to humans (zoonotic spillover) is critical for pandemic preparedness.

This tool analyzes the **VP8* domain of the VP4 spike protein**, which is responsible for receptor binding and host specificity. 

It uses a dual-model approach:
1. **ESM-2 Deep Learning Embeddings**: Extracts rich, structural 3D-aware features from the protein sequence using Meta's ESM-2 (35M) protein language model.
2. **K-mer Frequency Analysis**: Captures evolutionary motifs and amino acid composition.
3. **XGBoost Classifier**: An optimized gradient-boosting model trained on thousands of curated human and animal rotavirus strains that predicts zoonotic potential with high accuracy.

## Features
- 🧬 **Sequence Alignment**: Automatically aligns uploaded nucleotide or protein sequences against the human Wa P[8] reference using BLOSUM62, handling frameshifts and translating 6 frames to extract the VP8* domain.
- 📊 **SHAP Explainability**: Unboxes the AI! Shows exactly which ESM-2 dimensions and amino acid k-mers drove the model's prediction.
- 💻 **Modern Dashboard**: A beautiful, responsive glassmorphic React interface.

## Local Installation

### Prerequisites
- Python 3.10+
- Node.js & npm

### Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/BioInfoma/Rotavirus-Host-Predictor.git
   cd Rotavirus-Host-Predictor
   ```

2. Setup Python Backend:
   ```bash
   python -m venv vp4_env
   source vp4_env/bin/activate  # Or `vp4_env\Scripts\activate` on Windows
   pip install -r requirements.txt
   ```

3. Setup React Frontend:
   ```bash
   cd WebApp/frontend
   npm install
   npm run build
   ```

4. Run the App:
   ```bash
   cd ../..
   uvicorn WebApp.backend.main:app --host 0.0.0.1 --port 8000
   ```
   Navigate to `http://127.0.0.1:8000` in your browser.

## Deployment (Hugging Face Spaces)
This application is fully Dockerized for 1-click deployment on Hugging Face Spaces.
Live Demo: [https://huggingface.co/spaces/Bionforma/Rotavirus-Host-Predictor](https://huggingface.co/spaces/Bionforma/Rotavirus-Host-Predictor)

## Author
Developed by [@BioInfoma](https://github.com/BioInfoma).
