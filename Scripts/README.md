# Rotavirus A VP4 Preprocessing Pipeline

**Publication-Quality Preprocessing for Machine Learning**

A comprehensive preprocessing pipeline that prepares Rotavirus A VP4 sequences for machine learning prediction of human adaptation and spillover potential.

## Project Structure

```
rotavirus_vp4_pipeline/
│
├── main.py                          # Main pipeline entry point
│
├── config/
│   ├── __init__.py
│   └── settings.py                  # Configuration and constants
│
├── preprocessing/
│   ├── __init__.py
│   ├── data_loader.py               # Load training and evaluation data
│   ├── sequence_validator.py        # QC filtering and translation
│   ├── vp8_extractor.py            # VP8* extraction by alignment
│   ├── deduplicator.py             # Deduplication and clustering
│   ├── metadata_normalizer.py       # Host/genotype standardization
│   ├── conflict_detector.py         # Detect ambiguous combinations
│   └── report_generator.py          # QC and composition reports
│
├── data_vp4/
│   ├── Training_data/
│   │   ├── VP4_training_metadata.xlsx
│   │   └── VP4_training_dataset.fasta
│   └── Evaluation_dataset/
│       ├── Eval_metadata_combined.xlsx
│       └── Eval_dataset_comined.fasta
│
├── references/                      # VP4 reference sequences
│   ├── Wa_P8_VP4_protein.fasta
│   └── DS1_P4_VP4_protein.fasta
│
├── analysis_ready/                  # Final outputs
│   ├── cleaned_metadata_training.csv
│   ├── cleaned_metadata_evaluation.csv
│   ├── cleaned_nucleotide_training.fasta
│   ├── cleaned_nucleotide_evaluation.fasta
│   ├── cleaned_protein_training.fasta
│   ├── cleaned_protein_evaluation.fasta
│   ├── cleaned_vp8_training.fasta
│   ├── cleaned_vp8_evaluation.fasta
│   └── reports/
│       ├── composition/             # Distribution reports
│       ├── qc_summary.csv
│       ├── sequence_validation_*.csv
│       ├── vp8_completeness_*.csv
│       ├── host_genotype_conflicts_*.csv
│       └── integrity_summary_*.txt
│
└── logs/                            # Execution logs

```

## Pipeline Steps

### 1. Data Loading
- Loads metadata from Excel files
- Loads nucleotide sequences from FASTA
- Matches sequences to metadata by accession
- **Initial counts:**
  - Training: 427 sequences
  - Evaluation: 229 sequences

### 2. Sequence Validation (QC Filtering)
Removes low-quality sequences based on:
- **Length filtering:** 800-2600 bp (captures complete VP8 domains)
- **Ambiguity filtering:** Max 1% N bases
- **ORF detection:** Finds longest ORF without internal stop codons
- **Translation:** Converts to protein sequences

### 3. VP8* Extraction (Critical Step)
- **Alignment-based detection** (NOT fixed-length proxy)
- Aligns full VP4 proteins to Wa (P[8]) and DS-1 (P[4]) references
- Extracts VP8* domain (AA 1-272 in aligned coordinates)
- **Coverage QC thresholds:**
  - Keep: ≥90% coverage
  - Review: 80-89% coverage
  - Discard: <80% coverage
- **References used:**
  - Wa (P[8]): Accession JX406750
  - DS-1 (P[4]): Accession AB910901

### 4. Exact Deduplication
- Removes identical nucleotide sequences (MD5 hash)
- Keeps first occurrence
- Generates duplicate report

### 5. Redundancy Clustering
- **Goal:** Prevent model from memorizing outbreak clusters
- **Method:** CD-HIT at 99% nucleotide identity (or fallback clustering)
- **Output:** Representative sequences only

### 6. Metadata Normalization
Standardizes all metadata fields:

**Host Normalization:**
- Homo sapiens → Human
- Sus scrofa → Porcine
- Bos taurus → Bovine
- Equus caballus → Equine
- Chiroptera → Bat
- Gallus gallus/other → Avian

**Genotype Normalization:**
- P8, P[8] → P[8]
- P4, P[4] → P[4]
- etc.

**Collection Year Extraction:**
- Extracts 4-digit year from collection date field

### 7. Conflict Detection
Flags (but does NOT remove) ambiguous host-genotype combinations:
- Human P[9], P[14], P[19], P[25] (zoonotic)
- Porcine P[6] (atypical for swine)
- Bat P[6] (atypical for bats)

**Important:** These sequences are preserved as they provide valuable intermediate adaptation examples for evaluation.

### 8. Output Generation
Produces three sequence formats:
- **Nucleotide FASTA:** Original DNA sequences
- **Protein FASTA:** Translated VP4 proteins
- **VP8* FASTA:** Extracted VP8 domains for ML

### 9. Report Generation
Comprehensive QC reports:
- `qc_summary.csv` - Pipeline filtering statistics
- `sequence_validation_*.csv` - Filtering reasons
- `vp8_completeness_*.csv` - VP8* coverage distribution
- `host_genotype_conflicts_*.csv` - Flagged ambiguities
- `integrity_summary_*.txt` - Data completeness checks
- `composition/host_distribution_*.csv` - Host breakdown
- `composition/genotype_distribution_*.csv` - Genotype breakdown
- `composition/label_distribution_*.csv` - Label breakdown
- `composition/adaptation_distribution_*.csv` - Adaptation group breakdown

## Running the Pipeline

### Prerequisites

```bash
pip install biopython pandas openpyxl --break-system-packages
```

### Basic Execution

```bash
cd /path/to/rotavirus_vp4_pipeline
python3 main.py
```

### Output

The pipeline generates:
1. **Analysis-ready sequences** in `/analysis_ready/`
2. **Comprehensive reports** in `/analysis_ready/reports/`
3. **Execution log** in `/logs/`

## Dataset Design

### Training Dataset (Final)
**Positive Class (Human-Adapted Anchors):**
- Host: Homo sapiens
- Genotypes: P[8], P[4], P[6]
- Label: Positive
- Adaptation Group: Human_Anchor
- Expected: ~100-150 sequences after filtering

**Negative Class (Animal Reservoirs):**
- Hosts: Bat, Avian, Equine
- Label: Negative
- Adaptation Group: Animal_Anchor
- Expected: ~50-100 sequences after filtering

### Evaluation Dataset (Final - NEVER used for training)
**Porcine Strains:**
- Host: Sus scrofa
- Genotypes: P[6], P[7], P[13]
- Label: Intermediate
- Adaptation Group: Porcine

**Bovine Strains:**
- Host: Bos taurus
- Genotypes: P[5], P[11]
- Label: Intermediate
- Adaptation Group: Bovine

**Human Intermediate Genotypes:**
- Host: Homo sapiens
- Genotypes: P[9], P[14], P[19], P[25]
- Label: Intermediate
- Adaptation Group: Human_Intermediate

## Configuration

Key parameters in `config/settings.py`:

```python
# Length filtering
SEQUENCE_LENGTH_MIN = 800  # bp
SEQUENCE_LENGTH_MAX = 2600  # bp

# Ambiguity threshold
AMBIGUITY_MAX_PERCENT = 1.0  # %

# VP8* coverage thresholds
VP8_COVERAGE_KEEP = 0.90  # >= 90%
VP8_COVERAGE_REVIEW = 0.80  # 80-89%

# Clustering
CDHIT_IDENTITY_THRESHOLD = 0.99  # 99% identity
```

## Important Notes

### 1. VP8* Extraction is CRITICAL
- Uses **alignment-based detection**, NOT fixed-length coordinates
- Handles divergent sequences across hosts and genotypes
- Coverage-based QC ensures biological completeness
- DO NOT use length proxies

### 2. Train/Evaluation Separation is STRICT
- Evaluation sequences NEVER appear in analysis-ready training outputs
- Pig and bovine sequences are ALWAYS evaluation-only
- Human P[6] may appear in training; Pig P[6] remains evaluation-only

### 3. Conflict Flagging
- Ambiguous combinations are FLAGGED, not REMOVED
- Preserved for evaluation dataset composition
- Provides intermediate adaptation examples

### 4. Reproducibility
- All filtering steps are logged
- MD5 hashes enable exact duplicate detection
- Cluster membership tracked for traceability
- Random seeds set for deterministic results

## Output Schema

### cleaned_metadata_[dataset].csv
```
accession              - Sequence ID
full_name              - Full sequence description
host                   - Normalized host (Human, Porcine, etc.)
genotype               - Normalized genotype (P[8], P[4], etc.)
label                  - Classification (Positive, Negative, Intermediate)
adaptation_group       - Group assignment
country                - Collection country
collection_year        - Extracted year
dataset_split          - training or evaluation
selected_frame         - ORF frame (0, 1, or 2)
protein_length         - Translated protein length
vp8_reference          - Wa_P8 or DS1_P4
vp8_coverage           - Fraction aligned VP8 (0-1)
vp8_status             - PASS, REVIEW, or FAIL
sequence_md5           - Hash for duplicate detection
protein_md5            - Hash for protein deduplication
```

### Sequence FASTAs
Standard FASTA format with:
- Header: `>accession`
- Sequence: 80 chars per line
- No gaps or modifications

## Troubleshooting

### Missing reference sequences
The pipeline includes hardcoded Wa and DS-1 references. For real VP4 references, download from NCBI and place in `/references/`.

### CD-HIT not available
Pipeline automatically falls back to simple clustering (99% identity in sliding window).

### Slow execution on large datasets
- Clustering is the most computationally intensive step
- CD-HIT clusters can be run in parallel on large servers
- Simple fallback is adequate for datasets <500 sequences

## Next Steps (Machine Learning)

The cleaned outputs are ready for:
1. **Feature extraction**
   - ESM-2 protein embeddings
   - k-mer motif analysis
   - Manually curated HBGA-binding residues

2. **Model training**
   - XGBoost classification
   - Leave-one-P-genotype-out validation
   - SHAP explainability analysis

3. **Evaluation**
   - Porcine/Bovine as biological intermediates
   - Human intermediate genotypes for robustness
   - Confidence intervals via bootstrap ensembles

## Citation

If using this pipeline, cite:
- Technical specification document
- Reference genomes (Wa JX406750, DS-1 AB910901)
- CD-HIT for sequence clustering

## Support

For issues or questions:
1. Check the execution log in `/logs/`
2. Review report details in `/analysis_ready/reports/`
3. Verify input data matches expected format
4. Ensure Python packages installed correctly

---

**Pipeline Version:** 1.0  
**Last Updated:** 2024  
**Status:** Production-Ready
