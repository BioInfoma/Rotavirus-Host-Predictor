import React, { useState, useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import axios from 'axios';
import { UploadCloud, ShieldAlert, CheckCircle2, FileText, Search, Activity, GitBranch } from 'lucide-react';
import './index.css';

// Import the new tabs
import StructureTab from './components/StructureTab';

function App() {
  const API_URL = import.meta.env.DEV ? 'http://127.0.0.1:8000' : '';
  const [activeTab, setActiveTab] = useState('single'); // single, phylo, structure
  
  // Single sequence state
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState(null);
  const [error, setError] = useState('');
  const [pastedSeq, setPastedSeq] = useState('');
  const [showAdvanced, setShowAdvanced] = useState(false);

  const onDrop = useCallback(async (acceptedFiles) => {
    const file = acceptedFiles[0];
    if (!file) return;

    setLoading(true);
    setError('');
    setResult(null);

    const formData = new FormData();
    formData.append('file', file);

    try {
      const response = await axios.post(`${API_URL}/predict`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });
      setResult(response.data);
    } catch (err) {
      setError(err.response?.data?.detail || err.message || 'An error occurred during prediction.');
    } finally {
      setLoading(false);
    }
  }, []);

  const handlePasteSubmit = async () => {
    if (!pastedSeq.trim()) return;
    setLoading(true);
    setError('');
    setResult(null);

    const formData = new FormData();
    formData.append('raw_sequence', pastedSeq);

    try {
      const response = await axios.post(`${API_URL}/predict`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });
      setResult(response.data);
    } catch (err) {
      setError(err.response?.data?.detail || err.message || 'An error occurred during prediction.');
    } finally {
      setLoading(false);
    }
  };

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    multiple: false
  });

  return (
    <div>
      <h1>Rotavirus Host Predictor</h1>
      <p style={{ color: '#94a3b8', marginBottom: '2rem' }}>
        Bioinformatics pipeline combining ESM-2 deep learning and XGBoost to predict zoonotic potential.
      </p>

      {/* Tab Selector Bar */}
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        gap: '1rem',
        marginBottom: '2.5rem'
      }}>
        {[
          { id: 'single', label: 'Single Sequence Analysis', icon: <Activity size={18} /> },
          { id: 'structure', label: '3D Structure Visualizer', icon: <Search size={18} /> }
        ].map(tab => (
          <button
            key={tab.id}
            onClick={() => setActiveTab(tab.id)}
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem',
              padding: '0.75rem 1.5rem',
              background: activeTab === tab.id ? 'rgba(192, 132, 252, 0.2)' : 'rgba(30, 41, 59, 0.4)',
              border: activeTab === tab.id ? '1px solid #c084fc' : '1px solid rgba(255, 255, 255, 0.05)',
              borderRadius: '9999px',
              color: activeTab === tab.id ? '#f8fafc' : '#94a3b8',
              cursor: 'pointer',
              fontSize: '0.9rem',
              fontWeight: '600',
              backdropFilter: 'blur(12px)',
              transition: 'all 0.2s ease'
            }}
          >
            {tab.icon}
            {tab.label}
          </button>
        ))}
      </div>

      {/* Render Active Tab */}
      {activeTab === 'single' && (
        <div className="glass-panel" style={{ maxWidth: '800px', margin: '0 auto' }}>
          
          {!loading && !result && (
            <>
              <div {...getRootProps()} className={`dropzone ${isDragActive ? 'active' : ''}`}>
                <input {...getInputProps()} />
                <UploadCloud className="icon" />
                {isDragActive ? (
                  <h3>Drop the FASTA file here ...</h3>
                ) : (
                  <div>
                    <h3>Drag & Drop your FASTA sequence here</h3>
                    <p style={{ color: '#64748b' }}>or click to select file</p>
                  </div>
                )}
              </div>

              <div style={{ margin: '2rem 0', textAlign: 'center', color: '#64748b' }}>
                <span>— OR —</span>
              </div>

              <div className="glass-panel" style={{ padding: '1.5rem', background: 'rgba(30, 41, 59, 0.5)' }}>
                <h4 style={{ margin: '0 0 1rem 0', color: '#e2e8f0', textAlign: 'left' }}>Paste Sequence</h4>
                <textarea 
                  value={pastedSeq}
                  onChange={(e) => setPastedSeq(e.target.value)}
                  placeholder="Paste raw nucleotide or protein sequence here (FASTA header is optional)..."
                  style={{ width: '100%', minHeight: '120px', padding: '1rem', background: '#0f172a', color: '#e2e8f0', border: '1px solid #334155', borderRadius: '8px', fontFamily: 'monospace', resize: 'vertical' }}
                />
                <button 
                  onClick={handlePasteSubmit}
                  disabled={!pastedSeq.trim()}
                  style={{ marginTop: '1rem', width: '100%', padding: '0.75rem', background: pastedSeq.trim() ? '#3b82f6' : '#1e293b', color: pastedSeq.trim() ? '#fff' : '#64748b', border: 'none', borderRadius: '8px', cursor: pastedSeq.trim() ? 'pointer' : 'not-allowed', fontWeight: '600', transition: 'background 0.2s' }}
                >
                  Analyze Pasted Sequence
                </button>
              </div>
            </>
          )}

          {loading && (
            <div style={{ padding: '3rem 0', textAlign: 'center' }}>
              <div className="loader"></div>
              <h3 style={{ marginTop: '1.5rem', color: '#e2e8f0' }}>Analyzing Sequence...</h3>
              <p style={{ color: '#64748b' }}>Running VP8* alignment and ESM-2 deep learning extraction</p>
            </div>
          )}

          {error && (
            <div className="error-msg">
              <ShieldAlert size={24} style={{ marginBottom: '0.5rem' }} />
              <div>{error}</div>
              <button 
                onClick={() => setError('')} 
                style={{ marginTop: '1rem', padding: '0.5rem 1rem', background: 'transparent', border: '1px solid #f87171', color: '#f87171', borderRadius: '4px', cursor: 'pointer' }}
              >
                Try Again
              </button>
            </div>
          )}

          {result && (
            <div style={{ animation: 'fadeIn 0.5s ease' }}>
              <div 
                className="score-container" 
                style={{ '--score': `${result.zoonotic_potential}%` }}
              >
                <div className="score-inner">
                  <div className="score-value">{result.zoonotic_potential.toFixed(1)}<span style={{ fontSize: '1.5rem'}}>%</span></div>
                  <div className="score-label">Zoonotic Potential</div>
                </div>
              </div>

              <div className={`status-badge ${result.is_human_adapted ? 'status-human' : 'status-animal'}`}>
                {result.is_human_adapted ? (
                  <><CheckCircle2 size={16} style={{ display: 'inline', verticalAlign: 'text-bottom', marginRight: '4px' }}/> Predicted: Human-Adapted</>
                ) : (
                  <><ShieldAlert size={16} style={{ display: 'inline', verticalAlign: 'text-bottom', marginRight: '4px' }}/> Predicted: Animal-Adapted (Zoonotic Risk)</>
                )}
              </div>

              <div className="details-grid">
                <div className="detail-item">
                  <div className="detail-label">Accession</div>
                  <div className="detail-value" style={{ wordBreak: 'break-all', fontSize: '1rem' }}>{result.accession}</div>
                </div>
                <div className="detail-item">
                  <div className="detail-label">VP8* Extracted Length</div>
                  <div className="detail-value">{result.vp8_length} aa</div>
                </div>
                <div className="detail-item">
                  <div className="detail-label">Alignment Coordinates</div>
                  <div className="detail-value">{result.alignment_start} - {result.alignment_end}</div>
                </div>
              </div>

              {/* Biological Interpretation */}
              {result.interpretation && (
                <div className="glass-panel" style={{ marginTop: '2.5rem', padding: '2rem', background: 'linear-gradient(135deg, rgba(59, 130, 246, 0.2), rgba(139, 92, 246, 0.2))', textAlign: 'left', border: '1px solid #8b5cf6', boxShadow: '0 4px 20px rgba(139, 92, 246, 0.15)' }}>
                  <h3 style={{ margin: '0 0 1rem 0', color: '#c4b5fd', display: 'flex', alignItems: 'center', gap: '0.75rem', fontSize: '1.4rem' }}>
                    <FileText size={24} /> Interpretation
                  </h3>
                  <p style={{ color: '#f8fafc', fontSize: '1.1rem', lineHeight: '1.8', margin: 0 }}>
                    {result.interpretation}
                  </p>
                </div>
              )}

              <div style={{ textAlign: 'center', marginTop: '2rem' }}>
                <button 
                  onClick={() => setShowAdvanced(!showAdvanced)}
                  style={{ background: 'transparent', border: '1px solid #334155', color: '#94a3b8', padding: '0.5rem 1rem', borderRadius: '4px', cursor: 'pointer', transition: 'all 0.2s' }}
                >
                  {showAdvanced ? 'Hide Advanced Math Details' : 'View Advanced Math Details (SHAP)'}
                </button>
              </div>

              {/* SHAP Feature Importance */}
              {showAdvanced && (
              <div className="glass-panel" style={{ marginTop: '1.5rem', padding: '1.5rem', background: 'rgba(30, 41, 59, 0.5)', textAlign: 'left' }}>
                <h4 style={{ margin: '0 0 0.5rem 0', color: '#e2e8f0' }}>Key Amino Acid Mutations Driving the Result (SHAP)</h4>
                <p style={{ color: '#94a3b8', fontSize: '0.85rem', marginBottom: '1.5rem', lineHeight: '1.4' }}>
                  <strong>How to interpret:</strong> SHAP values show the mathematical impact of a specific feature on the model's decision. 
                  Higher positive values strongly push the prediction towards Human-adapted, while negative values push it towards Animal-adapted. 
                  "esm_dim" refers to structural deep learning embeddings, while amino acid letters represent specific K-mer motifs.
                </p>
                <div style={{ display: 'flex', gap: '2rem' }}>
                  <div style={{ flex: 1 }}>
                    <h5 style={{ color: '#f87171', margin: '0 0 0.5rem 0' }}>Animal-Adapted Signals</h5>
                    {result.top_animal_features.map((f, i) => (
                      <div key={i} style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.25rem', fontSize: '0.9rem' }}>
                        <span style={{ color: '#cbd5e1', fontFamily: 'monospace' }}>{f.feature}</span>
                        <span style={{ color: '#f87171' }}>{f.impact.toFixed(3)}</span>
                      </div>
                    ))}
                  </div>
                  <div style={{ flex: 1 }}>
                    <h5 style={{ color: '#34d399', margin: '0 0 0.5rem 0' }}>Human-Adapted Signals</h5>
                    {result.top_human_features.map((f, i) => (
                      <div key={i} style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '0.25rem', fontSize: '0.9rem' }}>
                        <span style={{ color: '#cbd5e1', fontFamily: 'monospace' }}>{f.feature}</span>
                        <span style={{ color: '#34d399' }}>+{f.impact.toFixed(3)}</span>
                      </div>
                    ))}
                  </div>
                </div>
              </div>
              )}

              {/* Sequence Alignment Visualization */}
              <div className="glass-panel" style={{ marginTop: '2rem', padding: '1.5rem', background: 'rgba(30, 41, 59, 0.5)', textAlign: 'left' }}>
                <h4 style={{ margin: '0 0 0.5rem 0', color: '#e2e8f0' }}>Sequence Alignment (vs Human Wa P[8])</h4>
                <p style={{ color: '#94a3b8', fontSize: '0.85rem', marginBottom: '1.5rem', lineHeight: '1.4' }}>
                  <strong>How to interpret:</strong> This shows your extracted sequence aligned against the prototypical Human Rotavirus Wa strain. 
                  <span style={{ color: '#f43f5e', fontWeight: 'bold' }}> Pink highlighted letters </span> 
                  indicate mutations where your sequence differs from the human reference, potentially contributing to its zoonotic risk.
                </p>
                <div style={{ overflowX: 'auto', background: '#0f172a', padding: '1rem', borderRadius: '8px', border: '1px solid #334155' }}>
                  <div style={{ fontFamily: 'monospace', whiteSpace: 'nowrap', fontSize: '0.9rem' }}>
                    <div style={{ display: 'flex', marginBottom: '0.25rem' }}>
                      <span style={{ width: '80px', color: '#94a3b8' }}>Reference</span>
                      <span style={{ color: '#94a3b8' }}>{result.ref_aligned}</span>
                    </div>
                    <div style={{ display: 'flex' }}>
                      <span style={{ width: '80px', color: '#60a5fa' }}>Extracted</span>
                      <span>
                        {result.query_aligned.split('').map((char, index) => {
                          const isMismatch = result.ref_aligned[index] !== char && char !== '-';
                          return (
                            <span key={index} style={{ color: isMismatch ? '#f43f5e' : '#e2e8f0', fontWeight: isMismatch ? 'bold' : 'normal' }}>
                              {char}
                            </span>
                          );
                        })}
                      </span>
                    </div>
                  </div>
                </div>
              </div>

              <button 
                onClick={() => setResult(null)} 
                style={{ marginTop: '2.5rem', padding: '0.75rem 2rem', background: '#c084fc', color: '#fff', border: 'none', borderRadius: '8px', cursor: 'pointer', fontSize: '1rem', fontWeight: '600', transition: 'background 0.2s' }}
                onMouseOver={(e) => e.target.style.background = '#a855f7'}
                onMouseOut={(e) => e.target.style.background = '#c084fc'}
              >
                Analyze Another Sequence
              </button>
            </div>
          )}

        </div>
      )}



      {activeTab === 'structure' && (
        <StructureTab />
      )}
    </div>
  );
}

export default App;
