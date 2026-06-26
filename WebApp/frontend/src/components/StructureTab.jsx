import React, { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import { ShieldAlert, Info, Layers, RefreshCw, FileText } from 'lucide-react';
import axios from 'axios';

function StructureTab() {
  const API_URL = import.meta.env.DEV ? 'http://127.0.0.1:8000' : '';
  
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [seqInput, setSeqInput] = useState('');
  const [data, setData] = useState(null);
  const [styleMode, setStyleMode] = useState('cartoon'); // cartoon, sphere, stick
  const [colorMode, setColorMode] = useState('shap'); // shap, mutation, uniform
  const [selectedResidue, setSelectedResidue] = useState(null);
  const [searchQuery, setSearchQuery] = useState('');
  
  const viewerContainerRef = useRef(null);
  const viewerRef = useRef(null);

  // Trigger analysis for structure mapping
  const handleAnalyze = async () => {
    if (!seqInput.trim()) return;
    setLoading(true);
    setError('');
    setData(null);
    setSelectedResidue(null);
    
    const formData = new FormData();
    formData.append('raw_sequence', seqInput);
    
    try {
      const response = await axios.post(`${API_URL}/predict/structure`, formData);
      setData(response.data);
    } catch (err) {
      setError(err.response?.data?.detail || err.message || 'Error loading structural mapping.');
    } finally {
      setLoading(false);
    }
  };

  // Re-draw or update colors on the 3D model
  useEffect(() => {
    if (!data || !viewerContainerRef.current) return;
    
    // Clear container
    viewerContainerRef.current.innerHTML = '';
    
    // Create viewer
    const viewer = $3Dmol.createViewer(viewerContainerRef.current, {
      backgroundColor: '#0f172a'
    });
    viewerRef.current = viewer;
    
    // Load PDB content
    viewer.addModel(data.pdb_content, 'pdb');
    
    // Map colors to PDB residues
    const residues = data.residues;
    
    // Get min/max SHAP for color scale normalization
    const shapVals = residues.map(r => r.shap_val);
    const maxShap = Math.max(...shapVals.map(Math.abs), 0.01);
    
    // Apply styling residue-by-residue
    residues.forEach(res => {
      let color = '#94a3b8'; // default grey
      
      if (colorMode === 'shap') {
        const val = res.shap_val;
        if (val > 0) {
          // Zoonotic/Animal adapting (Red gradient)
          const intensity = Math.min(Math.floor((val / maxShap) * 255), 255);
          const hex = intensity.toString(16).padStart(2, '0');
          color = `#ff${(255 - intensity).toString(16).padStart(2, '0')}${(255 - intensity).toString(16).padStart(2, '0')}`;
        } else if (val < 0) {
          // Human adapting (Blue gradient)
          const intensity = Math.min(Math.floor((Math.abs(val) / maxShap) * 255), 255);
          color = `#${(255 - intensity).toString(16).padStart(2, '0')}${(255 - intensity).toString(16).padStart(2, '0')}ff`;
        } else {
          color = '#e2e8f0';
        }
      } else if (colorMode === 'mutation') {
        color = res.is_mutation ? '#f43f5e' : '#10b981'; // pink-red for mutation, green for reference match
      } else {
        color = '#c084fc'; // default purple
      }
      
      // Set atom styles for this residue
      const sel = { chain: 'A', resi: res.pdb_res_num };
      
      if (styleMode === 'cartoon') {
        viewer.setStyle(sel, { cartoon: { color: color } });
      } else if (styleMode === 'sphere') {
        viewer.setStyle(sel, { sphere: { color: color, scale: 0.9 } });
      } else if (styleMode === 'stick') {
        viewer.setStyle(sel, { stick: { color: color, radius: 0.25 } });
      }
    });
    
    // Zoom and render
    viewer.zoomTo();
    viewer.render();
    
    // Resize handler
    const handleResize = () => viewer.resize();
    window.addEventListener('resize', handleResize);
    
    return () => {
      window.removeEventListener('resize', handleResize);
    };
  }, [data, styleMode, colorMode]);

  // Click handler on table rows to highlight specific residue in viewer
  const handleSelectResidue = (res) => {
    setSelectedResidue(res);
    if (!viewerRef.current) return;
    
    const viewer = viewerRef.current;
    
    // Clear any existing shapes/labels
    viewer.removeAllLabels();
    
    // Add label at residue position
    const sel = { chain: 'A', resi: res.pdb_res_num };
    
    viewer.addLabel(`Res ${res.pdb_res_num}: ${res.input_aa || 'Gap'} (SHAP: ${res.shap_val.toFixed(4)})`, {
      position: sel,
      background: 'rgba(15, 23, 42, 0.85)',
      fontColor: '#f8fafc',
      borderColor: '#c084fc',
      borderThickness: 1,
      fontSize: 12
    });
    
    // Zoom in on the residue slightly
    viewer.zoomTo(sel, 600);
  };

  const filteredResidues = data
    ? data.residues.filter(r => 
        r.pdb_res_num.toString().includes(searchQuery) ||
        (r.input_aa && r.input_aa.toLowerCase().includes(searchQuery.toLowerCase()))
      )
    : [];

  return (
    <div style={{ textAlign: 'left' }}>
      <h2 style={{ fontSize: '1.75rem', fontWeight: '700', marginBottom: '1rem', color: '#c084fc' }}>
        Interactive 3D Structure & Mutation Visualizer
      </h2>
      <p style={{ color: '#94a3b8', marginBottom: '2rem', fontSize: '0.95rem' }}>
        Map predictive SHAP scores and mutations directly onto the Rotavirus VP8* 3D structure (reference PDB: 2DWR).
        Red residues represent zoonotic adaptation drivers, blue represents human-adapted, and highlighted residues indicate mutations.
      </p>

      {/* Input Stage */}
      {!data && (
        <div className="glass-panel" style={{ padding: '2rem', marginBottom: '2rem' }}>
          <label style={{ display: 'block', fontSize: '0.9rem', fontWeight: '600', color: '#e2e8f0', marginBottom: '0.75rem' }}>
            Paste Rotavirus VP4 sequence (or VP8* domain) to map onto 3D structure:
          </label>
          <textarea
            placeholder="Paste FASTA or raw sequence here..."
            value={seqInput}
            onChange={(e) => setSeqInput(e.target.value)}
            style={{
              width: '100%',
              height: '150px',
              background: 'rgba(15, 23, 42, 0.6)',
              border: '1px solid rgba(255, 255, 255, 0.1)',
              borderRadius: '8px',
              color: '#f8fafc',
              padding: '1rem',
              fontFamily: 'monospace',
              fontSize: '0.9rem',
              resize: 'vertical',
              outline: 'none',
              boxSizing: 'border-box'
            }}
          />
          <button
            onClick={handleAnalyze}
            disabled={loading || !seqInput.trim()}
            style={{
              marginTop: '1.5rem',
              padding: '0.75rem 2rem',
              background: '#c084fc',
              color: '#fff',
              border: 'none',
              borderRadius: '8px',
              cursor: seqInput.trim() ? 'pointer' : 'not-allowed',
              fontSize: '1rem',
              fontWeight: '600',
              opacity: seqInput.trim() ? 1 : 0.6,
              transition: 'background 0.2s',
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem'
            }}
          >
            {loading && <RefreshCw style={{ animation: 'spin 1.5s linear infinite' }} size={18} />}
            {loading ? 'Aligning and Mapping...' : 'Visualize 3D Structure'}
          </button>
          
          {error && <div className="error-msg" style={{ marginTop: '1.5rem' }}>{error}</div>}
        </div>
      )}

      {/* Visualization Stage */}
      {data && (
        <div style={{ display: 'flex', flexDirection: 'column', gap: '2rem' }}>
          
          {/* Biological Interpretation */}
          {data.interpretation && (
            <div className="glass-panel" style={{ padding: '2rem', background: 'linear-gradient(135deg, rgba(59, 130, 246, 0.2), rgba(139, 92, 246, 0.2))', textAlign: 'left', border: '1px solid #8b5cf6', boxShadow: '0 4px 20px rgba(139, 92, 246, 0.15)' }}>
              <h3 style={{ margin: '0 0 1rem 0', color: '#c4b5fd', display: 'flex', alignItems: 'center', gap: '0.75rem', fontSize: '1.4rem' }}>
                <FileText size={24} /> Interpretation
              </h3>
              <p style={{ color: '#f8fafc', fontSize: '1.1rem', lineHeight: '1.8', margin: 0 }}>
                {data.interpretation}
              </p>
            </div>
          )}

          <div style={{ display: 'grid', gridTemplateColumns: '1.5fr 1fr', gap: '2rem', minHeight: '600px' }}>
          
          {/* Left panel: 3D viewer + controls */}
          <div>
            <div className="glass-panel" style={{ padding: '1rem', marginBottom: '1rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
              <div style={{ display: 'flex', gap: '1.5rem' }}>
                {/* Representation selector */}
                <div>
                  <span style={{ fontSize: '0.75rem', textTransform: 'uppercase', color: '#94a3b8', display: 'block', marginBottom: '0.25rem' }}>Style</span>
                  <div style={{ display: 'flex', background: 'rgba(15, 23, 42, 0.6)', padding: '2px', borderRadius: '6px', border: '1px solid rgba(255, 255, 255, 0.05)' }}>
                    {['cartoon', 'sphere', 'stick'].map(style => (
                      <button
                        key={style}
                        onClick={() => setStyleMode(style)}
                        style={{
                          padding: '0.25rem 0.75rem',
                          background: styleMode === style ? '#c084fc' : 'transparent',
                          color: styleMode === style ? '#fff' : '#94a3b8',
                          border: 'none',
                          borderRadius: '4px',
                          cursor: 'pointer',
                          fontSize: '0.8rem',
                          fontWeight: '600',
                          textTransform: 'capitalize'
                        }}
                      >
                        {style}
                      </button>
                    ))}
                  </div>
                </div>

                {/* Color scheme selector */}
                <div>
                  <span style={{ fontSize: '0.75rem', textTransform: 'uppercase', color: '#94a3b8', display: 'block', marginBottom: '0.25rem' }}>Coloring</span>
                  <div style={{ display: 'flex', background: 'rgba(15, 23, 42, 0.6)', padding: '2px', borderRadius: '6px', border: '1px solid rgba(255, 255, 255, 0.05)' }}>
                    {[
                      { id: 'shap', label: 'SHAP Risk' },
                      { id: 'mutation', label: 'Mutations' },
                      { id: 'uniform', label: 'Solid' }
                    ].map(col => (
                      <button
                        key={col.id}
                        onClick={() => setColorMode(col.id)}
                        style={{
                          padding: '0.25rem 0.75rem',
                          background: colorMode === col.id ? '#c084fc' : 'transparent',
                          color: colorMode === col.id ? '#fff' : '#94a3b8',
                          border: 'none',
                          borderRadius: '4px',
                          cursor: 'pointer',
                          fontSize: '0.8rem',
                          fontWeight: '600'
                        }}
                      >
                        {col.label}
                      </button>
                    ))}
                  </div>
                </div>
              </div>

              <button
                onClick={() => setData(null)}
                style={{
                  padding: '0.4rem 1rem',
                  background: 'rgba(244, 63, 94, 0.1)',
                  color: '#f43f5e',
                  border: '1px solid rgba(244, 63, 94, 0.2)',
                  borderRadius: '6px',
                  cursor: 'pointer',
                  fontSize: '0.8rem',
                  fontWeight: '600'
                }}
              >
                Reset
              </button>
            </div>

            {/* The WebGL Canvas */}
            <div 
              ref={viewerContainerRef} 
              className="glass-panel" 
              style={{ 
                height: '500px', 
                padding: 0, 
                overflow: 'hidden', 
                position: 'relative', 
                borderRadius: '16px',
                border: '1px solid rgba(255, 255, 255, 0.1)'
              }}
            />
            
            {/* Color Legend */}
            {colorMode === 'shap' && (
              <div style={{ marginTop: '1rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center', background: 'rgba(30, 41, 59, 0.3)', padding: '0.75rem 1rem', borderRadius: '8px', border: '1px solid rgba(255, 255, 255, 0.05)' }}>
                <span style={{ fontSize: '0.8rem', color: '#60a5fa', fontWeight: '600' }}>← Human-Adapted (Negative SHAP)</span>
                <div style={{ flexGrow: 1, height: '8px', margin: '0 1rem', borderRadius: '4px', background: 'linear-gradient(90deg, #60a5fa, #f8fafc, #f43f5e)' }} />
                <span style={{ fontSize: '0.8rem', color: '#f43f5e', fontWeight: '600' }}>Zoonotic Potential (Positive SHAP) →</span>
              </div>
            )}

            {colorMode === 'mutation' && (
              <div style={{ marginTop: '1rem', display: 'flex', gap: '1.5rem', background: 'rgba(30, 41, 59, 0.3)', padding: '0.75rem 1rem', borderRadius: '8px', border: '1px solid rgba(255, 255, 255, 0.05)' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                  <div style={{ width: '12px', height: '12px', borderRadius: '3px', background: '#f43f5e' }} />
                  <span style={{ fontSize: '0.8rem', color: '#f8fafc' }}>Mutation / Mismatch vs PDB Reference</span>
                </div>
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                  <div style={{ width: '12px', height: '12px', borderRadius: '3px', background: '#10b981' }} />
                  <span style={{ fontSize: '0.8rem', color: '#f8fafc' }}>Identical to Reference</span>
                </div>
              </div>
            )}
          </div>

          {/* Right panel: Sidebar table listing residues */}
          <div className="glass-panel" style={{ display: 'flex', flexDirection: 'column', padding: '1.5rem', maxHeight: '600px' }}>
            <h3 style={{ fontSize: '1.1rem', fontWeight: '600', color: '#f8fafc', marginBottom: '0.75rem' }}>
              Residue Mapping Details
            </h3>
            
            {/* Search filter */}
            <input
              type="text"
              placeholder="Search residue number or AA..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              style={{
                width: '100%',
                padding: '0.5rem 0.75rem',
                background: 'rgba(15, 23, 42, 0.6)',
                border: '1px solid rgba(255, 255, 255, 0.1)',
                borderRadius: '6px',
                color: '#f8fafc',
                fontSize: '0.85rem',
                marginBottom: '1rem',
                outline: 'none',
                boxSizing: 'border-box'
              }}
            />

            {/* Table wrapper */}
            <div style={{ flexGrow: 1, overflowY: 'auto', border: '1px solid rgba(255, 255, 255, 0.05)', borderRadius: '6px' }}>
              <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.85rem' }}>
                <thead>
                  <tr style={{ background: 'rgba(15, 23, 42, 0.8)', borderBottom: '1px solid rgba(255, 255, 255, 0.05)' }}>
                    <th style={{ padding: '0.5rem', textAlign: 'left', color: '#94a3b8' }}>PDB #</th>
                    <th style={{ padding: '0.5rem', textAlign: 'left', color: '#94a3b8' }}>Query</th>
                    <th style={{ padding: '0.5rem', textAlign: 'right', color: '#94a3b8' }}>SHAP</th>
                    <th style={{ padding: '0.5rem', textAlign: 'center', color: '#94a3b8' }}>Mut</th>
                  </tr>
                </thead>
                <tbody>
                  {filteredResidues.map((res, index) => {
                    const isSelected = selectedResidue?.pdb_res_num === res.pdb_res_num;
                    return (
                      <tr
                        key={index}
                        onClick={() => handleSelectResidue(res)}
                        style={{
                          borderBottom: '1px solid rgba(255, 255, 255, 0.03)',
                          cursor: 'pointer',
                          background: isSelected 
                            ? 'rgba(192, 132, 252, 0.15)' 
                            : res.is_mutation 
                              ? 'rgba(244, 63, 94, 0.03)' 
                              : 'transparent',
                          transition: 'background 0.2s'
                        }}
                        onMouseOver={(e) => { if (!isSelected) e.currentTarget.style.background = 'rgba(255, 255, 255, 0.02)'; }}
                        onMouseOut={(e) => { if (!isSelected) e.currentTarget.style.background = res.is_mutation ? 'rgba(244, 63, 94, 0.03)' : 'transparent'; }}
                      >
                        <td style={{ padding: '0.5rem', fontWeight: '600', color: '#e2e8f0' }}>{res.pdb_res_num}</td>
                        <td style={{ padding: '0.5rem', fontFamily: 'monospace' }}>
                          {res.input_aa ? (
                            <span>
                              {res.input_aa} <span style={{ color: '#94a3b8', fontSize: '0.75rem' }}>(#{res.input_res_num})</span>
                            </span>
                          ) : (
                            <span style={{ color: '#64748b', fontStyle: 'italic' }}>Gap</span>
                          )}
                        </td>
                        <td style={{ 
                          padding: '0.5rem', 
                          textAlign: 'right', 
                          fontWeight: '600',
                          color: res.shap_val > 0 ? '#f43f5e' : res.shap_val < 0 ? '#60a5fa' : '#94a3b8'
                        }}>
                          {res.shap_val.toFixed(4)}
                        </td>
                        <td style={{ padding: '0.5rem', textAlign: 'center' }}>
                          {res.is_mutation && (
                            <span style={{ color: '#f43f5e', fontSize: '0.75rem', fontWeight: 'bold' }}>Yes</span>
                          )}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>

            {/* Selected detail box */}
            {selectedResidue && (
              <div style={{ marginTop: '1rem', padding: '0.75rem', background: 'rgba(15, 23, 42, 0.6)', borderRadius: '6px', border: '1px solid rgba(192, 132, 252, 0.3)' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.25rem', marginBottom: '0.25rem' }}>
                  <Info size={14} style={{ color: '#c084fc' }} />
                  <span style={{ fontSize: '0.8rem', fontWeight: '600', color: '#c084fc' }}>Residue Details</span>
                </div>
                <span style={{ fontSize: '0.75rem', color: '#cbd5e1', display: 'block' }}>
                  PDB Residue: <strong>{selectedResidue.pdb_res_num}</strong>
                </span>
                <span style={{ fontSize: '0.75rem', color: '#cbd5e1', display: 'block' }}>
                  Query Residue: <strong>{selectedResidue.input_res_num ? `${selectedResidue.input_aa} (Position ${selectedResidue.input_res_num})` : 'Gap'}</strong>
                </span>
                <span style={{ fontSize: '0.75rem', color: '#cbd5e1', display: 'block' }}>
                  SHAP Contribution: <strong style={{ color: selectedResidue.shap_val > 0 ? '#f43f5e' : selectedResidue.shap_val < 0 ? '#60a5fa' : '#cbd5e1' }}>{selectedResidue.shap_val.toFixed(5)}</strong>
                </span>
                <span style={{ fontSize: '0.75rem', color: '#cbd5e1', display: 'block' }}>
                  Mutation status: <strong>{selectedResidue.is_mutation ? 'Mutated vs reference 2DWR' : 'Identical to reference'}</strong>
                </span>
              </div>
            )}
          </div>
          
        </div>
        </div>
      )}
    </div>
  );
}

export default StructureTab;
