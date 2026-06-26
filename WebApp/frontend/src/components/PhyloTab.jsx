import React, { useState, useEffect, useRef, useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import axios from 'axios';
import { UploadCloud, RefreshCw, ShieldAlert, CheckCircle2, FileText, Info } from 'lucide-react';
import * as d3 from 'd3';
import { phylotree } from 'phylotree';

// CSS styling for phylotree (injecting statically so D3 classes render beautifully)
const injectPhyloStyles = () => {
  const id = 'phylotree-styles';
  if (document.getElementById(id)) return;
  const style = document.createElement('style');
  style.id = id;
  style.innerHTML = `
    .phylotree-container svg {
      background: #0f172a;
      border-radius: 12px;
      border: 1px solid rgba(255, 255, 255, 0.05);
    }
    .phylotree-container .node {
      font-size: 11px;
      font-family: 'Inter', sans-serif;
      fill: #cbd5e1;
      cursor: pointer;
    }
    .phylotree-container .node:hover text {
      fill: #c084fc !important;
      font-weight: bold;
    }
    .phylotree-container .node.leaf text {
      font-weight: 600;
    }
    .phylotree-container .link {
      fill: none;
      stroke: #475569;
      stroke-width: 2px;
      transition: stroke 0.3s;
    }
    .phylotree-container .link:hover {
      stroke: #c084fc;
      stroke-width: 3px;
    }
  `;
  document.head.appendChild(style);
};

function PhyloTab() {
  const API_URL = import.meta.env.DEV ? 'http://127.0.0.1:8000' : '';
  
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [data, setData] = useState(null); // { newick, predictions, model_used, alignment_length, num_sequences }
  const [selectedLeaf, setSelectedLeaf] = useState(null);
  
  const treeContainerRef = useRef(null);

  useEffect(() => {
    injectPhyloStyles();
  }, []);

  const onDrop = useCallback(async (acceptedFiles) => {
    const file = acceptedFiles[0];
    if (!file) return;

    setLoading(true);
    setError('');
    setData(null);
    setSelectedLeaf(null);

    const formData = new FormData();
    formData.append('file', file);

    try {
      const response = await axios.post(`${API_URL}/predict/phylo`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });
      setData(response.data);
    } catch (err) {
      setError(err.response?.data?.detail || err.message || 'An error occurred during phylogenetic tree building.');
    } finally {
      setLoading(false);
    }
  }, []);

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: { 'text/plain': ['.fasta', '.fa', '.fna', '.faa', '.txt'] },
    multiple: false
  });

  // Render tree using phylotree
  useEffect(() => {
    if (!data || !treeContainerRef.current) return;
    
    // Clear container
    treeContainerRef.current.innerHTML = '';
    
    try {
      const width = treeContainerRef.current.clientWidth || 700;
      const height = Math.max(data.num_sequences * 25, 400); // Dynamic height based on leaves count
      
      const tree = new phylotree(data.newick);
      
      // Render
      tree.render({
        container: treeContainerRef.current,
        width: width,
        height: height,
        "align-tips": true,
        zoom: true,
        "show-scale": true,
        "show-bootstrap": true
      });
      
      // Apply custom leaf node colors and click handlers
      const svg = d3.select(treeContainerRef.current).select("svg");
      
      // Select all leaf nodes and color them based on predictions
      svg.selectAll(".node.leaf").each(function() {
        const node = d3.select(this);
        const nameText = node.select("text").text();
        
        // Match name in predictions
        const prediction = data.predictions[nameText];
        if (prediction) {
          const score = prediction.score;
          // Interpolate color: Human (Green, #10b981) to Animal/Zoonotic (Red, #f43f5e)
          const color = score >= 50 ? '#f43f5e' : '#10b981';
          
          // Style text and add circle marker
          node.select("text")
            .style("fill", color)
            .style("font-size", "12px");
            
          node.append("circle")
            .attr("r", 5)
            .attr("cx", -8)
            .attr("cy", -2)
            .style("fill", color)
            .style("stroke", "#1e293b")
            .style("stroke-width", "1.5px");
            
          // Add click event
          this.onclick = () => {
            setSelectedLeaf({
              id: nameText,
              ...prediction
            });
          };
        }
      });
      
      // Color internal branches (D3 styling)
      svg.selectAll(".link").style("stroke", "#475569");
      
    } catch (err) {
      console.error("Error rendering tree with phylotree.js:", err);
      // Fallback is handled in the render JSX
    }
  }, [data]);

  return (
    <div style={{ textAlign: 'left' }}>
      <h2 style={{ fontSize: '1.75rem', fontWeight: '700', marginBottom: '1rem', color: '#c084fc' }}>
        Phylogenetic Spillover Map
      </h2>
      <p style={{ color: '#94a3b8', marginBottom: '2rem', fontSize: '0.95rem' }}>
        Upload a multi-sequence FASTA file of Rotavirus VP4. The app aligns them using **MAFFT**, constructs a maximum-likelihood phylogenetic tree using **IQ-TREE**, runs predictions, and colors the branches by zoonotic potential.
      </p>

      {/* Upload Box */}
      {!data && (
        <div 
          {...getRootProps()} 
          className={`dropzone ${isDragActive ? 'active' : ''}`}
          style={{
            border: '2px dashed rgba(192, 132, 252, 0.4)',
            borderRadius: '16px',
            padding: '3rem 2rem',
            textAlign: 'center',
            background: isDragActive ? 'rgba(192, 132, 252, 0.05)' : 'rgba(30, 41, 59, 0.3)',
            cursor: 'pointer',
            transition: 'border-color 0.2s, background-color 0.2s'
          }}
        >
          <input {...getInputProps()} />
          <UploadCloud size={48} style={{ color: '#c084fc', marginBottom: '1rem', animation: loading ? 'bounce 1s infinite' : 'none' }} />
          {loading ? (
            <div>
              <p style={{ color: '#e2e8f0', fontSize: '1.1rem', fontWeight: '600', marginBottom: '0.5rem' }}>
                Running Phylogenetic Analysis...
              </p>
              <p style={{ color: '#94a3b8', fontSize: '0.85rem' }}>
                Performing MAFFT Alignment & IQ-TREE Maximum-Likelihood modeling. This can take up to 45 seconds.
              </p>
              <div style={{ marginTop: '1.5rem', display: 'flex', justifyContent: 'center' }}>
                <RefreshCw style={{ animation: 'spin 1.5s linear infinite', color: '#c084fc' }} size={24} />
              </div>
            </div>
          ) : (
            <div>
              <p style={{ color: '#e2e8f0', fontSize: '1.1rem', fontWeight: '600' }}>
                Drag & drop multi-sequence FASTA file here, or click to select
              </p>
              <p style={{ color: '#94a3b8', fontSize: '0.85rem', marginTop: '0.5rem' }}>
                Supports FASTA format (.fasta, .fa, .txt). Limit: 2 to 100 sequences.
              </p>
            </div>
          )}
        </div>
      )}

      {error && <div className="error-msg" style={{ marginTop: '1rem' }}>{error}</div>}

      {/* Results View */}
      {data && (
        <div style={{ display: 'grid', gridTemplateColumns: '1.5fr 1fr', gap: '2rem' }}>
          
          {/* Left Panel: The Tree rendering container */}
          <div>
            <div className="glass-panel" style={{ padding: '1rem', marginBottom: '1rem', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
              <div>
                <span style={{ fontSize: '0.8rem', color: '#94a3b8' }}>Substitution Model Selected: </span>
                <span style={{ fontSize: '0.9rem', fontWeight: '700', color: '#c084fc' }}>{data.model_used}</span>
                <span style={{ margin: '0 0.75rem', color: 'rgba(255,255,255,0.1)' }}>|</span>
                <span style={{ fontSize: '0.8rem', color: '#94a3b8' }}>MSA Length: </span>
                <span style={{ fontSize: '0.9rem', fontWeight: '700', color: '#60a5fa' }}>{data.alignment_length} bp</span>
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
                Upload New FASTA
              </button>
            </div>

            <div 
              className="phylotree-container"
              ref={treeContainerRef} 
              style={{ 
                minHeight: '450px', 
                background: '#0f172a',
                padding: '1rem',
                borderRadius: '16px',
                border: '1px solid rgba(255, 255, 255, 0.1)',
                boxShadow: '0 8px 32px 0 rgba(0, 0, 0, 0.3)'
              }}
            />
            
            {/* Legend */}
            <div style={{ marginTop: '1rem', display: 'flex', gap: '1.5rem', background: 'rgba(30, 41, 59, 0.3)', padding: '0.75rem 1rem', borderRadius: '8px', border: '1px solid rgba(255, 255, 255, 0.05)' }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <div style={{ width: '12px', height: '12px', borderRadius: '50%', background: '#10b981' }} />
                <span style={{ fontSize: '0.8rem', color: '#cbd5e1' }}>Human-Adapted (&lt;50% zoonotic potential)</span>
              </div>
              <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                <div style={{ width: '12px', height: '12px', borderRadius: '50%', background: '#f43f5e' }} />
                <span style={{ fontSize: '0.8rem', color: '#cbd5e1' }}>Zoonotic Potential / Animal-Adapted (&ge;50% potential)</span>
              </div>
              <div style={{ fontSize: '0.8rem', color: '#94a3b8', marginLeft: 'auto' }}>
                * Click on leaf nodes to view prediction analysis.
              </div>
            </div>
          </div>

          {/* Right Panel: Selected leaf details */}
          <div className="glass-panel" style={{ padding: '1.5rem', display: 'flex', flexDirection: 'column', height: 'fit-content' }}>
            {selectedLeaf ? (
              <div>
                <h3 style={{ fontSize: '1.25rem', fontWeight: '700', color: '#cbd5e1', marginBottom: '1.25rem', overflowWrap: 'break-word' }}>
                  {selectedLeaf.id}
                </h3>
                
                {/* Score badge */}
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.75rem', marginBottom: '1.5rem' }}>
                  <div style={{
                    padding: '0.75rem 1.25rem',
                    background: selectedLeaf.is_human_adapted ? 'rgba(16, 185, 129, 0.1)' : 'rgba(244, 63, 94, 0.1)',
                    border: selectedLeaf.is_human_adapted ? '1px solid rgba(16, 185, 129, 0.2)' : '1px solid rgba(244, 63, 94, 0.2)',
                    borderRadius: '8px',
                    color: selectedLeaf.is_human_adapted ? '#10b981' : '#f43f5e',
                    fontSize: '1.5rem',
                    fontWeight: '700'
                  }}>
                    {selectedLeaf.score.toFixed(1)}%
                  </div>
                  <div>
                    <span style={{ fontSize: '0.75rem', textTransform: 'uppercase', color: '#94a3b8', display: 'block' }}>Prediction</span>
                    <span style={{ fontSize: '0.95rem', fontWeight: '600', color: '#f8fafc' }}>
                      {selectedLeaf.is_human_adapted ? 'Human-Adapted' : 'Animal-Adapted / Zoonotic'}
                    </span>
                  </div>
                </div>

                <div style={{ display: 'flex', flexDirection: 'column', gap: '0.75rem', marginBottom: '1.5rem', fontSize: '0.85rem' }}>
                  <div style={{ display: 'flex', justifyContent: 'space-between', borderBottom: '1px solid rgba(255,255,255,0.03)', paddingBottom: '0.5rem' }}>
                    <span style={{ color: '#94a3b8' }}>VP8* Sequence Length:</span>
                    <span style={{ fontWeight: '600', color: '#f8fafc' }}>{selectedLeaf.vp8_length} aa</span>
                  </div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', borderBottom: '1px solid rgba(255,255,255,0.03)', paddingBottom: '0.5rem' }}>
                    <span style={{ color: '#94a3b8' }}>Alignment Range:</span>
                    <span style={{ fontWeight: '600', color: '#f8fafc' }}>Residues {selectedLeaf.alignment_start} – {selectedLeaf.alignment_end}</span>
                  </div>
                </div>

                {/* Biological Interpretation */}
                <div style={{
                  background: 'linear-gradient(135deg, rgba(192, 132, 252, 0.15) 0%, rgba(244, 72, 182, 0.05) 100%)',
                  padding: '1.25rem',
                  borderRadius: '12px',
                  border: '1px solid rgba(192, 132, 252, 0.2)'
                }}>
                  <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '0.75rem' }}>
                    <FileText size={18} style={{ color: '#c084fc' }} />
                    <span style={{ fontWeight: '600', color: '#c084fc', fontSize: '0.95rem' }}>Biological Interpretation</span>
                  </div>
                  <p style={{ fontSize: '0.85rem', lineHeight: '1.45', color: '#cbd5e1', margin: 0 }}>
                    {selectedLeaf.interpretation}
                  </p>
                </div>
              </div>
            ) : (
              <div style={{ textAlign: 'center', padding: '3rem 1rem', color: '#94a3b8' }}>
                <Info size={36} style={{ color: '#475569', marginBottom: '1rem' }} />
                <p style={{ fontSize: '0.9rem', margin: 0 }}>
                  Click on any colored leaf node in the phylogenetic tree to view its host adaptation prediction and molecular interpretation details.
                </p>
              </div>
            )}
          </div>
          
        </div>
      )}
    </div>
  );
}

export default PhyloTab;
