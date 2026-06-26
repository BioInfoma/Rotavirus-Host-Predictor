import sys

with open('WebApp/frontend/src/App.jsx', 'r') as f:
    content = f.read()

# 1. Add showAdvanced state
if 'const [showAdvanced, setShowAdvanced] = useState(false);' not in content:
    content = content.replace(
        "const [pastedSeq, setPastedSeq] = useState('');",
        "const [pastedSeq, setPastedSeq] = useState('');\n  const [showAdvanced, setShowAdvanced] = useState(false);"
    )

# 2. Update Biological Interpretation styling
old_interp = """              {/* Biological Interpretation */}
              {result.interpretation && (
                <div className="glass-panel" style={{ marginTop: '2rem', padding: '1.5rem', background: 'linear-gradient(135deg, rgba(59, 130, 246, 0.15), rgba(139, 92, 246, 0.15))', textAlign: 'left', borderLeft: '4px solid #8b5cf6' }}>
                  <h4 style={{ margin: '0 0 0.75rem 0', color: '#c4b5fd', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                    <FileText size={18} /> Biological Interpretation
                  </h4>
                  <p style={{ color: '#e2e8f0', fontSize: '0.95rem', lineHeight: '1.7', margin: 0 }}>
                    {result.interpretation}
                  </p>
                </div>
              )}"""

new_interp = """              {/* Biological Interpretation */}
              {result.interpretation && (
                <div className="glass-panel" style={{ marginTop: '2.5rem', padding: '2rem', background: 'linear-gradient(135deg, rgba(59, 130, 246, 0.2), rgba(139, 92, 246, 0.2))', textAlign: 'left', border: '1px solid #8b5cf6', boxShadow: '0 4px 20px rgba(139, 92, 246, 0.15)' }}>
                  <h3 style={{ margin: '0 0 1rem 0', color: '#c4b5fd', display: 'flex', alignItems: 'center', gap: '0.75rem', fontSize: '1.4rem' }}>
                    <FileText size={24} /> Plain English Breakdown
                  </h3>
                  <p style={{ color: '#f8fafc', fontSize: '1.1rem', lineHeight: '1.8', margin: 0 }}>
                    {result.interpretation}
                  </p>
                </div>
              )}"""

content = content.replace(old_interp, new_interp)

# 3. Add toggle button and wrap SHAP section
old_shap_header = """              {/* SHAP Feature Importance */}
              <div className="glass-panel" style={{ marginTop: '2rem', padding: '1.5rem', background: 'rgba(30, 41, 59, 0.5)', textAlign: 'left' }}>
                <h4 style={{ margin: '0 0 0.5rem 0', color: '#e2e8f0' }}>Top Factors Driving Prediction</h4>"""

new_shap_header = """              <div style={{ textAlign: 'center', marginTop: '2rem' }}>
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
                <h4 style={{ margin: '0 0 0.5rem 0', color: '#e2e8f0' }}>Key Amino Acid Mutations Driving the Result (SHAP)</h4>"""

content = content.replace(old_shap_header, new_shap_header)

# Wrap end of SHAP section
old_shap_footer = """                  </div>
                </div>
              </div>

              {/* Sequence Alignment Visualization */}"""

new_shap_footer = """                  </div>
                </div>
              </div>
              )}

              {/* Sequence Alignment Visualization */}"""

content = content.replace(old_shap_footer, new_shap_footer)

with open('WebApp/frontend/src/App.jsx', 'w') as f:
    f.write(content)
