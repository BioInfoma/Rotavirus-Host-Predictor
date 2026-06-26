import sys
with open('WebApp/backend/main.py', 'r') as f:
    lines = f.readlines()
out = []
skip = False
for l in lines:
    if l.startswith('class PhyloNodePrediction(BaseModel):'):
        skip = True
    elif l.startswith('class StructureResidue(BaseModel):'):
        skip = False
    
    if l.startswith('@app.post("/predict/phylo"'):
        skip = True
    elif l.startswith('@app.post("/predict/structure"'):
        skip = False
        
    if not skip:
        out.append(l)

with open('WebApp/backend/main.py', 'w') as f:
    f.writelines(out)
