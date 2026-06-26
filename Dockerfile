FROM python:3.12-slim

# Install Node.js, MAFFT, and wget
RUN apt-get update && apt-get install -y curl mafft wget && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash -
RUN apt-get install -y nodejs

# Install IQ-TREE 2 Linux binary
RUN wget -q https://github.com/iqtree/iqtree2/releases/download/v2.3.6/iqtree-2.3.6-Linux-intel.tar.gz \
    && tar -xzf iqtree-2.3.6-Linux-intel.tar.gz \
    && cp iqtree-2.3.6-Linux-intel/bin/iqtree2 /usr/local/bin/iqtree2 \
    && rm -rf iqtree-2.3.6-Linux-intel*

# Set up the working directory
WORKDIR /app

# Copy the requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

# Pre-download the static reference PDB structure
RUN mkdir -p /app/WebApp/backend/static \
    && wget -q -O /app/WebApp/backend/static/2dwr.pdb "https://files.rcsb.org/download/2DWR.pdb"

# Build the React frontend
WORKDIR /app/WebApp/frontend
RUN npm install
RUN npm run build

# Return to the main app directory
WORKDIR /app

# Expose port 7860 (Hugging Face Spaces default)
EXPOSE 7860

# Start the FastAPI server
CMD ["uvicorn", "WebApp.backend.main:app", "--host", "0.0.0.0", "--port", "7860"]
