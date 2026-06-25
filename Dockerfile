FROM python:3.12-slim

# Install Node.js
RUN apt-get update && apt-get install -y curl
RUN curl -fsSL https://deb.nodesource.com/setup_22.x | bash -
RUN apt-get install -y nodejs

# Set up the working directory
WORKDIR /app

# Copy the requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

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
