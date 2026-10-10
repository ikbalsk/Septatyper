# Septatype server - ciprofloxacin resistance prediction from FASTA
FROM python:3.12-slim

# BLAST+ (system package, Debian)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python deps
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# App files
COPY app.py .
COPY septatype.db .
COPY references/ ./references/
COPY static/ ./static/

ENV BLAST_BIN=/usr/bin
EXPOSE 8000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
