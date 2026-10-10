FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends ncbi-blast+ && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY app.py .
COPY septatype.db .
RUN mkdir -p references static
COPY gyrA_cds_ref.fasta parC_cds_ref.fasta recN_cds_ref.fasta qnrB_ref.fasta qnrS_ref.fasta references/
COPY index.html static/
ENV BLAST_BIN=/usr/bin
EXPOSE 8000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
