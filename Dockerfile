# SeptaTyper v2, flat repo layout: every file sits in the repo root, the image rebuilds the folders.
FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends ncbi-blast+ && rm -rf /var/lib/apt/lists/*
WORKDIR /srv
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN mkdir -p app/static scripts data refs && touch app/__init__.py
COPY main.py rules.py septatype.py grammar.py config.py app/
COPY index.html style.css app.js app/static/
COPY build_db.py fetch_bvbrc.py scripts/
COPY Table_S1.xlsx data/
COPY gyrA_cds_ref.fasta parC_cds_ref.fasta recN_cds_ref.fasta qnrB_ref.fasta qnrS_ref.fasta refs/
RUN python scripts/build_db.py
ENV PORT=8000 DATA_DIR=/srv/data
EXPOSE 8000
CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port ${PORT}"]
