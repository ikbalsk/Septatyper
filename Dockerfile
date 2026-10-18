# SeptaTyper v2, flat repo layout: every file sits in the repo root, the image rebuilds the folders.
FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends ncbi-blast+ && rm -rf /var/lib/apt/lists/*
WORKDIR /srv
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN mkdir -p app/static scripts data refs && touch app/__init__.py
COPY main.py rules.py septatype.py grammar.py config.py bvbrc.py mlst.py persist.py app/
COPY index.html style.css app.js app/static/
COPY build_db.py fetch_bvbrc.py scripts/
COPY Table_S1.xlsx data/
COPY gyrA.fasta parC.fasta recN.fasta qnrB2.fasta qnrS1.fasta refs/
COPY mlst_ecoli_achtman_4.fasta mlst_ecoli_achtman_4.txt mlst_senterica_achtman_2.fasta mlst_senterica_achtman_2.txt mlst_klebsiella.fasta mlst_klebsiella.txt refs/
RUN python scripts/build_db.py
ENV PORT=8000 DATA_DIR=/srv/data
EXPOSE 8000
CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port ${PORT}"]
