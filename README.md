# SeptaTyper Web Database

Ciprofloxacin resistance prediction plus a browsable repository of 2,734
Enterobacterales isolates. Two modules at NAR Database quality:

- **Module A (Predict):** upload a nucleotide FASTA, get the septatype and
  resistance call directly via BLAST+.
- **Module B (Browse):** searchable isolate repository with BV-BRC genome links,
  septatype frequency charts, geographic breakdown, REST API, and CSV/JSON downloads.

## What is inside

- `app.py` - FastAPI backend (prediction + REST API v1)
- `septatype.db` - SQLite with 2,734 isolates and 119 septatype rules
- `references/` - 5 verified NCBI reference sequences (gyrA, parC, recN, qnrB, qnrS)
- `static/index.html` - phone-friendly web interface (Predict, Browse, Charts, API, About)
- `Dockerfile`, `requirements.txt`
- `fetch_references.py` - re-downloads references from NCBI if needed
- `build_db.py` - rebuilds septatype.db from Table_S1.xlsx (needs the xlsx)

## REST API

- `GET /api/v1/isolates?species=E.%20coli&site=urine&phenotype=Resistant&page=1&per_page=50`
- `GET /api/v1/isolates/562.7710` (includes BV-BRC URL)
- `GET /api/v1/rules?species=K.%20pneumoniae&site=urine`
- `GET /api/v1/stats/septatypes`, `/api/v1/stats/geography`, `/api/v1/stats/accuracy`
- `GET /api/v1/download/csv`, `/api/v1/download/json`
- `POST /api/v1/predict` (form: species, site, fasta file)

## Deploy (phone-friendly, about 10 minutes)

You deploy once from a computer, then use it from your phone's browser.

### Option A: Render.com (free, no credit card)

1. Put this folder in a GitHub repo (drag and drop the files on github.com).
2. Go to render.com, sign up, click New, then Web Service.
3. Connect your repo. Settings:
   - Environment: Docker
   - Plan: Free
4. Click Create. Wait about 5 minutes for the first build.
5. Open the URL Render gives you (like `septatype.onrender.com`) on your phone.
   Bookmark it. Upload a FASTA and you get the call.

Note: the free plan sleeps after inactivity, so the first upload after a
while takes ~30 seconds to wake up.

### Option B: Run on your own computer (then open from phone)

1. Install Docker Desktop.
2. In this folder, run:
   ```
   docker build -t septatype .
   docker run -p 8000:8000 septatype
   ```
3. On your phone (same Wi-Fi), open `http://YOUR-COMPUTER-IP:8000`.
   Find your computer's IP in Wi-Fi settings.

### Option C: Run without Docker

Needs Python 3.10+ and BLAST+ installed (`apt install ncbi-blast+`).

```
pip install -r requirements.txt
BLAST_BIN=/usr/bin uvicorn app:app --host 0.0.0.0 --port 8000
```

## Using it from your phone

1. Open the server URL in your phone browser.
2. Pick the species and isolation site.
3. Tap "Choose file" and pick your FASTA file (or paste sequences).
4. Tap "Predict phenotype".
5. You get a big Resistant / Susceptible / Unknown badge, the septatype,
   supporting isolate counts, and details of which residues and genes were found.

Max upload is 20 Mb per request.

## How the sequence analysis works

- BLAST+ (`blastn`) searches your FASTA for each reference gene.
- gyrA/parC/recN: the best hit's alignment is walked codon by codon; the
  codons at GyrA 83/87, ParC 80/84, RecN 346 are translated to amino acids.
  A codon spanning a gap or containing N gives "?" and the call becomes Unknown.
- qnrB/qnrS: present if a hit reaches 80% identity AND 90% coverage of the
  reference gene; otherwise absent.
- The 7 values form the septatype, looked up in the SQLite rules table:
  majority phenotype wins, ties go to Susceptible, unseen septatypes return
  Unknown (never forced).

## Reference sequences (all from NCBI, none invented)

- gyrA: E. coli K-12 MG1655 NC_000913.3 complement(2336793..2339420),
  protein NP_416734.1; verified codons 83=TCG(S), 87=GAC(D)
- parC: E. coli K-12 MG1655 NC_000913.3 complement(3163715..3165973),
  protein NP_417491.1; verified codons 80=AGC(S), 84=GAA(E)
- recN: E. coli K-12 MG1655 NC_000913.3 2751795..2753456,
  protein YP_026172.1; verified codon 346=CTT(L)
- qnrB: DQ351241.1 (K. pneumoniae plasmid pMG298 qnrB1, 681 bp)
- qnrS: LC909223.1 (E. coli qnrS1, 655 bp)

To re-download: `python3 fetch_references.py`

## Tested

- Synthetic mutant genome (gyrA S83L+D87N, parC S80I, qnrS1 present):
  correctly read L-N-I-E-L-0-1 and called Resistant.
- Wild-type references: correctly read S-D-S-E-L-0-0 and called Susceptible.
- All 7 per-group accuracies recomputed from the SQLite rules match the
  published values (98.9, 98.1, 99.4, 99.8, 98.7, 96.8, 97.2).

Not tested: a full multi-megabase genome upload through the phone browser
(no test genome available here); the 20 Mb cap and BLAST+ timeouts are set
for typical use. Research use only, not for clinical decisions.
