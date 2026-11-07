# CiproTyper: Software Architecture, Database Schema, Algorithms, and REST API Reference

**Documentation:** Return to [`README.md`](README.md) for the biological overview, ciprotype rules, checkpoint explanations, and user guidelines.  
**Live Server:** [https://septatyper-is.onrender.com/](https://septatyper-is.onrender.com/) | **Swagger / OpenAPI UI:** [https://septatyper-is.onrender.com/docs](https://www.google.com/search?q=https://septatyper-is.onrender.com/docs)  
**GitHub Repository:** [https://github.com/ikbalsk/Septatyper](https://www.google.com/search?q=https://github.com/ikbalsk/Septatyper) | **Contact:** [ikbalsk.personal@gmail.com](mailto:ikbalsk.personal@gmail.com)

## 1. Software Stack and Repository Layout

CiproTyper is built as a lightweight, self-contained Python web service and SQLite repository packaged inside a Docker container.

* **Backend Runtime:** Python `3.12` (`python:3.12-slim` Debian base image), FastAPI `0.162.0`, Starlette `1.10.0`, Uvicorn `0.54.0` (`uvicorn[standard]`), `python-multipart==0.0.32`.
* **Bioinformatics Engine:** NCBI BLAST+ (`makeblastdb` and `blastn` installed via `apt-get install ncbi-blast+`), Biopython `1.88` (`Bio.SeqIO`, `Bio.Seq`, `Bio.Align.PairwiseAligner`), `openpyxl==3.1.5` (for parsing `Table_S1.xlsx`).
* **Frontend Interface:** Vanilla JavaScript (`app.js`, no frontend build step or external framework), semantic HTML5 (`index.html`), responsive light-blue and white CSS (`style.css`), and locally bundled `plotly.min.js` (`v2.35.2`, with CDN fallback protected by Subresource Integrity hash `sha384-cCVCZkAjYNxaYKbM8lsArLznDF/SvMFr1jcZrvOpSTCa0W40ZAdLzHCEulnUa5i7`).
* **Testing Suite:** `pytest==9.1.1` and `httpx==0.28.1` (`requirements-dev.txt`), comprising **54 automated test functions** covering pipeline accuracy reproduction, grammar checkpoints, reverse-strand and indel handling, gzip uploads, MLST typing, curation state transitions, and simulated BV-BRC and GitHub backup APIs.

### Flat Repository to Container Folder Mapping (`Dockerfile`)

So that the repository can be updated directly through the GitHub web upload interface from a laptop or mobile phone, all files sit in a flat layout in the repository root. During `docker build`, the `Dockerfile` organizes them into `/srv`:

| Container Path | Files Copied from Repository Root | Role in the System |
| :--- | :--- | :--- |
| `/srv/app/` | `main.py`, `rules.py`, `septatype.py`, `grammar.py`, `config.py`, `bvbrc.py`, `mlst.py`, `persist.py`, `countries.py`, `__init__.py` | FastAPI endpoints, BLAST/alignment pipeline, 3-checkpoint grammar, SQLite access, fast MLST engine, BV-BRC client, and GitHub cloud backup. |
| `/srv/app/static/` | `index.html`, `style.css`, `app.js`, `plotly.min.js` | Single-page application with hash-routed views (`#/`, `#/predict`, `#/browse`, `#/analytics`, `#/submit`, `#/theory`, `#/data`, `#/about`, `#/curate`). |
| `/srv/scripts/` | `build_db.py`, `fetch_bvbrc.py` | `build_db.py` parses `Table_S1.xlsx` at Docker build time to create `/srv/data/seed.db`; `fetch_bvbrc.py` is the CLI BV-BRC updater. |
| `/srv/data/` | `Table_S1.xlsx`, `seed.db` (and live `septatype.db` in `DATA_DIR`) | Founding 2,734-isolate Excel dataset, read-only baked `seed.db`, live SQLite database (`septatype.db`), and rotated `.bak` backups. |
| `/srv/refs/` | `gyrA.fasta`, `parC.fasta`, `recN.fasta`, `qnrB2.fasta`, `qnrS1.fasta`, `mlst_ecoli_achtman_4.{fasta,txt}`, `mlst_senterica_achtman_2.{fasta,txt}`, `mlst_klebsiella.{fasta,txt}` | Species-specific chromosomal CDS references, plasmid `qnr` alleles, and PubMLST allele sequences and profile tables (snapshot `2026-03-11`). |

## 2. Internal Algorithms and Implementation Details

### 2.1 Sequence Alignment, Codon Mapping, and `qnr` Screening (`app/septatype.py`)

1. **Starlette Form Parser Patch (`_raise_form_field_limit`):** By default, Starlette rejects ordinary form fields larger than 1 MB (`1,048,576` bytes), which blocks pasted bacterial genomes. On startup, `main.py` patches `MultiPartParser` and `FormParser` to allow form fields up to `(MAX_PASTE_MB + 1) MB` (`11 MB`).
2. **Concurrency Semaphore and In-Memory SHA-256 Cache:** To run reliably on small free-tier cloud instances, `main.py` limits concurrent BLAST analyses using `BLAST_SLOTS = threading.BoundedSemaphore(2)` (returning HTTP `429` if busy) and caches up to `_CACHE_MAX = 32` recent results in `_CACHE` keyed by `(sha256(text), species, qnrB2_override, qnrS1_override)` under a thread lock (`_CACHE_LOCK`).
3. **BLAST+ Execution (`_blast`):** Builds a temporary nucleotide database (`makeblastdb -dbtype nucl -parse_seqids`) inside `tempfile.TemporaryDirectory()` and runs `blastn -task blastn -evalue 1e-20 -dust no -max_target_seqs 50` with a timeout of `BLAST_TIMEOUT_S = 180` seconds.
4. **Contig Extraction and Biopython Global Alignment (`_extract_region`, `_ref_to_region_map`):** For the highest-bitscore hit of `gyrA`, `parC`, and `recN` with `pident >= 75.0%` (`MIN_ID_CHROM`), the matching contig slice is extended by the unaligned reference ends plus `FLANK = 30` bp on each side (and reverse-complemented if `sstart > send`). It is then globally aligned to the species reference CDS using `Bio.Align.PairwiseAligner` (`mode = "global"`, `match_score = 2`, `mismatch_score = -3`, `open_gap_score = -8`, `extend_gap_score = -1`, `end_gap_score = 0`).
5. **Exact 3-Base Codon Extraction (`_codon`):** From the alignment coordinate array `m`, for target 1-based amino-acid position `pos`, nucleotide indices `i = (pos - 1) * 3`, `i + 1`, and `i + 2` must map to valid, strictly consecutive indices (`b == a + 1` and `c == b + 1`). Any gap returns `None`; any base outside `ACGT` translates to `"X"`; stop codons translate to `"*"`.

### 2.2 Fast Locate-Then-Lookup MLST Engine and Species Gate (`app/mlst.py`, `app/septatype.py`)

BLASTing all ~13,000 PubMLST alleles against a 5 Mb assembly takes ~7.6 seconds. `app/mlst.py` reduces this to **~0.1 seconds** (0.6 seconds total pipeline time):

1. **In-Memory Caching (`@lru_cache(maxsize=4)`):** Loads the 7-locus profile table (`_profiles`) and builds an exact sequence-to-allele hash map (`lookup[locus][seq] = allele`) plus a list of **4 evenly spaced representative alleles per locus** (`N_REPS = 4`, 28 sequences total in `mlst_reps.fa`).
2. **Pass 1 (`mlstfast`):** BLASTs only the 28 representative alleles against the genome. For each locus, the best hit's coordinates define the exact full-length allele window (`_window`) on the contig, which is sliced out and looked up in `lookup[locus]` in $O(1)$ time.
3. **Pass 2 (`mlstslow` fallback):** Only loci whose extracted window did not match a known allele in `lookup[locus]` are BLASTed against the full allele FASTA for those specific loci (`_call_slow`), identifying either an exact match (`pident >= 99.999%`, `qstart == 1`, `qend == qlen`), a novel nearest allele (`~allele` at `qcov >= 0.90` and `pident >= 90.0%`), or a missing locus (`"-"`).
4. **Species Verification Gate (`MLST_SPECIES_MIN_ID = 90.0%`):** For genome assemblies (`>= 1,000,000` bp), all 7 MLST housekeeping genes must match the selected species scheme at **`>= 90.0%` nucleotide identity** (reduced from `95.0%`). When any MLST locus falls below `90.0%` identity, the pipeline displays *"Please check the species selected and resubmit"*, **does not assign a Sequence Type (ST), does not apply the 3 checkpoints, and returns no phenotype prediction**.

### 2.3 Three-Checkpoint Grammar Evaluation (`app/grammar.py`)

When a valid isolate passes the `90.0%` species and chromosomal identity gates (`CP1_MIN_ID = 90.0%`) and carries a novel ciprotype unseen in the Tier 1 `septatypes` reference table, `app/grammar.py` evaluates the three checkpoints:

| Amino acid Position | Susceptible | Resistant |
| :--- | :--- | :--- |
| **GyrA-83** | **S**,T,**Y** — Polar with -OH group | A,**L**,**I**,V,**F**,W,M,C,G — Hydrophobic |
| **GyrA-87** | **D**,E,**N**,Q — Negatively Charged or its derivative | **D**,E,**N**,Q,S,**Y**,T,**G**, **A**,L,I,V,F,W,M,C — Negatively charged or uncharged |
| **ParC-80** | **S**,T,**Y** — Polar with -OH group | G,A,L,**I**,V,**F**,W,M,C — Hydrophobic |
| **ParC-84** | **D**,**E** — Negatively Charged | D,**E**,**G**,A,L,I,**V**,F,W,M,C — Negatively Charged or hydrophobic |

*\*Residues observed in our dataset are marked in bold and remaining are possible residues based on physicochemical features of the observed one.*

* **Checkpoint 1 Single-Column Rule:** To assign a phenotype on Checkpoint 1, **all 4 positions of the sample's tetratype (`GyrA-83`, `GyrA-87`, `ParC-80`, `ParC-84`) must fall within one single column** (either all 4 are in the `Resistant` column -> **Resistant**, or all 4 are in the `Susceptible` column -> **Susceptible**; shared residues such as `D,E,N,Q` at `GyrA-87` and `D,E` at `ParC-84` appear in both columns, so one column must still contain all 4 residues of the tetratype). If the 4 residues do **not** all fall within one single column (mixed across columns so neither column carries all 4 residues), **Checkpoint 1 is not applied** (no hydrophobicity `<H>` value or tie-breaker is used), and the phenotype is decided solely by **Checkpoint 2 and Checkpoint 3**.
* **Checkpoint 2 (`RecN-346`):** Calls **Resistant** if `RecN-346` is `F` or `L` encoded by `ctg`; calls **Susceptible** if `RecN-346` is `L` encoded by `cta` (or other/unknown codon).
* **Checkpoint 3 (`qnrB2` and `qnrS1`):** Calls **Resistant** if either `qnrB2` or `qnrS1` (or both) is present (`1`); calls **Susceptible** if both are absent (`0`).
* **Combined Verdict:** Returns **Resistant** if **any evaluated checkpoint** is Resistant (for example, when Checkpoint 1 is not applied, if Checkpoint 2 is Resistant and Checkpoint 3 is Susceptible, the phenotype is Resistant), and returns **Susceptible** **only when all evaluated checkpoints are Susceptible**.

## 3. SQLite Database Schema and Persistence Architecture

### 3.1 Seven SQLite Tables (`seed.db` and `septatype.db`)

| Table Name | Primary Key | Columns & Purpose |
| :--- | :--- | :--- |
| **`isolates`** | `genome_id` | `31` columns (`ISO_COLS` + `year_note`): `genome_id`, `genome_name`, `isolate_name`, `species`, `habitat`, `lab_phenotype`, `year_raw`, `year_start`, `year_end`, `year_note`, `country`, `iso3`, `st`, `gyrA_83`, `gyrA_87`, `parC_80`, `parC_84`, `recN_346`, `qnrB2`, `qnrS1`, `septatype`, `tetratype`, `septatyper_call`, `bvbrc_call`, `bvbrc_method`, `source` (`'founding'` or `'community'`), `date_added`, `submitter`, `submitted_date`, `affiliation`, `publication`. Indexed on `(species, habitat, septatype)`, `(st)`, and `(country)`. |
| **`septatypes`** | `(species, habitat, septatype)` | Tier 1 lookup rule table: `species`, `habitat`, `septatype`, `tetratype`, `n_isolates`, `n_resistant`, `n_susceptible`, `predicted` (`'R'` or `'S'`), `confidence` (`max(R, S)/n_isolates`). Stores both source-specific rows and pooled `habitat = 'ALL'` rows. |
| **`species_stats`** | `(species, habitat)` | Cohort summary metrics: `species`, `habitat`, `n_isolates`, `n_resistant`, `n_septatypes`, `accuracy`, `grammar_accuracy`. |
| **`aa_properties`** | `aa` | The 20 standard amino acids: `aa`, `name`, `class` (polar with -OH group, negatively charged, hydrophobic, etc.). |
| **`submissions`** | `job_id` | Community submission queue (`26` columns): `job_id` (12-char hex), `created`, `status` (`needs_review`, `pending`, `approved`, `rejected`), `species`, `habitat`, `genome_id`, `genome_name`, `country`, `year`, `source_text`, `lab_phenotype`, `email`, `submitter`, `affiliation`, `publication`, `septatype`, `tier`, `call`, `novel`, `matches`, `result_json`, `curator_note`, `decided`, `decided_by`, `consent_version`, `consent_at`. |
| **`audit_log`** | `id` (AUTOINCREMENT) | Curator decision trail: `id`, `at` (UTC timestamp), `who` (curator name), `ip`, `action` (`approve` or `reject`), `job_id`, `detail`. |
| **`meta`** | `key` | Build and state metadata: `built`, `seed_built`, `source` (`Table_S1.xlsx`), `n_isolates` (`2734`), `bvbrc_loaded` (`0` or `1`). |

### 3.2 Startup Migration, Seed Refresh, and Automatic GitHub Cloud Backup

1. **Idempotent Schema Migration (`_migrate` in `rules.py`):** Runs on startup to add any missing columns to older `isolates` or `submissions` tables, create `audit_log`, and populate founding isolate attribution (`FOUNDING_SUBMITTER = "Ghosh et al"`, `FOUNDING_DATE = "2026-10-16"`, `FOUNDING_AFFILIATION = "JISIASR, Kolkata"`).
2. **Non-Destructive Seed Refresh (`_refresh_from_seed` in `rules.py`):** When `REFRESH_DB_ON_START=1` and the baked `seed.db` has a newer `built` timestamp than `seed_built` in the live `septatype.db`, the server creates a timestamped `.bak` copy in `DATA_DIR/backups/` (keeping the last 3 backups), replaces `source='founding'` rows and base tables from `seed.db` inside a transaction, and replays all `source='community'` isolates into `septatypes`.
3. **8-Second Debounced GitHub Cloud Backup (`app/persist.py`):** When `BACKUP_GITHUB_REPO` and `BACKUP_GITHUB_TOKEN` are set, any new submission or curator decision triggers `persist.schedule_save()`, which waits `DELAY_S = 8` seconds and pushes the `submissions` table as base64-encoded JSON to the private GitHub repository via `PUT /repos/{owner}/{repo}/contents/{path}`. On startup, if `submissions` is empty, `restore_if_configured()` downloads the JSON backup, inserts missing submissions, and replays `rules.apply_approval()` for all `approved` records.

## 4. Complete REST API Specification (`/api/v1`)

Base URL: `https://septatyper-is.onrender.com/api/v1` (Interactive OpenAPI documentation at `/docs`).

| Method & Endpoint | Auth / Rate Limit | Parameters / Request Payload | Response Summary |
| :--- | :--- | :--- | :--- |
| `GET /api/v1/meta` | Public | None | Supported species, habitats, 12 sources, 249 countries, isolate counts, reference file status, upload/paste limits (`60`/`10` MB), `consent_version`, and filter facets. |
| `POST /api/v1/predict` | Public (`RATE_PREDICT_PER_HOUR`, max 2 concurrent BLASTs) | Form data: `species` (req), `habitat` (`"ALL"`), `file` or `sequence`, `qnrB2` (`"auto"`\|`"present"`\|`"absent"`), `qnrS1` (`"auto"`\|`"present"`\|`"absent"`). | Full sequence analysis (`genes`, `positions`, `qnr`, `checkpoint1`, `recN_codon`, `mlst`) and phenotype `prediction` (`call`, `tier`, `confidence`, `reference_sts`, `criteria`, `mechanism`). Halts without ST or phenotype if MLST species check is `< 90.0%`. |
| `POST /api/v1/predict/bvbrc` | Public (`limit_predict`) | JSON: `{"genome_id": "562.7710", "species": "auto", "fallback_species": "", "habitat": "ALL"}`. | Fetches assembly FASTA and metadata from BV-BRC API, auto-detects species, runs pipeline, flags if in founding set, and returns prediction + `bvbrc` block. |
| `POST /api/v1/predict/manual` | Public | JSON: `{"species": "ecoli", "habitat": "ALL", "gyrA_83": "L", "gyrA_87": "N", "parC_80": "I", "parC_84": "E", "recN_346": "L", "qnrB2": false, "qnrS1": false, "recN_codon": "CTA"}`. | Evaluates 5-residue profile without BLAST; returns Tier 1 lookup and Tier 2 `criteria` compared against wild-type `SDSE`. |
| `GET /api/v1/isolates` | Public | Query: `species`, `habitat`, `st`, `country`, `septatype`, `tetratype`, `phenotype`, `year_from`, `year_to`, `wins`, `source`, `q`, `sort`, `desc`, `page`, `per_page` (`1..1000`). | Paginated isolates with metadata, `date_display`, `star` (CiproTyper win over BV-BRC), and `bvbrc_url`. |
| `GET /api/v1/isolates/{genome_id}` | Public | Path: `genome_id` | Single isolate record (or HTTP `404`). |
| `GET /api/v1/rules` | Public | Query: `species` (req), `habitat` (`"ALL"`). | Returns `septatypes` rule table rows and `species_stats` for the species and source. |
| `GET /api/v1/grammar` | Public | Query: `residues` (4 letters), `species` (`ecoli`), `reference` (`SDSE`), `recn` (`L`), `codon` (`""`), `qnrB` (`false`), `qnrS` (`false`). | Evaluates the 3 checkpoints (single-column Checkpoint 1 table, Checkpoint 2 `RecN-346`, and Checkpoint 3 `qnrB2`/`qnrS1`) and returns `positions` and `checkpoints` blocks. |
| `GET /api/v1/stats/accuracy` | Public | None | `species_stats` rows, `meta`, `benchmark` (BV-BRC head-to-head stats), and `edgecases` (`wildtype_tetratype_by_qnr`, `senterica_FDSE`). |
| `GET /api/v1/stats/geography` | Public | Query: `species`, `st`, `tetratype`, `phenotype`, `septatype`. | Country counts (`iso3`, `country`, `n`, `res`) plus `no_country` count. |
| `GET /api/v1/stats/septatypes` | Public | Query: `species` (req), `top` (`1..40`, default `12`). | Top septatypes for the species broken down by `habitat`. |
| `GET /api/v1/aminoacids` | Public | None | 20 amino acids and their physicochemical classes. |
| `GET /api/v1/download/{kind}` | Public | Path: `kind` (`csv` or `json`) + isolate filter params. | Streams `ciprotyper_isolates.csv` (formula-injection sanitized via `_csv_safe`) or `ciprotyper_isolates.json`. |
| `GET /api/v1/bvbrc/{genome_id}` | Public | Path: `genome_id` (`^\d{1,8}\.\d{1,6}$`). | Looks up official `genome_name` from BV-BRC API. |
| `GET /api/v1/bvbrc-status` | Public | None | Background BV-BRC batch annotation status (`loaded`, `running`, `done`, `total`, `error`, `names`, `calls`). |
| `POST /api/v1/submissions` | Public (`20`/hr/IP) | Form: `species`, `habitat` (not `ALL`), `submitter`, `email`, `consent` (`=="2026-11-1"`), `file` or `sequence`, plus optional metadata. | Analyses and deletes genome, inserts derived record into `submissions`, triggers GitHub backup, returns `job_id` and `status`. |
| `GET /api/v1/admin/queue` | `X-Admin-Token` | Query: `status` (`needs_review`, `pending`, `approved`, `rejected`, or `""`). | Lists submissions and queue counts across all 5 tabs. |
| `POST /api/v1/admin/decide` | `X-Admin-Token` | JSON: `{"job_id": "...", "action": "approve"\|"reject", "note": "...", "curator": "..."}`. | Approves, rejects, revokes an approved isolate, or re-approves a rejected isolate; records in `audit_log`; triggers GitHub backup. |
| `GET /api/v1/admin/backup` | `X-Admin-Token` | None | Downloads `ciprotyper_backup.json`. |
| `POST /api/v1/admin/restore` | `X-Admin-Token` | JSON backup payload | Restores missing submissions and replays approved isolates into `isolates` and `septatypes`. |
| `GET /api/v1/admin/persist-status` | `X-Admin-Token` | None | Returns GitHub auto-backup configuration and last saved/restored timestamps. |
| `POST /api/v1/admin/bvbrc/load` | `X-Admin-Token` | None | Starts background BV-BRC batch loader thread. |
| `GET /healthz` | Public | None | Container readiness check (`{"ok": true, "blast": true}`). |

### HTTP Status Codes

* `200 OK`: Request succeeded.
* `400 Bad Request`: Invalid input (protein sequence uploaded, unknown species/source, invalid BV-BRC ID, missing consent or submitter fields).
* `401 Unauthorized`: Wrong `X-Admin-Token` header (compared via constant-time `hmac.compare_digest`).
* `404 Not Found`: Genome ID or submission `job_id` not found.
* `413 Payload Too Large`: Upload exceeds `60 MB` (or `240 MB` decompressed) or pasted sequence exceeds `10 MB`.
* `429 Too Many Requests`: More than 2 concurrent BLAST analyses running, IP rate limit exceeded (`RATE_SUBMIT_PER_HOUR = 20`), or IP locked out by `FailGate` (`10` wrong curator tokens -> `900` s lockout).
* `502 Bad Gateway`: Upstream BV-BRC API unreachable during `/predict/bvbrc`.
* `503 Service Unavailable`: Reference FASTA files missing or `ADMIN_TOKEN` not configured on server.

## 5. Local Installation, Docker Deployment, and Environment Variables

### 5.1 Running Locally or with Docker

```bash
# Option A: Run with Docker (recommended, includes NCBI BLAST+ and builds seed.db automatically)
git clone https://github.com/ikbalsk/Septatyper.git
cd Septatyper
docker build -t ciprotyper .
docker run --rm -p 8000:8000 -e ADMIN_TOKEN="your-secret-token" ciprotyper

# Option B: Run directly with Python 3.12 and system BLAST+
sudo apt-get update && sudo apt-get install -y ncbi-blast+
pip install -r requirements-dev.txt
python scripts/build_db.py
uvicorn app.main:app --host 0.0.0.0 --port 8000
python -m pytest tests
```

### 5.2 Environment Variables (`app/config.py`)

| Variable | Default | Description |
| :--- | :--- | :--- |
| `ADMIN_TOKEN` | `""` (disabled) | Secret token required in `X-Admin-Token` header for curator endpoints. |
| `ADMIN_MAX_FAILS` / `ADMIN_LOCK_S` | `10` / `900` | Locks out an IP for `900` seconds (15 min) after `10` failed curator token attempts. |
| `RATE_SUBMIT_PER_HOUR` | `20` | Max community submissions per IP per hour (`0` = unlimited). |
| `RATE_PREDICT_PER_HOUR` | `0` | Max genome predictions per IP per hour (`0` = unlimited). |
| `CP1_MIN_ID` | `90` (`90.0%`) | Minimum identity across `gyrA`, `parC`, and `recN` required to apply Checkpoint 1. |
| `MLST_SPECIES_MIN_ID` | `90` (`90.0%`) | Minimum identity across all 7 MLST loci (on genomes `>= 1 Mb`) to pass species check. Below `90.0%`, prompts *"Please check the species selected and resubmit"*, assigns no ST, does not apply the 3 checkpoints, and gives no phenotype result. |
| `MAX_UPLOAD_MB` / `MAX_PASTE_MB` | `60` / `10` | Upload size limit (`60` MB) and pasted sequence limit (`10` MB). |
| `BLAST_TIMEOUT_S` | `180` | Subprocess timeout in seconds for `makeblastdb` and `blastn`. |
| `REFRESH_DB_ON_START` | `1` | Refreshes founding rows from `seed.db` on startup while preserving community submissions. |
| `BVBRC_AUTOLOAD` | `1` | Fetches BV-BRC genome names and computational predictions in background at startup. |
| `BACKUP_GITHUB_REPO` | Unset | Private GitHub repository (`owner/repo`) for automatic JSON backups. |
| `BACKUP_GITHUB_TOKEN` | Unset | Fine-grained GitHub token with `Contents: Read and write` on `BACKUP_GITHUB_REPO`. |
| `BACKUP_GITHUB_BRANCH` / `_PATH` | Default / `septatyper_backup.json` | Optional target branch and filename in the backup repository. |
| `DATA_DIR` / `PORT` | `/srv/data` / `8000` | Directory for live `septatype.db` and HTTP listening port. |
