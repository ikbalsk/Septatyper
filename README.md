# Ciprotyper

A prediction tool and a curated repository for ciprofloxacin resistance in Enterobacterales. It predicts resistance from five codons and two plasmid-borne genes, and it holds a reference repository of 2,734 genomes.

---

## 1. Overview

Ciprotyper is three things in one service: a web server, a REST API and a curated repository. As a web server and API it predicts the ciprofloxacin phenotype (Resistant or Susceptible) of a bacterial isolate from assembled nucleotide sequence. It supports *Escherichia coli*, *Shigella sonnei*, *Salmonella enterica* subsp. I and *Klebsiella pneumoniae*.

The method is based on the septatype, a seven-field genotype described in Ghosh et al., Scientific Reports, doi:10.1038/s41598-026-75690-y. The first four fields are the residues at GyrA-83, GyrA-87, ParC-80 and ParC-84. The fifth is the residue at RecN-346. The last two record the presence (1) or absence (0) of qnrB2 and qnrS1. An example is `L-N-I-E-L-0-0`.

Every call is explained as three checkpoints. The server also reports the MLST sequence type of the sample.

**Repository.** The service keeps a curated repository of reference isolates. Each record holds the genome identifier (linked to BV-BRC), species, isolation source, country, year, sequence type, the five residues and qnr status, the septatype, the laboratory phenotype, the Ciprotyper call, the BV-BRC computational call, and the submitter, date, affiliation and publication. Users can filter and browse the records, view them by country and septatype, and download them as CSV or JSON. Anyone can submit a new isolate. A curator checks each submission, and only approved isolates join the repository and update the prediction rules. Only derived genotype data are stored, never the genome.

**Intended use:** research only. The tool is not a diagnostic device and has not been validated for clinical use.

| Item | Value |
|---|---|
| Reference isolates | 2,734 (E. coli 400, S. sonnei 334, S. enterica 510, K. pneumoniae 1,490) |
| Resistant isolates | 1,682 (E. coli 193, S. sonnei 151, S. enterica 6, K. pneumoniae 1,319) |
| Isolation sources | Stool 1,105; urine 1,231; respiratory 398 |
| Distinct septatypes | 64 |
| Countries | 28 |
| Accuracy of the majority rule on the reference set | 98.1 percent (resubstitution, see section 7) |
| Implementation | Python 3.12, FastAPI, SQLite, NCBI BLAST+, Biopython; plain JavaScript front end |

---

## 2. Methods

### 2.1 Sequence analysis

1. **Input.** Nucleotide FASTA, plain or gzip-compressed (`.fa .fasta .fna .fas .txt .gz`), up to 60 MB by default. Other input modes are pasted text, a batch of up to ten files, a BV-BRC genome identifier, and manual entry of the five residues.
2. **Validation.** The file must parse as FASTA. If fewer than 90 percent of bases are A, C, G, T or N, the input is rejected as protein.
3. **Gene location.** `makeblastdb` and `blastn` (e-value 1e-20) locate gyrA, parC and recN using one species-specific reference record per gene. A hit is accepted at 75 percent identity or higher. A warning is issued below 92 percent identity or below 90 percent coverage.
4. **Codon reading.** The matched region is extracted (reverse-complemented if on the minus strand) and globally aligned to the reference with the Biopython pairwise aligner (match 2, mismatch -3, gap open -8, gap extend -1). The codons at the five positions are read through the alignment. A position with a gap, a stop or an ambiguous base is not called and produces a warning. Reading codons rather than residues preserves synonymous differences, which matter at RecN-346.
5. **qnr detection.** qnrB2 and qnrS1 reference alleles are searched with blastn. A gene is scored present at 80 percent identity or higher and 90 percent coverage or higher. For gene-only input the user may set presence manually, because absence cannot be shown from a partial file.
6. **MLST.** Schemes from PubMLST: Achtman for E. coli, S. sonnei and S. enterica, and Diancourt for K. pneumoniae. Four spread-out alleles per locus locate the locus. The matching window is cut out and looked up exactly in a dictionary of known alleles. Only a locus without an exact hit is searched again against its full allele set. A sequence type is reported only if all seven alleles match exactly and the profile exists in the scheme. A failure in typing never prevents the prediction.
7. **Septatype assembly.** The seven values are joined with hyphens. If any of the five residues cannot be read, no septatype and no call are produced.

### 2.2 Phenotype call

**Tier 1: majority rule.** If the septatype has been observed for the same species and isolation source, the most frequent laboratory phenotype among those reference isolates is the call. Ties are resolved as Susceptible. If the septatype has not been observed for that source but has for the species, isolates from all sources are used and the result is flagged.

**Tier 2: grammar.** Used for septatypes not seen before. Three checkpoints are evaluated, and they are also displayed beside every Tier 1 call.

| Checkpoint | Rule |
|---|---|
| 1. Residue table | Susceptible if GyrA-83, GyrA-87, ParC-80 and ParC-84 all lie in the Susceptible column (S,Y,T / D,N,E / S,Y,T / D,E). Resistant if all four lie in the Resistant column (F,L,I,W / N,D,G,Y,A,E / I,F,W / E,V,D). One exception pair, S at ParC-80 with K at ParC-84, is Resistant. |
| 2. RecN-346 | Susceptible if leucine encoded by `cta`. Resistant if leucine encoded by `ctg`, or phenylalanine. |
| 3. qnr | Resistant if qnrB2 or qnrS1 is present. |

The combined call is Susceptible only if all three checkpoints are Susceptible, and Resistant if any one is Resistant.

**Rules added by the implementation.** Checkpoint 1 does not cover samples whose four residues do not fall in a single column (for example L-D-S-E). In that case the server counts the positions outside the Susceptible column. Two or more gives Resistant and none gives Susceptible. Exactly one is decided by the mean hydrophobicity of the QRDR helices (H) against a cutoff of 0.334. In addition, S. sonnei carries one rule derived from the data: a single change at GyrA-83 or GyrA-87 is Resistant.

H is calculated on the Eisenberg normalised consensus scale. The GyrA helix covers residues 82 to 93 (12 residues) and the ParC helix covers residues 80 to 90 (11 residues), with the sample's residues placed at positions 83, 87, 80 and 84. H is the mean of the GyrA helix mean and the ParC helix mean. The wild-type helices (S-D-S-E) give 0.270, the resistant reference (L-N-I-E) gives 0.397, and the cutoff is their midpoint, 0.334. The cutoff was fixed from these two reference helices and was not tuned on the isolates.

### 2.3 Submission and curation

1. A user submits a genome with species, isolation source, submitter name and contact email (all required). Optional fields are a BV-BRC identifier, country, year, laboratory phenotype, affiliation and publication (DOI or link).
2. The genome is analysed and deleted. Only the derived genotype, the entered metadata and the result are stored. The submission date is set by the server.
3. A submission is marked "needs review" if its septatype is new, if its laboratory phenotype disagrees with the prediction, or if no call was possible. Otherwise it is "pending".
4. A curator, authenticated by a shared token in the `X-Admin-Token` header, approves or rejects. Approval requires a laboratory phenotype and a complete septatype. It adds the isolate to the repository and updates the rule table, so a new septatype becomes known. An approved submission can later be rejected, which removes the isolate and its vote. A rejected submission can be approved again.
5. Only approved submissions appear in Browse, the analytics and the map.

### 2.4 BV-BRC integration

- **Prediction by BV-BRC identifier.** The server retrieves the genome record and assembly from the BV-BRC API, determines the species from the record, analyses the assembly and discards it. A genome that BV-BRC assigns to another organism is refused and is not guessed. A genome that is one of the reference isolates is flagged, because its result is not an independent test. Only identifiers of the form `digits.digits` are accepted, and only the fixed BV-BRC host is contacted.
- **Bulk annotation.** At start-up, unless `BVBRC_AUTOLOAD=0`, a background task retrieves official genome names and BV-BRC computational ciprofloxacin predictions for the reference isolates. This enables the BV-BRC column in Browse, the marker for isolates where Ciprotyper agrees with the laboratory result and BV-BRC does not, and the head-to-head table. The task fails softly.

### 2.5 Persistence

Free hosting discards local files on restart. Two mechanisms preserve submissions and curation decisions:

- Manual download and restore of a JSON backup from the curation page.
- Automatic backup to a private GitHub repository (`BACKUP_GITHUB_REPO`, `BACKUP_GITHUB_TOKEN`). The file is written about 8 seconds after each change and restored at start-up when no submissions are present.

The backup contains derived genotypes, metadata, curator decisions and any contact emails, and never the genome. The repository must therefore be private.

---

## 3. Architecture

### 3.1 Components

| File | Role |
|---|---|
| `main.py` | HTTP endpoints, upload handling, result cache, concurrency limit, admin authentication, static page |
| `septatype.py` | Pipeline: FASTA parsing, BLAST, alignment, codon reading, qnr detection, septatype |
| `grammar.py` | Residue tables, checkpoints, hydrophobicity, explanations |
| `rules.py` | Database access, majority rule, statistics, browse queries, curation |
| `mlst.py` | Sequence typing |
| `bvbrc.py` | BV-BRC client |
| `persist.py` | Backup and restore |
| `config.py` | Paths, environment variables, species and source lists, founding-isolate metadata |
| `countries.py` | ISO 3166 country names and alpha-3 codes |
| `build_db.py` | Builds the seed database from `Table_S1.xlsx` |
| `fetch_bvbrc.py` | Command-line version of the BV-BRC annotation |
| `index.html`, `app.js`, `style.css` | Complete front end (ten pages, hash routing, no framework) |

The Dockerfile places the application in `app/`, the front end in `app/static/`, the build scripts in `scripts/`, the table in `data/` and the reference sequences in `refs/`. The seed database is built at image build time.

### 3.2 Database

SQLite with the tables `isolates`, `septatypes` (the rule table, one row per species, source and septatype), `species_stats`, `aa_properties`, `submissions` and `meta`. A read-only seed database is built into the image. A live copy is created in `DATA_DIR` on first start and holds submissions.

### 3.3 Reference sequences

One reference record per species for each of gyrA, parC and recN (header format `gene|Species|accession`), and single qnrB2 and qnrS1 alleles (CP032523.1). Chromosomal accessions are CP058222.1 (E. coli), AAVEOC010000165.1, AAVEOC010000011.1 and AAVEOC010000127.1 (S. sonnei gyrA, parC and recN), CP163542.1 (S. enterica) and CP089383.1 (K. pneumoniae). The residues at the five positions in all records are S, D, S, E and L. MLST alleles and profiles were taken from the tseemann/mlst repository (download dated 2026-03-11).

### 3.4 Web interface

Home, Predict, Browse, Analytics, Submit, Theory, Calc (worked hydrophobicity calculation), API and data, About, and Curate (token protected).

---

## 4. REST API

A REST API is a set of web addresses that a program can call to get or send data over HTTP, without using the web page. Each address (endpoint) does one job and returns structured data (JSON). It lets users run predictions or download the repository from scripts. Base path `/api/v1`. Interactive documentation is served at `/docs`.

| Endpoint | Purpose |
|---|---|
| `GET /meta` | Species, sources, counts, reference status |
| `POST /predict` | Form fields `species`, `habitat`, `file` or `sequence`, optional `qnrB2` and `qnrS1` (`auto`, `present`, `absent`) |
| `POST /predict/bvbrc` | JSON with `genome_id`; one genome per call |
| `POST /predict/manual` | JSON with five residues, qnr flags and optional `recN_codon` |
| `GET /isolates`, `GET /isolates/{id}` | Filtered and paged isolates (up to 1000 per page) with `h_gyr`, `h_parc`, `h_qrdr` |
| `GET /rules` | Rule table for a species and source |
| `GET /grammar` | Checkpoints for supplied residues |
| `GET /stats/accuracy`, `/stats/geography`, `/stats/septatypes`, `/stats/cutoff` | Summary statistics |
| `GET /aminoacids` | Hydrophobicity scale |
| `GET /download/csv`, `/download/json` | Export, optionally filtered |
| `GET /bvbrc/{id}`, `GET /bvbrc-status` | Name lookup and loader status |
| `POST /submissions` | Community submission |
| `GET /admin/queue`, `POST /admin/decide`, `GET /admin/backup`, `POST /admin/restore`, `GET /admin/persist-status`, `POST /admin/bvbrc/load` | Curator functions (token required) |
| `GET /healthz` | Health check |

Status codes: 400 invalid input; 401 wrong token; 413 file too large; 429 server busy; 502 BV-BRC unreachable; 503 reference files missing or curation disabled.

---

## 5. Installation and operation

### 5.1 Local installation

Requirements: Python 3.12, NCBI BLAST+, and the packages in `requirements.txt` (fastapi, uvicorn, python-multipart, biopython, openpyxl, httpx, pytest).

```
pip install -r requirements.txt
sudo apt install ncbi-blast+
python scripts/build_db.py
uvicorn app.main:app --reload
python -m pytest tests
```

The application is served at http://localhost:8000.

### 5.2 Deployment

The application is distributed as a Docker image and runs on any container host. On Render:

1. Place the project files in a private GitHub repository.
2. Create a Web Service from the repository with Environment set to Docker and health check path `/healthz`.
3. Set `ADMIN_TOKEN`. Optionally set the backup variables.
4. Attach a persistent disk at `/srv/data` on a paid plan. On the free plan the service sleeps when idle and may take about a minute to resume, and local files are lost on restart.

The same image runs as a Hugging Face Docker Space with `PORT=7860` and `DATA_DIR=/tmp/septatype-data`. Storage there is also ephemeral, so the GitHub backup is advised.

### 5.3 Configuration

| Variable | Default | Meaning |
|---|---|---|
| `ADMIN_TOKEN` | empty | Curator token. An empty value disables curation. |
| `DATA_DIR` | `/srv/data` in Docker | Location of the live database |
| `DB_PATH`, `SEED_DB`, `REFS_DIR` | derived | Path overrides |
| `PORT` | 8000 | Listening port |
| `MAX_UPLOAD_MB` | 60 | Upload limit |
| `BLAST_TIMEOUT_S` | 180 | Time limit per BLAST call |
| `BVBRC_AUTOLOAD` | 1 | Set to 0 to disable the start-up BV-BRC task |
| `BACKUP_GITHUB_REPO`, `_TOKEN`, `_BRANCH`, `_PATH` | none | Automatic backup. The token needs Contents read and write access to that one repository. |

Parameters set in code: identity and coverage thresholds (`septatype.py`); `LOW_SUPPORT = 5` and `MIXED = 0.80` (`rules.py`); `H_CUTOFF`, `SUSCEPTIBLE_SET`, `RESISTANT_SET` and `RESISTANT_EXCEPTION` (`grammar.py`).

### 5.4 Updating a running service

Push the changed files to the repository and let the host rebuild. Confirm that every file named in the Dockerfile is present.

The live database is copied from the seed only when it does not exist. On a persistent disk, later changes to `build_db.py` or `Table_S1.xlsx` do not reach an existing live database, because the start-up migration only adds missing columns. Back up, remove the live database and let it rebuild, or add a versioned migration.

### 5.5 Refreshing reference data

MLST files can be replaced from the tseemann/mlst repository when PubMLST adds alleles. Reference genes follow the header convention in section 3.3. After any change to Table S1, run `build_db.py` again.

---

## 6. Testing

The suite (`tests/test_pipeline.py`, `tests/conftest.py`) contains 54 test functions. They cover reproduction of the published accuracy and counts, grammar checkpoints and boundary cases, wild-type and resistant synthetic genomes (including reverse-strand genes and indels), qnr detection, gzip equivalence, MLST with real alleles, submission, approval, rejection and re-approval, backup and restore with a simulated GitHub endpoint, BV-BRC with a simulated API, per-species references, hydrophobicity arithmetic and upgrade of an older database.

The tests embed real reference genes in synthetic genomes. They do not use real assemblies, and the BV-BRC service is simulated. They show that the software behaves as designed. They do not measure predictive performance.

---

## 7. Validation status and limitations

1. **No independent validation.** The reported 98.1 percent accuracy is resubstitution accuracy. The rule table and the score use the same 2,734 isolates, and the grammar was written while examining them. Only the 23 in-house Shigella isolates described in the source article are independent. Performance on unseen genomes is not yet known.
2. **Real assemblies.** The pipeline has not yet been run on real assemblies. Genomes from BV-BRC should be passed through `/api/v1/predict` and compared with Table S1.
3. **BV-BRC comparison.** The BV-BRC integration has not been tested against the live service, and no claim is made about the head-to-head comparison.
4. **Assumed parameters.** The cutoff of 0.334, the classification of glycine, proline and cysteine as hydrophobic, and the treatment of an acidic ParC-84 as resistant only next to a hydrophobic ParC-80 are assumptions of the implementation. They are not part of the published checkpoints.
5. **qnr coverage.** Detection uses one qnrB2 and one qnrS1 allele. Family members below 80 percent identity to these alleles are missed. Presence of a gene does not demonstrate expression, and some qnr-positive isolates in the data are susceptible.
6. **RecN-346 codon.** Table S1 contains no codons, so the ctg rule of checkpoint 2 could not be scored on the reference isolates.
7. **Scope.** Four species and one drug. The user selects the species, and E. coli and S. sonnei cannot be separated with these genes. The classes are unbalanced: S. enterica has 6 resistant isolates among 510.
8. **Input.** Assembled nucleotide sequence only. Raw reads are not supported.
9. **Provenance.** The laboratory breakpoint standard behind the phenotypes, and the MLST scheme that produced the sequence types in Table S1, are not documented. Achtman and Diancourt were assumed.

---

## 8. Known issues and maintenance notes

| # | Severity | Issue | Recommended action |
|---|---|---|---|
| 1 | Medium | CSV export writes user-entered fields (submitter, affiliation, genome name, publication) without neutralising a leading `=`, `+`, `-` or `@`. This allows formula injection when the file is opened in a spreadsheet. | Prefix such cells with an apostrophe on export. |
| 2 | Medium | A live database on a persistent disk is not refreshed from a new seed (section 5.4). | Add a versioned migration. |
| 3 | Medium | `scripts/validate_genomes.py` selects a column `phenotype`, but the column is `lab_phenotype`. The script will fail. | Correct and test it. |
| 4 | Medium | `scripts/fetch_references.py` writes `refs/<species>/<gene>.fasta` and qnr family files, which the application does not read. | Update or remove. |
| 5 | Low | The project README states 13 tests, and the suite has 54. | Correct the text. |
| 6 | Low | Version labels are inconsistent (README v2, FastAPI 2.0, User-Agent 2.0 and 2.1). | Use a single version string. |
| 7 | Low | `@app.on_event("startup")` is deprecated in current FastAPI. | Move to lifespan handlers. |
| 8 | Low | The in-memory result cache is modified by several threads without a lock. | Add a lock. |
| 9 | Low | `/admin/bvbrc/load` returns an unclear `started` value. | Return the actual result. |
| 10 | Low | Dependencies are unpinned, and pytest is installed in the production image. | Pin versions and separate development requirements. |

**Security controls in place:** parameterised SQL with allow-listed column and sort names; strict identifier validation and a fixed host for BV-BRC requests; HTML escaping of all server text in the page; constant-time token comparison; upload size and decompression limits; deletion of uploaded genomes; a limit of two concurrent analyses; subprocess calls without a shell; no use of browser storage.

**Security gaps:** no rate limiting or spam control on `/submissions` and `/predict`; a single shared curator token without lockout or an audit trail; Google Fonts and the Plotly library loaded from third-party CDNs without integrity checks; submissions lost on restart when the backup variables are not set.

---

## 9. Licensing, data governance and compliance

No licence file is included yet. Without one, the code remains under default copyright and cannot be legally reused. A licence must be chosen before release.

### 9.1 Code

- Add a `LICENSE` file. MIT or Apache-2.0 permit broad reuse. GPL-3.0 requires derivatives to remain open.
- Confirm the copyright holder under the institute's policy on software written by staff or students.
- No third-party code is copied into the repository.

### 9.2 Dependencies

| Component | Licence |
|---|---|
| FastAPI, openpyxl, pytest, Plotly.js | MIT |
| Uvicorn, httpx | BSD-3-Clause |
| python-multipart | Apache-2.0 |
| Biopython | Biopython License (permissive) |
| NCBI BLAST+ | US government work, public domain |
| JetBrains Mono, Public Sans (loaded from Google Fonts) | SIL Open Font License |
| `python:3.12-slim` base image | PSF and Debian package licences |

### 9.3 Data

| Data | Source | Requirement |
|---|---|---|
| Table S1 (2,734 isolates) | Ghosh et al., Scientific Reports | Follow the licence of the source article and attribute it. Founding records carry the submitter, date, affiliation and DOI. |
| MLST alleles and profiles | PubMLST (Jolley et al. 2018) via tseemann/mlst | Follow the PubMLST data terms and cite PubMLST. Only the data files are used, not the tseemann/mlst code. |
| Reference genes and qnr alleles | GenBank accessions listed in section 3.3 | Open sequence data. Cite the accessions. |
| BV-BRC records and genomes | BV-BRC API | Follow the BV-BRC terms of use and request limits. Genomes are processed and discarded, not redistributed. |
| Hydrophobicity scale | Eisenberg et al. 1984 | Cite. |

A licence for the published repository data (CSV and JSON) must be chosen, for example CC BY 4.0, and stated on the site.

### 9.4 Personal data and submissions

- The service stores submitter names and contact emails. A privacy notice stating purpose, retention and a contact for deletion should be published. Depending on user location, the GDPR (EU and UK) or the DPDP Act (India) may apply.
- The Submit page should state the licence under which submitted data are published, and should ask submitters to confirm that they may share the data and that it contains no patient-identifying information.
- The backup repository contains emails and must remain private.
- The statement that genomes are not stored must remain true.

### 9.5 Regulatory and publication points

- The tool must keep its "research use only" label. Clinical use would fall under in-vitro-diagnostic regulation in most jurisdictions.
- Where a journal requires it, the use of AI-assisted code must be declared in the manuscript.
- Required citations: Ghosh et al.; BLAST+; Biopython; PubMLST; Eisenberg et al. 1984; BV-BRC.
- Code identifiers and the backup file name still use the earlier name "septatyper". Confirm that the product name is free of conflicts before wider promotion.

---

## 10. Release checklist

1. Choose and add a code licence and a data licence.
2. Publish the privacy notice and add submission consent text.
3. Validate on independent, real genomes with laboratory phenotypes, and report the result. Repair `validate_genomes.py` first.
4. Test BV-BRC identifier mode and the bulk annotation with a known identifier (for example 562.7710).
5. Set `ADMIN_TOKEN` and the private GitHub backup, and confirm the last-saved time on the curation page.
6. Correct the CSV export. Consider rate limiting and self-hosting or integrity-checking the CDN assets.
7. Confirm the cutoff, residue classes and exception pair against the source article.
8. Pin dependency versions. Correct the test count and the version label.
9. Add all citations.

---

## 11. Glossary

**Septatype:** seven-field genotype string. **Tetratype:** its first four residues. **QRDR:** quinolone resistance-determining region. **H:** mean hydrophobicity of the QRDR helices. **ST:** MLST sequence type. **Founding isolates:** the 2,734 isolates of Table S1. **Resubstitution accuracy:** accuracy measured on the data used to build the rules. **Curation queue:** submissions awaiting a curator decision.

---

## References

1. Ghosh et al. Scientific Reports. doi:10.1038/s41598-026-75690-y.
2. Camacho C, et al. BLAST+: architecture and applications. BMC Bioinformatics 2009;10:421.
3. Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009;25:1422.
4. Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res 2018;3:124.
5. Eisenberg D, Schwarz E, Komaromy M, Wall R. Analysis of membrane and surface protein sequences with the hydrophobic moment plot. J Mol Biol 1984;179:125.
6. Olson RD, et al. Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC). Nucleic Acids Res 2023;51:D678.
