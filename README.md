# Ciprotyper v3.5

Ciprofloxacin resistance prediction and a curated repository of 2,734 Enterobacterales genomes
(*E. coli*, *S. sonnei*, *S. enterica* subsp. I, *K. pneumoniae*).

| Module | What it does |
|---|---|
| A Predict | Upload or paste FASTA (or enter residues). BLAST finds gyrA, parC, recN and the qnr genes, codons are read, the septatype is built and the phenotype called. Shows the sample's MLST sequence type, the STs of the reference isolates behind the majority-rule call, the RecN-346 codon (cta/ctg flag), combined mean hydrophobicity <H>, and a plain-language reason. |
| B Browse | Filter by species, site, ST, septatype, country, year. Columns: Genome ID (BV-BRC link, isolate name beneath), species and source, country, septatype, Ciprotyper predicted phenotype (star where it matches the laboratory result but the BV-BRC computational call does not), BV-BRC Phenotype Info (Lab, Computationally predicted), submitter, date, affiliation and publication link. 100 rows per page by default, 25 to 500 or a custom number up to 1000. Small text, a header that stays put while the rows scroll inside the table box, an always-visible horizontal scroll bar, and a Hide filters button to give the table room. The CSV follows the table columns. |
| C Submit | Community form with BV-BRC name lookup. Submitter and contact email are required, the date is set automatically, affiliation and publication (a DOI or link) are optional. There is no MIC field. Submissions wait in a curation queue. New septatypes and phenotype/prediction disagreements are flagged. A curator approves or rejects. |
| D Theory | The context-dependent grammar, an interactive schematic of the Mg2+ bridge, validation tables, edge cases. |

REST API under `/api/v1` (interactive docs at `/docs`).

## How a call is made

1. **Database majority rule** (the paper's rule). If the septatype was seen for that species and site,
   the majority phenotype of those isolates is the call. Ties go to Susceptible. Reproduces the paper's accuracies.
2. **Three checkpoints** (fallback only, also shown next to every majority call). Checkpoint 1 follows the syntax rules table: GyrA-83,
   GyrA-87, ParC-80 and ParC-84 are Susceptible if all four residues are in the Susceptible column (S,Y,T | D,N,E | S,Y,T | D,E) and Resistant
   if all four are in the Resistant column (F,L,I,W | N,D,G,Y,A,E | I,F,W | E,V,D, with the pair S at ParC-80 plus K at ParC-84 as the one exception, as in Fig. 6b of the paper). The residues follow the paper's Table S6 and Fig. 6. The table
   is `SUSCEPTIBLE_SET`, `RESISTANT_SET` and `RESISTANT_EXCEPTION` (a pair) in `app/grammar.py`. Checkpoint 2: RecN-346 is Susceptible if L encoded by cta,
   Resistant if L encoded by ctg or F. Checkpoint 3: Resistant if qnrB2 or qnrS1 is present. The isolate is Susceptible if all checkpoints are
   Susceptible and Resistant if any is. Where the four residues are not all in one column the written rules are silent, so the site counts
   positions outside the Susceptible column (two or more: Resistant, none: Susceptible) and settles one position with the QRDR helix `<H>`.
   *S. sonnei* keeps one extra rule from the data (a single GyrA-83 or GyrA-87 change is Resistant). The result is flagged as new and sent to curation.

## Run locally

```bash
pip install -r requirements-dev.txt         # the runtime packages plus pytest; the production image installs requirements.txt only
sudo apt install ncbi-blast+                 # or: conda install -c bioconda blast
python scripts/build_db.py
uvicorn app.main:app --reload                # http://localhost:8000
python -m pytest tests                       # the full test suite
```

## Deploy on Render (free)

1. Push the files to a **private** GitHub repository.
2. Render, New, Web Service, pick the repo, Environment: Docker, Plan: Free.
3. Environment variables: `ADMIN_TOKEN` = a long random string (turns on the curator page at `#/curate`).
4. For permanent submissions, attach a disk mounted at `/srv/data` (paid plan) or host on an institute server:
   `docker run -p 8000:8000 -v septa-data:/srv/data -e ADMIN_TOKEN=... septatyper`.
   Without a disk, submissions and approved isolates are lost whenever the free instance restarts.

## Fill the BV-BRC columns (once, needs internet)

```bash
python scripts/fetch_bvbrc.py --dry-run      # look at what BV-BRC returns first
python scripts/fetch_bvbrc.py
```
This fills official genome names and BV-BRC computational ciprofloxacin predictions, which switches on the
BV-BRC column, the gold star and the head-to-head table on the Theory page. On a hosted site, run it
against the live database (`--db /srv/data/septatype.db`) from the server shell or a one-off job.

## Honest limits

- **Not run on real genomes.** Tests use real reference genes in synthetic genomes (both strands, indels, qnr).
  Run real BV-BRC genomes through `/api/v1/predict` and compare with Table S1 before release.
- **`fetch_bvbrc.py` is untested** (BV-BRC was unreachable from the build environment). Use `--dry-run` first.
- **Accuracy is resubstitution accuracy.** Rules and scores use the same 2,734 isolates. The grammar was written
  while looking at the same data. Only the 23 in-house *Shigella* isolates in the paper are independent.
- **Assumptions to confirm:** `<H>` is the mean Eisenberg hydrophobicity of the whole QRDR helices: GyrA 82-93 `D(S)AVY(D)TIVRMA` (sum over 12) and ParC 80-90 `(S)ACY(E)AMVLMA` (sum over 11), with the sample's residues at GyrA-83/87 and ParC-80/84, and `<H>` = (GyrA mean + ParC mean) / 2. G, P, C count as hydrophobic. The cutoff is 0.334, the midpoint between the wild-type helices S-D-S-E (0.270) and the resistant reference L-N-I-E (0.397). That rule was chosen after seeing that the wild-type value alone (0.270) turns *S. enterica* F-D-S-E resistant whatever the codon, so confirm the cutoff against the paper. References are one gene record per species (refs/gyrA.fasta, parC.fasta, recN.fasta) plus single qnrB2 and qnrS1 alleles, so qnr family members below 80% identity to those alleles would be missed.
- **Not built:** raw-read input, Mol* 3D structures (a labelled schematic is used instead), user accounts, Celery job queue.
- The Theory schematic is a reading of the grammar, not a structure. Align its wording with Figure 6.
- Fonts and the map library load from CDNs. The site works without them (fallback fonts, table instead of map).
- Declare AI-assisted code in the manuscript as required by the journal, and add a data licence.

## Sequence types (MLST)

`app/mlst.py` types the uploaded sequence with PubMLST schemes: Achtman for *E. coli*, *S. sonnei* and *S. enterica*,
Diancourt for *K. pneumoniae*. Alleles and profiles are in `refs/mlst_<scheme>.fasta` and `.txt` (taken from the
tseemann/mlst database, downloaded 2026-03-11; cite PubMLST, Jolley et al. 2018, Wellcome Open Res 3:124).
An ST is reported only when all seven alleles match exactly and the profile exists. Otherwise the nearest alleles
are shown. Gene-only uploads cannot be typed. The reference STs shown for a majority-rule call are the `st`
values of Table S1 for the isolates that carry the same septatype. Refresh the scheme files from the tseemann/mlst
repository when PubMLST adds new alleles.

## Speed and explanation panel (v2.2)

- A 5 Mb genome is analysed in under a second on one CPU (was about 8 s). MLST now locates each locus with four
  spread-out alleles and looks the exact window up in a dictionary, and only falls back to BLASTing a locus's
  full allele set when it is not an exact hit. Browsers compress large uploads before sending, and the server
  keeps the last 32 results so the same file is never analysed twice.
- **Predict from a BV-BRC genome ID.** The Predict page has a "BV-BRC ID" tab. Type one ID (for example 562.7710) or up to 10 separated by
  spaces, commas or new lines. For each ID the server asks BV-BRC for the genome record (`/api/v1/predict/bvbrc`, one genome per call, the page
  loops through the list), reads the species from the record, downloads the assembly (FASTA output of the `genome_sequence` collection, with a JSON
  fallback that rebuilds the FASTA), analyses it and discards it. One ID opens the full result, several give the batch table with CSV download.
  A genome BV-BRC names as another organism is refused, never guessed. A genome that is one of the 2,734 reference isolates is flagged, because
  its call is not an independent test. Only IDs of the form `digits.digits` are ever put in a URL, and only on the fixed BV-BRC host. The tests
  fake the two BV-BRC calls: the real download could not be tried from the build environment, so test it once from the live site with a known ID.
- The result explains each call in three criteria, laid out like the manuscript's residue-class table:
  checkpoint 1 as a two-column residue table (Susceptible, Resistant) in which the sample's residue is lit and blinks, with the QRDR helix <H> worked out when the table leaves it open, checkpoint 2 RecN-346 and checkpoint 3 qnr, then the combined call.
  The criterion that decides the sample is highlighted and pulses.
- `#/calc` shows the worked calculation of the QRDR helix `<H>` for any four residues (every helix residue, the GyrA mean over 12, the ParC mean over 11, the average of the two), where the 0.334 cutoff comes from, and a data check of it against the reference isolates (`/api/v1/stats/cutoff`, which also scores the whole grammar at each cutoff). `<H>` alone is weak when exactly one position is in the Resistant column, so it is only a last tie-breaker. Browse and the CSV/JSON downloads give `h_gyr`, `h_parc` and `h_qrdr` for every isolate.

## Curation lists

The curator page (`#/curate`) has Needs review, Pending, Approved, Rejected and All lists with counts. An approved
submission can be rejected later: its isolate leaves the repository and its vote leaves the rules table. A rejected
submission can be approved again. Contact for suggestions and problems: ikbalsk.personal@gmail.com.

## Keeping submissions on a free plan

Free hosting forgets its files when the server restarts. Two free ways to keep community submissions and curation decisions:

1. **Manual:** on the curation page, tap **Download backup** now and then, and **Restore from backup file** after a restart.
2. **Automatic, through a private GitHub repository:**
   1. Create a new **private** repository, for example `septatyper-data`.
   2. In GitHub, Settings, Developer settings, Fine-grained tokens: create a token for that one repository with
      **Contents: Read and write**.
   3. In Render, open the service, Environment, and add `BACKUP_GITHUB_REPO` = `your-name/septatyper-data` and
      `BACKUP_GITHUB_TOKEN` = the token. Optional: `BACKUP_GITHUB_BRANCH` and `BACKUP_GITHUB_PATH`.
   4. The server saves a few seconds after every submission or decision and restores at start-up when it has no
      submissions. The curation page shows the backup status.

The backup holds the derived genotype, metadata, decisions and any e-mail addresses people gave, never the genome, so
keep the repository private. Approved submissions are replayed on restore. A backup problem never blocks the site.

## Isolation sources and countries

The isolation source is a drop-down in Predict, Browse and Submit: stool (fecal), urine, respiratory, blood, wound or
skin, cerebrospinal fluid, abdominal, genital tract, environment, food, animal and other (plus "unknown" in Predict).
Only the first three have reference isolates. For any other source the call uses the isolates of that species from all
sources and says so, until curators approve submissions from that source. Submitters pick their country from a
searchable list of all countries (`app/countries.py`, ISO 3166 with alpha-3 codes), so approved submissions from any
country appear on the map. Approved submissions count in the analytics; pending and rejected ones do not.


## Name and field names

The tool is called Ciprotyper. The word *septatype* (the typing string such as S-D-S-E-L-0-0) is unchanged. Some internal names still carry the old
spelling so that live databases and backups keep working: the database column and JSON field `septatyper_call`, the accuracy fields
`septatyper_*` and the default backup file name `septatyper_backup.json` (set `BACKUP_GITHUB_PATH` to change it).

## Submitter, date, affiliation and publication

Every isolate has these four fields. The 2,734 founding isolates carry the values in `app/config.py` (`FOUNDING_*`): Ghosh et al, 16 October 2026,
JISIASR, Kolkata, https://doi.org/10.1038/s41598-026-75690-y. A community submission stores what the submitter typed, with the date set to the
day of submission. A database made by an older version is upgraded on start-up (columns added, founding isolates filled in).


## Pasted sequence size

A pasted sequence can be up to 10 MB (`MAX_PASTE_MB`), enough for a whole genome; above that the page and the API say so and ask for a file upload
(up to `MAX_UPLOAD_MB`, 60 MB). Starlette's form parsers stop at 1 MB for ordinary form fields, which is where pasted text arrives, so
`app/main.py` raises that limit (`_raise_form_field_limit`) for both the multipart parser a browser uses and the urlencoded one a script may use.

## Credit

Ciprotyper and its database are based on:

Ghosh, A., Deb, P., Bose, P. *et al.* Generic context-dependence governs clonal emergence of ciprofloxacin resistance mutations. *Sci Rep* (2026).
https://doi.org/10.1038/s41598-026-75690-y

The citation is in the page footer, on the home page, in a "Credit and how to cite" box on the About page (with a copy button), on the Theory page,
at the bottom of every result and in the printed batch report. Please cite it if you use Ciprotyper or its data.


## Safety settings and what changed in v3.5

All are environment variables on the server; the defaults are meant to leave normal use untouched.

| Setting | Default | What it does |
|---|---|---|
| `RATE_SUBMIT_PER_HOUR` | 20 | Submissions allowed per visitor per hour (`0` = no limit). Attempts count even if the form was incomplete. |
| `RATE_PREDICT_PER_HOUR` | 0 (off) | Genome predictions, including BV-BRC ID lookups, per visitor per hour. Off by default; set it if the site is abused. |
| `ADMIN_MAX_FAILS`, `ADMIN_LOCK_S` | 10, 900 | After this many wrong curator tokens from one address, that address waits 15 minutes. |
| `REFRESH_DB_ON_START` | 1 | When the server is updated with a new seed, bring the founding data of an already saved database up to date (see below). `0` turns it off. |
| `MAX_PASTE_MB` | 10 | Largest pasted sequence. |

A visitor is identified by the first address in `X-Forwarded-For`, else the connecting address. The limits are in memory and reset when the server restarts. They are spam control, not a security boundary.

- **Saved databases.** On a server with a persistent disk, the saved database used to keep the founding data it was first copied with. Now, when the seed changes (a new `Table_S1.xlsx`, new rules or scores), the first start after the update saves a copy of the old database in `backups/` next to it (the last 3 are kept), replaces the founding isolates and the rule tables, and keeps everything the community added, including its votes. If anything goes wrong nothing changes and the site starts as before. On a server without a persistent disk the database is copied fresh from the seed each time, so nothing is needed.
- **Who decided.** Each approval or rejection records the time and the address, plus the name the curator typed in the optional "Your name" box (`audit_log` table, and shown in the queue). There is still one shared token, so the name is what the curator chooses to enter.
- **Spreadsheets.** In the CSV downloads, text that starts with `=`, `+`, `-` or `@` gets a leading apostrophe so Excel cannot run it as a formula.
- **Map library.** Plotly 2.35.2 is now served with the site (`static/plotly.min.js`, MIT licence). If that file is missing, the page falls back to the jsDelivr copy only with a fixed integrity hash, so a changed file is refused. Google Fonts still load from Google (a stylesheet cannot be given an integrity hash; it cannot run code).
- **Backups.** If the GitHub backup is not set, the Curation page now shows a clear warning that submissions and approvals are lost when the server restarts (unless it has a persistent disk).
- **Pinned packages.** `requirements.txt` holds the exact versions the tests ran against; `pytest` and `httpx` moved to `requirements-dev.txt`.
- **Scripts.** `scripts/validate_genomes.py` uses the right column (`lab_phenotype`). `scripts/fetch_references.py` now writes the layout the app reads (`refs/gyrA.fasta` with one record per species, `qnrB2.fasta`, `qnrS1.fasta`), one reference record per gene as the shipped files were cut, and never touches `refs/` unless run with `--install`.
