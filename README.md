# SeptaTyper v2

Ciprofloxacin resistance prediction and a curated repository of 2,734 Enterobacterales genomes
(*E. coli*, *S. sonnei*, *S. enterica* subsp. I, *K. pneumoniae*).

| Module | What it does |
|---|---|
| A Predict | Upload or paste FASTA (or enter residues). BLAST finds gyrA, parC, recN and the qnr genes, codons are read, the septatype is built and the phenotype called. Shows the sample's MLST sequence type, the STs of the reference isolates behind the majority-rule call, the RecN-346 codon (cta/ctg flag), combined mean hydrophobicity <H>, and a plain-language reason. |
| B Browse | Filter by species, site, ST, septatype, country, year. Dual-display Genome ID (BV-BRC link) with isolate name beneath. Lab, BV-BRC and SeptaTyper columns with a star where SeptaTyper matches the laboratory phenotype and the BV-BRC computational prediction does not. CSV/JSON download. Septatype-by-site bars and a country map. |
| C Submit | Community form with BV-BRC name lookup. Submissions wait in a curation queue. New septatypes and phenotype/prediction disagreements are flagged. A curator approves or rejects. |
| D Theory | The context-dependent grammar, an interactive schematic of the Mg2+ bridge, validation tables, edge cases. |

REST API under `/api/v1` (interactive docs at `/docs`).

## How a call is made

1. **Database majority rule** (the paper's rule). If the septatype was seen for that species and site,
   the majority phenotype of those isolates is the call. Ties go to Susceptible. Reproduces the paper's accuracies.
2. **Physicochemical grammar** (fallback only). For a septatype never seen for that species, four positions vote
   (hydroxyl group lost, negative charge lost), the QRDR helix `<H>` settles single votes, then the qnr, *S. sonnei* and *S. enterica* ctg overrides apply.
   The result is flagged as new and sent to curation. Rules and residue classes are in `app/grammar.py`.

## Run locally

```bash
pip install -r requirements.txt
sudo apt install ncbi-blast+                 # or: conda install -c bioconda blast
python scripts/build_db.py
uvicorn app.main:app --reload                # http://localhost:8000
python -m pytest tests                       # 13 tests
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
- The result explains each call in three criteria, laid out like the manuscript's residue-class table:
  1 residue class at each position, 2 mean hydrophobicity <H> against the cutoff, 3 septatype overrides.
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
