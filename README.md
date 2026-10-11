# SeptaTyper v2

Ciprofloxacin resistance prediction and a curated repository of 2,734 Enterobacterales genomes
(*E. coli*, *S. sonnei*, *S. enterica* subsp. I, *K. pneumoniae*).

| Module | What it does |
|---|---|
| A Predict | Upload or paste FASTA (or enter residues). BLAST finds gyrA, parC, recN and the qnr genes, codons are read, the septatype is built and the phenotype called. Shows RecN-346 codon (cta/ctg flag), combined mean hydrophobicity <H>, and a plain-language reason. |
| B Browse | Filter by species, site, ST, septatype, country, year. Dual-display Genome ID (BV-BRC link) with isolate name beneath. Lab, BV-BRC and SeptaTyper columns with a star where SeptaTyper is right and BV-BRC is wrong. CSV/JSON download. Septatype-by-site bars and a country map. |
| C Submit | Community form with BV-BRC name lookup. Submissions wait in a curation queue. New septatypes and phenotype/prediction disagreements are flagged. A curator approves or rejects. |
| D Theory | The context-dependent grammar, an interactive schematic of the Mg2+ bridge, validation tables, edge cases. |

REST API under `/api/v1` (interactive docs at `/docs`).

## How a call is made

1. **Database majority rule** (the paper's rule). If the septatype was seen for that species and site,
   the majority phenotype of those isolates is the call. Ties go to Susceptible. Reproduces the paper's accuracies.
2. **Physicochemical grammar** (fallback only). For a septatype never seen for that species, four positions vote
   (hydroxyl lost, negative charge lost), `<H>` settles single votes, then the qnr, *S. sonnei* and *S. enterica* ctg overrides apply.
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
- **Assumptions to confirm:** `<H>` is the mean Eisenberg value of the four GyrA/ParC residues; G, P, C count as
  hydrophobic; the cutoff is 0.0; one shared *E. coli* K-12 reference set is used for all species.
- **Not built:** raw-read input, Mol* 3D structures (a labelled schematic is used instead), user accounts, Celery job queue.
- The Theory schematic is a reading of the grammar, not a structure. Align its wording with Figure 6.
- Fonts and the map library load from CDNs. The site works without them (fallback fonts, table instead of map).
- Declare AI-assisted code in the manuscript as required by the journal, and add a data licence.
