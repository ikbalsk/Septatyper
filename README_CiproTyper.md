# CiproTyper

**Septatype-based prediction of ciprofloxacin susceptibility in selected Enterobacterales**

CiproTyper is a web server, REST API and curated reference-isolate browser for exploring ciprofloxacin susceptibility in four supported bacterial groups. It analyses assembled nucleotide sequences, identifies selected chromosomal quinolone-resistance-associated positions and plasmid-mediated quinolone-resistance markers, and combines these features into a seven-field genotype called a **septatype**.

- **Web application:** https://septatyper-is.onrender.com/
- **Source code:** https://github.com/ikbalsk/Septatyper
- **Interactive API documentation:** https://septatyper-is.onrender.com/docs
- **Support:** ikbalsk.personal@gmail.com

> **Research use only.** CiproTyper is a research and surveillance resource. It is not a diagnostic test and must not be used to make treatment decisions. Predictions should be interpreted alongside laboratory antimicrobial susceptibility testing (AST).

## Overview

Ciprofloxacin resistance in Enterobacterales is associated with changes in the quinolone-resistance-determining regions (QRDRs) of DNA gyrase and topoisomerase IV, as well as plasmid-mediated mechanisms. The effect of an individual substitution is not necessarily the same in every genetic background. The septatype framework therefore records a combination of QRDR residues, a position in the DNA-repair protein RecN, and two plasmid-mediated quinolone-resistance markers. The biological rationale and the founding reference dataset are described by Ghosh *et al.* (2026) [1].

CiproTyper currently supports:

- *Escherichia coli*
- *Klebsiella pneumoniae*
- *Salmonella enterica* subsp. I
- *Shigella sonnei*

The available isolation-source cohorts are those represented in the reference dataset: stool, urine and respiratory specimens. Not every species is represented in every source category; predictions should be made only for the species and source combinations offered by the application.

## Reference dataset

The founding reference collection contains **2,734 genome records** representing **28 countries** and **64 distinct septatypes** [1]. The species distribution is shown below.

| Organism | Reference isolates |
|---|---:|
| *E. coli* | 400 |
| *K. pneumoniae* | 1,490 |
| *S. enterica* subsp. I | 510 |
| *S. sonnei* | 334 |
| **Total** | **2,734** |

The source distribution in the founding collection is 1,105 stool isolates, 1,231 urine isolates and 398 respiratory isolates. The live repository may change when additional submissions have been reviewed and approved; the counts above refer to the founding collection described in the project documentation and associated study.

### Interpreting the reported accuracy

The repository documentation reports **98.1% accuracy for the majority-rule predictions on the founding reference set** [1]. This is a resubstitution result: the same reference isolates were used to define the reference rules and to calculate this value. It is not an independent external validation estimate, and it should not be presented as evidence that the method will achieve the same accuracy on new, unrelated clinical isolates. Accuracy alone can also conceal poor performance for a less frequent phenotype. Sensitivity, specificity, class balance, confidence intervals and independent validation should be reported when evaluating performance in new datasets.

## Septatype definition

A septatype is a hyphen-delimited string with seven fields:

```text
[GyrA-83]-[GyrA-87]-[ParC-80]-[ParC-84]-[RecN-346]-[qnrB2]-[qnrS1]
```

| Field | Feature | Value recorded |
|---|---|---|
| 1 | GyrA position 83 | Amino-acid residue |
| 2 | GyrA position 87 | Amino-acid residue |
| 3 | ParC position 80 | Amino-acid residue |
| 4 | ParC position 84 | Amino-acid residue |
| 5 | RecN position 346 | Amino-acid residue |
| 6 | `qnrB2` | `1` if detected; otherwise `0` |
| 7 | `qnrS1` | `1` if detected; otherwise `0` |

For example, `L-N-I-E-L-0-0` indicates leucine at GyrA-83, asparagine at GyrA-87, isoleucine at ParC-80, glutamate at ParC-84, leucine at RecN-346, and no detected `qnrB2` or `qnrS1` marker.

**Important RecN detail:** the fifth field records the amino acid, not the nucleotide triplet. The prediction logic can also use the underlying RecN-346 codon. In the described rule, leucine encoded by `cta` is treated differently from leucine encoded by `ctg`; phenylalanine is also a resistance-associated state. Consequently, the septatype string alone does not preserve every nucleotide-level distinction at this position. When the nucleotide sequence is available, CiproTyper reads the codon from the alignment; manual entry can include the optional RecN codon where supported.

## Analysis workflow

### 1. Input and sequence validation

Users can submit assembled nucleotide sequences in FASTA format, paste FASTA text, analyse several assemblies in one batch, provide a BV-BRC [7] genome identifier, or enter the relevant residues manually. The web interface accepts `.fa`, `.fasta`, `.fna` and gzip-compressed FASTA files, with a default upload limit of 60 MB. Pasted sequence input is limited to 10 MB. Batch input supports up to 10 genome files.

Raw FASTQ reads are not accepted. Reads must first be assembled into contigs or a genome using an appropriate assembly workflow. Protein sequences are not suitable input for the nucleotide-based gene-identification stage.

### 2. Identification of chromosomal loci

NCBI BLAST+ [2] is used to locate `gyrA`, `parC` and `recN` in the submitted nucleotide sequence. The implementation uses a species-specific reference for each chromosomal locus. The reported search settings include an e-value threshold of `1e-20`; hits are accepted from 75% identity, while lower identity or coverage triggers a warning. The matching region is oriented as required and aligned to its reference using Biopython's pairwise alignment functionality [3]. The five target codons are then read from the alignment.

A position that cannot be resolved because of a gap, stop codon or ambiguous base is not assigned. If one of the five required positions cannot be called, a complete septatype and phenotype prediction are not returned.

### 3. Detection of plasmid-mediated markers

The sequence is also screened against the `qnrB2` and `qnrS1` reference alleles. A hit is recorded as present when it meets both thresholds:

- at least 80% nucleotide identity; and
- at least 90% query coverage.

For gene-only input, absence cannot be inferred from the sequence provided. The interface therefore allows users to specify the qnr status manually in the relevant sequence-input workflow.

### 4. MLST assignment

CiproTyper also reports a seven-locus multilocus sequence type (MLST), when the required allele profile can be matched. The project uses the Achtman scheme for *E. coli*, *S. sonnei* and *S. enterica* subsp. I, and the Diancourt scheme for *K. pneumoniae* [4,6]. A sequence type is reported only when all seven alleles match a recognised profile. Failure to assign an MLST type does not, by itself, prevent septatype analysis.

### 5. Septatype construction and phenotype prediction

After locus extraction, the five amino-acid states and two qnr flags are joined to form the septatype. The prediction procedure is organised into two tiers.

**Tier 1 – reference majority rule.** If the septatype is present for the selected species and isolation source, the most frequent laboratory phenotype among matching reference isolates is used. Ties are resolved as susceptible. If the septatype is not represented for that source but is present elsewhere for the same species, the species-level reference records are used and the result is flagged accordingly.

**Tier 2 – context-dependent grammar.** If the septatype has not been observed in the reference data, the implementation evaluates the residue and marker rules. These include:

1. **QRDR residue combination.** GyrA-83, GyrA-87, ParC-80 and ParC-84 are interpreted together rather than as independent additive substitutions. One documented resistant exception is the ParC-80/ParC-84 pair S-K. For combinations that do not fall wholly within the susceptible or resistant residue groups, the implementation uses the number of residues outside the susceptible group and, for the one-residue borderline case, the QRDR helix hydrophobicity calculation described below. A *S. sonnei*-specific rule is also implemented for single changes at GyrA-83 or GyrA-87.
2. **RecN-346.** Leucine encoded by `cta` is assigned to the susceptible state; leucine encoded by `ctg`, or phenylalanine, is assigned to the resistance-associated state.
3. **qnr markers.** Detection of either `qnrB2` or `qnrS1` is treated as a resistance-associated result in the grammar.

The grammar returns susceptible only when the applicable checkpoints support susceptibility; a resistance-associated checkpoint produces a resistant call. For a previously observed septatype, the Tier 1 majority result remains the primary call. The checkpoint explanations displayed by the interface should therefore not be mistaken for independent validation of the Tier 1 call.

### QRDR helix hydrophobicity

For the specified borderline combinations, CiproTyper calculates a mean hydrophobicity value, *H*, using the normalized consensus scale of Eisenberg *et al.* [5] The GyrA segment spans residues 82–93, and the ParC segment spans residues 80–90. The implementation calculates the mean for each segment and then takes the mean of the two segment means. The documented wild-type reference (S-D-S-E) has *H* = 0.270; the resistant reference (L-N-I-E) has *H* = 0.397. Their midpoint, **0.334**, is used as the decision threshold. This is a sequence-based descriptor used by the rule set; it is not a direct molecular-dynamics or three-dimensional binding-energy calculation.

## Web application

The web interface provides the following main functions:

- **Predict:** analyse a genome assembly, pasted FASTA sequence, batch of assemblies, BV-BRC genome identifier or manually entered residue profile.
- **Browse:** search and filter the reference isolates by organism, isolation source, MLST sequence type, country, septatype, phenotype and collection year. Filtered records can be exported as CSV or JSON.
- **Analytics:** inspect septatype distributions by isolation source and explore the geographic distribution of records.
- **Theory:** review the biological rationale and explore the residue-based grammar and QRDR illustration.
- **Submit:** submit an isolate and its laboratory phenotype for curator review.
- **API and data:** access API information and repository data downloads.

The site also provides project information and supporting explanatory pages. The exact layout may change as the application is updated.

## REST API

The API is intended for reproducible queries and integration into analysis scripts. The base path is:

```text
https://septatyper-is.onrender.com/api/v1
```

Interactive endpoint schemas are available at [`/docs`](https://septatyper-is.onrender.com/docs). The endpoint set documented in the repository includes:

| Method and endpoint | Purpose |
|---|---|
| `GET /meta` | Retrieve supported species, sources, record counts and reference metadata. |
| `POST /predict` | Predict from an uploaded FASTA file or supplied sequence; accepts species and habitat fields and optional qnr overrides. |
| `POST /predict/bvbrc` | Retrieve and analyse an assembly using a BV-BRC genome identifier. |
| `POST /predict/manual` | Predict from the five residue states and qnr flags, with an optional RecN codon. |
| `GET /isolates` | Retrieve filtered and paginated reference-isolate records. |
| `GET /isolates/{id}` | Retrieve one isolate record. |
| `GET /rules` | Retrieve the rule table for a species and isolation source. |
| `GET /grammar` | Retrieve the checkpoint interpretation for supplied residue states. |
| `GET /stats/accuracy` | Retrieve summary accuracy statistics. |
| `GET /stats/geography` | Retrieve country-level counts. |
| `GET /stats/septatypes` | Retrieve septatype frequency summaries. |
| `GET /stats/cutoff` | Retrieve the hydrophobicity cutoff information. |
| `GET /aminoacids` | Retrieve the amino-acid hydrophobicity scale used by the implementation. |
| `GET /download/csv` and `GET /download/json` | Download repository records, optionally filtered. |
| `GET /bvbrc/{id}` and `GET /bvbrc-status` | Retrieve BV-BRC lookup or loader status information. |
| `POST /submissions` | Add an isolate to the curation queue. |
| `GET /healthz` | Check service health. |

Administrative endpoints are available for curation and persistence operations and require the configured curator token. Use the interactive API documentation for the current request schemas, accepted parameter values and response fields; do not infer them from examples written for earlier API versions.

### Example: query the reference collection

```bash
curl -G "https://septatyper-is.onrender.com/api/v1/isolates" \
  --data-urlencode "species=Escherichia coli" \
  --data-urlencode "habitat=Urine" \
  --data-urlencode "septatype=L-N-I-E-L-0-0" \
  --data-urlencode "limit=5"
```

The hosted service may take longer to respond after a period of inactivity because the hosting platform can suspend idle services. A delayed first response should not automatically be interpreted as a failed analysis.

## Local deployment

The repository includes a Dockerfile and the application dependencies. The documented runtime is based on Python 3.12, FastAPI, Uvicorn, SQLite, NCBI BLAST+ and Biopython. Docker is the recommended starting point for reproducing the deployed environment because it can provide the required system-level BLAST+ installation alongside the Python application.

```bash
git clone https://github.com/ikbalsk/Septatyper.git
cd Septatyper
docker build -t ciprotyper .
docker run --rm -p 8000:8000 ciprotyper
```

After the container starts, open `http://localhost:8000`. The API documentation should be available at `http://localhost:8000/docs`.

For a deployment that enables curator functions, configure a strong `ADMIN_TOKEN` using the deployment platform's secret-management facility. Do not commit curator tokens, database credentials or backup tokens to the repository. Refer to the Dockerfile and the configuration module for the environment variables supported by the checked-out version. The repository documentation lists `DATA_DIR`, `DB_PATH`, `SEED_DB`, `REFS_DIR`, `PORT`, `MAX_UPLOAD_MB`, `BLAST_TIMEOUT_S`, `BVBRC_AUTOLOAD` and the `BACKUP_GITHUB_*` settings.

## Reference-data curation and data handling

The Submit page accepts new isolate records for review. A submission includes the sequence file and required submitter details, with optional isolate metadata such as BV-BRC identifier, country, collection year, affiliation and publication information. A laboratory phenotype is required for an isolate to be considered for approval. A curator reviews the submission before it is incorporated into the public reference collection.

For routine prediction, uploaded sequence data are processed for analysis and then deleted; the service documentation states that genome files are not retained. For community submissions, derived genotype results and the metadata supplied with the submission are retained for curation and, if approved, repository display. Submitters should not upload patient-identifiable information or data they are not authorised to share. The service's backup configuration can include submitter contact information and must be protected accordingly.

## Limitations

- **Restricted organism scope:** the rules are implemented for the four supported bacterial groups only. They should not be extrapolated to other Enterobacterales without validation.
- **Restricted resistance markers:** the septatype includes selected QRDR positions, the RecN-346 state and two qnr reference markers. It does not represent every known mechanism of fluoroquinolone resistance.
- **Assembly input:** raw FASTQ reads are not accepted. Sequence quality, assembly completeness and contamination can affect locus detection.
- **Reference dependence:** calls depend on the curated reference isolates, their laboratory phenotypes and the selected isolation source. Sparse or imbalanced cohorts can limit the reliability of phenotype estimates.
- **Resubstitution performance:** the reported 98.1% reference-set accuracy is not an independent test-set result. Prospective or external validation is needed before making claims about performance on unseen isolates.
- **No clinical validation:** the output is intended for research and surveillance, not patient-level diagnosis or antibiotic selection.

## Availability, licensing and maintenance

The source code is publicly accessible at the GitHub repository and the web service is available at the URL above. However, **a `LICENSE` file was not present in the repository when this README was reviewed**. Until the rights holder selects and adds an explicit software licence, users should not assume that the code may be reused, modified or redistributed. The data and reference sequences may be subject to separate terms from the software; in particular, do not assign a CC BY 4.0 licence to the reference dataset without confirming that the source publication and data providers permit it.

For a *Bioinformatics* Application Note submission, the authors should resolve the software and data licensing terms, archive the exact source release and test data in a persistent repository (for example, Zenodo) and provide the DOI. The journal's current software guidance also asks authors to ensure long-term availability and maintenance of the server, with support requirements extending beyond the date of publication.

## How to cite

If you use the septatype framework or its founding reference dataset, cite the study describing the method:

Ghosh A, Deb P, Bose P, Sk I, Biswas S, Mukhopadhyay AK, Ganguly D, Chattopadhyay S. Generic context-dependence governs clonal emergence of ciprofloxacin resistance mutations. *Scientific Reports*. 2026. https://doi.org/10.1038/s41598-026-75690-y.

Please also cite CiproTyper as software after a versioned software release and persistent archive DOI are available. Until then, cite the GitHub repository and record the access date in accordance with the target journal's software-citation guidance.

## References

1. Ghosh A, Deb P, Bose P, Sk I, Biswas S, Mukhopadhyay AK, Ganguly D, Chattopadhyay S. Generic context-dependence governs clonal emergence of ciprofloxacin resistance mutations. *Scientific Reports*. 2026. doi: [10.1038/s41598-026-75690-y](https://doi.org/10.1038/s41598-026-75690-y).
2. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. doi: [10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421).
3. Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–1423. doi: [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163).
4. Diancourt L, Passet V, Verhoef J, Grimont PAD, Brisse S. Multilocus sequence typing of *Klebsiella pneumoniae* nosocomial isolates. *Journal of Clinical Microbiology*. 2005;43(8):4178–4182. doi: [10.1128/JCM.43.8.4178-4182.2005](https://doi.org/10.1128/JCM.43.8.4178-4182.2005).
5. Eisenberg D, Schwarz E, Komaromy M, Wall R. Analysis of membrane and surface protein sequences with the hydrophobic moment plot. *Journal of Molecular Biology*. 1984;179(1):125–142. doi: [10.1016/0022-2836(84)90309-7](https://doi.org/10.1016/0022-2836(84)90309-7).
6. Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. *Wellcome Open Research*. 2018;3:124. doi: [10.12688/wellcomeopenres.14826.1](https://doi.org/10.12688/wellcomeopenres.14826.1).
7. Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. *Nucleic Acids Research*. 2023;51(D1):D678–D689. doi: [10.1093/nar/gkac1003](https://doi.org/10.1093/nar/gkac1003).

## Contact

For questions about predictions, the reference collection or software maintenance, contact **ikbalsk.personal@gmail.com** or open an issue in the [GitHub repository](https://github.com/ikbalsk/Septatyper).
