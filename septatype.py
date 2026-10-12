"""Core analysis: find gyrA/parC/recN and the qnr genes in an uploaded nucleotide
FASTA, read the five septatype codons, and assemble the septatype string.

Works for whole-genome assemblies and for gene-level FASTA files alike, because
genes are located by BLAST rather than by contig name.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from .config import BLAST_TIMEOUT_S, POSITIONS, REFS_DIR, WILD

# thresholds. QNR values are the ones stated in the manuscript Methods
# (>=80% identity, >=90% gene-length coverage). The chromosomal ones are
# practical guards, not taken from the paper.
MIN_ID_CHROM = 75.0
WARN_ID_CHROM = 92.0
MIN_COV = 0.90
QNR_MIN_ID = 80.0
QNR_MIN_COV = 0.90
FLANK = 30


class InputError(ValueError):
    """The user's file is unusable."""


class RefsMissing(RuntimeError):
    """Reference sequences have not been installed yet."""


# ---------------------------------------------------------------- references
# refs/gyrA.fasta, parC.fasta, recN.fasta hold one record per species (header: gene|Species name|accession).
# refs/qnrB2.fasta and qnrS1.fasta hold the qnr reference alleles.
REF_SPECIES = {"ecoli": "Escherichia coli", "ssonnei": "Shigella sonnei", "senterica": "Salmonella enterica",
               "kpneumoniae": "Klebsiella pneumoniae"}


def ref_files(species: str = "") -> dict[str, Path]:
    out = {g: REFS_DIR / f"{g}.fasta" for g in ("gyrA", "parC", "recN")}
    out["qnrB"], out["qnrS"] = REFS_DIR / "qnrB2.fasta", REFS_DIR / "qnrS1.fasta"
    return out


def _species_record(path: Path, species: str):
    want = REF_SPECIES[species].lower()
    for rec in SeqIO.parse(path, "fasta"):
        if want in rec.description.lower():
            return rec
    return None


def refs_status() -> dict:
    from .config import SPECIES
    f, out = ref_files(), {}
    for sp in SPECIES:
        missing = [p.name for p in f.values() if not p.exists()]
        if not missing:
            missing = [f"{g}.fasta has no {REF_SPECIES[sp]} record" for g in ("gyrA", "parC", "recN") if _species_record(f[g], sp) is None]
        out[sp] = {"ready": not missing, "missing": missing}
    return out


def _load_refs(species: str):
    f = ref_files()
    missing = [p.name for p in f.values() if not p.exists()]
    if missing:
        raise RefsMissing("Reference sequences not installed: " + ", ".join(missing))
    chrom = {}
    for g in ("gyrA", "parC", "recN"):
        rec = _species_record(f[g], species)
        if rec is None:
            raise RefsMissing(f"{g}.fasta has no record for {REF_SPECIES[species]}.")
        chrom[g] = rec
    return chrom, f["qnrB"], f["qnrS"]


# --------------------------------------------------------------------- input
def read_fasta(path: Path):
    """Return ({id: seq}, {id: original name}, total_bp). Accepts .gz."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, names = {}, {}
    with opener(path, "rt") as fh:
        try:
            for i, rec in enumerate(SeqIO.parse(fh, "fasta")):
                s = str(rec.seq).upper().replace("U", "T")
                if s:
                    seqs[f"c{i}"], names[f"c{i}"] = s, rec.id
        except (UnicodeDecodeError, OSError, ValueError):
            raise InputError("Could not read the file as FASTA. It must be plain text starting with a '>' header line (or .gz).")
    if not seqs:
        raise InputError("No FASTA records found. The file must start with a '>' header line.")
    total = sum(len(s) for s in seqs.values())
    good = sum(s.count(c) for s in seqs.values() for c in "ACGTN")
    if good / total < 0.9:
        raise InputError("This looks like protein sequence. Upload nucleotide FASTA (genome or gene sequences).")
    return seqs, names, total


# --------------------------------------------------------------------- BLAST
def _blast(query: Path, asm: Path, workdir: Path, tag: str) -> list[dict]:
    db = workdir / f"db_{tag}"
    try:
        subprocess.run(["makeblastdb", "-in", str(asm), "-dbtype", "nucl", "-parse_seqids", "-out", str(db)],
                       check=True, capture_output=True, timeout=BLAST_TIMEOUT_S)
        fmt = "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore"
        res = subprocess.run(
            ["blastn", "-task", "blastn", "-query", str(query), "-db", str(db), "-outfmt", fmt,
             "-evalue", "1e-20", "-dust", "no", "-max_target_seqs", "50"],
            check=True, capture_output=True, text=True, timeout=BLAST_TIMEOUT_S)
    except FileNotFoundError:
        raise RuntimeError("BLAST+ is not installed on the server (makeblastdb/blastn not found).")
    except subprocess.TimeoutExpired:
        raise InputError("Analysis timed out. Try a smaller file or a single assembly.")
    except subprocess.CalledProcessError as e:
        raise InputError("BLAST failed on this input: " + (e.stderr.decode()[:200] if isinstance(e.stderr, bytes) else str(e.stderr)[:200]))
    keys = ["qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    hits = []
    for line in res.stdout.splitlines():
        v = line.split("\t")
        h = dict(zip(keys, v))
        for k in ("pident", "bitscore"):
            h[k] = float(h[k])
        for k in ("length", "qlen", "qstart", "qend", "sstart", "send"):
            h[k] = int(h[k])
        h["qcov"] = (h["qend"] - h["qstart"] + 1) / h["qlen"]
        hits.append(h)
    return hits


# ------------------------------------------------------------ codon calling
def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -8, -1
    a.end_gap_score = 0   # flanks of the extracted region are free
    return a


def _extract_region(seq: str, hit: dict) -> str:
    """Cut the contig so it spans the whole reference gene (plus a small flank)."""
    s, e = hit["sstart"], hit["send"]
    fwd = s <= e
    lo, hi = (s, e) if fwd else (e, s)
    left, right = hit["qstart"] - 1, hit["qlen"] - hit["qend"]
    if fwd:
        lo, hi = lo - left - FLANK, hi + right + FLANK
    else:
        lo, hi = lo - right - FLANK, hi + left + FLANK
    lo, hi = max(1, lo), min(len(seq), hi)
    region = seq[lo - 1:hi]
    return region if fwd else str(Seq(region).reverse_complement())


def _ref_to_region_map(region: str, ref: str) -> list[int]:
    aln = _aligner().align(region, ref)[0]
    m = [-1] * len(ref)
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(t1 - t0):
            m[q0 + k] = t0 + k
    return m


def _codon(m: list[int], region: str, pos: int):
    i = (pos - 1) * 3
    if i + 2 >= len(m):
        return None
    a, b, c = m[i], m[i + 1], m[i + 2]
    if min(a, b, c) < 0 or b != a + 1 or c != b + 1:
        return None
    return region[a:a + 3]


def _aa(codon: str | None) -> str | None:
    if codon is None:
        return None
    if any(c not in "ACGT" for c in codon):
        return "X"
    return str(Seq(codon).translate())


# -------------------------------------------------------------------- main
def septatype_string(aas: list[str], qnrB2: bool, qnrS1: bool) -> str:
    return "-".join(list(aas) + [str(int(qnrB2)), str(int(qnrS1))])


def analyse(fasta: Path, species: str, qnr_override: dict | None = None) -> dict:
    """qnr_override: {'qnrB2': 'present'|'absent'|None, 'qnrS1': ...} for gene-level input."""
    qnr_override = qnr_override or {}
    chrom, qnrB_fa, qnrS_fa = _load_refs(species)
    seqs, names, total = read_fasta(fasta)
    warnings: list[str] = []

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        asm = td / "asm.fa"
        with open(asm, "w") as fh:
            for k, s in seqs.items():
                fh.write(f">{k}\n{s}\n")
        q1 = td / "chrom.fa"
        with open(q1, "w") as fh:
            for g, rec in chrom.items():
                fh.write(f">{g}\n{str(rec.seq).upper()}\n")
        hits = _blast(q1, asm, td, "chrom")

        genes, positions = {}, {}
        for g, rec in chrom.items():
            ref = str(rec.seq).upper()
            cand = [h for h in hits if h["qseqid"] == g]
            best = max(cand, key=lambda h: h["bitscore"]) if cand else None
            info = {"found": False}
            if best and best["pident"] >= MIN_ID_CHROM:
                info = {"found": True, "identity": round(best["pident"], 1), "coverage": round(100 * best["qcov"], 1),
                        "contig": names[best["sseqid"]], "strand": "+" if best["sstart"] <= best["send"] else "-"}
                if best["qcov"] < MIN_COV:
                    warnings.append(f"{g}: only {info['coverage']}% of the reference gene was covered (partial gene or contig break).")
                if best["pident"] < WARN_ID_CHROM:
                    warnings.append(f"{g}: identity to the {species} reference is {info['identity']}%. "
                                    f"identity to the {REF_SPECIES[species]} reference is lower than expected. Check the selected species and the assembly.")
                region = _extract_region(seqs[best["sseqid"]], best)
                gmap = _ref_to_region_map(region, ref)
            else:
                warnings.append(f"{g}: not found (best identity {best['pident'] if best else 0}%). "
                                "Wrong species, or the gene is missing from the input.")
                gmap = region = None
            genes[g] = info
            for gg, pos in POSITIONS:
                if gg != g:
                    continue
                key = f"{g}{pos}"
                p = {"gene": g, "pos": pos, "codon": None, "aa": None}
                ref_codon = ref[(pos - 1) * 3:(pos - 1) * 3 + 3]
                p["ref_aa"] = _aa(ref_codon) if len(ref_codon) == 3 else None
                if gmap is not None:
                    p["codon"] = _codon(gmap, region, pos)
                    p["aa"] = _aa(p["codon"])
                if p["aa"] is None and gmap is not None:
                    warnings.append(f"{g}-{pos}: position not covered or contains an indel.")
                if p["aa"] in ("X", "*"):
                    warnings.append(f"{g}-{pos}: codon {p['codon']} is ambiguous or a stop; cannot be called.")
                if p["ref_aa"] and p["ref_aa"] != WILD[(g, pos)]:
                    warnings.append(f"{g}-{pos}: the installed reference has {p['ref_aa']} here, expected {WILD[(g, pos)]}. "
                                    "Check reference numbering.")
                positions[key] = p

        # qnr families: presence/absence
        qnr = {}
        for fam, fa, label in (("qnrB", qnrB_fa, "qnrB2"), ("qnrS", qnrS_fa, "qnrS1")):
            h = _blast(fa, asm, td, fam)
            ok = [x for x in h if x["pident"] >= QNR_MIN_ID and x["qcov"] >= QNR_MIN_COV]
            best = max(ok, key=lambda x: x["bitscore"]) if ok else None
            present = best is not None
            src = "detected"
            ov = qnr_override.get(label)
            if ov in ("present", "absent"):
                present, src = ov == "present", "set by user"
            qnr[label] = {"present": present, "source": src,
                          "best_hit": (f"{best['qseqid']} ({best['pident']:.1f}% id, {100*best['qcov']:.0f}% cov)" if best else None)}
        if total < 1_000_000 and not any(qnr_override.get(k) for k in ("qnrB2", "qnrS1")):
            warnings.append("The input is smaller than a bacterial genome. qnr 'absent' only means no qnr "
                            "sequence was in the file. Set qnrB/qnrS manually if you know them.")

    aas = [positions[f"{g}{p}"]["aa"] for g, p in POSITIONS]
    complete = all(a and a not in ("X", "*") for a in aas)
    sep = septatype_string(aas, qnr["qnrB2"]["present"], qnr["qnrS1"]["present"]) if complete else None
    return {"genes": genes, "positions": positions, "qnr": qnr, "septatype": sep,
            "recN_codon": positions["recN346"]["codon"],
            "recN_syn": ({"CTA": "cta", "CTG": "ctg"}.get(positions["recN346"]["codon"] or "", "other")
                         if positions["recN346"]["aa"] == "L" else None),
            "warnings": warnings,
            "input": {"records": len(seqs), "total_bp": total}}
