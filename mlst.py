"""MLST typing of the uploaded sequence with PubMLST schemes (database from the tseemann/mlst repository).

Schemes: E. coli and S. sonnei use the Achtman scheme, S. enterica the Achtman scheme, K. pneumoniae the
Diancourt scheme. An ST is only reported when all seven loci match a known allele exactly and the profile
exists in the scheme. Otherwise the nearest alleles are shown and the ST is left open.

Speed: instead of BLASTing every allele (13,000 sequences) against the genome, a few spread-out alleles per
locus locate the locus, the matching window is cut out and looked up in a dictionary of known alleles.
Only a locus that is not an exact hit falls back to BLASTing that locus's full allele set.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from .config import REFS_DIR

SCHEMES = {"ecoli": "ecoli_achtman_4", "ssonnei": "ecoli_achtman_4", "senterica": "senterica_achtman_2",
           "kpneumoniae": "klebsiella"}
NAMES = {"ecoli_achtman_4": "Achtman, 7 loci", "senterica_achtman_2": "Achtman, 7 loci", "klebsiella": "Diancourt, 7 loci"}
N_REPS = 4


def query_path(species: str):
    p = REFS_DIR / f"mlst_{SCHEMES[species]}.fasta"
    return p if p.exists() else None


@lru_cache(maxsize=4)
def _profiles(scheme: str):
    with open(REFS_DIR / f"mlst_{scheme}.txt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        loci = [c for c in rd.fieldnames[1:] if c != "clonal_complex"]
        prof = {tuple(r[l] for l in loci): r["ST"] for r in rd}
    return loci, prof


@lru_cache(maxsize=4)
def _alleles(scheme: str):
    lookup, full = {}, {}
    for rec in SeqIO.parse(REFS_DIR / f"mlst_{scheme}.fasta", "fasta"):
        locus, allele = rec.id.rsplit("_", 1)
        seq = str(rec.seq).upper()
        full.setdefault(locus, []).append((allele, seq))
        lk = lookup.setdefault(locus, {})
        if seq not in lk or int(allele) < int(lk[seq]):
            lk[seq] = allele
    reps = {}
    for locus, lst in full.items():
        idx = sorted({round(i * (len(lst) - 1) / (N_REPS - 1)) for i in range(N_REPS)})
        reps[locus] = [lst[i] for i in idx]
    return lookup, reps, full


def _window(seq: str, h: dict) -> str:
    """The stretch of the contig that would hold a full-length allele, given one representative hit."""
    s, e = h["sstart"], h["send"]
    fwd = s <= e
    lo, hi = (s, e) if fwd else (e, s)
    left, right = h["qstart"] - 1, h["qlen"] - h["qend"]
    lo, hi = (lo - left, hi + right) if fwd else (lo - right, hi + left)
    if lo < 1 or hi > len(seq):
        return ""
    w = seq[lo - 1:hi]
    return (w if fwd else str(Seq(w).reverse_complement())).upper()


def _call_slow(hits: list[dict]):
    """Exact allele, nearest allele, or missing, from a locus's full allele set BLASTed against the genome."""
    exact, near = set(), None
    for h in hits:
        allele = h["qseqid"].rsplit("_", 1)[1]
        if h["pident"] >= 99.999 and h["qstart"] == 1 and h["qend"] == h["qlen"]:
            exact.add(allele)
        elif h["qcov"] >= 0.9 and h["pident"] >= 90 and (near is None or h["bitscore"] > near[0]):
            near = (h["bitscore"], allele)
    if exact:
        return sorted(exact, key=int)[0], "exact"
    if near:
        return f"~{near[1]}", "novel"
    return "-", "missing"


def _by_locus(hits, locus):
    return [h for h in hits if h["qseqid"].rsplit("_", 1)[0] == locus]


def type_genome(species: str, seqs: dict, blast, workdir) -> dict:
    """blast(query_fasta_path, tag) -> list of hit dicts against the uploaded sequence."""
    scheme = SCHEMES[species]
    loci, prof = _profiles(scheme)
    lookup, reps, full = _alleles(scheme)
    wd = Path(workdir)
    q = wd / "mlst_reps.fa"
    q.write_text("".join(f">{l}_{a}\n{s}\n" for l in loci for a, s in reps[l]))
    hits = blast(q, "mlstfast")
    calls, state, slow = {}, {}, []
    for l in loci:
        hl = _by_locus(hits, l)
        if not hl:
            calls[l], state[l] = "-", "missing"
            continue
        best = max(hl, key=lambda x: x["bitscore"])
        allele = lookup[l].get(_window(seqs[best["sseqid"]], best))
        if allele:
            calls[l], state[l] = allele, "exact"
        else:
            slow.append(l)
    if slow:
        q2 = wd / "mlst_slow.fa"
        q2.write_text("".join(f">{l}_{a}\n{s}\n" for l in slow for a, s in full[l]))
        hits2 = blast(q2, "mlstslow")
        for l in slow:
            calls[l], state[l] = _call_slow(_by_locus(hits2, l))
    states = [state[l] for l in loci]
    out = {"scheme": scheme, "scheme_name": NAMES[scheme], "alleles": calls, "st": None,
           "loci_total": len(loci), "loci_exact": states.count("exact")}
    if all(s == "exact" for s in states):
        st = prof.get(tuple(calls[l] for l in loci))
        if st:
            out.update(st="ST" + st, status="typed", note="")
        else:
            out.update(status="new_combination", note="All seven alleles are known but this combination has no ST in the scheme yet.")
    elif "missing" in states:
        out.update(status="incomplete", note=f"Only {states.count('exact') + states.count('novel')} of {len(loci)} MLST loci were found, so no ST can be given. A gene-only upload cannot be typed.")
    else:
        out.update(status="novel_allele", note="At least one allele differs from every known allele. The nearest alleles are shown, marked ~.")
    return out
