"""MLST typing of the uploaded sequence with PubMLST schemes (database from the tseemann/mlst repository).

Schemes: E. coli and S. sonnei use the Achtman scheme, S. enterica the Achtman scheme, K. pneumoniae the
Diancourt scheme. The alleles of a scheme are BLASTed against the uploaded sequence. A locus is an exact
match at 100% identity over the full allele length. An ST is only reported when all seven loci are exact
and the profile exists in the scheme. Otherwise the nearest alleles are shown and the ST is left open.
"""
from __future__ import annotations

import csv
from functools import lru_cache

from .config import REFS_DIR

SCHEMES = {"ecoli": "ecoli_achtman_4", "ssonnei": "ecoli_achtman_4", "senterica": "senterica_achtman_2",
           "kpneumoniae": "klebsiella"}
NAMES = {"ecoli_achtman_4": "Achtman, 7 loci", "senterica_achtman_2": "Achtman, 7 loci", "klebsiella": "Diancourt, 7 loci"}


def query_path(species: str):
    p = REFS_DIR / f"mlst_{SCHEMES[species]}.fasta"
    return p if p.exists() else None


@lru_cache(maxsize=4)
def _profiles(scheme: str):
    with open(REFS_DIR / f"mlst_{scheme}.txt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        loci = [c for c in rd.fieldnames[1:] if c != "clonal_complex"]
        prof = {tuple(r[l] for l in loci): (r["ST"], (r.get("clonal_complex") or "").strip()) for r in rd}
    return loci, prof


def type_hits(species: str, hits: list[dict]) -> dict:
    scheme = SCHEMES[species]
    loci, prof = _profiles(scheme)
    exact, near = {}, {}
    for h in hits:
        locus, allele = h["qseqid"].rsplit("_", 1)
        full = h["pident"] >= 99.999 and h["qstart"] == 1 and h["qend"] == h["qlen"]
        if full:
            exact.setdefault(locus, set()).add(allele)
        elif h["qcov"] >= 0.9 and h["pident"] >= 90:
            if locus not in near or h["bitscore"] > near[locus][0]:
                near[locus] = (h["bitscore"], allele, round(h["pident"], 1))
    calls, state = {}, []
    for l in loci:
        if l in exact:
            a = sorted(exact[l], key=int)[0]
            calls[l] = a
            state.append("exact")
        elif l in near:
            calls[l] = f"~{near[l][1]}"
            state.append("novel")
        else:
            calls[l] = "-"
            state.append("missing")
    out = {"scheme": scheme, "scheme_name": NAMES[scheme], "alleles": calls, "st": None, "cc": None,
           "loci_total": len(loci), "loci_exact": state.count("exact")}
    if all(s == "exact" for s in state):
        hit = prof.get(tuple(calls[l] for l in loci))
        if hit:
            out.update(st="ST" + hit[0], cc=hit[1] or None, status="typed", note="")
        else:
            out.update(status="new_combination", note="All seven alleles are known but this combination has no ST in the scheme yet.")
    elif "missing" in state:
        out.update(status="incomplete", note=f"Only {state.count('exact') + state.count('novel')} of {len(loci)} MLST loci were found, so no ST can be given. A gene-only upload cannot be typed.")
    else:
        out.update(status="novel_allele", note="At least one allele differs from every known allele. The nearest alleles are shown, marked ~.")
    return out
