"""Context-dependent grammar (physicochemical rules) for the GyrA/ParC tetratype.

Used only as a fallback when a septatype was never seen in the reference isolates.
Residue classes and voting rules follow the project's written rules; the places where
the rules were silent are marked ASSUMPTION and are easy to change here.

  GyrA-83, ParC-80 : residue with an -OH group (S, T, Y) keeps the Mg2+-water bridge:
                     no vote. Anything else (hydrophobic or no hydroxyl) votes Resistant.
  GyrA-87          : acidic residue (D, E) keeps the negative charge: no vote. Anything else votes Resistant.
  ParC-84          : a non-acidic residue votes Resistant. An acidic residue (E, D) votes Resistant
                     only when ParC-80 is hydrophobic (context dependence).
  Two or more votes = Resistant, zero = Susceptible. A single vote is settled by the combined
  mean hydrophobicity <H> of the four residues against 0.0 (the mean of the Eisenberg scale).
  S. sonnei with a single vote at GyrA-83 or GyrA-87 -> Resistant (data show L-D-S-E and S-Y-S-E are resistant).
  Septatype overrides (applied to a Susceptible baseline):
    qnrB or qnrS present            -> Resistant
    S. sonnei with RecN-346 F       -> Resistant
    S. enterica F-D-S-E, RecN ctg   -> Resistant (cta stays Susceptible)
"""
from __future__ import annotations

# Eisenberg normalised consensus hydrophobicity (Eisenberg et al. 1984, J Mol Biol 179:125)
EISENBERG = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29, "Q": -0.85, "E": -0.74, "G": 0.48,
             "H": -0.40, "I": 1.38, "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12, "S": -0.18,
             "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}
NAMES = {"A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate", "C": "Cysteine", "Q": "Glutamine",
         "E": "Glutamate", "G": "Glycine", "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
         "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine", "T": "Threonine",
         "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine"}
POLAR_OH = set("STY")
ACIDIC = set("DE")
HYDROPHOBIC = set("AVILMFWPGC")  # ASSUMPTION: non-polar side chains incl. G, P, C
BASIC = set("KRH")
H_CUTOFF = 0.0


def aa_class(a: str) -> str:
    if a in POLAR_OH:
        return "polar, hydroxyl"
    if a in ACIDIC:
        return "acidic (negative)"
    if a in BASIC:
        return "basic (positive)"
    if a in HYDROPHOBIC:
        return "hydrophobic"
    return "polar, uncharged"


def mean_h(aas4) -> float:
    return round(sum(EISENBERG[a] for a in aas4) / 4, 3)


def votes(g83, g87, p80, p84) -> list[dict]:
    p80_hydrophobic = p80 in HYDROPHOBIC
    v = [
        {"position": "GyrA-83", "residue": g83, "vote": g83 not in POLAR_OH,
         "why": "hydroxyl lost" if g83 not in POLAR_OH else "hydroxyl kept"},
        {"position": "GyrA-87", "residue": g87, "vote": g87 not in ACIDIC,
         "why": "negative charge lost" if g87 not in ACIDIC else "negative charge kept"},
        {"position": "ParC-80", "residue": p80, "vote": p80 not in POLAR_OH,
         "why": "hydroxyl lost" if p80 not in POLAR_OH else "hydroxyl kept"},
    ]
    if p84 in ACIDIC:
        v.append({"position": "ParC-84", "residue": p84, "vote": p80_hydrophobic,
                  "why": "acidic residue next to hydrophobic ParC-80" if p80_hydrophobic else "acidic residue in a polar context"})
    else:
        v.append({"position": "ParC-84", "residue": p84, "vote": True, "why": "negative charge lost"})
    return v


def predict(species: str, aas, qnrB: bool, qnrS: bool, recn_codon: str | None = None) -> dict:
    """aas = [g83, g87, p80, p84, r346]. Returns call, rule fired, votes, <H>."""
    g83, g87, p80, p84, r346 = aas
    vs = votes(g83, g87, p80, p84)
    n = sum(1 for x in vs if x["vote"])
    h = mean_h([g83, g87, p80, p84])
    if n >= 2:
        call, rule = "Resistant", f"{n} of 4 positions vote Resistant"
    elif n == 0:
        call, rule = "Susceptible", "no position votes Resistant"
    elif species == "ssonnei" and next(x for x in vs if x["vote"])["position"] in ("GyrA-83", "GyrA-87"):
        call, rule = "Resistant", "S. sonnei single GyrA-83 or GyrA-87 change (L-D-S-E, S-Y-S-E pattern)"
    else:
        call = "Resistant" if h > H_CUTOFF else "Susceptible"
        rule = f"1 vote; mean hydrophobicity {h} is {'above' if h > H_CUTOFF else 'at or below'} {H_CUTOFF}"
    if call == "Susceptible":
        if qnrB or qnrS:
            call, rule = "Resistant", "qnr gene present on a susceptible baseline (PMQR override)"
        elif species == "ssonnei" and r346 == "F":
            call, rule = "Resistant", "S. sonnei with RecN-346 F on a susceptible baseline"
        elif species == "senterica" and [g83, g87, p80, p84] == list("FDSE") and (recn_codon or "").upper() == "CTG":
            call, rule = "Resistant", "S. enterica F-D-S-E with RecN-346 codon ctg"
    return {"call": call, "rule": rule, "votes": vs, "n_votes": n, "mean_h": h}


def mechanism(aas) -> list[dict]:
    """Plain-language reading of each position, for the Theory page and result panel.
    A schematic reading of the grammar, not a structure calculation."""
    g83, g87, p80, p84, _ = aas
    out = []
    for pos, res in (("GyrA-83", g83), ("ParC-80", p80)):
        ok = res in POLAR_OH
        out.append({"position": pos, "residue": res, "effect": "intact" if ok else "lost", "kind": "H-bond",
                    "text": f"{res} keeps a hydroxyl that can hydrogen-bond the water-Mg2+ bridge to ciprofloxacin." if ok
                    else f"{res} has no hydroxyl, so the hydrogen bond to the water-Mg2+ bridge cannot form."})
    ok87 = g87 in ACIDIC
    out.append({"position": "GyrA-87", "residue": g87, "effect": "intact" if ok87 else "lost", "kind": "charge",
                "text": f"{g87} keeps the negative charge that holds the Mg2+ bridge." if ok87
                else f"{g87} removes the negative charge, so the Mg2+ bridge is held more weakly."})
    if p84 in ACIDIC and p80 in HYDROPHOBIC:
        out.append({"position": "ParC-84", "residue": p84, "effect": "context", "kind": "charge",
                    "text": f"{p84} is still acidic, but ParC-80 is hydrophobic, so the pocket around it is less polar (context-dependent)."})
    elif p84 in ACIDIC:
        out.append({"position": "ParC-84", "residue": p84, "effect": "intact", "kind": "charge",
                    "text": f"{p84} keeps the negative charge next to the Mg2+ bridge."})
    else:
        out.append({"position": "ParC-84", "residue": p84, "effect": "lost", "kind": "charge",
                    "text": f"{p84} is not acidic, so the negative charge near the Mg2+ bridge is lost."})
    return out
