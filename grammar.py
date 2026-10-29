"""Context-dependent grammar (physicochemical rules) for the GyrA/ParC tetratype.

Used only as a fallback when a septatype was never seen in the reference isolates.
Residue classes and voting rules follow the project's written rules; the places where
the rules were silent are marked ASSUMPTION and are easy to change here.

  GyrA-83, ParC-80 : residue with an -OH group (S, T, Y) keeps the Mg2+-water bridge:
                     no vote. Anything else (hydrophobic/no hydroxyl group) votes Resistant.
  GyrA-87          : acidic residue (D, E) keeps the negative charge: no vote. Anything else votes Resistant.
  ParC-84          : a non-acidic residue votes Resistant. An acidic residue (E, D) votes Resistant
                     only when ParC-80 is hydrophobic (context dependence).
  Two or more votes = Resistant, zero = Susceptible. A single vote is settled by the QRDR helix mean
  hydrophobicity <H> against the cutoff (see below).
  <H> is read over the whole QRDR helix of both enzymes, with the sample's residues put in at the four
  tetratype positions:
      GyrA helix 82-93  D(S)AVY(D)TIVRMA   12 residues, mean = sum / 12   (S = GyrA-83, D = GyrA-87)
      ParC helix 80-90  (S)ACY(E)AMVLMA    11 residues, mean = sum / 11   (S = ParC-80, E = ParC-84)
      <H> = (GyrA mean + ParC mean) / 2
  Cutoff: the midpoint between the <H> of the susceptible wild-type helices (S-D-S-E, 0.270) and the <H> of the
  resistant reference helices (L-N-I-E, 0.397), which is 0.334. It is set from these two reference helices,
  not by searching the isolates for the best score. A cutoff at the wild-type value alone (0.270) would turn
  single-vote septatypes such as S. enterica F-D-S-E (0.327) resistant whatever the codon.
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
POLAR_UNCHARGED = set("STYNQ")     # polar side chain, no charge
HYDROPHOBIC = set("AVILMFWPGC")  # ASSUMPTION: non-polar side chains incl. G, P, C
BASIC = set("KRH")


def aa_class(a: str) -> str:
    if a in POLAR_OH:
        return "polar, hydroxyl group"
    if a in ACIDIC:
        return "acidic (negative)"
    if a in BASIC:
        return "basic (positive)"
    if a in HYDROPHOBIC:
        return "hydrophobic"
    return "polar, uncharged"


# QRDR helices (E. coli numbering). The four tetratype positions are replaced by the sample's residues.
GYR_HELIX, GYR_START = "DSAVYDTIVRMA", 82      # GyrA 82-93; GyrA-83 and GyrA-87 sit at helix places 2 and 6
PARC_HELIX, PARC_START = "SACYEAMVLMA", 80     # ParC 80-90; ParC-80 and ParC-84 sit at helix places 1 and 5
WILD4 = "SDSE"                                 # wild-type, susceptible residues at the four positions
RES4 = "LNIE"                                  # resistant reference: the commonest fully changed tetratype in the reference isolates


def _helices(aas4):
    g83, g87, p80, p84 = aas4
    g, p = list(GYR_HELIX), list(PARC_HELIX)
    g[83 - GYR_START], g[87 - GYR_START] = g83, g87
    p[80 - PARC_START], p[84 - PARC_START] = p80, p84
    return g, p


def _helix_means(aas4):
    g, p = _helices(aas4)
    return sum(EISENBERG[a] for a in g) / len(g), sum(EISENBERG[a] for a in p) / len(p)


def mean_h(aas4) -> float:
    """QRDR helix mean hydrophobicity: (GyrA helix mean over 12 residues + ParC helix mean over 11 residues) / 2."""
    hg, hp = _helix_means(aas4)
    return round((hg + hp) / 2, 3)


def _unrounded(aas4) -> float:
    hg, hp = _helix_means(aas4)
    return (hg + hp) / 2


H_WILD, H_RES = mean_h(list(WILD4)), mean_h(list(RES4))                        # 0.270 and 0.397
H_CUTOFF = round((_unrounded(list(WILD4)) + _unrounded(list(RES4))) / 2, 3)   # 0.334, the midpoint of the two reference helices


def votes(g83, g87, p80, p84) -> list[dict]:
    p80_hydrophobic = p80 in HYDROPHOBIC
    v = [
        {"position": "GyrA-83", "residue": g83, "vote": g83 not in POLAR_OH,
         "why": "hydroxyl group lost" if g83 not in POLAR_OH else "hydroxyl group kept"},
        {"position": "GyrA-87", "residue": g87, "vote": g87 not in ACIDIC,
         "why": "negative charge lost" if g87 not in ACIDIC else "negative charge kept"},
        {"position": "ParC-80", "residue": p80, "vote": p80 not in POLAR_OH,
         "why": "hydroxyl group lost" if p80 not in POLAR_OH else "hydroxyl group kept"},
    ]
    if p84 in ACIDIC:
        v.append({"position": "ParC-84", "residue": p84, "vote": p80_hydrophobic,
                  "why": "acidic residue next to hydrophobic ParC-80" if p80_hydrophobic else "acidic residue in a polar context"})
    else:
        v.append({"position": "ParC-84", "residue": p84, "vote": True, "why": "negative charge lost"})
    return v


POS = ("GyrA-83", "GyrA-87", "ParC-80", "ParC-84")
# Checkpoint 1, the syntax rules table: the residues that make each position Susceptible or Resistant.
# A residue can sit in both columns (D at GyrA-87, E at ParC-84), so the four positions are read together.
# The one exception counts as Resistant and is shown only when the sample carries it: the ParC pair S-K (Fig. 6b of the paper),
# that is S at ParC-80 together with K at ParC-84. S alone (the susceptible S-E state) or K alone is not an exception.
SUSCEPTIBLE_SET = {"GyrA-83": "SYT", "GyrA-87": "DNE", "ParC-80": "SYT", "ParC-84": "DE"}
RESISTANT_SET = {"GyrA-83": "FLIW", "GyrA-87": "NDGYAE", "ParC-80": "IFW", "ParC-84": "EVD"}
RESISTANT_EXCEPTION = {"ParC-80": "S", "ParC-84": "K"}      # one pair: both must be present
ROWS = [(p, " ".join(SUSCEPTIBLE_SET[p]), " ".join(RESISTANT_SET[p])) for p in POS]


def _cell(pos: str, res: str, pair: bool = False) -> dict:
    """The two columns of one row of the table, each residue flagged when it is the sample's. An exception residue
    is added to the Resistant column only when the sample carries the whole S-K pair (pair=True)."""
    exc = RESISTANT_EXCEPTION.get(pos, "") if pair and res == RESISTANT_EXCEPTION.get(pos) else ""
    return {"S": [{"residue": a, "hit": a == res, "exception": False} for a in SUSCEPTIBLE_SET[pos]],
            "R": [{"residue": a, "hit": a == res, "exception": False} for a in RESISTANT_SET[pos]]
                 + ([{"residue": exc, "hit": True, "exception": True}] if exc else []),
            "fits_s": res in SUSCEPTIBLE_SET[pos], "fits_r": res in RESISTANT_SET[pos] or bool(exc), "exception": bool(exc)}


def checkpoint1(species: str, aas4) -> dict:
    """Checkpoint 1, the context-sensitive grammar of GyrA-83, GyrA-87, ParC-80 and ParC-84, from the syntax rules table.
    Susceptible when all four residues are in the Susceptible column, Resistant when all four are in the Resistant column
    (the exception pair S-K included). Where the four fit neither, the written rules are silent: the residue votes and <H>
    settle it (ASSUMPTION), and S. sonnei keeps its extra single-change rule."""
    pair = aas4[2] == RESISTANT_EXCEPTION["ParC-80"] and aas4[3] == RESISTANT_EXCEPTION["ParC-84"]
    cells = []
    for pos, res in zip(POS, aas4):
        c = _cell(pos, res, pair)
        cells.append({"position": pos, "residue": res, **c})
    pattern = "susceptible" if all(c["fits_s"] for c in cells) else "resistant" if all(c["fits_r"] for c in cells) else None
    vs = votes(*aas4)
    n = sum(1 for x in vs if x["vote"])
    h = mean_h(list(aas4))
    voter = next((x["position"] for x in vs if x["vote"]), None)
    used_exc = [c["position"] for c in cells if c["exception"]] if pattern == "resistant" else []
    override = None
    if pattern == "resistant":
        call, by = "Resistant", "pattern"
        rule = "all four residues are in the Resistant column" + (" (exception: the ParC pair S-K)" if used_exc else "")
    elif species == "ssonnei" and n == 1 and voter in ("GyrA-83", "GyrA-87"):
        call, by, override = "Resistant", "override", "sonnei_gyr"
        rule = "S. sonnei single GyrA-83 or GyrA-87 change"
    elif pattern == "susceptible":
        call, by, rule = "Susceptible", "pattern", "all four residues are in the Susceptible column"
    elif n >= 2:
        call, by, rule = "Resistant", "votes", f"not all four in one column; {n} of 4 positions vote Resistant"
    elif n == 0:
        call, by, rule = "Susceptible", "votes", "not all four in one column; no position votes Resistant"
    else:
        cut = H_CUTOFF if _CUT[0] is None else _CUT[0]
        call = "Resistant" if h > cut else "Susceptible"
        by = "hydrophobicity"
        rule = f"not all four in one column; 1 vote, so the QRDR helix mean hydrophobicity {h} is {'above' if h > cut else 'at or below'} {cut}"
    return {"call": call, "pattern": pattern, "resolved_by": by, "rule": rule, "cells": cells, "votes": vs, "n_votes": n,
            "mean_h": h, "override": override, "exceptions": used_exc}


def checkpoint2(r346: str, codon: str | None = None) -> dict:
    """Checkpoint 2, RecN-346. Susceptible if leucine encoded by cta; Resistant if leucine encoded by ctg, or phenylalanine."""
    cd = (codon or "").upper() or None
    if r346 == "F":
        return {"call": "Resistant", "residue": r346, "codon": cd, "text": "RecN-346 is phenylalanine (F)"}
    if r346 == "L" and cd == "CTG":
        return {"call": "Resistant", "residue": r346, "codon": cd, "text": "RecN-346 is leucine encoded by ctg"}
    if r346 == "L" and cd == "CTA":
        return {"call": "Susceptible", "residue": r346, "codon": cd, "text": "RecN-346 is leucine encoded by cta"}
    if r346 == "L":
        why = f"leucine encoded by {cd.lower()}, which is neither cta nor ctg" if cd else "leucine, codon not given, so ctg cannot be checked"
        return {"call": "Susceptible", "residue": r346, "codon": cd, "text": f"RecN-346 is {why}: no resistant signal"}
    return {"call": "Susceptible", "residue": r346, "codon": cd, "text": f"RecN-346 is {r346}, which the rules do not cover: no resistant signal"}


def checkpoint3(qnrB: bool, qnrS: bool) -> dict:
    """Checkpoint 3, PMQR genes. Resistant if qnrB2 or qnrS1 (or both) is present, Susceptible if both are absent."""
    on = [n for n, v in (("qnrB2", qnrB), ("qnrS1", qnrS)) if v]
    return {"call": "Resistant" if on else "Susceptible", "qnrB2": bool(qnrB), "qnrS1": bool(qnrS),
            "text": (" and ".join(on) + " present") if on else "qnrB2 and qnrS1 both absent"}


_CUT = [None]      # lets the calculation page score other <H> cutoffs without changing the one the site uses


def predict(species: str, aas, qnrB: bool, qnrS: bool, recn_codon: str | None = None, cutoff: float | None = None) -> dict:
    """aas = [g83, g87, p80, p84, r346]. The three checkpoints, then: Susceptible if all say Susceptible,
    Resistant if any says Resistant."""
    g83, g87, p80, p84, r346 = aas
    _CUT[0] = cutoff
    try:
        c1 = checkpoint1(species, [g83, g87, p80, p84])
    finally:
        _CUT[0] = None
    c2, c3 = checkpoint2(r346, recn_codon), checkpoint3(qnrB, qnrS)
    if c1["call"] == "Resistant":
        call, rule, decided = "Resistant", c1["rule"], c1["resolved_by"]
    elif c2["call"] == "Resistant":
        call, rule, decided = "Resistant", f"checkpoint 2: {c2['text']}", "checkpoint2"
    elif c3["call"] == "Resistant":
        call, rule, decided = "Resistant", f"checkpoint 3: {c3['text']}", "checkpoint3"
    else:
        call, rule, decided = "Susceptible", c1["rule"], c1["resolved_by"]
    return {"call": call, "rule": rule, "votes": c1["votes"], "n_votes": c1["n_votes"], "mean_h": c1["mean_h"],
            "override": c1["override"], "decided_by": decided,
            "checkpoints": {"residues": c1, "recn": c2, "qnr": c3}}


def h_block(aas4) -> dict:
    """QRDR helix hydrophobicity worked out for four tetratype residues: every residue's H for the GyrA helix
    (12 residues) and the ParC helix (11 residues), each helix mean, and the combined mean."""
    g, p = _helices(aas4)

    def blk(res, start, label, subs):
        vals = [{"position": f"{label}-{start + i}", "residue": a, "h": EISENBERG[a], "sub": (start + i) in subs}
                for i, a in enumerate(res)]
        tot = sum(v["h"] for v in vals)
        return {"gene": label, "start": start, "end": start + len(res) - 1, "n": len(res), "values": vals,
                "sum": round(tot, 3), "mean": round(tot / len(res), 3)}
    return {"residues": "".join(aas4), "gyr": blk(g, GYR_START, "GyrA", (83, 87)),
            "parc": blk(p, PARC_START, "ParC", (80, 84)), "mean": mean_h(aas4)}


def explain(species: str, aas, g: dict, ref4=None, qnr=(False, False), codon=None) -> dict:
    """The three checkpoints behind a call, laid out like the authors' residue-class table, plus the <H> worked out
    for the cases where the residue classes alone do not decide checkpoint 1."""
    c1, c2, c3 = g["checkpoints"]["residues"], g["checkpoints"]["recn"], g["checkpoints"]["qnr"]
    pos = []
    for cell in c1["cells"]:
        fs, fr = cell["fits_s"], cell["fits_r"]
        column = "both" if fs and fr else "susceptible" if fs else "resistant" if fr else "neither"
        note = {"both": "in both columns", "neither": "in neither column"}.get(column, "")
        if cell["exception"]:
            note = "exception residue"
        pos.append({"position": cell["position"], "residue": cell["residue"], "class": aa_class(cell["residue"]), "column": column,
                    "fits_susceptible": fs, "fits_resistant": fr, "exception": cell["exception"],
                    "susceptible": cell["S"], "resistant": cell["R"], "note": note})
    sample = h_block([p["residue"] for p in pos])
    ref = h_block(list(ref4) if ref4 else list(WILD4))
    return {"decided_by": g["decided_by"], "n_resistant_positions": g["n_votes"], "positions": pos,
            "checkpoints": {"residues": {"call": c1["call"], "pattern": c1["pattern"], "resolved_by": c1["resolved_by"], "rule": c1["rule"],
                                         "n_votes": c1["n_votes"], "override": c1["override"], "exceptions": c1["exceptions"]},
                            "recn": c2, "qnr": c3},
            "hydrophobicity": {"gyr": sample["gyr"], "parc": sample["parc"], "mean": sample["mean"], "cutoff": H_CUTOFF,
                               "side": "above" if sample["mean"] > H_CUTOFF else "at_or_below"},
            "reference": ref, "grammar_call": g["call"], "rule": g["rule"]}


def mechanism(aas) -> list[dict]:
    """Plain-language reading of each position, for the result panel and the Theory page.
    A schematic reading of the grammar, not a structure calculation."""
    g83, g87, p80, p84, _ = aas
    out = []
    for pos, res in (("GyrA-83", g83), ("ParC-80", p80)):
        ok = res in POLAR_OH
        out.append({"position": pos, "residue": res, "effect": "intact" if ok else "lost", "kind": "H-bond",
                    "text": f"{res} keeps its hydroxyl group, which can hydrogen-bond the water-Mg2+ bridge that holds ciprofloxacin." if ok
                    else f"{res} has no hydroxyl group, so it cannot hydrogen-bond the water-Mg2+ bridge."})
    ok87 = g87 in ACIDIC
    out.append({"position": "GyrA-87", "residue": g87, "effect": "intact" if ok87 else "lost", "kind": "charge",
                "text": f"{g87} keeps the negative charge that holds the Mg2+ bridge." if ok87
                else f"{g87} has no negative charge, so the Mg2+ bridge is held more weakly."})
    if p84 in ACIDIC and p80 in HYDROPHOBIC:
        out.append({"position": "ParC-84", "residue": p84, "effect": "context", "kind": "charge",
                    "text": f"{p84} is still acidic, but ParC-80 is hydrophobic, so the pocket around it is less polar."})
    elif p84 in ACIDIC:
        out.append({"position": "ParC-84", "residue": p84, "effect": "intact", "kind": "charge",
                    "text": f"{p84} keeps the negative charge next to the Mg2+ bridge."})
    else:
        out.append({"position": "ParC-84", "residue": p84, "effect": "lost", "kind": "charge",
                    "text": f"{p84} is not acidic, so the negative charge near the Mg2+ bridge is lost."})
    return out
