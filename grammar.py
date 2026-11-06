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
# Checkpoint 1, the syntax rules table. Each position has a Susceptible and a Resistant column. A residue can sit in both
# columns (D and E at GyrA-87, E and D at ParC-84), so the four positions are read together.
# The residue lists below are the DATA behind the table (they decide which column a residue is marked in). The table on the
# site no longer prints them: it shows the amino-acid PROPERTIES of each column (PROP_TAGS) and lights up the sample's residue
# and the properties that residue has. Edit a list here to move a residue between columns.
#   GyrA-83 Resistant : hydrophobic residues without an -OH group (F L I W, and A V M)
#   GyrA-87 Susceptible: negatively charged residue or its amide derivative (D E, and the amides N Q)
#   GyrA-87 Resistant  : hydrophobic, or negatively charged residue or its amide derivative
#   ParC-80 Resistant  : hydrophobic residues without an -OH group (I F W, and L A V M), plus the exception S
#   ParC-84 Resistant  : hydrophobic or negatively charged residue, plus the exception K
# The one exception counts as Resistant and is marked only when the sample carries it: the ParC pair S-K (Fig. 6b of the paper),
# that is S at ParC-80 together with K at ParC-84. S alone (the susceptible S-E state) or K alone is not an exception.
SUSCEPTIBLE_SET = {"GyrA-83": "SYT", "GyrA-87": "DNEQ", "ParC-80": "SYT", "ParC-84": "DE"}
RESISTANT_SET = {"GyrA-83": "FLIWAVM", "GyrA-87": "NDGYAEQ", "ParC-80": "IFWLAVM", "ParC-84": "EVD"}
RESISTANT_EXCEPTION = {"ParC-80": "S", "ParC-84": "K"}      # one pair: both must be present
ROWS = [(p, " ".join(SUSCEPTIBLE_SET[p]), " ".join(RESISTANT_SET[p])) for p in POS]
# The wording of each column, shown beside the residues in the API: (Susceptible sentence, Resistant sentence)
CLASS_TEXT = {"GyrA-83": ("Polar residue with -OH group", "Hydrophobic residue without -OH group"),
              "GyrA-87": ("Negatively charged residue or its amide derivative", "Hydrophobic residue, or negatively charged residue or its amide derivative"),
              "ParC-80": ("Polar residue with -OH group", "Hydrophobic residue without -OH group"),
              "ParC-84": ("Negatively charged residue", "Hydrophobic or negatively charged residue")}

# Amino-acid properties used by the table. A property lights up for the sample only in a column the residue is marked in.
#   hydrophobic : positive on the Eisenberg scale (A V I L M F W P G C and Y; ASSUMPTION for display, Y is +0.26)
#   oh          : has an -OH group (S T Y)          no_oh : has none
#   polar       : polar side chain (S T Y N Q)
#   neg         : negatively charged (D E)          neg_deriv : amide derivative of a negatively charged residue (N from D, Q from E)
RES_PROPS = {a: {k for k, on in (("hydrophobic", EISENBERG[a] > 0), ("oh", a in POLAR_OH), ("no_oh", a not in POLAR_OH),
                                 ("polar", a in POLAR_UNCHARGED), ("neg", a in ACIDIC), ("neg_deriv", a in "NQ")) if on}
             for a in EISENBERG}
_P = {"polar": "Polar", "oh": "-OH group", "no_oh": "No -OH group", "hydrophobic": "Hydrophobic",
      "neg": "Negatively charged", "neg_deriv": "Amide derivative of a negatively charged residue"}
# the property tags printed in each column of the table, in order: (key, label)
PROP_TAGS = {"GyrA-83": {"S": ["polar", "oh"], "R": ["hydrophobic", "no_oh"]},
             "GyrA-87": {"S": ["neg", "neg_deriv"], "R": ["hydrophobic", "neg", "neg_deriv"]},
             "ParC-80": {"S": ["polar", "oh"], "R": ["hydrophobic", "no_oh"]},
             "ParC-84": {"S": ["neg"], "R": ["hydrophobic", "neg"]}}
EXCEPTION_TAG = {"ParC-80": "Exception: S together with K at ParC-84", "ParC-84": "Exception: K together with S at ParC-80"}


def _tags(pos: str, side: str, res: str, member: bool, exc: bool) -> list[dict]:
    """The property tags of one column. A tag is lit (hit) when the sample residue is marked in this column and has the property."""
    out = [{"key": k, "label": _P[k], "hit": bool(member and k in RES_PROPS.get(res, ())), "exception": False} for k in PROP_TAGS[pos][side]]
    if side == "R" and pos in EXCEPTION_TAG:
        out.append({"key": "exception", "label": EXCEPTION_TAG[pos], "hit": exc, "exception": True})
    return out


def _cell(pos: str, res: str, pair: bool = False) -> dict:
    """The two columns of one row of the table. A residue is marked in a column when it belongs to it; an exception residue
    is marked in the Resistant column only when the sample carries the whole S-K pair (pair=True)."""
    exc = RESISTANT_EXCEPTION.get(pos, "") if pair and res == RESISTANT_EXCEPTION.get(pos) else ""
    in_s, in_r = res in SUSCEPTIBLE_SET[pos], res in RESISTANT_SET[pos] or bool(exc)
    return {"S": [{"residue": a, "hit": a == res, "exception": False} for a in SUSCEPTIBLE_SET[pos]],
            "R": [{"residue": a, "hit": a == res, "exception": False} for a in RESISTANT_SET[pos]]
                 + ([{"residue": exc, "hit": True, "exception": True}] if exc else []),
            "S_props": _tags(pos, "S", res, in_s, False), "R_props": _tags(pos, "R", res, in_r, bool(exc)),
            "fits_s": in_s, "fits_r": in_r, "exception": bool(exc)}


def checkpoint1(species: str, aas4) -> dict:
    """Checkpoint 1, the context-sensitive grammar of GyrA-83, GyrA-87, ParC-80 and ParC-84, from the syntax rules table.
    Each of the four sample residues is marked in the Susceptible column, the Resistant column, or both (D at GyrA-87 and
    E or D at ParC-84 sit in both). The column with MORE marks decides: Susceptible when more positions are marked in the
    Susceptible column, Resistant when more are marked in the Resistant column (all four in one column is the clearest case).
    Only when the two columns have the same number of marks does the table give no answer; then the residue votes and <H>
    settle it (ASSUMPTION). S. sonnei keeps its extra single-change rule."""
    pair = aas4[2] == RESISTANT_EXCEPTION["ParC-80"] and aas4[3] == RESISTANT_EXCEPTION["ParC-84"]
    cells = []
    for pos, res in zip(POS, aas4):
        c = _cell(pos, res, pair)
        cells.append({"position": pos, "residue": res, **c})
    n_s, n_r = sum(1 for c in cells if c["fits_s"]), sum(1 for c in cells if c["fits_r"])
    pattern = "susceptible" if n_s == 4 else "resistant" if n_r == 4 else None
    vs = votes(*aas4)
    n = sum(1 for x in vs if x["vote"])
    h = mean_h(list(aas4))
    voter = next((x["position"] for x in vs if x["vote"]), None)
    used_exc = [c["position"] for c in cells if c["exception"]] if n_r > n_s else []
    override = None
    if pattern == "resistant":
        call, by = "Resistant", "pattern"
        rule = "all four residues are in the Resistant column" + (" (exception: the ParC pair S-K)" if used_exc else "")
    elif species == "ssonnei" and n == 1 and voter in ("GyrA-83", "GyrA-87"):
        call, by, override = "Resistant", "override", "sonnei_gyr"
        rule = "S. sonnei single GyrA-83 or GyrA-87 change"
    elif pattern == "susceptible":
        call, by, rule = "Susceptible", "pattern", "all four residues are in the Susceptible column"
    elif n_s > n_r:
        call, by = "Susceptible", "marks"
        rule = f"{n_s} positions are marked in the Susceptible column and {n_r} in the Resistant column: the column with more marks decides"
    elif n_r > n_s:
        call, by = "Resistant", "marks"
        rule = (f"{n_r} positions are marked in the Resistant column and {n_s} in the Susceptible column: the column with more marks decides"
                + (" (exception: the ParC pair S-K)" if used_exc else ""))
    elif n >= 2:
        call, by, rule = "Resistant", "votes", f"both columns have {n_s} marks; {n} of 4 positions vote Resistant"
    elif n == 0:
        call, by, rule = "Susceptible", "votes", f"both columns have {n_s} marks; no position votes Resistant"
    else:
        cut = H_CUTOFF if _CUT[0] is None else _CUT[0]
        call = "Resistant" if h > cut else "Susceptible"
        by = "hydrophobicity"
        rule = f"both columns have {n_s} marks; 1 vote, so the QRDR helix mean hydrophobicity {h} is {'above' if h > cut else 'at or below'} {cut}"
    return {"call": call, "pattern": pattern, "resolved_by": by, "rule": rule, "cells": cells, "votes": vs, "n_votes": n,
            "mean_h": h, "override": override, "exceptions": used_exc, "marks": {"susceptible": n_s, "resistant": n_r}}


def checkpoint2(r346: str, codon: str | None = None) -> dict:
    """Checkpoint 2, RecN-346. Susceptible if leucine encoded by cta; Resistant if leucine encoded by ctg, or phenylalanine."""
    cd = (codon or "").upper() or None
    if not r346:
        return {"call": "Not determined", "residue": None, "codon": cd, "text": "RecN-346 could not be read from this sequence: no signal from checkpoint 2"}
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


def predict_without_checkpoint1(r346, recn_codon, qnrB: bool, qnrS: bool, why: str) -> dict:
    """The call when checkpoint 1 cannot be applied (a gene matches the species reference at too low an identity):
    checkpoints 2 and 3 only. Resistant if either says Resistant, else Susceptible."""
    c2, c3 = checkpoint2(r346, recn_codon), checkpoint3(qnrB, qnrS)
    c1 = {"call": None, "applied": False, "pattern": None, "resolved_by": "not_applied", "rule": why, "cells": [], "votes": [],
          "n_votes": 0, "mean_h": None, "override": None, "exceptions": []}
    if c2["call"] == "Resistant":
        call, decided, rule = "Resistant", "checkpoint2", f"checkpoint 2: {c2['text']} (checkpoint 1 not applied)"
    elif c3["call"] == "Resistant":
        call, decided, rule = "Resistant", "checkpoint3", f"checkpoint 3: {c3['text']} (checkpoint 1 not applied)"
    else:
        call, decided, rule = "Susceptible", "checkpoints_2_3", "checkpoints 2 and 3 are both without a resistant signal (checkpoint 1 not applied)"
    return {"call": call, "rule": rule, "votes": [], "n_votes": 0, "mean_h": None, "override": None, "decided_by": decided,
            "checkpoints": {"residues": c1, "recn": c2, "qnr": c3}}


def explain_not_applied(g: dict, why: str) -> dict:
    """The criteria block for a call made without checkpoint 1."""
    c1, c2, c3 = g["checkpoints"]["residues"], g["checkpoints"]["recn"], g["checkpoints"]["qnr"]
    return {"decided_by": g["decided_by"], "n_resistant_positions": 0, "positions": [],
            "checkpoints": {"residues": {"call": None, "applied": False, "pattern": None, "resolved_by": "not_applied", "rule": why,
                                         "n_votes": 0, "override": None, "exceptions": []}, "recn": c2, "qnr": c3},
            "hydrophobicity": None, "reference": None, "grammar_call": g["call"], "rule": g["rule"]}


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
                    "susceptible": cell["S"], "resistant": cell["R"], "susceptible_props": cell["S_props"],
                    "resistant_props": cell["R_props"], "note": note,
                    "susceptible_label": CLASS_TEXT[cell["position"]][0], "resistant_label": CLASS_TEXT[cell["position"]][1]})
    sample = h_block([p["residue"] for p in pos])
    ref = h_block(list(ref4) if ref4 else list(WILD4))
    return {"decided_by": g["decided_by"], "n_resistant_positions": g["n_votes"], "positions": pos,
            "checkpoints": {"residues": {"call": c1["call"], "pattern": c1["pattern"], "resolved_by": c1["resolved_by"], "rule": c1["rule"],
                                         "n_votes": c1["n_votes"], "override": c1["override"], "exceptions": c1["exceptions"],
                                         "marks": c1["marks"]},
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
