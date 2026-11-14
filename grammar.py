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
  Checkpoint 1 is read from the residue table (all four residues in one column), checkpoint 2 from RecN-346 (cta Susceptible,
  ctg or F Resistant) and checkpoint 3 from qnrB2/qnrS1. The votes and <H> above are only worked out for information.
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
# ---------------------------------------------------------------------------------------------------------------------
# Checkpoint 1, the syntax table (the authors' "Checkpoint_1" table). Every position has a Susceptible and a Resistant column.
#   residues : the residues of the column, in the table's order; each string is one printed line
#   observed : the residues seen in the reference dataset (printed in bold); the others are possible residues, inferred from the
#              physicochemical features of the observed ones
#   props    : the property text of the column, line by line; each piece is (text, residues it describes), None = connecting word
# A residue can sit in both columns (D, E and N, Q at GyrA-87, D and E at ParC-84), so the four positions are read together.
# Edit this table to change the rule; nothing else needs to change.
TABLE = {
    "GyrA-83": {"S": {"residues": ["STY"], "observed": "SY", "props": [[("Polar with -OH group", "STY")]]},
                "R": {"residues": ["ALIVFWMCG"], "observed": "LIF", "props": [[("Hydrophobic", "ALIVFWMCG")]]}},
    "GyrA-87": {"S": {"residues": ["DENQ"], "observed": "DN",
                      "props": [[("Negatively Charged", "DE")], [("or its derivative", "NQ")]]},
                "R": {"residues": ["DENQSYTG", "ALIVFWMC"], "observed": "DNYGA",
                      "props": [[("Negatively charged", "DE")], [("or uncharged", "NQSYTGALIVFWMC")]]}},
    "ParC-80": {"S": {"residues": ["STY"], "observed": "S", "props": [[("Polar with -OH group", "STY")]]},
                "R": {"residues": ["GALIVFWMC"], "observed": "I", "props": [[("Hydrophobic", "GALIVFWMC")]]}},
    "ParC-84": {"S": {"residues": ["DE"], "observed": "E", "props": [[("Negatively Charged", "DE")]]},
                "R": {"residues": ["DEGALIVFWMC"], "observed": "EV",
                      "props": [[("Negatively Charged", "DE"), (" or ", None), ("hydrophobic", "GALIVFWMC")]]}},
}
TABLE_NOTE = ("*Residues observed in our dataset are marked in bold and remaining are possible residues based on "
              "physicochemical features of the observed one")
SUSCEPTIBLE_SET = {p: "".join(TABLE[p]["S"]["residues"]) for p in POS}
RESISTANT_SET = {p: "".join(TABLE[p]["R"]["residues"]) for p in POS}
ROWS = [(p, " ".join(SUSCEPTIBLE_SET[p]), " ".join(RESISTANT_SET[p])) for p in POS]
# The property text of each column as one sentence (Susceptible, Resistant)
CLASS_TEXT = {p: tuple(" ".join(" ".join("".join(t for t, _ in line) for line in TABLE[p][k]["props"]).split()) for k in ("S", "R"))
              for p in POS}
# The exception: the PAIR S at ParC-80 together with K at ParC-84 (the S-K pair, F-N-S-K and similar tetratypes) is counted as
# Resistant. S alone is the susceptible state and K alone is not in the table, so the exception applies only when both are present.
RESISTANT_EXCEPTION = {"ParC-80": "S", "ParC-84": "K"}
EXCEPTION_TEXT = "Exception: S at ParC-80 together with K at ParC-84 is Resistant"


def _side(pos: str, k: str, res: str, pair: bool = False) -> tuple[list[dict], list[dict], bool]:
    """One column of one row: its residues (a residue of the sample is a hit), its property text (a piece lights up when the
    sample residue is in the column and has that property), and whether the sample residue is in the column. In the Resistant
    column of ParC-80 and ParC-84 the exception residue (S, K) is added when the sample carries the whole S-K pair."""
    col = TABLE[pos][k]
    member = res in "".join(col["residues"])
    chips = [{"residue": a, "line": i, "observed": a in col["observed"], "hit": a == res, "exception": False}
             for i, line in enumerate(col["residues"]) for a in line]
    props = [{"text": t, "line": i, "hit": bool(member and rs and res in rs), "exception": False}
             for i, line in enumerate(col["props"]) for t, rs in line]
    if k == "R" and pair and RESISTANT_EXCEPTION.get(pos) == res:
        member = True
        chips.append({"residue": res, "line": len(col["residues"]) - 1, "observed": True, "hit": True, "exception": True})
        props.append({"text": EXCEPTION_TEXT, "line": len(col["props"]), "hit": True, "exception": True})
    return chips, props, member


def _cell(pos: str, res: str, pair: bool = False) -> dict:
    """The two columns of one row of the table, with the sample's residue marked."""
    s_chips, s_props, in_s = _side(pos, "S", res, pair)
    r_chips, r_props, in_r = _side(pos, "R", res, pair)
    return {"S": s_chips, "R": r_chips, "S_props": s_props, "R_props": r_props, "fits_s": in_s, "fits_r": in_r,
            "exception": bool(pair and RESISTANT_EXCEPTION.get(pos) == res)}


def checkpoint1(species: str, aas4) -> dict:
    """Checkpoint 1, the syntax table of GyrA-83, GyrA-87, ParC-80 and ParC-84. Each sample residue is marked in the column or
    columns of the table it belongs to. To give a call ALL FOUR residues must be in one column: all four in the Susceptible
    column is Susceptible, all four in the Resistant column is Resistant. If they are not all in one column (some residues
    may be in both), checkpoint 1 is INCONCLUSIVE and the phenotype comes from checkpoints 2 and 3 alone. The hydrophobicity
    <H> plays no part in this decision (it is only worked out for information). The S-K pair at ParC-80/ParC-84 counts as
    Resistant (RESISTANT_EXCEPTION)."""
    pair = aas4[2] == RESISTANT_EXCEPTION["ParC-80"] and aas4[3] == RESISTANT_EXCEPTION["ParC-84"]
    cells = []
    for pos, res in zip(POS, aas4):
        c = _cell(pos, res, pair)
        cells.append({"position": pos, "residue": res, **c})
    all_s, all_r = all(c["fits_s"] for c in cells), all(c["fits_r"] for c in cells)
    pattern = "susceptible" if all_s else "resistant" if all_r else None
    vs = votes(*aas4)
    n = sum(1 for x in vs if x["vote"])
    override = None                  # kept in the reply for compatibility; there is no override any more
    applied = True
    used_exc = [c["position"] for c in cells if c["exception"]] if all_r else []
    if all_r:
        call, by, rule = "Resistant", "pattern", "all four residues are in the Resistant column" + (" (exception: the S-K pair at ParC-80 and ParC-84)" if used_exc else "")
    elif all_s:
        call, by, rule = "Susceptible", "pattern", "all four residues are in the Susceptible column"
    else:
        call, by, applied = None, "not_applied", False
        rule = "the four residues are not all in one column of the table, so checkpoint 1 is inconclusive"
    return {"call": call, "applied": applied, "pattern": pattern, "resolved_by": by, "rule": rule, "cells": cells, "votes": vs,
            "n_votes": n, "mean_h": mean_h(list(aas4)), "override": override, "exceptions": used_exc,
            "marks": {"susceptible": sum(1 for c in cells if c["fits_s"]), "resistant": sum(1 for c in cells if c["fits_r"])}}


def checkpoint2(r346: str, codon: str | None = None) -> dict:
    """Checkpoint 2, RecN-346. Susceptible if leucine (L) encoded by cta. Resistant if leucine encoded by ctg, or phenylalanine (F).
    Anything else (leucine with another codon or no codon given, any other amino acid, or RecN-346 not read) is NOT considered:
    the call is "Not considered" and checkpoint 2 takes no part in the phenotype."""
    cd = (codon or "").upper() or None
    out = lambda call, text: {"call": call, "applied": call != "Not considered", "residue": r346 or None, "codon": cd, "text": text}
    if r346 == "F":
        return out("Resistant", "RecN-346 is phenylalanine (F)")
    if r346 == "L" and cd == "CTG":
        return out("Resistant", "RecN-346 is leucine encoded by ctg")
    if r346 == "L" and cd == "CTA":
        return out("Susceptible", "RecN-346 is leucine encoded by cta")
    if not r346:
        return out("Not considered", "RecN-346 could not be read from this sequence, so checkpoint 2 is not considered")
    if r346 == "L":
        why = f"leucine encoded by {cd.lower()}, which is neither cta nor ctg" if cd else "leucine, but the codon is not known (cta and ctg cannot be told apart)"
        return out("Not considered", f"RecN-346 is {why}, so checkpoint 2 is not considered")
    return out("Not considered", f"RecN-346 is {r346}, which is not leucine or phenylalanine, so checkpoint 2 is not considered")


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
        call, decided, rule = "Resistant", "checkpoint2", f"checkpoint 2: {c2['text']} (checkpoint 1 inconclusive)"
    elif c3["call"] == "Resistant":
        call, decided, rule = "Resistant", "checkpoint3", f"checkpoint 3: {c3['text']} (checkpoint 1 inconclusive)"
    else:
        call, decided, rule = "Susceptible", "checkpoints_2_3", "checkpoints 2 and 3 are both without a resistant signal (checkpoint 1 inconclusive)"
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
    elif c1["call"] is None:        # checkpoint 1 inconclusive: Susceptible because no checkpoint that was considered says Resistant
        call, rule, decided = "Susceptible", "checkpoints 2 and 3 are both without a resistant signal (checkpoint 1 inconclusive)", "checkpoints_2_3"
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
        pos.append({"position": cell["position"], "residue": cell["residue"], "class": aa_class(cell["residue"]), "column": column,
                    "fits_susceptible": fs, "fits_resistant": fr, "exception": cell["exception"],
                    "susceptible": cell["S"], "resistant": cell["R"], "susceptible_props": cell["S_props"],
                    "resistant_props": cell["R_props"], "note": note,
                    "susceptible_label": CLASS_TEXT[cell["position"]][0], "resistant_label": CLASS_TEXT[cell["position"]][1]})
    sample = h_block([p["residue"] for p in pos])
    ref = h_block(list(ref4) if ref4 else list(WILD4))
    return {"decided_by": g["decided_by"], "n_resistant_positions": g["n_votes"], "positions": pos,
            "checkpoints": {"residues": {"call": c1["call"], "pattern": c1["pattern"], "resolved_by": c1["resolved_by"], "rule": c1["rule"],
                                         "applied": c1["applied"], "n_votes": c1["n_votes"], "override": c1["override"], "exceptions": c1["exceptions"],
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
