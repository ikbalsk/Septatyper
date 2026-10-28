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


def predict(species: str, aas, qnrB: bool, qnrS: bool, recn_codon: str | None = None, cutoff: float | None = None) -> dict:
    """aas = [g83, g87, p80, p84, r346]. Returns call, rule fired, votes, <H>."""
    g83, g87, p80, p84, r346 = aas
    vs = votes(g83, g87, p80, p84)
    n = sum(1 for x in vs if x["vote"])
    h = mean_h([g83, g87, p80, p84])
    cut = H_CUTOFF if cutoff is None else cutoff       # only the calculation page tries other cutoffs
    override, decided = None, "votes"
    if n >= 2:
        call, rule = "Resistant", f"{n} of 4 positions vote Resistant"
    elif n == 0:
        call, rule = "Susceptible", "no position votes Resistant"
    elif species == "ssonnei" and next(x for x in vs if x["vote"])["position"] in ("GyrA-83", "GyrA-87"):
        call, rule, override, decided = "Resistant", "S. sonnei single GyrA-83 or GyrA-87 change (L-D-S-E, S-Y-S-E pattern)", "sonnei_gyr", "override"
    else:
        call = "Resistant" if h > cut else "Susceptible"
        rule, decided = f"1 vote; QRDR helix mean hydrophobicity {h} is {'above' if h > cut else 'at or below'} {cut}", "hydrophobicity"
    if call == "Susceptible":
        if qnrB or qnrS:
            call, rule, override, decided = "Resistant", "qnr gene present on a susceptible baseline (PMQR override)", "qnr", "override"
        elif species == "ssonnei" and r346 == "F":
            call, rule, override, decided = "Resistant", "S. sonnei with RecN-346 F on a susceptible baseline", "sonnei_recn", "override"
        elif species == "senterica" and [g83, g87, p80, p84] == list("FDSE") and (recn_codon or "").upper() == "CTG":
            call, rule, override, decided = "Resistant", "S. enterica F-D-S-E with RecN-346 codon ctg", "senterica_ctg", "override"
    return {"call": call, "rule": rule, "votes": vs, "n_votes": n, "mean_h": h, "override": override, "decided_by": decided}


OVERRIDES = [("qnr", None, "qnr gene present"),
             ("sonnei_gyr", "ssonnei", "One change at GyrA-83 or GyrA-87"),
             ("sonnei_recn", "ssonnei", "RecN-346 is F"),
             ("senterica_ctg", "senterica", "F-D-S-E with the ctg codon at RecN-346")]
ROWS = [("GyrA-83", "Polar residue with -OH group", "Hydrophobic residue"),
        ("GyrA-87", "Negatively charged residue", "Polar or uncharged residue"),
        ("ParC-80", "Polar residue with -OH group", "Hydrophobic residue"),
        ("ParC-84", "Negatively charged residue", "Uncharged or negatively charged residue")]


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


def _terms(pos: str, res: str, vote: bool):
    """The words of the authors' residue-class table, each flagged when it describes this residue.
    Returns (susceptible terms, resistant terms, susceptible suffix, resistant suffix)."""
    charged = ACIDIC | BASIC
    if pos in ("GyrA-83", "ParC-80"):
        S, R, sS, sR = [("Polar residue with -OH group", res in POLAR_OH)], [("Hydrophobic", res in HYDROPHOBIC)], "", " residue"
    elif pos == "GyrA-87":
        S = [("Negatively charged", res in ACIDIC)]
        R, sS, sR = [("Polar", res in POLAR_UNCHARGED), ("uncharged", res not in charged)], " residue", " residue"
    else:
        S = [("Negatively charged", res in ACIDIC and not vote)]
        R, sS, sR = [("Uncharged", res not in charged), ("negatively charged", res in ACIDIC and vote)], " residue", " residue"
    active = R if vote else S
    if not any(hit for _, hit in active):      # no single word fits, so the whole cell is the match
        active[:] = [(t, True) for t, _ in active]
    return S, R, sS, sR


def explain(species: str, aas, g: dict, ref4=None, qnr=(False, False), codon=None) -> dict:
    """Three criteria behind a call, laid out like the authors' residue-class table.
    1 residue class per position, 2 QRDR helix mean hydrophobicity against the cutoff, 3 septatype overrides."""
    g83, g87, p80, p84, _ = aas
    pos = []
    for v, (name, s_txt, r_txt), res in zip(g["votes"], ROWS, (g83, g87, p80, p84)):
        note = ""
        if name in ("GyrA-83", "ParC-80") and res not in POLAR_OH and res not in HYDROPHOBIC:
            note = "no -OH group, counted as resistant"
        elif name == "GyrA-87" and res in BASIC:
            note = "positively charged, counted as resistant"
        elif name == "ParC-84" and res in BASIC:
            note = "positively charged, counted as resistant"
        elif name == "ParC-84" and res in ACIDIC:
            note = ("acidic, but ParC-80 is hydrophobic: counted as resistant" if v["vote"]
                    else "acidic, and ParC-80 is not hydrophobic: stays susceptible")
        S, R, sS, sR = _terms(name, res, v["vote"])
        fmt = lambda terms, suf: [{"text": t, "hit": h} for t, h in terms]  # noqa: E731
        pos.append({"position": name, "residue": res, "class": aa_class(res), "column": "resistant" if v["vote"] else "susceptible",
                    "susceptible": s_txt, "resistant": r_txt,
                    "susceptible_parts": fmt(S, sS), "susceptible_suffix": sS,
                    "resistant_parts": fmt(R, sR), "resistant_suffix": sR, "note": note})
    sample = h_block([p["residue"] for p in pos])
    ref = h_block(list(ref4) if ref4 else list(WILD4))
    voters = [v["position"] for v in g["votes"] if v["vote"]]
    met = {"qnr": bool(qnr[0] or qnr[1]),
           "sonnei_gyr": species == "ssonnei" and len(voters) == 1 and voters[0] in ("GyrA-83", "GyrA-87"),
           "sonnei_recn": species == "ssonnei" and aas[4] == "F",
           "senterica_ctg": species == "senterica" and [g83, g87, p80, p84] == list("FDSE") and (codon or "").upper() == "CTG"}
    ov = [{"key": k, "text": t, "fired": g["override"] == k, "met": met[k]} for k, sp, t in OVERRIDES if sp in (None, species)]
    return {"decided_by": g["decided_by"], "n_resistant_positions": g["n_votes"], "positions": pos,
            "hydrophobicity": {"gyr": sample["gyr"], "parc": sample["parc"], "mean": sample["mean"], "cutoff": H_CUTOFF,
                               "side": "above" if sample["mean"] > H_CUTOFF else "at_or_below"},
            "reference": ref, "overrides": ov, "grammar_call": g["call"], "rule": g["rule"]}


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
