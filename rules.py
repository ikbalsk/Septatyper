"""Database access, prediction (database lookup first, grammar fallback second), statistics, curation."""
from __future__ import annotations

import datetime
import json
import shutil
import sqlite3
import uuid

from . import grammar
from . import countries
from .config import DB_PATH, HABITATS, SEED_DB, SPECIES

LOW_SUPPORT = 5       # fewer reference isolates than this -> flag
MIXED = 0.80          # majority share below this -> flag
ISO_COLS = ("genome_id", "genome_name", "isolate_name", "species", "habitat", "lab_phenotype", "year_raw",
            "year_start", "year_end", "country", "iso3", "st", "gyrA_83", "gyrA_87", "parC_80", "parC_84",
            "recN_346", "qnrB2", "qnrS1", "septatype", "tetratype", "septatyper_call", "bvbrc_call",
            "bvbrc_method", "source", "date_added")


def ensure_db():
    """Live DB lives in DATA_DIR (put it on a persistent disk). First start copies the baked seed."""
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not DB_PATH.exists():
        if not SEED_DB.exists():
            raise RuntimeError("Seed database missing. Run: python scripts/build_db.py")
        shutil.copy(SEED_DB, DB_PATH)


def conn() -> sqlite3.Connection:
    ensure_db()
    c = sqlite3.connect(DB_PATH, check_same_thread=False)
    c.row_factory = sqlite3.Row
    return c


def _persist():
    """Ask the backup module to save soon. Does nothing unless automatic backup is configured."""
    try:
        from . import persist
        persist.schedule_save()
    except Exception:  # noqa: BLE001  a backup problem must never break a submission or a decision
        pass


def _now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def _label(p: str) -> str:
    return "Resistant" if p == "R" else "Susceptible"


def reference_sts(db, species: str, habitat: str, sep: str, top: int = 10) -> list[dict]:
    """Sequence types of the reference isolates that carry this septatype (the isolates behind a majority-rule call)."""
    where, args = "species=? AND septatype=?", [species, sep]
    if habitat != "ALL":
        where += " AND habitat=?"
        args.append(habitat)
    rows = db.execute(f"""SELECT COALESCE(st,'ST not recorded') st, COUNT(*) n, SUM(lab_phenotype='Resistant') r
                          FROM isolates WHERE {where} GROUP BY 1 ORDER BY n DESC, st""", args).fetchall()
    out = [{"st": r["st"], "n": r["n"], "n_resistant": r["r"], "n_susceptible": r["n"] - r["r"]} for r in rows[:top]]
    rest = rows[top:]
    if rest:
        n, r = sum(x["n"] for x in rest), sum(x["r"] for x in rest)
        out.append({"st": f"{len(rest)} other STs", "n": n, "n_resistant": r, "n_susceptible": n - r})
    return out


def reference_sites(db, species: str, habitat: str, sep: str) -> list[dict]:
    """Isolation sites (fecal, urine, ...) of the reference isolates that carry this septatype."""
    where, args = "species=? AND septatype=?", [species, sep]
    if habitat != "ALL":
        where += " AND habitat=?"
        args.append(habitat)
    rows = db.execute(f"""SELECT habitat, COUNT(*) n, SUM(lab_phenotype='Resistant') r FROM isolates
                          WHERE {where} GROUP BY habitat ORDER BY n DESC, habitat""", args).fetchall()
    return [{"habitat": r["habitat"], "n": r["n"], "n_resistant": r["r"], "n_susceptible": r["n"] - r["r"]} for r in rows]


def reference_countries(db, species: str, habitat: str, sep: str) -> list[dict]:
    """Countries of the reference isolates that carry this septatype, with the sites they came from."""
    where, args = "species=? AND septatype=?", [species, sep]
    if habitat != "ALL":
        where += " AND habitat=?"
        args.append(habitat)
    rows = db.execute(f"""SELECT COALESCE(country,'not recorded') country, habitat, COUNT(*) n,
                          SUM(lab_phenotype='Resistant') r FROM isolates WHERE {where}
                          GROUP BY 1, habitat""", args).fetchall()
    by = {}
    for r in rows:
        d = by.setdefault(r["country"], {"country": r["country"], "n": 0, "n_resistant": 0, "n_susceptible": 0, "sites": {}})
        d["n"] += r["n"]
        d["n_resistant"] += r["r"]
        d["n_susceptible"] += r["n"] - r["r"]
        d["sites"][r["habitat"]] = r["n"]
    return sorted(by.values(), key=lambda d: (-d["n"], d["country"]))


# ---------------------------------------------------------------- prediction
def _rule(db, species, habitat, sep):
    return db.execute("SELECT * FROM septatypes WHERE species=? AND habitat=? AND septatype=?",
                      (species, habitat, sep)).fetchone()


def predict(species: str, habitat: str, sep: str, recn_codon: str | None = None) -> dict:
    if species not in SPECIES:
        raise ValueError("Unknown species")
    if habitat not in HABITATS[species]:
        raise ValueError(f"Habitat {habitat} not available for {SPECIES[species]}")
    parts = sep.split("-")
    aas, qb, qs = parts[:5], parts[5] == "1", parts[6] == "1"
    db = conn()
    try:
        rule = _rule(db, species, habitat, sep)
        flags: list[str] = []
        out = {"septatype": sep, "species": species, "habitat": habitat, "mean_h": grammar.mean_h(aas[:4]),
               "mechanism": grammar.mechanism(aas), "novel": False}
        g = grammar.predict(species, aas, qb, qs, recn_codon)
        out["grammar"] = {"call": g["call"], "rule": g["rule"], "votes": g["votes"]}
        out["criteria"] = grammar.explain(species, aas, g, None, (qb, qs), recn_codon)
        pooled = None if rule or habitat == "ALL" else _rule(db, species, "ALL", sep)
        use = rule or pooled
        if use:
            n, r, s = use["n_isolates"], use["n_resistant"], use["n_susceptible"]
            out.update(call=_label(use["predicted"]), tier="database", known=True, n_isolates=n, n_resistant=r,
                       n_susceptible=s, confidence=round(use["confidence"], 3))
            out["reference_sts"] = reference_sts(db, species, "ALL" if pooled else habitat, sep)
            out["reference_sites"] = reference_sites(db, species, "ALL" if pooled else habitat, sep)
            out["reference_countries"] = reference_countries(db, species, "ALL" if pooled else habitat, sep)
            if pooled:
                out["tier"] = "database_other_site"
                flags.append("Not seen for this isolation source, so the call uses isolates of this species from all sources.")
            if n < LOW_SUPPORT:
                flags.append(f"Only {n} reference isolate(s) carry this septatype, so treat the call with care.")
            if use["confidence"] < MIXED:
                flags.append(f"The reference isolates with this septatype are mixed: {r} resistant, {s} susceptible.")
            if r == s:
                flags.append("The reference isolates are split evenly, so the paper's rule calls it susceptible.")
            if g["call"] != out["call"]:
                flags.append(f"The grammar alone would call this {g['call']}. The database majority is used.")
        else:
            out.update(call=g["call"], tier="grammar", known=False, novel=True, n_isolates=0, n_resistant=0,
                       n_susceptible=0, confidence=None, reference_sts=[], reference_sites=[], reference_countries=[])
            flags.append("New septatype: not seen in the reference isolates for this species, "
                         "so the call comes from the grammar and is less certain. It is flagged for curation.")
            rows = db.execute("SELECT septatype,n_isolates,n_resistant,n_susceptible,predicted FROM septatypes "
                              "WHERE species=? AND habitat=? AND tetratype=? ORDER BY n_isolates DESC",
                              (species, "ALL", "-".join(parts[:4]))).fetchall()
            out["similar"] = [dict(x) for x in rows]
        if habitat == "ALL":
            calls = {}
            for h in HABITATS[species]:
                if h != "ALL":
                    rr = _rule(db, species, h, sep)
                    if rr:
                        calls[h] = _label(rr["predicted"])
            out["by_habitat"] = calls
            if len(set(calls.values())) > 1:
                flags.append("The call differs by source (" + ", ".join(f"{k}: {v}" for k, v in calls.items())
                             + "). Pick the isolate's source for a firmer answer.")
        out["flags"] = flags
        return out
    finally:
        db.close()


# ---------------------------------------------------------------- statistics
def accuracy_stats() -> dict:
    db = conn()
    try:
        rows = [dict(r) for r in db.execute("SELECT * FROM species_stats ORDER BY species, habitat")]
        meta = {r["key"]: r["value"] for r in db.execute("SELECT * FROM meta")}
    finally:
        db.close()
    return {"rows": rows, "meta": meta}


def edgecases() -> dict:
    """Resistance among isolates with the wild-type S-D-S-E tetratype, split by qnr and species."""
    db = conn()
    try:
        q = db.execute("""SELECT species, qnrB2, qnrS1, COUNT(*) n,
                          SUM(lab_phenotype='Resistant') res FROM isolates
                          WHERE tetratype='S-D-S-E' AND source='founding' GROUP BY species, qnrB2, qnrS1
                          ORDER BY species, qnrB2, qnrS1""").fetchall()
        s = db.execute("""SELECT habitat, septatype, COUNT(*) n, SUM(lab_phenotype='Resistant') res FROM isolates
                          WHERE species='senterica' AND tetratype='F-D-S-E' AND source='founding'
                          GROUP BY habitat, septatype ORDER BY n DESC""").fetchall()
    finally:
        db.close()
    return {"wildtype_tetratype_by_qnr": [dict(x) for x in q], "senterica_FDSE": [dict(x) for x in s]}


def benchmark() -> dict:
    db = conn()
    try:
        loaded = db.execute("SELECT value FROM meta WHERE key='bvbrc_loaded'").fetchone()
        n_all = db.execute("SELECT COUNT(*) FROM isolates WHERE source='founding'").fetchone()[0]
        sc = db.execute("SELECT SUM(septatyper_call=lab_phenotype) FROM isolates WHERE source='founding'").fetchone()[0]
        res = {"loaded": bool(loaded and loaded[0] == "1"), "n_isolates": n_all,
               "septatyper_correct": sc, "septatyper_accuracy": round(100 * sc / n_all, 1)}
        if res["loaded"]:
            r = db.execute("""SELECT COUNT(*) n, SUM(septatyper_call=lab_phenotype) s_ok, SUM(bvbrc_call=lab_phenotype) b_ok,
                              SUM(septatyper_call=lab_phenotype AND bvbrc_call<>lab_phenotype) wins,
                              SUM(septatyper_call<>lab_phenotype AND bvbrc_call=lab_phenotype) losses
                              FROM isolates WHERE bvbrc_call IS NOT NULL AND source='founding'""").fetchone()
            res.update(compared=r["n"], septatyper_ok=r["s_ok"], bvbrc_ok=r["b_ok"], wins=r["wins"], losses=r["losses"])
            res["examples"] = [dict(x) for x in db.execute(
                """SELECT genome_id, isolate_name, species, septatype, lab_phenotype, bvbrc_call, septatyper_call
                   FROM isolates WHERE septatyper_call=lab_phenotype AND bvbrc_call<>lab_phenotype
                   AND bvbrc_call IS NOT NULL LIMIT 25""")]
    finally:
        db.close()
    return res


def geography(species="", st="", tetratype="", phenotype="", missing=False, septatype=""):
    """Counts per country. With missing=True, returns how many matching isolates have no country."""
    where, args = ["iso3 IS NULL" if missing else "iso3 IS NOT NULL"], []
    for col, val in (("species", species), ("st", st), ("tetratype", tetratype), ("septatype", septatype), ("lab_phenotype", phenotype)):
        if val:
            where.append(f"{col}=?")
            args.append(val)
    db = conn()
    try:
        if missing:
            return db.execute(f"SELECT COUNT(*) FROM isolates WHERE {' AND '.join(where)}", args).fetchone()[0]
        rows = db.execute(f"""SELECT iso3, country, COUNT(*) n, SUM(lab_phenotype='Resistant') res FROM isolates
                              WHERE {' AND '.join(where)} GROUP BY iso3, country ORDER BY n DESC""", args).fetchall()
    finally:
        db.close()
    return [dict(r) for r in rows]


def septatype_freq(species: str, top: int = 12) -> dict:
    db = conn()
    try:
        rows = db.execute("""SELECT septatype, habitat, COUNT(*) n, SUM(lab_phenotype='Resistant') res FROM isolates
                             WHERE species=? GROUP BY septatype, habitat""", (species,)).fetchall()
    finally:
        db.close()
    tot, by = {}, {}
    for r in rows:
        tot[r["septatype"]] = tot.get(r["septatype"], 0) + r["n"]
        by.setdefault(r["septatype"], {})[r["habitat"]] = r["n"]
    keys = sorted(tot, key=lambda k: -tot[k])[:top]
    return {"species": species, "total": sum(tot.values()),
            "rows": [{"septatype": k, "total": tot[k], "by_habitat": by[k]} for k in keys]}


def facets() -> dict:
    """Every value with its isolate count, most frequent first, for the filter lists."""
    db = conn()
    try:
        q = lambda sql: [[r[0], r[1]] for r in db.execute(sql)]  # noqa: E731
        f = {"st": q("SELECT st, COUNT(*) FROM isolates WHERE st IS NOT NULL GROUP BY st ORDER BY COUNT(*) DESC, st"),
             "country": q("SELECT country, COUNT(*) FROM isolates WHERE country IS NOT NULL GROUP BY country ORDER BY COUNT(*) DESC, country"),
             "tetratype": q("SELECT tetratype, COUNT(*) FROM isolates GROUP BY tetratype ORDER BY COUNT(*) DESC, tetratype"),
             "septatype": q("SELECT septatype, COUNT(*) FROM isolates GROUP BY septatype ORDER BY COUNT(*) DESC, septatype"),
             "habitat": q("SELECT habitat, COUNT(*) FROM isolates GROUP BY habitat ORDER BY COUNT(*) DESC, habitat"),
             "year": [r[0] for r in db.execute("SELECT MIN(year_start), MAX(COALESCE(year_end,year_start)) FROM isolates")]}
    finally:
        db.close()
    return f


def _row_h(d: dict):
    """QRDR helix hydrophobicity of an isolate, from its four tetratype residues: GyrA helix mean, ParC helix mean, combined."""
    try:
        a4 = [d["gyrA_83"], d["gyrA_87"], d["parC_80"], d["parC_84"]]
        hg, hp = grammar._helix_means(a4)
    except (KeyError, TypeError):
        return None, None, None
    return round(hg, 3), round(hp, 3), grammar.mean_h(a4)


SORTS = {"genome_id": "genome_id", "species": "species, habitat", "st": "st", "year": "year_start",
         "country": "country", "septatype": "septatype", "lab": "lab_phenotype"}


def isolates_query(f: dict, sort="genome_id", desc=False, page=1, per_page=25, limit=True):
    where, args = [], []
    for col in ("species", "habitat", "st", "country", "septatype", "tetratype"):
        if f.get(col):
            where.append(f"{col}=?")
            args.append(f[col])
    if f.get("phenotype"):
        where.append("lab_phenotype=?")
        args.append(f["phenotype"])
    if f.get("year_from"):
        where.append("COALESCE(year_end, year_start) >= ?")
        args.append(int(f["year_from"]))
    if f.get("year_to"):
        where.append("COALESCE(year_start, year_end) <= ?")
        args.append(int(f["year_to"]))
    if f.get("wins"):
        where.append("septatyper_call=lab_phenotype AND bvbrc_call IS NOT NULL AND bvbrc_call<>lab_phenotype")
    if f.get("source"):
        where.append("source=?")
        args.append(f["source"])
    if f.get("q"):
        where.append("(genome_id LIKE ? OR isolate_name LIKE ? OR genome_name LIKE ?)")
        args += [f"%{f['q']}%"] * 3
    w = (" WHERE " + " AND ".join(where)) if where else ""
    order = SORTS.get(sort, "genome_id") + (" DESC" if desc else "")
    db = conn()
    try:
        total = db.execute("SELECT COUNT(*) FROM isolates" + w, args).fetchone()[0]
        sql = f"SELECT {','.join(ISO_COLS)} FROM isolates{w} ORDER BY {order}, genome_id"
        if limit:
            rows = db.execute(sql + " LIMIT ? OFFSET ?", args + [per_page, (page - 1) * per_page]).fetchall()
        else:
            rows = db.execute(sql, args).fetchall()
    finally:
        db.close()
    out = []
    for r in rows:
        d = dict(r)
        d["h_gyr"], d["h_parc"], d["h_qrdr"] = _row_h(d)
        d["star"] = bool(d["bvbrc_call"] and d["septatyper_call"] == d["lab_phenotype"] and d["bvbrc_call"] != d["lab_phenotype"])
        d["bvbrc_url"] = f"https://www.bv-brc.org/view/Genome/{d['genome_id']}" if not d["genome_id"].startswith("SUB-") else None
        out.append(d)
    return {"total": total, "rows": out}


# ------------------------------------------------------------------ curation
def create_submission(d: dict, result: dict) -> dict:
    pred = result.get("prediction") or {}
    call = pred.get("call")
    novel = bool(pred.get("novel"))
    lab = d.get("lab_phenotype") or None
    matches = None if not (lab and call) else int(lab == call)
    status = "needs_review" if (novel or matches == 0 or not call) else "pending"
    job = uuid.uuid4().hex[:12]
    db = conn()
    try:
        db.execute("INSERT INTO submissions VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                   (job, _now(), status, d["species"], d["habitat"], d.get("genome_id"), d.get("genome_name"),
                    d.get("country"), d.get("year"), d.get("source_text"), d.get("mic"), lab, d.get("email"),
                    result.get("septatype"), pred.get("tier"), call, int(novel), matches,
                    json.dumps(result), None, None))
        db.commit()
    finally:
        db.close()
    _persist()
    return {"job_id": job, "status": status, "novel": novel, "matches": matches}


def list_queue(status: str = "") -> list[dict]:
    db = conn()
    try:
        sql = ("SELECT job_id,created,status,species,habitat,genome_id,genome_name,country,year,source_text,mic,"
               "lab_phenotype,septatype,tier,call,novel,matches,curator_note,decided,result_json FROM submissions")
        rows = db.execute(sql + (" WHERE status=?" if status else "") + " ORDER BY created DESC",
                          (status,) if status else ()).fetchall()
    finally:
        db.close()
    out = []
    for r in rows:
        d = dict(r)
        try:
            d["st"] = ((json.loads(d.pop("result_json") or "{}").get("analysis") or {}).get("mlst") or {}).get("st")
        except (TypeError, ValueError):
            d["st"] = None
        out.append(d)
    return out


def queue_counts() -> dict:
    db = conn()
    try:
        c = {r[0]: r[1] for r in db.execute("SELECT status, COUNT(*) FROM submissions GROUP BY status")}
    finally:
        db.close()
    return {k: c.get(k, 0) for k in ("needs_review", "pending", "approved", "rejected")} | {"all": sum(c.values())}


def _revoke(db, s) -> None:
    """Undo an approval: take the isolate out of the database and its vote out of the rules."""
    gid = s["genome_id"] or f"SUB-{s['job_id']}"
    db.execute("DELETE FROM isolates WHERE genome_id=? AND source='community'", (gid,))
    for hab in {s["habitat"], "ALL"}:
        r = db.execute("SELECT * FROM septatypes WHERE species=? AND habitat=? AND septatype=?",
                       (s["species"], hab, s["septatype"])).fetchone()
        if not r:
            continue
        res_ = r["n_resistant"] - (1 if s["lab_phenotype"] == "Resistant" else 0)
        sus_ = r["n_susceptible"] - (1 if s["lab_phenotype"] == "Susceptible" else 0)
        if res_ + sus_ <= 0:
            db.execute("DELETE FROM septatypes WHERE species=? AND habitat=? AND septatype=?", (s["species"], hab, s["septatype"]))
        else:
            db.execute("UPDATE septatypes SET n_isolates=?, n_resistant=?, n_susceptible=?, predicted=?, confidence=? "
                       "WHERE species=? AND habitat=? AND septatype=?",
                       (res_ + sus_, res_, sus_, "R" if res_ > sus_ else "S", max(res_, sus_) / (res_ + sus_),
                        s["species"], hab, s["septatype"]))


def apply_approval(db, s) -> str:
    """Add an approved submission's isolate to the repository and its vote to the rules. Returns the genome ID."""
    aas = s["septatype"].split("-")
    gid = s["genome_id"] or f"SUB-{s['job_id']}"
    if db.execute("SELECT 1 FROM isolates WHERE genome_id=?", (gid,)).fetchone():
        raise ValueError(f"Genome ID {gid} is already in the database.")
    country = countries.canonical(s["country"]) or s["country"]
    iso3 = countries.iso3(country)
    if country and not iso3:
        r = db.execute("SELECT iso3 FROM isolates WHERE country=? AND iso3 IS NOT NULL LIMIT 1", (country,)).fetchone()
        iso3 = r[0] if r else None
    ys = int(s["year"]) if s["year"] and str(s["year"]).isdigit() else None
    try:
        mst = ((json.loads(s["result_json"]).get("analysis") or {}).get("mlst") or {}).get("st")
    except (TypeError, ValueError):
        mst = None
    db.execute(f"INSERT INTO isolates ({','.join(ISO_COLS)}) VALUES ({','.join('?'*len(ISO_COLS))})",
               (gid, s["genome_name"], s["genome_name"] or gid, s["species"], s["habitat"], s["lab_phenotype"],
                s["year"], ys, ys, country, iso3, mst, aas[0], aas[1], aas[2], aas[3], aas[4], int(aas[5]),
                int(aas[6]), s["septatype"], "-".join(aas[:4]), s["call"], None, None, "community", _now()))
    # the approved isolate joins the rules table, so a new septatype becomes a known one
    for hab in {s["habitat"], "ALL"}:
        r = db.execute("SELECT * FROM septatypes WHERE species=? AND habitat=? AND septatype=?",
                       (s["species"], hab, s["septatype"])).fetchone()
        res_, sus_ = (r["n_resistant"], r["n_susceptible"]) if r else (0, 0)
        if s["lab_phenotype"] == "Resistant":
            res_ += 1
        else:
            sus_ += 1
        db.execute("INSERT OR REPLACE INTO septatypes VALUES (?,?,?,?,?,?,?,?,?,?)",
                   (s["species"], hab, s["septatype"], "-".join(aas[:4]), res_ + sus_, res_, sus_,
                    "R" if res_ > sus_ else "S", max(res_, sus_) / (res_ + sus_), grammar.mean_h(aas[:4])))
    return gid


def decide(job_id: str, action: str, note: str = "") -> dict:
    """approve: pending, needs review or rejected -> approved. reject: pending, needs review or approved -> rejected.
    Rejecting an approved submission removes its isolate again and takes its vote out of the rules."""
    db = conn()
    try:
        s = db.execute("SELECT * FROM submissions WHERE job_id=?", (job_id,)).fetchone()
        if not s:
            raise KeyError("Unknown submission")
        status = s["status"]
        if action == "reject":
            if status == "rejected":
                raise ValueError("Already rejected")
            if status == "approved":
                _revoke(db, s)
            db.execute("UPDATE submissions SET status='rejected', curator_note=?, decided=? WHERE job_id=?",
                       (note or ("Rejected after approval" if status == "approved" else ""), _now(), job_id))
            db.commit()
            _persist()
            return {"status": "rejected", "revoked": status == "approved"}
        if action != "approve":
            raise ValueError("action must be approve or reject")
        if status == "approved":
            raise ValueError("Already approved")
        if s["lab_phenotype"] not in ("Resistant", "Susceptible") or not s["septatype"]:
            raise ValueError("Cannot approve without a laboratory phenotype and a complete septatype.")
        gid = apply_approval(db, s)
        db.execute("UPDATE submissions SET status='approved', curator_note=?, decided=? WHERE job_id=?",
                   (note, _now(), job_id))
        db.commit()
        _persist()
        return {"status": "approved", "genome_id": gid}
    finally:
        db.close()


def cutoff_stats() -> dict:
    """How well the QRDR helix mean hydrophobicity separates resistant from susceptible isolates, from the founding data."""
    db = conn()
    try:
        rows = db.execute("""SELECT tetratype, SUM(lab_phenotype='Resistant') r, SUM(lab_phenotype='Susceptible') s
                             FROM isolates WHERE source='founding' GROUP BY tetratype""").fetchall()
    finally:
        db.close()
    groups = []
    for t in rows:
        a = t["tetratype"].split("-")
        groups.append((grammar.mean_h(a), sum(1 for v in grammar.votes(*a) if v["vote"]), t["r"], t["s"]))
    n_r, n_s = sum(g[2] for g in groups), sum(g[3] for g in groups)
    cut = grammar.H_CUTOFF
    grid = sorted({round(0.20 + 0.02 * i, 2) for i in range(16)} | {cut, grammar.H_WILD})

    def table(sub):
        tot = sum(g[2] + g[3] for g in sub)
        out = []
        for c in grid:
            ok = sum((g[2] if g[0] > c else g[3]) for g in sub)
            out.append({"cutoff": c, "accuracy": round(100 * ok / tot, 1) if tot else None})
        return out, tot

    allt, tot_all = table(groups)
    one = [g for g in groups if g[1] == 1]
    onet, tot_one = table(one)
    best = lambda t: max(t, key=lambda x: (x["accuracy"], -abs(x["cutoff"] - cut)))  # noqa: E731
    hs = [g[0] for g in groups]
    # the whole grammar (votes, hydrophobicity, overrides) scored over all founding isolates at each cutoff, as in the
    # stored grammar accuracy: this is the figure that matters for choosing a cutoff
    db = conn()
    try:
        sep_rows = db.execute("""SELECT species, septatype, SUM(lab_phenotype='Resistant') r, SUM(lab_phenotype='Susceptible') s
                                 FROM isolates WHERE source='founding' GROUP BY species, septatype""").fetchall()
    finally:
        db.close()
    gram = []
    for c in grid:
        ok = 0
        for t in sep_rows:
            a = t["septatype"].split("-")
            call = grammar.predict(t["species"], a[:5], a[5] == "1", a[6] == "1", None, c)["call"]
            ok += t["r"] if call == "Resistant" else t["s"]
        gram.append({"cutoff": c, "accuracy": round(100 * ok / (n_r + n_s), 1)})
    return {"n_isolates": n_r + n_s, "cutoff": cut, "wild_type": grammar.h_block(list(grammar.WILD4)), "resistant_reference": grammar.h_block(list(grammar.RES4)),
            "helix": {"gyrA": {"start": grammar.GYR_START, "sequence": grammar.GYR_HELIX},
                      "parC": {"start": grammar.PARC_START, "sequence": grammar.PARC_HELIX}},
            "range": [min(hs), max(hs)] if hs else None,
            "mean_resistant": round(sum(g[0] * g[2] for g in groups) / n_r, 3),
            "mean_susceptible": round(sum(g[0] * g[3] for g in groups) / n_s, 3),
            "all": allt, "best_all": best(allt), "n_single_vote": tot_one,
            "single_vote": onet, "best_single_vote": best(onet) if one else None,
            "grammar_all": gram, "best_grammar_all": best(gram),
            "single_vote_resistant": sum(g[2] for g in one)}
