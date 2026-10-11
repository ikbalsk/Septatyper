"""Database access, prediction (database lookup first, grammar fallback second), statistics, curation."""
from __future__ import annotations

import datetime
import json
import shutil
import sqlite3
import uuid

from . import grammar
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


def _now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def _label(p: str) -> str:
    return "Resistant" if p == "R" else "Susceptible"


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
        pooled = None if rule or habitat == "ALL" else _rule(db, species, "ALL", sep)
        use = rule or pooled
        if use:
            n, r, s = use["n_isolates"], use["n_resistant"], use["n_susceptible"]
            out.update(call=_label(use["predicted"]), tier="database", known=True, n_isolates=n, n_resistant=r,
                       n_susceptible=s, confidence=round(use["confidence"], 3))
            if pooled:
                out["tier"] = "database_other_site"
                flags.append("Not seen at this isolation site. The call uses isolates of this species from all sites.")
            if n < LOW_SUPPORT:
                flags.append(f"Low support: only {n} reference isolate(s) carry this septatype.")
            if use["confidence"] < MIXED:
                flags.append(f"Mixed phenotype: {r} resistant and {s} susceptible among reference isolates.")
            if r == s:
                flags.append("Tie in the reference data (called Susceptible by the paper's rule).")
            if g["call"] != out["call"]:
                flags.append(f"The physicochemical grammar would call this {g['call']}. The database majority is used.")
        else:
            out.update(call=g["call"], tier="grammar", known=False, novel=True, n_isolates=0, n_resistant=0,
                       n_susceptible=0, confidence=None)
            flags.append("New septatype: never seen in the reference isolates for this species. "
                         "The call comes from the physicochemical grammar and is less certain. It is flagged for curation.")
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
                flags.append("The call differs between sites (" + ", ".join(f"{k}: {v}" for k, v in calls.items())
                             + "). Choose the isolate's site for a reliable result.")
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


def geography(species="", st="", tetratype="", phenotype="") -> list[dict]:
    where, args = ["iso3 IS NOT NULL"], []
    for col, val in (("species", species), ("st", st), ("tetratype", tetratype), ("lab_phenotype", phenotype)):
        if val:
            where.append(f"{col}=?")
            args.append(val)
    db = conn()
    try:
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
    db = conn()
    try:
        f = {"st": [r[0] for r in db.execute("SELECT st FROM isolates WHERE st IS NOT NULL GROUP BY st ORDER BY COUNT(*) DESC")],
             "country": [r[0] for r in db.execute("SELECT country FROM isolates WHERE country IS NOT NULL GROUP BY country ORDER BY COUNT(*) DESC")],
             "tetratype": [r[0] for r in db.execute("SELECT tetratype FROM isolates GROUP BY tetratype ORDER BY COUNT(*) DESC")],
             "year": [r[0] for r in db.execute("SELECT MIN(year_start), MAX(COALESCE(year_end,year_start)) FROM isolates")]}
    finally:
        db.close()
    return f


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
    return {"job_id": job, "status": status, "novel": novel, "matches": matches}


def list_queue(status: str = "") -> list[dict]:
    db = conn()
    try:
        sql = ("SELECT job_id,created,status,species,habitat,genome_id,genome_name,country,year,source_text,mic,"
               "lab_phenotype,septatype,tier,call,novel,matches,curator_note,decided FROM submissions")
        rows = db.execute(sql + (" WHERE status=?" if status else "") + " ORDER BY created DESC",
                          (status,) if status else ()).fetchall()
    finally:
        db.close()
    return [dict(r) for r in rows]


def decide(job_id: str, action: str, note: str = "") -> dict:
    db = conn()
    try:
        s = db.execute("SELECT * FROM submissions WHERE job_id=?", (job_id,)).fetchone()
        if not s:
            raise KeyError("Unknown submission")
        if s["status"] in ("approved", "rejected"):
            raise ValueError(f"Already {s['status']}")
        if action == "reject":
            db.execute("UPDATE submissions SET status='rejected', curator_note=?, decided=? WHERE job_id=?",
                       (note, _now(), job_id))
            db.commit()
            return {"status": "rejected"}
        if action != "approve":
            raise ValueError("action must be approve or reject")
        if s["lab_phenotype"] not in ("Resistant", "Susceptible") or not s["septatype"]:
            raise ValueError("Cannot approve without a laboratory phenotype and a complete septatype.")
        aas = s["septatype"].split("-")
        gid = s["genome_id"] or f"SUB-{job_id}"
        if db.execute("SELECT 1 FROM isolates WHERE genome_id=?", (gid,)).fetchone():
            raise ValueError(f"Genome ID {gid} is already in the database.")
        country = s["country"]
        iso3 = None
        if country:
            r = db.execute("SELECT iso3 FROM isolates WHERE country=? AND iso3 IS NOT NULL LIMIT 1", (country,)).fetchone()
            iso3 = r[0] if r else None
        ys = int(s["year"]) if s["year"] and str(s["year"]).isdigit() else None
        db.execute(f"INSERT INTO isolates ({','.join(ISO_COLS)}) VALUES ({','.join('?'*len(ISO_COLS))})",
                   (gid, s["genome_name"], s["genome_name"] or gid, s["species"], s["habitat"], s["lab_phenotype"],
                    s["year"], ys, ys, country, iso3, None, aas[0], aas[1], aas[2], aas[3], aas[4], int(aas[5]),
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
            pred = "R" if res_ > sus_ else "S"
            db.execute("INSERT OR REPLACE INTO septatypes VALUES (?,?,?,?,?,?,?,?,?,?)",
                       (s["species"], hab, s["septatype"], "-".join(aas[:4]), res_ + sus_, res_, sus_, pred,
                        max(res_, sus_) / (res_ + sus_), grammar.mean_h(aas[:4])))
        db.execute("UPDATE submissions SET status='approved', curator_note=?, decided=? WHERE job_id=?",
                   (note, _now(), job_id))
        db.commit()
        return {"status": "approved", "genome_id": gid}
    finally:
        db.close()
