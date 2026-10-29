"""Ciprotyper API (v1) and web interface.

Modules: A predict | B browse + benchmark | C community submission + curation | D theory and validation.
All SQL is parameterised. Uploaded sequences are analysed in a temporary folder and deleted.
"""
from __future__ import annotations

import csv
import hmac
import io
import json
import os
import re
import tempfile
import threading
import urllib.parse
import urllib.request
from pathlib import Path

from fastapi import Body, Depends, FastAPI, File, Form, Header, HTTPException, Query, UploadFile
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import bvbrc, countries, grammar, persist, rules, septatype
from .config import ADMIN_TOKEN, HABITATS, MAX_UPLOAD_MB, SOURCES, SPECIES

app = FastAPI(title="Ciprotyper API", version="2.0",
              description="Ciprofloxacin resistance prediction and isolate repository for Enterobacterales.")
STATIC = Path(__file__).parent / "static"
AA = set("ACDEFGHIKLMNPQRSTVWY")
BLAST_SLOTS = threading.BoundedSemaphore(2)       # free instances are small: at most 2 analyses at once
GENOME_ID = re.compile(r"^\d{1,8}\.\d{1,6}$")


# ------------------------------------------------------------------ helpers
def _check(species: str, habitat: str):
    if species not in SPECIES:
        raise HTTPException(400, "Unknown species.")
    if habitat not in HABITATS[species]:
        raise HTTPException(400, f"Isolation source '{habitat}' is not available.")


_CACHE: dict = {}          # (sha256, species, qnr overrides) -> analysis; the same file is never analysed twice
_CACHE_MAX = 32


def _analyse_text(text: str, species: str, qnr: dict, name: str = "pasted"):
    """Write FASTA text or bytes to a temp file, run the pipeline, delete the file."""
    import hashlib
    key = (hashlib.sha256(text.encode()).hexdigest(), species, qnr.get("qnrB2"), qnr.get("qnrS1"))
    if key in _CACHE:
        return _CACHE[key]
    if not BLAST_SLOTS.acquire(blocking=False):
        raise HTTPException(429, "The server is busy with other analyses. Try again in a minute.")
    try:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "input.fa"
            path.write_text(text)
            res = septatype.analyse(path, species, qnr)
            if len(_CACHE) >= _CACHE_MAX:
                _CACHE.pop(next(iter(_CACHE)))
            _CACHE[key] = res
            return res
    except septatype.InputError as e:
        raise HTTPException(400, str(e))
    except septatype.RefsMissing as e:
        raise HTTPException(503, str(e))
    except RuntimeError as e:
        raise HTTPException(500, str(e))
    finally:
        BLAST_SLOTS.release()


def _read_upload(file: UploadFile) -> str | bytes:
    limit = MAX_UPLOAD_MB * 1024 * 1024
    buf, size = io.BytesIO(), 0
    while chunk := file.file.read(1 << 20):
        size += len(chunk)
        if size > limit:
            raise HTTPException(413, f"File is larger than {MAX_UPLOAD_MB} MB.")
        buf.write(chunk)
    data = buf.getvalue()
    if (file.filename or "").endswith(".gz") or data[:2] == b"\x1f\x8b":
        import gzip
        try:
            data = gzip.decompress(data)
        except OSError:
            raise HTTPException(400, "Could not decompress the .gz file.")
        if len(data) > 4 * limit:
            raise HTTPException(413, "Decompressed file is too large.")
    try:
        return data.decode()
    except UnicodeDecodeError:
        raise HTTPException(400, "The file is not plain-text FASTA.")


def _result(species, habitat, analysis, sep):
    pred = rules.predict(species, habitat, sep, (analysis or {}).get("recN_codon")) if sep else None
    if pred and pred.get("criteria"):
        pos = (analysis or {}).get("positions") or {}
        ref4 = [(pos.get(k) or {}).get("ref_aa") for k in ("gyrA83", "gyrA87", "parC80", "parC84")]
        ref4 = ref4 if all(ref4) else list("SDSE")      # manual entry has no sample genome: compare with the wild type
        pred["criteria"]["reference"] = grammar.h_block(ref4)
    warnings = (analysis or {}).get("warnings", []) + ((pred or {}).get("flags", []))
    return {"status": "ok" if sep else "incomplete", "species": species, "species_name": SPECIES[species],
            "habitat": habitat, "septatype": sep, "analysis": analysis, "prediction": pred, "warnings": warnings}


def admin(x_admin_token: str = Header(default="")):
    if not ADMIN_TOKEN:
        raise HTTPException(503, "Curation is disabled: set the ADMIN_TOKEN environment variable on the server.")
    if not hmac.compare_digest(x_admin_token.strip().encode("utf-8"), ADMIN_TOKEN.encode("utf-8")):
        raise HTTPException(401, "Wrong curator token. Check it matches ADMIN_TOKEN in Render exactly.")


@app.on_event("startup")
def _startup():
    rules.ensure_db()
    if os.environ.get("BVBRC_AUTOLOAD", "1") != "0" and not bvbrc.is_loaded():
        bvbrc.start()
    if persist.configured():
        threading.Thread(target=persist.restore_if_configured, daemon=True).start()


@app.get("/api/v1/bvbrc-status")
def bvbrc_status():
    return bvbrc.status()


@app.post("/api/v1/admin/bvbrc/load", dependencies=[Depends(admin)])
def bvbrc_load():
    started = bvbrc.start() if not bvbrc.STATE["running"] else None
    return bvbrc.status() | {"started": not bvbrc.STATE["running"] or started is None}


# --------------------------------------------------------------------- meta
@app.get("/api/v1/meta")
def meta():
    s = rules.accuracy_stats()
    db = rules.conn()
    try:
        n = db.execute("SELECT COUNT(*) FROM isolates").fetchone()[0]
        nsp = db.execute("SELECT COUNT(DISTINCT species) FROM isolates").fetchone()[0]
        nsep = db.execute("SELECT COUNT(DISTINCT septatype) FROM isolates").fetchone()[0]
        ncountry = db.execute("SELECT COUNT(DISTINCT country) FROM isolates").fetchone()[0]
        ncomm = db.execute("SELECT COUNT(*) FROM isolates WHERE source='community'").fetchone()[0]
    finally:
        db.close()
    return {"species": SPECIES, "habitats": HABITATS, "stats": s["rows"], "build": s["meta"],
            "sources": [{"key": k, "label": l, "short": sh} for k, l, sh in SOURCES], "countries": sorted(countries.COUNTRIES),
            "counts": {"genomes": n, "species": nsp, "septatypes": nsep, "countries": ncountry, "community": ncomm},
            "references": septatype.refs_status(), "max_upload_mb": MAX_UPLOAD_MB,
            "curation_enabled": bool(ADMIN_TOKEN), "facets": rules.facets()}


# ------------------------------------------------------------------ predict
@app.post("/api/v1/predict")
def predict_upload(file: UploadFile | None = File(None), sequence: str = Form(""), species: str = Form(...),
                   habitat: str = Form("ALL"), qnrB2: str = Form("auto"), qnrS1: str = Form("auto")):
    _check(species, habitat)
    if file is not None and file.filename:
        text = _read_upload(file)
    elif sequence.strip():
        text = sequence if sequence.lstrip().startswith(">") else ">pasted\n" + sequence
        if len(text) > MAX_UPLOAD_MB * 1024 * 1024:
            raise HTTPException(413, "Pasted text is too large.")
    else:
        raise HTTPException(400, "Provide a FASTA file or paste sequence.")
    analysis = _analyse_text(text, species, {"qnrB2": qnrB2, "qnrS1": qnrS1})
    return _result(species, habitat, analysis, analysis["septatype"])


class BvbrcReq(BaseModel):
    genome_id: str
    species: str = "auto"            # "auto" reads the species from the BV-BRC record
    fallback_species: str = ""       # used only when auto cannot decide
    habitat: str = "ALL"


@app.post("/api/v1/predict/bvbrc")
def predict_bvbrc(m: BvbrcReq):
    """Predict from a BV-BRC genome ID: the assembly is fetched from BV-BRC, analysed and discarded. One genome per call."""
    gid = m.genome_id.strip()
    if not bvbrc.valid_id(gid):
        raise HTTPException(400, f"'{gid[:30]}' is not a BV-BRC genome ID. An ID looks like 562.7710.")
    try:
        meta = bvbrc.genome_meta(gid)
    except Exception as e:
        raise HTTPException(502, f"Could not reach BV-BRC ({type(e).__name__}). Try again in a minute.")
    if not meta:
        raise HTTPException(404, f"BV-BRC has no genome with ID {gid}.")
    name, detected, notes = meta.get("genome_name") or gid, bvbrc.detect_species(meta), []
    species = m.species
    if species == "auto":
        named = bool((meta.get("species") or "").strip() or (meta.get("genome_name") or "").strip())
        if detected:
            species = detected
        elif named:      # the record names another organism: never guess
            raise HTTPException(400, f"BV-BRC calls this genome '{name}', which is not one of the four species Ciprotyper covers.")
        else:
            species = m.fallback_species
            if species not in SPECIES:
                raise HTTPException(400, "BV-BRC does not name the species of this genome. Choose the species and try again.")
            notes.append(f"BV-BRC did not name the species, so {SPECIES[species]} (your selection) was used.")
    elif species in SPECIES and detected and detected != species:
        notes.append(f"BV-BRC calls this genome '{name}', but {SPECIES[species]} was chosen.")
    if species not in SPECIES:
        raise HTTPException(400, "Unknown species.")
    habitat = m.habitat
    if habitat not in HABITATS[species]:
        notes.append(f"Isolation source '{habitat}' is not available for {SPECIES[species]}, so all sources were used.")
        habitat = "ALL"
    try:
        text = bvbrc.genome_fasta(gid, MAX_UPLOAD_MB * 1024 * 1024)
    except bvbrc.TooLarge as e:
        raise HTTPException(413, str(e))
    except Exception as e:
        raise HTTPException(502, f"Could not download the assembly from BV-BRC ({type(e).__name__}: {str(e)[:120]}).")
    analysis = _analyse_text(text, species, {"qnrB2": "auto", "qnrS1": "auto"})
    res = _result(species, habitat, analysis, analysis["septatype"])
    db = rules.conn()
    try:
        ref = db.execute("SELECT lab_phenotype FROM isolates WHERE genome_id=? AND source='founding'", (gid,)).fetchone()
    finally:
        db.close()
    if ref:
        notes.append(f"This genome is one of the reference isolates (laboratory phenotype {ref[0]}), so its call is not an independent test.")
    res["bvbrc"] = {"genome_id": gid, "genome_name": name, "species_detected": detected, "isolation_source": meta.get("isolation_source"),
                    "year": meta.get("collection_year"), "country": meta.get("isolation_country"),
                    "url": f"https://www.bv-brc.org/view/Genome/{gid}", "in_reference": bool(ref), "lab_phenotype": ref[0] if ref else None}
    res["warnings"] = list(res.get("warnings") or []) + notes
    return res


class Manual(BaseModel):
    species: str
    habitat: str = "ALL"
    gyrA_83: str
    gyrA_87: str
    parC_80: str
    parC_84: str
    recN_346: str
    qnrB2: bool = False
    qnrS1: bool = False
    recN_codon: str = ""


@app.post("/api/v1/predict/manual")
def predict_manual(m: Manual):
    _check(m.species, m.habitat)
    aas = [x.strip().upper() for x in (m.gyrA_83, m.gyrA_87, m.parC_80, m.parC_84, m.recN_346)]
    if any(len(a) != 1 or a not in AA for a in aas):
        raise HTTPException(400, "Each position needs one amino-acid letter (for example S, D, I, E, L).")
    codon = m.recN_codon.strip().upper()
    if codon and not re.fullmatch(r"[ACGT]{3}", codon):
        raise HTTPException(400, "RecN codon must be three letters from A, C, G, T.")
    sep = septatype.septatype_string(aas, m.qnrB2, m.qnrS1)
    analysis = {"recN_codon": codon or None, "warnings": []}
    return _result(m.species, m.habitat, analysis, sep) | {"analysis": None, "recn_codon": codon or None}


# ---------------------------------------------------------------- database
def _filters(species="", habitat="", st="", country="", septatype_="", tetratype="", phenotype="", year_from="",
             year_to="", wins=False, source="", q=""):
    for v in (year_from, year_to):
        if v and not str(v).isdigit():
            raise HTTPException(400, "Years must be numbers.")
    return dict(species=species, habitat=habitat, st=st, country=country, septatype=septatype_, tetratype=tetratype,
                phenotype=phenotype, year_from=year_from, year_to=year_to, wins=wins, source=source, q=q[:80])


@app.get("/api/v1/isolates")
def isolates(species: str = "", habitat: str = "", st: str = "", country: str = "", septatype: str = "",
             tetratype: str = "", phenotype: str = "", year_from: str = "", year_to: str = "", wins: bool = False,
             source: str = "", q: str = "", sort: str = "genome_id", desc: bool = False,
             page: int = Query(1, ge=1), per_page: int = Query(25, ge=1, le=1000)):
    f = _filters(species, habitat, st, country, septatype, tetratype, phenotype, year_from, year_to, wins, source, q)
    r = rules.isolates_query(f, sort, desc, page, per_page)
    return {"total": r["total"], "page": page, "per_page": per_page, "rows": r["rows"]}


@app.get("/api/v1/isolates/{genome_id}")
def isolate(genome_id: str):
    r = rules.isolates_query({"q": genome_id}, limit=False)
    for row in r["rows"]:
        if row["genome_id"] == genome_id:
            return row
    raise HTTPException(404, "Genome ID not found.")


@app.get("/api/v1/rules")
def rules_table(species: str, habitat: str = "ALL"):
    _check(species, habitat)
    db = rules.conn()
    try:
        rows = db.execute("SELECT septatype,tetratype,n_isolates,n_resistant,n_susceptible,predicted,confidence,mean_h "
                          "FROM septatypes WHERE species=? AND habitat=? ORDER BY n_isolates DESC, septatype",
                          (species, habitat)).fetchall()
        st = db.execute("SELECT * FROM species_stats WHERE species=? AND habitat=?", (species, habitat)).fetchone()
    finally:
        db.close()
    out = [dict(r) for r in rows]
    for r in out:                      # recomputed, so a live database built before the helix change still shows current values
        r["mean_h"] = grammar.mean_h(r["tetratype"].split("-"))
    return {"stats": dict(st) if st else None, "rows": out}


@app.get("/api/v1/stats/accuracy")
def stats_accuracy():
    return rules.accuracy_stats() | {"benchmark": rules.benchmark(), "edgecases": rules.edgecases()}


@app.get("/api/v1/stats/geography")
def stats_geo(species: str = "", st: str = "", tetratype: str = "", phenotype: str = "", septatype: str = ""):
    return {"rows": rules.geography(species, st, tetratype, phenotype, septatype=septatype),
            "no_country": rules.geography(species, st, tetratype, phenotype, missing=True, septatype=septatype)}


@app.get("/api/v1/stats/septatypes")
def stats_sep(species: str, top: int = Query(12, ge=1, le=40)):
    if species not in SPECIES:
        raise HTTPException(400, "Unknown species.")
    return rules.septatype_freq(species, top)


@app.get("/api/v1/stats/cutoff")
def stats_cutoff():
    return rules.cutoff_stats()


@app.get("/api/v1/grammar")
def grammar_explain(residues: str, species: str = "ecoli", reference: str = "SDSE", recn: str = "L", codon: str = "",
                    qnrB: bool = False, qnrS: bool = False):
    """The three checkpoints worked for four residues (GyrA-83, GyrA-87, ParC-80, ParC-84) next to a reference.
    Optional RecN-346 residue and codon and qnr presence. Used by the calculation page."""
    r, ref, rn, cd = residues.strip().upper(), reference.strip().upper(), recn.strip().upper(), codon.strip().upper()
    if species not in SPECIES or any(len(x) != 4 or any(a not in AA for a in x) for x in (r, ref)) or rn not in AA:
        raise HTTPException(400, "Give four amino-acid letters for the sample and for the reference, for example LNIE and SDSE.")
    aas = list(r) + [rn]
    g = grammar.predict(species, aas, qnrB, qnrS, cd or None)
    return grammar.explain(species, aas, g, list(ref), (qnrB, qnrS), cd or None)


@app.get("/api/v1/aminoacids")
def aminoacids():
    return {"scale": "Eisenberg normalised consensus (1984)",
            "rows": [{"aa": a, "name": grammar.NAMES[a], "hydrophobicity": h, "class": grammar.aa_class(a)}
                     for a, h in grammar.EISENBERG.items()]}


# ----------------------------------------------------------------- downloads
SITE_NAMES = {k: label for k, label, _ in SOURCES} | {"ALL": "Source not stated"}
CSV_COLUMNS = [("genome_id", "Genome ID"), ("genome_name", "Genome name"), ("species_name", "Species"), ("source_name", "Isolation source"),
               ("country", "Country"), ("septatype", "Septatype"), ("septatyper_call", "Ciprotyper predicted phenotype"),
               ("lab_phenotype", "BV-BRC Phenotype Info: Lab"), ("bvbrc_call", "BV-BRC Phenotype Info: Computationally predicted"),
               ("submitter", "Submitter"), ("date_display", "Date"), ("affiliation", "Affiliation"), ("publication", "Publication"),
               ("st", "ST"), ("year_raw", "Year"), ("community", "Community submission")]


def _download(kind: str, f: dict):
    rows = rules.isolates_query(f, limit=False)["rows"]
    if kind == "json":
        return JSONResponse(rows, headers={"Content-Disposition": "attachment; filename=ciprotyper_isolates.json"})
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow([h for _, h in CSV_COLUMNS])
    for r in rows:
        r = dict(r, species_name=SPECIES.get(r["species"], r["species"]), source_name=SITE_NAMES.get(r["habitat"], r["habitat"]),
                 community="yes" if r.get("source") == "community" else "no")
        w.writerow(["" if r.get(k) is None else r.get(k) for k, _ in CSV_COLUMNS])
    return StreamingResponse(iter([buf.getvalue()]), media_type="text/csv",
                             headers={"Content-Disposition": "attachment; filename=ciprotyper_isolates.csv"})


@app.get("/api/v1/download/{kind}")
def download(kind: str, species: str = "", habitat: str = "", st: str = "", country: str = "", septatype: str = "",
             tetratype: str = "", phenotype: str = "", year_from: str = "", year_to: str = "", wins: bool = False,
             source: str = "", q: str = ""):
    if kind not in ("csv", "json"):
        raise HTTPException(404, "Use csv or json.")
    return _download(kind, _filters(species, habitat, st, country, septatype, tetratype, phenotype, year_from, year_to, wins, source, q))


# ----------------------------------------------------------- BV-BRC lookup
@app.get("/api/v1/bvbrc/{genome_id}")
def bvbrc_lookup(genome_id: str):
    """Ask BV-BRC for the official genome name. Fails soft: submissions still work without it."""
    if not GENOME_ID.match(genome_id):
        raise HTTPException(400, "A BV-BRC genome ID looks like 562.7710.")
    url = ("https://www.bv-brc.org/api/genome/?" + urllib.parse.quote(f"eq(genome_id,{genome_id})&select(genome_id,genome_name)", safe="=&(),")
           + "&http_accept=application/json")
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Ciprotyper/2.0"}), timeout=8) as r:
            data = json.loads(r.read().decode())
    except Exception as e:  # network, timeout, bad JSON
        return {"found": False, "genome_id": genome_id, "reason": f"BV-BRC could not be reached ({type(e).__name__})."}
    if data:
        return {"found": True, "genome_id": genome_id, "genome_name": data[0].get("genome_name")}
    return {"found": False, "genome_id": genome_id, "reason": "No genome with this ID at BV-BRC."}


# --------------------------------------------------------------- submissions
@app.post("/api/v1/submissions")
def submit(species: str = Form(...), habitat: str = Form(...), genome_id: str = Form(""), genome_name: str = Form(""),
           country: str = Form(""), year: str = Form(""), source_text: str = Form(""),
           lab_phenotype: str = Form(""), email: str = Form(""), submitter: str = Form(""), affiliation: str = Form(""),
           publication: str = Form(""), file: UploadFile | None = File(None),
           sequence: str = Form(""), qnrB2: str = Form("auto"), qnrS1: str = Form("auto")):
    _check(species, habitat)
    if habitat == "ALL":
        raise HTTPException(400, "Choose the isolation source.")
    country_name = None
    if country.strip():
        country_name = countries.canonical(country)
        if not country_name:
            raise HTTPException(400, "Pick the country from the list.")
    if genome_id and not GENOME_ID.match(genome_id.strip()):
        raise HTTPException(400, "A BV-BRC genome ID looks like 562.7710. Leave it empty if the genome is not at BV-BRC.")
    if lab_phenotype and lab_phenotype not in ("Resistant", "Susceptible"):
        raise HTTPException(400, "Phenotype must be Resistant or Susceptible.")
    if year and not re.fullmatch(r"\d{4}", year.strip()):
        raise HTTPException(400, "Year must be four digits.")
    if not submitter.strip():
        raise HTTPException(400, "Enter the submitter's name.")
    if not email.strip():
        raise HTTPException(400, "Enter a contact email. Curators use it if they have a question.")
    if not re.fullmatch(r"[^@\s]+@[^@\s]+\.[^@\s]+", email.strip()):
        raise HTTPException(400, "Email looks wrong.")
    pub = publication.strip()
    if pub:
        if re.fullmatch(r"(?i)(doi:\s*)?10\.\d{4,9}/\S+", pub):
            pub = "https://doi.org/" + re.sub(r"(?i)^doi:\s*", "", pub)
        elif not re.fullmatch(r"https?://[^\s<>\"']{4,290}", pub):
            raise HTTPException(400, "Publication should be a DOI (10.1038/...) or a web link starting with http.")
    if file is not None and file.filename:
        text = _read_upload(file)
    elif sequence.strip():
        text = sequence if sequence.lstrip().startswith(">") else ">pasted\n" + sequence
    else:
        raise HTTPException(400, "Provide a genome FASTA (file or pasted).")
    analysis = _analyse_text(text, species, {"qnrB2": qnrB2, "qnrS1": qnrS1})
    result = _result(species, habitat, analysis, analysis["septatype"])
    # keep only derived genotype data. The genome itself is never stored.
    d = dict(species=species, habitat=habitat, genome_id=genome_id.strip() or None, genome_name=genome_name.strip()[:200] or None,
             country=country_name, year=year.strip() or None, source_text=None,
             lab_phenotype=lab_phenotype or None, email=email.strip(), submitter=submitter.strip()[:120],
             affiliation=affiliation.strip()[:200] or None, publication=pub or None)
    sub = rules.create_submission(d, result)
    return sub | {"result": result}


@app.get("/api/v1/admin/queue", dependencies=[Depends(admin)])
def queue(status: str = ""):
    return {"rows": rules.list_queue(status), "counts": rules.queue_counts()}


class Decision(BaseModel):
    job_id: str
    action: str
    note: str = ""


@app.post("/api/v1/admin/decide", dependencies=[Depends(admin)])
def decide(d: Decision):
    try:
        return rules.decide(d.job_id, d.action, d.note[:300])
    except KeyError as e:
        raise HTTPException(404, str(e))
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.get("/api/v1/admin/backup", dependencies=[Depends(admin)])
def backup_download():
    return JSONResponse(persist.export_backup(), headers={"Content-Disposition": "attachment; filename=ciprotyper_backup.json"})


@app.post("/api/v1/admin/restore", dependencies=[Depends(admin)])
def backup_restore(data: dict = Body(...)):
    try:
        res = persist.import_backup(data)
    except ValueError as e:
        raise HTTPException(400, str(e))
    persist.schedule_save()
    return res


@app.get("/api/v1/admin/persist-status", dependencies=[Depends(admin)])
def persist_status():
    return persist.status()


@app.get("/healthz")
def healthz():
    import shutil
    return {"ok": True, "blast": bool(shutil.which("blastn"))}


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


app.mount("/static", StaticFiles(directory=STATIC), name="static")
