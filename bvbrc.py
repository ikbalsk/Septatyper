"""Fetch BV-BRC genome names and computational ciprofloxacin predictions from inside the server.

Runs in a background thread at start-up (unless BVBRC_AUTOLOAD=0) and on demand from the curator page.
Fails soft: if BV-BRC cannot be reached the site works normally and the Browse page says so.
It is built to survive a slow BV-BRC: every request is retried with a pause, a batch that still fails is skipped and tried
again later, what was fetched is saved after every batch (so a stop or restart loses nothing and the next run resumes), and a
supervisor keeps trying in the background (BVBRC_RETRY_S apart, up to MAX_ATTEMPTS times) until every genome is done.
The free Render plan forgets the database on restart, so the fetch simply repeats after each start.
"""
from __future__ import annotations

import json
import os
import re
import socket
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter

from . import rules
from .config import VERSION

API = "https://www.bv-brc.org/api"
UA = {"User-Agent": f"Ciprotyper/{VERSION} (research)", "Accept": "application/json"}
# A saved copy of the BV-BRC answers (bvbrc_cache.json in the repo root). When it exists, every start fills all genomes from it in
# seconds, so the comparison is there at once after a restart or a redeploy. Without it (or if GitHub cannot be reached) the
# loader asks BV-BRC itself. Make the file once from /api/v1/bvbrc-cache after a full load and upload it to the repo.
CACHE_URL = os.environ.get("BVBRC_CACHE_URL", "https://raw.githubusercontent.com/ikbalsk/Septatyper/main/bvbrc_cache.json")
BATCH = 25                                    # genomes per request pair (smaller replies, fewer timeouts)
TIMEOUT_S = 90                                # read timeout of one loader request
RETRY_PAUSES = (3, 8, 20)                     # seconds to wait before the 2nd, 3rd and 4th try of one request
RETRY_S = int(os.environ.get("BVBRC_RETRY_S", "600"))    # pause between background passes while genomes are still missing
MAX_ATTEMPTS = int(os.environ.get("BVBRC_MAX_ATTEMPTS", "12"))
STATE = {"running": False, "done": 0, "total": 0, "error": None, "names": 0, "calls": 0, "evidence": {},
         "failed": 0, "attempt": 0, "retrying": False, "source": None}
_LOCK = threading.Lock()


def _retryable(e: Exception) -> bool:
    """Timeouts, dropped connections, server errors (5xx) and unreadable replies are worth another try; a 4xx is not."""
    if isinstance(e, urllib.error.HTTPError):
        return e.code >= 500 or e.code == 429
    return isinstance(e, (TimeoutError, socket.timeout, urllib.error.URLError, ConnectionError, json.JSONDecodeError, OSError))


def _get(endpoint: str, rql: str, timeout: int = 60, retries: int = 0):
    """One BV-BRC query. The loader asks for retries (with pauses); a visitor's own request (a single genome) does not wait."""
    url = f"{API}/{endpoint}/?{urllib.parse.quote(rql, safe='=&(),.')}&http_accept=application/json"
    pauses = (0,) + tuple(RETRY_PAUSES[:retries])
    last = None
    for pause in pauses:
        if pause:
            time.sleep(pause)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            last = e
            if not _retryable(e):
                raise
    raise last


def is_loaded() -> bool:
    db = rules.conn()
    try:
        r = db.execute("SELECT value FROM meta WHERE key='bvbrc_loaded'").fetchone()
    finally:
        db.close()
    return bool(r and r[0] == "1")


def status() -> dict:
    return {"loaded": is_loaded(), **STATE}


def export_cache() -> dict:
    """The BV-BRC names and calls of the founding isolates as a small JSON document (what bvbrc_cache.json holds)."""
    db = rules.conn()
    try:
        rows = db.execute("SELECT genome_id, genome_name, bvbrc_call, bvbrc_method FROM isolates WHERE source='founding' "
                          "AND (genome_name IS NOT NULL OR bvbrc_call IS NOT NULL) ORDER BY genome_id").fetchall()
    finally:
        db.close()
    return {"source": "BV-BRC", "made": time.strftime("%Y-%m-%d", time.gmtime()), "genomes": {r[0]: [r[1], r[2], r[3]] for r in rows}}


def load_cache(db) -> bool:
    """Fast path: fill every genome from the saved copy. True when it was found, read and held computational calls."""
    try:
        req = urllib.request.Request(CACHE_URL, headers=UA)
        with urllib.request.urlopen(req, timeout=30) as r:
            data = json.loads(r.read().decode())
        genomes = data["genomes"]
        if not isinstance(genomes, dict) or not genomes:
            return False
    except Exception:
        return False
    ids = [r[0] for r in db.execute("SELECT genome_id FROM isolates WHERE source='founding'")]
    n_calls = 0
    for gid in ids:
        nm, call, meth = (genomes.get(gid) or [None, None, None])[:3]
        db.execute("UPDATE isolates SET genome_name=COALESCE(?,genome_name), bvbrc_call=?, bvbrc_method=? WHERE genome_id=?", (nm, call, meth, gid))
        db.execute("INSERT OR IGNORE INTO bvbrc_done VALUES (?)", (gid,))
        n_calls += 1 if call else 0
    if not n_calls:
        db.rollback()
        return False
    db.execute("INSERT OR REPLACE INTO meta VALUES ('bvbrc_loaded','1')")
    db.commit()
    STATE.update(total=len(ids), done=len(ids), calls=n_calls, names=sum(1 for g in ids if (genomes.get(g) or [None])[0]),
                 failed=0, error=None, source="cache")
    return True


def run() -> bool:
    """One pass: fetch every genome that is not done yet, saving after each batch. Returns False when a pass is already running."""
    with _LOCK:
        if STATE["running"]:
            return False
        STATE.update(running=True, error=None, failed=0, retrying=False)
    db = None
    try:
        db = rules.conn()
        db.execute("CREATE TABLE IF NOT EXISTS bvbrc_done (genome_id TEXT PRIMARY KEY)")      # genomes already fetched (resume point)
        ids = [r[0] for r in db.execute("SELECT genome_id FROM isolates WHERE source='founding' ORDER BY genome_id")]
        done_ids = {r[0] for r in db.execute("SELECT genome_id FROM bvbrc_done")}
        if not done_ids and load_cache(db):         # nothing fetched yet: try the saved copy first (seconds instead of minutes)
            return True
        pending = [g for g in ids if g not in done_ids]
        STATE.update(total=len(ids), done=len(ids) - len(pending), source="BV-BRC")
        evid, failed, last_err = Counter(STATE.get("evidence") or {}), [], None
        for i in range(0, len(pending), BATCH):
            chunk = pending[i:i + BATCH]
            lst = ",".join(chunk)
            try:
                names = {g["genome_id"]: g.get("genome_name") for g in
                         _get("genome", f"in(genome_id,({lst}))&select(genome_id,genome_name)&limit({BATCH})", TIMEOUT_S, 3)}
                amr = _get("genome_amr", f"and(in(genome_id,({lst})),eq(antibiotic,ciprofloxacin))"
                                         f"&select(genome_id,resistant_phenotype,evidence,computational_method)&limit(1000)", TIMEOUT_S, 3)
            except Exception as e:              # this batch could not be fetched: skip it now, it is tried again in the next pass
                failed += chunk
                last_err = f"{type(e).__name__}: {str(e)[:100]}"
                STATE.update(failed=len(failed), error=last_err)
                time.sleep(2)
                continue
            calls = {}
            for r in amr:
                evid[r.get("evidence")] += 1
                if "omputational" in (r.get("evidence") or "") and r.get("resistant_phenotype"):
                    calls[r["genome_id"]] = (r["resistant_phenotype"].strip().capitalize(), r.get("computational_method"))
            for gid in chunk:
                ph, method = calls.get(gid, (None, None))
                db.execute("UPDATE isolates SET genome_name=COALESCE(?,genome_name), bvbrc_call=?, bvbrc_method=? WHERE genome_id=?",
                           (names.get(gid), ph, method, gid))
                db.execute("INSERT OR IGNORE INTO bvbrc_done VALUES (?)", (gid,))
            db.commit()                         # saved: a stop or restart from here on loses nothing
            n_done = db.execute("SELECT COUNT(*) FROM bvbrc_done").fetchone()[0]
            n_calls = db.execute("SELECT COUNT(*) FROM isolates WHERE source='founding' AND bvbrc_call IS NOT NULL").fetchone()[0]
            n_names = db.execute("SELECT COUNT(*) FROM isolates WHERE source='founding' AND genome_name IS NOT NULL").fetchone()[0]
            STATE.update(done=n_done, names=n_names, calls=n_calls, evidence=dict(evid))
            time.sleep(0.3)
        n_done = db.execute("SELECT COUNT(*) FROM bvbrc_done").fetchone()[0]
        n_calls = db.execute("SELECT COUNT(*) FROM isolates WHERE source='founding' AND bvbrc_call IS NOT NULL").fetchone()[0]
        STATE.update(done=n_done, calls=n_calls, failed=len(failed))
        if failed:
            STATE["error"] = f"{last_err} ({len(failed)} of {len(ids)} genomes still to fetch)"
        elif n_calls:
            db.execute("INSERT OR REPLACE INTO meta VALUES ('bvbrc_loaded','1')")
            db.commit()
            STATE["error"] = None
        else:
            STATE["error"] = "BV-BRC returned names but no computational ciprofloxacin predictions"
    except Exception as e:  # database error, schema change
        STATE["error"] = f"{type(e).__name__}: {str(e)[:120]}"
    finally:
        if db is not None:
            db.close()
        STATE["running"] = False
    return True


_SUP = {"on": False}
_WAKE = threading.Event()                       # lets the curator page start the next pass without waiting


def _supervise() -> None:
    """Keep running passes until every genome is fetched, pausing RETRY_S between passes (at most MAX_ATTEMPTS passes)."""
    try:
        while STATE["attempt"] < MAX_ATTEMPTS:
            STATE["attempt"] += 1
            run()
            if is_loaded():
                break
            if STATE["attempt"] >= MAX_ATTEMPTS:
                break
            STATE["retrying"] = True
            _WAKE.wait(RETRY_S)
            _WAKE.clear()
    finally:
        STATE["retrying"] = False
        with _LOCK:
            _SUP["on"] = False


def start() -> None:
    """Start (once) the background loader; a second call while it is already working does nothing."""
    with _LOCK:
        STATE["attempt"] = 0
        if _SUP["on"]:
            _WAKE.set()                         # it is waiting between passes: start the next pass now
            return
        _SUP["on"] = True
    threading.Thread(target=_supervise, daemon=True).start()


# ----------------------------------------------------------------------------------------------
# One genome on demand: the Predict page's "BV-BRC ID" mode. Only a validated ID is ever put into a
# URL, and only on the fixed BV-BRC host, so nothing the visitor types can point the server elsewhere.
GENOME_ID = re.compile(r"^\d{1,9}\.\d{1,6}$")
SPECIES_NAMES = (("escherichia coli", "ecoli"), ("shigella sonnei", "ssonnei"),
                 ("salmonella enterica", "senterica"), ("klebsiella pneumoniae", "kpneumoniae"))


class TooLarge(Exception):
    pass


def valid_id(s: str) -> bool:
    return bool(GENOME_ID.match(s or ""))


def genome_meta(gid: str):
    """The BV-BRC record of one genome (name, species, source, year, country), or None when the ID does not exist."""
    rows = _get("genome", f"eq(genome_id,{gid})&select(genome_id,genome_name,species,isolation_source,collection_year,"
                          f"isolation_country,host_name,mlst)&limit(1)")
    return rows[0] if rows else None


def detect_species(meta: dict):
    """Our species key from the BV-BRC species field or genome name, or None when it is not one of the four."""
    text = f"{meta.get('species') or ''} {meta.get('genome_name') or ''}".lower()
    for needle, key in SPECIES_NAMES:
        if needle in text:
            return key
    return None


def _download(url: str, accept: str, max_bytes: int) -> bytes:
    req = urllib.request.Request(url, headers={**UA, "Accept": accept})
    buf = bytearray()
    with urllib.request.urlopen(req, timeout=120) as r:
        while chunk := r.read(1 << 20):
            buf += chunk
            if len(buf) > max_bytes:
                raise TooLarge(f"The genome is larger than {max_bytes // 1048576} MB.")
    return bytes(buf)


def genome_fasta(gid: str, max_bytes: int) -> str:
    """The assembly of one genome as FASTA text. First the API's FASTA output for the genome_sequence collection;
    if that does not come back as FASTA, the contigs are fetched as JSON and written out here."""
    base = f"{API}/genome_sequence/?eq(genome_id,{gid})"
    try:
        text = _download(base + "&limit(25000)&http_accept=application/dna%2Bfasta", "application/dna+fasta", max_bytes).decode("utf-8", "replace")
        if text.lstrip().startswith(">"):
            return text
    except TooLarge:
        raise
    except Exception:
        pass                                   # fall through to the JSON route
    rows = json.loads(_download(base + "&select(sequence_id,accession,description,sequence)&limit(25000)&http_accept=application/json",
                                "application/json", max_bytes * 2).decode())
    out = []
    for r in rows:
        if r.get("sequence"):
            out.append(f">{r.get('accession') or r.get('sequence_id') or gid} {r.get('description') or ''}".rstrip() + "\n" + r["sequence"])
    if not out:
        raise RuntimeError("BV-BRC returned no sequence for this genome ID")
    return "\n".join(out) + "\n"
