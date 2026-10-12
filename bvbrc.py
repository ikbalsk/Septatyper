"""Fetch BV-BRC genome names and computational ciprofloxacin predictions from inside the server.

Runs in a background thread at start-up (unless BVBRC_AUTOLOAD=0) and on demand from the curator page.
Fails soft: if BV-BRC cannot be reached the site works normally and the Browse page says so.
The free Render plan forgets the database on restart, so the fetch simply repeats after each start.
"""
from __future__ import annotations

import json
import threading
import time
import urllib.parse
import urllib.request
from collections import Counter

from . import rules

API = "https://www.bv-brc.org/api"
UA = {"User-Agent": "SeptaTyper/2.1 (research)", "Accept": "application/json"}
BATCH = 50
STATE = {"running": False, "done": 0, "total": 0, "error": None, "names": 0, "calls": 0, "evidence": {}}
_LOCK = threading.Lock()


def _get(endpoint: str, rql: str):
    url = f"{API}/{endpoint}/?{urllib.parse.quote(rql, safe='=&(),.')}&http_accept=application/json"
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return json.loads(r.read().decode())


def is_loaded() -> bool:
    db = rules.conn()
    try:
        r = db.execute("SELECT value FROM meta WHERE key='bvbrc_loaded'").fetchone()
    finally:
        db.close()
    return bool(r and r[0] == "1")


def status() -> dict:
    return {"loaded": is_loaded(), **STATE}


def run() -> bool:
    with _LOCK:
        if STATE["running"]:
            return False
        STATE.update(running=True, done=0, total=0, error=None, names=0, calls=0, evidence={})
    try:
        db = rules.conn()
        ids = [r[0] for r in db.execute("SELECT genome_id FROM isolates WHERE source='founding' ORDER BY genome_id")]
        STATE["total"] = len(ids)
        names, calls, evid = {}, {}, Counter()
        for i in range(0, len(ids), BATCH):
            lst = ",".join(ids[i:i + BATCH])
            for g in _get("genome", f"in(genome_id,({lst}))&select(genome_id,genome_name)&limit({BATCH})"):
                names[g["genome_id"]] = g.get("genome_name")
            amr = _get("genome_amr", f"and(in(genome_id,({lst})),eq(antibiotic,ciprofloxacin))"
                                     f"&select(genome_id,resistant_phenotype,evidence,computational_method)&limit(1000)")
            for r in amr:
                evid[r.get("evidence")] += 1
                if "omputational" in (r.get("evidence") or "") and r.get("resistant_phenotype"):
                    calls[r["genome_id"]] = (r["resistant_phenotype"].strip().capitalize(), r.get("computational_method"))
            STATE.update(done=min(i + BATCH, len(ids)), names=len(names), calls=len(calls), evidence=dict(evid))
            time.sleep(0.3)
        if not names and not calls:
            raise RuntimeError("BV-BRC returned no data for these genome IDs")
        for gid in ids:
            ph, method = calls.get(gid, (None, None))
            db.execute("UPDATE isolates SET genome_name=COALESCE(?,genome_name), bvbrc_call=?, bvbrc_method=? WHERE genome_id=?",
                       (names.get(gid), ph, method, gid))
        if calls:
            db.execute("INSERT OR REPLACE INTO meta VALUES ('bvbrc_loaded','1')")
        else:
            STATE["error"] = "BV-BRC returned names but no computational ciprofloxacin predictions"
        db.commit()
        db.close()
    except Exception as e:  # network, bad JSON, schema change
        STATE["error"] = f"{type(e).__name__}: {str(e)[:120]}"
    finally:
        STATE["running"] = False
    return True


def start() -> None:
    threading.Thread(target=run, daemon=True).start()
