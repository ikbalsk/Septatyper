"""Keep community submissions and curation decisions across restarts, for free.

Free hosting forgets its files when the server restarts. This module saves what a curator needs to a small JSON
backup and can bring it back:

  * Download backup / Restore backup: a curator does it by hand from the curation page.
  * Automatic copy to a private GitHub repository: set BACKUP_GITHUB_REPO (owner/name) and BACKUP_GITHUB_TOKEN
    (a fine-grained token with Contents read and write on that one repository). BACKUP_GITHUB_BRANCH and
    BACKUP_GITHUB_PATH are optional. The server saves a few seconds after every submission or decision and
    restores at start-up when it has no submissions yet.

The backup holds the submissions table (derived genotype, metadata and decisions, never the genome) and
e-mail addresses if people gave them, so keep the repository private. Approved submissions are replayed on
restore: their isolates and their votes in the rules table are added again.
"""
from __future__ import annotations

import base64
import json
import os
import threading
import urllib.error
import urllib.request

from . import rules
from .config import HABITATS, SPECIES

API = "https://api.github.com"
STATE = {"last_saved": None, "last_restored": None, "error": None}
_LOCK = threading.Lock()
_PENDING = {"on": False}
DELAY_S = 8


def configured() -> bool:
    return bool(os.environ.get("BACKUP_GITHUB_REPO") and os.environ.get("BACKUP_GITHUB_TOKEN"))


def _columns(db) -> list[str]:
    return [r[1] for r in db.execute("PRAGMA table_info(submissions)")]


def export_backup() -> dict:
    db = rules.conn()
    try:
        rows = [dict(r) for r in db.execute("SELECT * FROM submissions ORDER BY created, job_id")]
    finally:
        db.close()
    return {"version": 1, "saved": rules._now(), "submissions": rows}


def import_backup(data: dict) -> dict:
    """Re-insert submissions that are not there yet and re-apply the approved ones."""
    if not isinstance(data, dict) or not isinstance(data.get("submissions"), list):
        raise ValueError("This does not look like a Ciprotyper backup.")
    db = rules.conn()
    added = applied = skipped = 0
    problems: list[str] = []
    try:
        cols = _columns(db)
        for r in sorted((x for x in data["submissions"] if isinstance(x, dict)), key=lambda x: (str(x.get("created") or ""), str(x.get("job_id")))):
            jid = str(r.get("job_id") or "")
            if not jid or r.get("species") not in SPECIES or r.get("habitat") not in HABITATS[r["species"]]:
                problems.append(f"{jid or '?'}: unusable row")
                continue
            if db.execute("SELECT 1 FROM submissions WHERE job_id=?", (jid,)).fetchone():
                skipped += 1
                continue
            vals = [r.get(c) if isinstance(r.get(c), (str, int, float, type(None))) else None for c in cols]
            db.execute(f"INSERT INTO submissions ({','.join(cols)}) VALUES ({','.join('?' * len(cols))})", vals)
            added += 1
            if r.get("status") == "approved":
                try:
                    rules.apply_approval(db, r)
                    applied += 1
                except (ValueError, KeyError, TypeError) as e:
                    problems.append(f"{jid}: {e}")
        db.commit()
    finally:
        db.close()
    return {"added": added, "approved_reapplied": applied, "already_present": skipped, "problems": problems[:20]}


# ------------------------------------------------------------------ GitHub
def _http(method: str, url: str, headers: dict | None = None, data: bytes | None = None):
    req = urllib.request.Request(url, data=data, method=method, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def _gh():
    repo, token = os.environ["BACKUP_GITHUB_REPO"], os.environ["BACKUP_GITHUB_TOKEN"]
    path = os.environ.get("BACKUP_GITHUB_PATH", "septatyper_backup.json")
    branch = os.environ.get("BACKUP_GITHUB_BRANCH", "")
    url = f"{API}/repos/{repo}/contents/{path}" + (f"?ref={branch}" if branch else "")
    h = {"Authorization": f"Bearer {token}", "Accept": "application/vnd.github+json", "X-GitHub-Api-Version": "2022-11-28",
         "User-Agent": "Ciprotyper-backup"}
    return url, h, branch


def _github_read():
    """Returns (sha, raw bytes) of the backup file, or (None, None) when it does not exist yet."""
    url, h, _ = _gh()
    st, body = _http("GET", url, h)
    if st == 404:
        return None, None
    if st != 200:
        raise RuntimeError(f"GitHub returned {st} when reading the backup")
    sha = json.loads(body)["sha"]
    st2, raw = _http("GET", url, {**h, "Accept": "application/vnd.github.raw+json"})
    if st2 != 200:
        raise RuntimeError(f"GitHub returned {st2} when downloading the backup")
    return sha, raw


def save_now() -> None:
    if not configured():
        return
    with _LOCK:
        _PENDING["on"] = False
        try:
            payload = json.dumps(export_backup(), ensure_ascii=False, indent=1).encode()
            url, h, branch = _gh()
            for _ in range(2):                       # one retry if the file changed under us
                sha, _raw = _github_read()
                body = {"message": "Ciprotyper backup", "content": base64.b64encode(payload).decode()}
                if sha:
                    body["sha"] = sha
                if branch:
                    body["branch"] = branch
                st, resp = _http("PUT", url.split("?")[0], {**h, "Content-Type": "application/json"}, json.dumps(body).encode())
                if st in (200, 201):
                    STATE.update(last_saved=rules._now(), error=None)
                    return
                if st != 409:
                    raise RuntimeError(f"GitHub returned {st} when saving the backup")
            raise RuntimeError("GitHub reported a conflict twice")
        except Exception as e:  # noqa: BLE001
            STATE["error"] = f"{type(e).__name__}: {str(e)[:150]}"


def schedule_save() -> None:
    """Save a few seconds from now, once for any number of quick changes."""
    if not configured() or _PENDING["on"]:
        return
    _PENDING["on"] = True
    t = threading.Timer(DELAY_S, save_now)
    t.daemon = True
    t.start()


def restore_if_configured() -> dict | None:
    """At start-up: when the server has no submissions yet, bring the saved ones back."""
    if not configured():
        return None
    try:
        db = rules.conn()
        try:
            n = db.execute("SELECT COUNT(*) FROM submissions").fetchone()[0]
        finally:
            db.close()
        if n:
            return None
        _sha, raw = _github_read()
        if not raw:
            return None
        res = import_backup(json.loads(raw))
        STATE.update(last_restored=rules._now(), error=None)
        return res
    except Exception as e:  # noqa: BLE001
        STATE["error"] = f"{type(e).__name__}: {str(e)[:150]}"
        return None


def status() -> dict:
    db = rules.conn()
    try:
        n = db.execute("SELECT COUNT(*) FROM submissions").fetchone()[0]
    finally:
        db.close()
    return {"github_backup": configured(), "submissions": n, **STATE}
