"""Fill the BV-BRC columns for the 2,734 founding isolates (run once, needs internet).

  python scripts/fetch_bvbrc.py            # writes official genome names + ciprofloxacin AMR predictions
  python scripts/fetch_bvbrc.py --dry-run  # only shows what BV-BRC returns for the first batch

For every genome ID it stores
  genome_name   the official BV-BRC genome name
  bvbrc_call    BV-BRC's computational ciprofloxacin phenotype (Resistant / Susceptible / Intermediate)
  bvbrc_method  the computational method BV-BRC reports

After it runs, the Browse table shows the BV-BRC column and the star for isolates where the
septatype call matches the lab result and BV-BRC's does not.

NOTE: written against BV-BRC's public RQL API (https://www.bv-brc.org/api). It could not be
run from the build environment, so check the --dry-run output first. If a field name differs,
the two constants below are the only things to change.
"""
import argparse, json, sqlite3, sys, time, urllib.parse, urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from app.config import DB_PATH, SEED_DB  # noqa: E402

API = "https://www.bv-brc.org/api"
UA = {"User-Agent": "Ceprotyper/2.0 (research)", "Accept": "application/json"}
BATCH = 50


def get(endpoint, rql):
    url = f"{API}/{endpoint}/?{urllib.parse.quote(rql, safe='=&(),.')}&http_accept=application/json"
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
        return json.loads(r.read().decode())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(DB_PATH if DB_PATH.exists() else SEED_DB))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    db = sqlite3.connect(a.db)
    ids = [r[0] for r in db.execute("SELECT genome_id FROM isolates WHERE source='founding' ORDER BY genome_id")]
    print(f"{len(ids)} genomes in {a.db}")
    names, calls, evid = {}, {}, Counter()
    for i in range(0, len(ids), BATCH):
        batch = ids[i:i + BATCH]
        lst = ",".join(batch)
        for g in get("genome", f"in(genome_id,({lst}))&select(genome_id,genome_name)&limit({BATCH})"):
            names[g["genome_id"]] = g.get("genome_name")
        amr = get("genome_amr", f"and(in(genome_id,({lst})),eq(antibiotic,ciprofloxacin))"
                                f"&select(genome_id,resistant_phenotype,evidence,computational_method)&limit(1000)")
        for r in amr:
            evid[r.get("evidence")] += 1
            if "omputational" in (r.get("evidence") or ""):
                calls[r["genome_id"]] = (r.get("resistant_phenotype"), r.get("computational_method"))
        print(f"  {min(i+BATCH, len(ids))}/{len(ids)}  names {len(names)}  computational calls {len(calls)}")
        if a.dry_run:
            print("first names:", list(names.items())[:3]); print("first calls:", list(calls.items())[:3])
            print("evidence values seen:", dict(evid)); return
        time.sleep(0.5)
    for gid in ids:
        ph, method = calls.get(gid, (None, None))
        ph = ph.capitalize() if ph else None
        db.execute("UPDATE isolates SET genome_name=?, bvbrc_call=?, bvbrc_method=? WHERE genome_id=?",
                   (names.get(gid), ph, method, gid))
    db.execute("INSERT OR REPLACE INTO meta VALUES ('bvbrc_loaded','1')")
    db.commit()
    print(f"Done. {len(names)} names, {len(calls)} computational predictions. Evidence values: {dict(evid)}")


if __name__ == "__main__":
    main()
