"""Build data/seed.db from Table_S1.xlsx.

Tables
  isolates       2,734 founding isolates (BV-BRC genome ID is the primary key)
  septatypes     rules table: majority phenotype per septatype for each species and site
  species_stats  isolate counts and reference accuracy per species and site
  aa_properties  Eisenberg hydrophobicity and residue class for the 20 amino acids
  submissions    community submissions and the curation queue
  meta           build info, grammar-engine scores

Rule: majority phenotype per septatype within each species/site. Ties go to Susceptible.
"""
import datetime, json, os, sqlite3, sys
from collections import defaultdict
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from app import grammar  # noqa: E402
from app.config import ROOT, SECTION_TO_SPECIES, SEED_DB  # noqa: E402

XLSX = Path(os.environ.get("TABLE_S1", ROOT / "data" / "Table_S1.xlsx"))

ISO3 = {"USA": "USA", "Australia": "AUS", "United Kingdom": "GBR", "India": "IND", "Thailand": "THA", "Kenya": "KEN",
        "Nepal": "NPL", "Mali": "MLI", "Jersey": "JEY", "Pakistan": "PAK", "Greece": "GRC", "Mozambique": "MOZ",
        "Switzerland": "CHE", "Gambia": "GMB", "Brazil": "BRA", "China": "CHN", "Tanzania": "TZA", "Italy": "ITA",
        "Germany": "DEU", "Bangladesh": "BGD", "Sweden": "SWE", "Taiwan": "TWN", "Canada": "CAN", "Honduras": "HND",
        "Norway": "NOR", "United Arab Emirates": "ARE", "Cambodia": "KHM", "Senegal": "SEN"}

SCHEMA = """
DROP TABLE IF EXISTS isolates; DROP TABLE IF EXISTS septatypes; DROP TABLE IF EXISTS species_stats;
DROP TABLE IF EXISTS aa_properties; DROP TABLE IF EXISTS submissions; DROP TABLE IF EXISTS meta;
CREATE TABLE isolates (
  genome_id TEXT PRIMARY KEY, genome_name TEXT, isolate_name TEXT, species TEXT, habitat TEXT,
  lab_phenotype TEXT, year_raw TEXT, year_start INTEGER, year_end INTEGER, year_note TEXT,
  country TEXT, iso3 TEXT, st TEXT,
  gyrA_83 TEXT, gyrA_87 TEXT, parC_80 TEXT, parC_84 TEXT, recN_346 TEXT, qnrB2 INTEGER, qnrS1 INTEGER,
  septatype TEXT, tetratype TEXT, septatyper_call TEXT,
  bvbrc_call TEXT, bvbrc_method TEXT, source TEXT DEFAULT 'founding', date_added TEXT);
CREATE INDEX ix_iso_a ON isolates(species, habitat, septatype);
CREATE INDEX ix_iso_b ON isolates(st); CREATE INDEX ix_iso_c ON isolates(country);
CREATE TABLE septatypes (
  species TEXT, habitat TEXT, septatype TEXT, tetratype TEXT, n_isolates INTEGER, n_resistant INTEGER,
  n_susceptible INTEGER, predicted TEXT, confidence REAL, mean_h REAL,
  PRIMARY KEY (species, habitat, septatype));
CREATE TABLE species_stats (species TEXT, habitat TEXT, n_isolates INTEGER, n_resistant INTEGER,
  n_septatypes INTEGER, accuracy REAL, grammar_accuracy REAL, PRIMARY KEY (species, habitat));
CREATE TABLE aa_properties (aa TEXT PRIMARY KEY, name TEXT, hydrophobicity REAL, class TEXT, scale TEXT);
CREATE TABLE submissions (
  job_id TEXT PRIMARY KEY, created TEXT, status TEXT, species TEXT, habitat TEXT,
  genome_id TEXT, genome_name TEXT, country TEXT, year TEXT, source_text TEXT, mic TEXT, lab_phenotype TEXT,
  email TEXT, septatype TEXT, tier TEXT, call TEXT, novel INTEGER, matches INTEGER, result_json TEXT,
  curator_note TEXT, decided TEXT);
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
"""


def parse():
    ws = openpyxl.load_workbook(XLSX, read_only=True).active
    section, out = None, []
    for r in list(ws.iter_rows(values_only=True))[3:]:
        if r[0] is None:
            continue
        if r[1] is None and r[2] is None:
            section = str(r[0]).strip()
            continue
        hab = section.split()[-1].upper()
        sp = SECTION_TO_SPECIES[section[: -len(hab)].strip()]
        g83, g87, p80, p84, r346 = [str(x).strip() for x in r[6:11]]
        b2, s1 = int(r[11]), int(r[12])
        raw = str(r[3]).strip() if r[3] else None
        ys = ye = None
        note = None
        if raw:
            parts = raw.replace(" ", "").split("-")
            try:
                nums = [int(x) for x in parts if x]
                ys, ye = nums[0], nums[-1]
                if ys < 1950:
                    note = f"Year range '{raw}' in Table S1 looks like a typo; start year left empty."
                    ys = None
            except ValueError:
                note = f"Unparsed year '{raw}'"
        country = str(r[4]).strip() if r[4] else None
        mlst = str(r[5]).strip() if r[5] and str(r[5]).strip() not in ("", "None", "nan") else None
        sep = "-".join([g83, g87, p80, p84, r346, str(b2), str(s1)])
        out.append(dict(
            genome_id=str(r[0]).strip(), genome_name=None, isolate_name=str(r[1]).strip(), species=sp, habitat=hab,
            lab_phenotype=str(r[2]).strip(), year_raw=raw, year_start=ys, year_end=ye, year_note=note,
            country=country, iso3=ISO3.get(country), st=mlst,
            gyrA_83=g83, gyrA_87=g87, parC_80=p80, parC_84=p84, recN_346=r346, qnrB2=b2, qnrS1=s1,
            septatype=sep, tetratype="-".join([g83, g87, p80, p84]), septatyper_call=None,
            bvbrc_call=None, bvbrc_method=None, source="founding", date_added=None))
    return out


def main():
    rows = parse()
    SEED_DB.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(SEED_DB)
    db.executescript(SCHEMA)

    groups = defaultdict(list)
    for r in rows:
        groups[(r["species"], r["habitat"])].append(r)
        groups[(r["species"], "ALL")].append(r)
    calls = {}
    for (sp, hab), items in groups.items():
        tab = defaultdict(lambda: [0, 0])
        for i in items:
            tab[i["septatype"]][0 if i["lab_phenotype"] == "Resistant" else 1] += 1
        correct = g_correct = 0
        for sep, (res, sus) in tab.items():
            pred = "R" if res > sus else "S"
            correct += res if pred == "R" else sus
            aas = sep.split("-")
            h = grammar.mean_h(aas[:4])
            db.execute("INSERT INTO septatypes VALUES (?,?,?,?,?,?,?,?,?,?)",
                       (sp, hab, sep, "-".join(aas[:4]), res + sus, res, sus, pred, max(res, sus) / (res + sus), h))
            if hab != "ALL":
                calls[(sp, hab, sep)] = "Resistant" if pred == "R" else "Susceptible"
        for i in items:
            aas = i["septatype"].split("-")
            g = grammar.predict(sp, aas[:5], aas[5] == "1", aas[6] == "1")
            g_correct += g["call"] == i["lab_phenotype"]
        db.execute("INSERT INTO species_stats VALUES (?,?,?,?,?,?,?)",
                   (sp, hab, len(items), sum(v[0] for v in tab.values()), len(tab),
                    round(100 * correct / len(items), 1), round(100 * g_correct / len(items), 1)))
    for r in rows:
        r["septatyper_call"] = calls[(r["species"], r["habitat"], r["septatype"])]
    cols = list(rows[0].keys())
    db.executemany(f"INSERT INTO isolates ({','.join(cols)}) VALUES ({','.join(':'+c for c in cols)})", rows)
    for a, h in grammar.EISENBERG.items():
        db.execute("INSERT INTO aa_properties VALUES (?,?,?,?,?)",
                   (a, grammar.NAMES[a], h, grammar.aa_class(a), "Eisenberg normalised consensus (1984)"))
    db.executemany("INSERT INTO meta VALUES (?,?)", [
        ("built", datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
        ("source", XLSX.name), ("n_isolates", str(len(rows))), ("bvbrc_loaded", "0")])
    db.commit()
    print(f"Built {SEED_DB} with {len(rows)} isolates")
    for r in db.execute("SELECT species,habitat,n_isolates,n_septatypes,accuracy,grammar_accuracy FROM species_stats ORDER BY 1,2"):
        print("  ", tuple(r))


if __name__ == "__main__":
    main()
