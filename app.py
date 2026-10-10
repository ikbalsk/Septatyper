#!/usr/bin/env python3
"""
SeptaTyper web database: ciprofloxacin resistance prediction + 2,734-isolate repository.
Modules A (prediction engine) + B (browse/search repository) at NAR Database quality.
"""
import os, re, json, sqlite3, subprocess, tempfile, shutil, csv, io
from pathlib import Path
from fastapi import FastAPI, UploadFile, File, Form, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

BASE = Path(__file__).parent
DB_PATH = BASE / "septatype.db"
REF_DIR = BASE / "references"
BLAST_BIN = os.environ.get("BLAST_BIN", str(BASE / "blast" / "bin"))

SPECIES = ["E. coli", "S. sonnei", "S. enterica", "K. pneumoniae"]
SITES = ["fecal", "urine", "respiratory"]
SITES_BY_SP = {"E. coli": ["fecal", "urine"], "S. sonnei": ["fecal"],
               "S. enterica": ["fecal"],
               "K. pneumoniae": ["fecal", "urine", "respiratory"]}
BV_BRC_URL = "https://www.bv-brc.org/view/Genome/"

REF_GENES = {
    "gyrA": ("gyrA_cds_ref.fasta", {83: "GyrA-83", 87: "GyrA-87"}),
    "parC": ("parC_cds_ref.fasta", {80: "ParC-80", 84: "ParC-84"}),
    "recN": ("recN_cds_ref.fasta", {346: "RecN-346"}),
}
QNR_REFS = {"qnrB": "qnrB_ref.fasta", "qnrS": "qnrS_ref.fasta"}

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

app = FastAPI(title="SeptaTyper", version="1.0",
              description="Ciprofloxacin resistance prediction and isolate repository for Enterobacterales")

def get_db():
    con = sqlite3.connect(DB_PATH)
    con.row_factory = sqlite3.Row
    return con

# ---------------------------------------------------------------- BLAST pipeline
def read_fasta(path):
    lines = Path(path).read_text().strip().split("\n")
    return lines[0], "".join(lines[1:]).upper().replace(" ", "").replace("\r", "")

def parse_fasta_multi(text):
    records, header, seq_parts = [], None, []
    for line in text.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts).upper()))
            header, seq_parts = line[1:].split()[0][:60], []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return [(h, s) for h, s in records if s and re.fullmatch(r"[ACGTNacgtn]+", s)]

def run_blastn(query_fasta, subject_fasta, outfmt="6"):
    cmd = [os.path.join(BLAST_BIN, "blastn"), "-query", query_fasta,
           "-subject", subject_fasta, "-outfmt", outfmt,
           "-evalue", "1e-10", "-max_target_seqs", "5", "-max_hsps", "3"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        raise RuntimeError(f"blastn failed: {r.stderr[:300]}")
    return r.stdout

def map_codon_to_query(ref_aln, qry_aln, codon_pos):
    want = {(codon_pos-1)*3, (codon_pos-1)*3+1, (codon_pos-1)*3+2}
    ref_idx, q_nts = 0, {}
    for rc, qc in zip(ref_aln, qry_aln):
        if rc != "-":
            if ref_idx in want:
                q_nts[ref_idx] = qc
            ref_idx += 1
    if len(q_nts) != 3:
        return None
    codon = "".join(q_nts[(codon_pos-1)*3 + i] for i in range(3))
    return None if ("-" in codon or "N" in codon) else codon

def analyze_genes(workdir, query_records):
    subject_fa = os.path.join(workdir, "subject.fasta")
    with open(subject_fa, "w") as f:
        for h, s in query_records:
            f.write(f">{h}\n")
            for i in range(0, len(s), 70):
                f.write(s[i:i+70] + "\n")
    result = {"residues": {}, "qnr": {}, "details": {}}
    for gene, (reffile, positions) in REF_GENES.items():
        _, ref_seq = read_fasta(REF_DIR / reffile)
        qf = os.path.join(workdir, f"{gene}_q.fasta")
        with open(qf, "w") as f:
            f.write(f">{gene}\n{ref_seq}\n")
        out = run_blastn(qf, subject_fa,
            "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq")
        best = None
        for line in out.strip().split("\n"):
            if not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 12:
                continue
            bs = float(p[9])
            if best is None or bs > best["bitscore"]:
                best = {"pident": float(p[2]), "length": int(p[3]), "bitscore": bs,
                        "sseqid": p[1], "qseq": p[10], "sseq": p[11]}
        det = {"found": best is not None}
        if best:
            det.update({"hit": best["sseqid"], "pident": round(best["pident"], 1),
                        "aln_len": best["length"]})
            for pos, label in positions.items():
                codon = map_codon_to_query(best["qseq"], best["sseq"], pos)
                aa = CODON_TABLE.get(codon, "?") if codon else "?"
                result["residues"][label] = aa
                det[label] = {"codon": codon or "?", "aa": aa}
        else:
            for label in positions.values():
                result["residues"][label] = "?"
        result["details"][gene] = det
    for qname, reffile in QNR_REFS.items():
        _, ref_seq = read_fasta(REF_DIR / reffile)
        ref_len = len(ref_seq)
        qf = os.path.join(workdir, f"{qname}_q.fasta")
        with open(qf, "w") as f:
            f.write(f">{qname}\n{ref_seq}\n")
        out = run_blastn(qf, subject_fa, "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore")
        bp, bc = 0, 0
        for line in out.strip().split("\n"):
            if not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 10:
                continue
            bp = max(bp, float(p[2])); bc = max(bc, int(p[3]) / ref_len * 100)
        present = bp >= 80 and bc >= 90
        result["qnr"][qname] = 1 if present else 0
        result["details"][qname] = {"found": present, "pident": round(bp, 1),
                                    "coverage": round(bc, 1)}
    return result

def predict_from_septatype(species, site, septatype):
    con = get_db()
    try:
        if site == "unknown":
            sites = SITES_BY_SP[species]
            per_site = {}
            for s in sites:
                r = con.execute(
                    "SELECT n, n_res, n_sus, call FROM rules WHERE species=? AND site=? AND septatype=?",
                    (species, s, septatype)).fetchone()
                if r:
                    per_site[s] = dict(r)
            if not per_site:
                combos = set()
                for s in sites:
                    row = con.execute("SELECT combos FROM four_combos WHERE species=? AND site=?",
                                      (species, s)).fetchone()
                    if row:
                        combos.update(row["combos"].split(","))
                return {"call": "Unknown", "septatype": septatype,
                        "four_combos": sorted(combos), "per_site": {}}
            nr = sum(v["n"] for v in per_site.values())
            nres = sum(v["n_res"] for v in per_site.values())
            call = "R" if nres > nr - nres else "S"
            return {"call": "Resistant" if call == "R" else "Susceptible",
                    "septatype": septatype, "n": nr, "n_res": nres, "n_sus": nr - nres,
                    "confidence": round(max(nres, nr-nres)/nr*100, 1),
                    "per_site": {s: ("Resistant" if v["call"] == "R" else "Susceptible")
                                 for s, v in per_site.items()},
                    "site_disagreement": len({v["call"] for v in per_site.values()}) > 1,
                    "low_support": nr < 5, "mixed": max(nres, nr-nres)/nr*100 < 80}
        r = con.execute(
            "SELECT n, n_res, n_sus, call FROM rules WHERE species=? AND site=? AND septatype=?",
            (species, site, septatype)).fetchone()
        if not r:
            row = con.execute("SELECT combos FROM four_combos WHERE species=? AND site=?",
                              (species, site)).fetchone()
            return {"call": "Unknown", "septatype": septatype,
                    "four_combos": sorted(row["combos"].split(",")) if row else []}
        nr, nres, nsus, call = r["n"], r["n_res"], r["n_sus"], r["call"]
        conf = round(max(nres, nsus)/nr*100, 1)
        return {"call": "Resistant" if call == "R" else "Susceptible",
                "septatype": septatype, "n": nr, "n_res": nres, "n_sus": nsus,
                "confidence": conf, "low_support": nr < 5, "mixed": conf < 80}
    finally:
        con.close()

# ---------------------------------------------------------------- frontend
@app.get("/", response_class=HTMLResponse)
def index():
    return (BASE / "static" / "index.html").read_text()

# ---------------------------------------------------------------- REST API v1
@app.get("/api/v1/species")
def v1_species():
    return {"species": [{"name": s, "sites": SITES_BY_SP[s]} for s in SPECIES]}

@app.get("/api/v1/isolates")
def v1_isolates(species: str = Query(None), site: str = Query(None),
                phenotype: str = Query(None), country: str = Query(None),
                mlst: str = Query(None), septatype: str = Query(None),
                q: str = Query(None), page: int = Query(1, ge=1),
                per_page: int = Query(50, ge=1, le=500)):
    con = get_db()
    try:
        where, params = [], []
        if species: where.append("species=?"); params.append(species)
        if site: where.append("site=?"); params.append(site)
        if phenotype: where.append("phenotype=?"); params.append(phenotype)
        if country: where.append("country LIKE ?"); params.append(f"%{country}%")
        if mlst: where.append("mlst=?"); params.append(mlst)
        if septatype: where.append("septatype=?"); params.append(septatype)
        if q:
            where.append("(genome_id LIKE ? OR isolate_name LIKE ? OR country LIKE ? OR mlst LIKE ?)")
            params += [f"%{q}%"]*4
        w = ("WHERE " + " AND ".join(where)) if where else ""
        total = con.execute(f"SELECT COUNT(*) c FROM isolates {w}", params).fetchone()["c"]
        rows = con.execute(
            f"""SELECT genome_id, isolate_name, phenotype, septatype, species, site,
                        year, country, mlst, gyra83, gyra87, parc80, parc84, recn346, qnrb, qnrs
                 FROM isolates {w} ORDER BY genome_id LIMIT ? OFFSET ?""",
            params + [per_page, (page-1)*per_page]).fetchall()
        items = []
        for r in rows:
            d = dict(r)
            d["bvbrc_url"] = BV_BRC_URL + str(d["genome_id"])
            items.append(d)
        return {"total": total, "page": page, "per_page": per_page,
                "pages": (total + per_page - 1)//per_page, "isolates": items}
    finally:
        con.close()

@app.get("/api/v1/isolates/{genome_id}")
def v1_isolate(genome_id: str):
    con = get_db()
    try:
        r = con.execute("SELECT * FROM isolates WHERE genome_id=?", (genome_id,)).fetchone()
        if not r:
            raise HTTPException(404, "Isolate not found")
        d = dict(r); d["bvbrc_url"] = BV_BRC_URL + genome_id
        return d
    finally:
        con.close()

@app.get("/api/v1/rules")
def v1_rules(species: str = Query(None), site: str = Query(None)):
    con = get_db()
    try:
        w, p = [], []
        if species: w.append("species=?"); p.append(species)
        if site: w.append("site=?"); p.append(site)
        wh = ("WHERE " + " AND ".join(w)) if w else ""
        rows = con.execute(
            f"SELECT species, site, septatype, n, n_res, n_sus, call FROM rules {wh} ORDER BY species, site, septatype", p)
        return {"rules": [dict(r) for r in rows]}
    finally:
        con.close()

@app.get("/api/v1/stats/septatypes")
def v1_stats_septatypes(species: str = Query(None), site: str = Query(None), top: int = Query(20, ge=1, le=100)):
    """Septatype frequency distribution for bar charts."""
    con = get_db()
    try:
        w, p = [], []
        if species: w.append("species=?"); p.append(species)
        if site: w.append("site=?"); p.append(site)
        wh = ("WHERE " + " AND ".join(w)) if w else ""
        rows = con.execute(
            f"""SELECT septatype, SUM(n) n, SUM(n_res) n_res, SUM(n_sus) n_sus
                FROM rules {wh} GROUP BY septatype ORDER BY n DESC LIMIT ?""", p + [top])
        return {"septatypes": [dict(r) for r in rows]}
    finally:
        con.close()

@app.get("/api/v1/stats/geography")
def v1_stats_geo():
    """Isolate counts per country for maps."""
    con = get_db()
    try:
        rows = con.execute(
            "SELECT country, COUNT(*) n FROM isolates GROUP BY country ORDER BY n DESC")
        return {"countries": [dict(r) for r in rows]}
    finally:
        con.close()

@app.get("/api/v1/stats/accuracy")
def v1_accuracy():
    targets = [("E. coli","fecal"),("E. coli","urine"),("S. sonnei","fecal"),
               ("S. enterica","fecal"),("K. pneumoniae","fecal"),
               ("K. pneumoniae","urine"),("K. pneumoniae","respiratory")]
    con = get_db()
    try:
        out = []
        for sp, site in targets:
            rows = con.execute("SELECT septatype, phenotype FROM isolates WHERE species=? AND site=?",
                               (sp, site)).fetchall()
            rules = {r["septatype"]: r["call"] for r in con.execute(
                "SELECT septatype, call FROM rules WHERE species=? AND site=?", (sp, site))}
            ok = sum(1 for r in rows if (r["phenotype"]=="Resistant") == (rules[r["septatype"]]=="R"))
            out.append({"group": f"{sp} ({site})", "n": len(rows),
                        "accuracy": round(ok/len(rows)*100, 1)})
        return {"accuracy": out,
                "note": "Measured on the same isolates that define the rules; not independent validation."}
    finally:
        con.close()

@app.get("/api/v1/download/csv")
def v1_csv():
    con = get_db()
    try:
        rows = con.execute("SELECT * FROM isolates ORDER BY genome_id").fetchall()
        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(["genome_id","isolate_name","phenotype","septatype","species","site",
                    "year","country","mlst","gyra83","gyra87","parc80","parc84",
                    "recn346","qnrb","qnrs","bvbrc_url"])
        for r in rows:
            d = dict(r)
            w.writerow([d["genome_id"],d["isolate_name"],d["phenotype"],d["septatype"],
                        d["species"],d["site"],d["year"],d["country"],d["mlst"],
                        d["gyra83"],d["gyra87"],d["parc80"],d["parc84"],d["recn346"],
                        d["qnrb"],d["qnrs"], BV_BRC_URL+str(d["genome_id"])])
        buf.seek(0)
        return StreamingResponse(iter([buf.getvalue()]), media_type="text/csv",
               headers={"Content-Disposition": "attachment; filename=septatyper_isolates.csv"})
    finally:
        con.close()

@app.get("/api/v1/download/json")
def v1_json():
    con = get_db()
    try:
        rows = con.execute("SELECT * FROM isolates ORDER BY genome_id").fetchall()
        data = []
        for r in rows:
            d = dict(r); d["bvbrc_url"] = BV_BRC_URL + str(d["genome_id"]); data.append(d)
        return JSONResponse(data, headers={"Content-Disposition": "attachment; filename=septatyper_isolates.json"})
    finally:
        con.close()

@app.post("/api/v1/predict")
async def v1_predict(species: str = Form(...), site: str = Form(...),
                     fasta: UploadFile = File(...)):
    if species not in SPECIES:
        raise HTTPException(400, "Unknown species")
    if site != "unknown" and site not in SITES_BY_SP[species]:
        raise HTTPException(400, "Invalid site for this species")
    content = (await fasta.read()).decode("utf-8", errors="ignore")
    records = parse_fasta_multi(content)
    if not records:
        raise HTTPException(400, "No valid nucleotide sequences found in FASTA")
    if sum(len(s) for _, s in records) > 20_000_000:
        raise HTTPException(400, "FASTA too large (max 20 Mb)")
    workdir = tempfile.mkdtemp(prefix="septatyper_")
    try:
        analysis = analyze_genes(workdir, records)
        res = analysis["residues"]
        order = ["GyrA-83", "GyrA-87", "ParC-80", "ParC-84", "RecN-346"]
        if any(res.get(k, "?") in ("?", "*") for k in order):
            missing = [k for k in order if res.get(k, "?") in ("?", "*")]
            return {"call": "Unknown",
                    "reason": f"Could not determine residues: {', '.join(missing)}",
                    "analysis": analysis}
        st = "-".join([res[k] for k in order] +
                      [str(analysis["qnr"]["qnrB"]), str(analysis["qnr"]["qnrS"])])
        pred = predict_from_septatype(species, site, st)
        pred["analysis"] = analysis; pred["species"] = species; pred["site"] = site
        return JSONResponse(pred)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
