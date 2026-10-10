#!/usr/bin/env python3
"""
Rebuild septatype.db from Table_S1.xlsx.
Run: python3 build_db.py /path/to/Table_S1.xlsx
Header row is at row 3 (1-based). Section rows split species/site.
"""
import pandas as pd, sqlite3, sys, os

XLSX = sys.argv[1] if len(sys.argv) > 1 else "Table_S1.xlsx"
DB = "septatype.db"

df = pd.read_excel(XLSX, sheet_name="Table S1", header=None)
data = df[df.iloc[:, 2].isin(["Resistant", "Susceptible"])].copy()
data.columns = ["genome_id", "isolate", "pheno", "year", "country", "mlst",
                "gyra83", "gyra87", "parc80", "parc84", "recn346", "qnrb", "qnrs"]
sections = [(3, "E. coli", "fecal"), (188, "E. coli", "urine"),
            (405, "S. sonnei", "fecal"), (740, "S. enterica", "fecal"),
            (1251, "K. pneumoniae", "fecal"), (1329, "K. pneumoniae", "urine"),
            (2345, "K. pneumoniae", "respiratory")]

def get_group(idx):
    for i, (srow, sp, site) in enumerate(sections):
        erow = sections[i+1][0] if i+1 < len(sections) else 10**9
        if srow < idx < erow:
            return sp, site
    return None, None

data["species"] = data.index.map(lambda i: get_group(i)[0])
data["site"] = data.index.map(lambda i: get_group(i)[1])
data["septatype"] = (data["gyra83"].astype(str) + "-" + data["gyra87"].astype(str) + "-" +
                     data["parc80"].astype(str) + "-" + data["parc84"].astype(str) + "-" +
                     data["recn346"].astype(str) + "-" + data["qnrb"].astype(str) + "-" +
                     data["qnrs"].astype(str))
data["year"] = data["year"].astype(str)
data["isolate"] = data["isolate"].astype(str).str.strip()

if os.path.exists(DB):
    os.remove(DB)
con = sqlite3.connect(DB)
iso = data[["genome_id", "isolate", "pheno", "year", "country", "mlst",
            "gyra83", "gyra87", "parc80", "parc84", "recn346", "qnrb", "qnrs",
            "species", "site", "septatype"]].copy()
iso.columns = ["genome_id", "isolate_name", "phenotype", "year", "country", "mlst",
               "gyra83", "gyra87", "parc80", "parc84", "recn346", "qnrb", "qnrs",
               "species", "site", "septatype"]
iso.to_sql("isolates", con, index=False)

rules = data.groupby(["species", "site", "septatype"]).agg(
    n=("pheno", "size"),
    n_res=("pheno", lambda s: (s == "Resistant").sum())).reset_index()
rules["n_sus"] = rules["n"] - rules["n_res"]
rules["call"] = (rules["n_res"] > rules["n_sus"]).map({True: "R", False: "S"})
rules[["species", "site", "septatype", "n", "n_res", "n_sus", "call"]].to_sql(
    "rules", con, index=False)

data["four"] = (data["gyra83"].astype(str) + "-" + data["gyra87"].astype(str) + "-" +
                data["parc80"].astype(str) + "-" + data["parc84"].astype(str))
fours = data.groupby(["species", "site"])["four"].unique().reset_index()
fours["combos"] = fours["four"].apply(lambda x: ",".join(sorted(x)))
fours[["species", "site", "combos"]].to_sql("four_combos", con, index=False)
con.execute("CREATE INDEX idx_rules ON rules(species, site, septatype)")
con.execute("CREATE INDEX idx_iso ON isolates(species, site)")
con.commit()
print("isolates:", con.execute("SELECT COUNT(*) FROM isolates").fetchone()[0])
print("rules:", con.execute("SELECT COUNT(*) FROM rules").fetchone()[0])
con.close()
