#!/usr/bin/env python3
"""
Fetch the 5 Septatype reference sequences from NCBI.
Run: python3 fetch_references.py
Output: references/*.fasta + references/PROVENANCE.txt

All sequences are downloaded from NCBI Entrez; nothing is hardcoded.
If NCBI is unreachable, the script fails with a clear message.
"""
import urllib.request, urllib.parse, json, time, sys
from pathlib import Path

UA = {"User-Agent": "SeptatypeTool/1.0 (research)"}
REF_DIR = Path(__file__).parent / "references"

def get(url, params, timeout=90):
    full = url + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(full, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout).read().decode()

def esearch(db, term, retmax=5):
    d = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                       {"db": db, "term": term, "retmode": "json", "retmax": retmax}))
    return d["esearchresult"]["idlist"]

def efetch(db, id_, rettype, extra=None):
    params = {"db": db, "id": id_, "rettype": rettype, "retmode": "text"}
    if extra:
        params.update(extra)
    return get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params)

def main():
    REF_DIR.mkdir(exist_ok=True)
    prov = []

    # 1-3: E. coli K-12 MG1655 gyrA, parC, recN CDS from NC_000913.3
    # Coordinates come from the proteins' /coded_by qualifiers (not hardcoded guesses)
    genes = [("gyrA", "16130166"), ("parC", "16130915"), ("recN", "49176247")]
    coords = {}
    for name, pid in genes:
        gb = efetch("protein", pid, "gb")
        import re
        m = re.search(r'/coded_by="([^"]+)"', gb)
        if not m:
            sys.exit(f"ERROR: could not find coded_by for {name}")
        coords[name] = m.group(1)
        print(f"{name}: {m.group(1)}")
        time.sleep(0.5)

    import re
    for name, loc in coords.items():
        m = re.match(r"(complement\()?NC_000913\.3:(\d+)\.\.(\d+)\)?", loc)
        comp = bool(m.group(1))
        s, e = int(m.group(2)), int(m.group(3))
        length = e - s + 1
        # fetch from s to end (efetch ignores seq_end), then take the right slice
        fasta = efetch("nuccore", "NC_000913.3", "fasta",
                       {"seq_start": s, "strand": 2 if comp else 1})
        seq = "".join(fasta.strip().split("\n")[1:])
        cds = seq[len(seq)-length:]  # last `length` bp = RC(forward[s..e])
        out = REF_DIR / f"{name}_cds_ref.fasta"
        with open(out, "w") as f:
            f.write(f">E_coli_K12_MG1655_{name} NC_000913.3 {loc}\n")
            for i in range(0, len(cds), 70):
                f.write(cds[i:i+70] + "\n")
        print(f"saved {out} ({len(cds)} bp)")
        prov.append(f"{name}_cds_ref.fasta\n  E. coli K-12 MG1655, NC_000913.3 {loc}")
        time.sleep(0.5)

    # 4-5: qnrB1 and qnrS1 gene sequences
    for gene_label, acc, fname in [
        ("qnrB1", "DQ351241.1", "qnrB_ref.fasta"),
        ("qnrS1", "LC909223.1", "qnrS_ref.fasta"),
    ]:
        fasta = efetch("nuccore", acc, "fasta")
        out = REF_DIR / fname
        with open(out, "w") as f:
            f.write(fasta if fasta.startswith(">") else f">{acc}\n{fasta}")
        seq = "".join(open(out).read().strip().split("\n")[1:])
        print(f"saved {out} ({len(seq)} bp)")
        prov.append(f"{fname}\n  {acc} ({len(seq)} bp)")

    with open(REF_DIR / "PROVENANCE.txt", "w") as f:
        f.write("Reference sequences for Septatype server (downloaded from NCBI)\n")
        f.write("-" * 70 + "\n")
        f.write("\n\n".join(prov) + "\n")
    print("\nDone. Verify codon positions match the expected wild type (S-D-S-E-L).")

if __name__ == "__main__":
    main()
