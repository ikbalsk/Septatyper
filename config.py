"""Shared constants and paths."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEED_DB = Path(os.environ.get("SEED_DB", ROOT / "data" / "seed.db"))      # read-only founding data, baked into the image
DATA_DIR = Path(os.environ.get("DATA_DIR", ROOT / "data"))
DB_PATH = Path(os.environ.get("DB_PATH", DATA_DIR / "septatype.db"))      # live copy: submissions + approved isolates
ADMIN_TOKEN = os.environ.get("ADMIN_TOKEN", "").strip()   # phone keyboards and dashboards often add stray spaces
REFS_DIR = Path(os.environ.get("REFS_DIR", ROOT / "refs"))
MAX_UPLOAD_MB = int(os.environ.get("MAX_UPLOAD_MB", "60"))
BLAST_TIMEOUT_S = int(os.environ.get("BLAST_TIMEOUT_S", "180"))

# internal key -> display name
SPECIES = {
    "ecoli": "Escherichia coli",
    "ssonnei": "Shigella sonnei",
    "senterica": "Salmonella enterica subsp. I",
    "kpneumoniae": "Klebsiella pneumoniae",
}
# habitats with data per species (ALL = pooled across habitats)
HABITATS = {
    "ecoli": ["FECAL", "URINE", "ALL"],
    "ssonnei": ["FECAL"],
    "senterica": ["FECAL"],
    "kpneumoniae": ["FECAL", "URINE", "RESPIRATORY", "ALL"],
}
# Table S1 section name -> species key
SECTION_TO_SPECIES = {
    "Escherichia coli": "ecoli",
    "Shigella sonnei": "ssonnei",
    "Salmonella enterica subspecies I": "senterica",
    "Klebsiella pneumoniae": "kpneumoniae",
}

# the five protein positions of the septatype (E. coli numbering) and the
# residue expected in a wild-type, susceptible reference
POSITIONS = (("gyrA", 83), ("gyrA", 87), ("parC", 80), ("parC", 84), ("recN", 346))
WILD = {("gyrA", 83): "S", ("gyrA", 87): "D", ("parC", 80): "S", ("parC", 84): "E", ("recN", 346): "L"}
