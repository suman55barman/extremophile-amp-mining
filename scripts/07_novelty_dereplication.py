import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import requests

BASE_DIR = Path("/home/suman_barman/projects/extremophile_amp_mining")
INPUT_FASTA = BASE_DIR / "data" / "03_safety" / "safety_filtered_amps.faa"
DB_DIR = BASE_DIR / "databases" / "benchmarks"
DRAMP_FASTA = DB_DIR / "dramp_amps.faa"
DB_PREFIX = DB_DIR / "dramp_db"

OUTPUT_DIR = BASE_DIR / "data" / "04_novelty"
OUTPUT_REPORT = OUTPUT_DIR / "dereplication_report.tsv"
OUTPUT_NOVEL_FASTA = OUTPUT_DIR / "novel_safe_amps.faa"

DRAMP_URL = "https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta"

STD_AMINO_ACID = set("ACDEFGHIKLMNPQRSTVWY")

# Threshold
EVALUE_CUTOFF = 1e-3
HOMOLOG_ID_THRESHOLD = 80.0
NOVEL_ID_THRESHOLD = 40.0

def sanitize_dramp_fasta(raw_path: Path, clean_path: Path) -> int:
    """Filters out synthetic peptidomimetics, non-canonical residues, and broken lines."""
    kept = 0
    with open(raw_path, "r", encoding="utf-8", errors="ignore") as f_in, \
         open (clean_path, "w", encoding="utf-8") as f_out:
        curr_header = None
        seq_parts = []
        for line in f_in:
            clean = line.strip().replace("\r", "")
            if not clean:
                continue
            if clean.startswith(">"):
                if curr_header and seq_parts:
                    seq = "".join(seq_parts).upper()
                    if seq and set(seq).issubset(STD_AMINO_ACID):
                        f_out.write(f"{curr_header}\n{seq}\n")
                        kept += 1
                curr_header = clean
                seq_parts = []
            else:
                seq_parts.append(clean)
        if curr_header and seq_parts:
            seq = "".join(seq_parts).upper()
            if seq and set(seq).issubset(STD_AMINO_ACID):
                f_out.write(f"{curr_header}\n{seq}\n")
                kept += 1
    return kept

def ensure_blast_database() -> None:
    """Verifies BLAST db files; downloads DRAMP and formats index if missing."""
    DB_DIR.mkdir(parents=True, exist_ok=True)

    # Check if BLAST db index files exist (.pin, .phr, .psq)
    blast_index_files = list(DB_DIR.glob(f"{DB_PREFIX.name}.p*"))
    if blast_index_files:
        print(f"Found existing BLAST database: {DB_PREFIX}")
        return

    print("BLAST database not found. Checking FASTA source...")
    if not DRAMP_FASTA.exists() or DRAMP_FASTA.stat().st_size == 0:
        print(f"Downloading DRAMP dataset from:\n    {DRAMP_URL}")
        try:
            header = {"User-Agent": "Mozilla/5.0 (Bioinfromatics Harvester; Linux x86_64)"}
            r = requests.get(DRAMP_URL, headers=header, timeout=60, allow_redirects=True)
            r.raise_for_status()
            with open(DRAMP_FASTA, "w", endcoding="utf-8") as f:
                f.write(r.text)
            print(f"Download complete: {DRAMP_FASTA.name} ({DRAMP_FASTA.stat().st_size // 1024} KB)")
        except requests.exceptions.RequestException as e:
            print(f"[ERROR] Failed to download DRAMP: {e}", file=sys.stderr)
            sys.exit(1)

    # Sanitize FASTA to remove synthetic/modified entries
    print("Sanitize DRAMP FASTA")
    clean_fasta = DB_DIR / "dramp_amps_clean.faa"
    kept_count = sanitize_dramp_fasta(DRAMP_FASTA, clean_fasta)
    clean_fasta.replace(DRAMP_FASTA)
    print(f"Sanitize refernce set: {kept_count} natural canonical AMPs retained.")
    print(f"Formatting BLAST database with makeblastdb...")

    # Remove any broken partial database files
    for stale in DB_DIR.glob(f"{DB_PREFIX.name}.p*"):
        stale.unlink(missing_ok=True)

    print("Building BLAST database via makeblastdb...")
    cmd = [
        "makeblastdb",
        "-in", str(DRAMP_FASTA),
        "-dbtype", "prot",
        "-out", str(DB_PREFIX)
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"BLAST database created at: {DB_PREFIX}\n")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"[ERROR] makeblastdb failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def parse_fasta(file_path: Path) -> List[Tuple[str, str]]:
    records = []
    if not file_path.exists():
        return records
    with open(file_path, "r", encoding="utf-8") as f:
        header, seq_parts = None, []
        for line in f:
            clean = line.strip().replace("\r", "")
            if not clean:
                continue
            if clean.startswith(">"):
                if header and seq_parts:
                    records.append((header, "".join(seq_parts).upper()))
                header = clean[1:]
                seq_parts = []
            else:
                seq_parts.append(clean)
        if header and seq_parts:
            records.append((header, "".join(seq_parts).upper()))
    return records

def run_blastp(query_fasta: Path, db_path: Path) -> Dict[str, dict]:
    """Runs local BLASTp and maps best target hit per query ID."""
    cmd = [
        "blastp",
        "-query", str(query_fasta),
        "-db", str(db_path),
        "-evalue", str(EVALUE_CUTOFF),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
        "-max_target_seqs", "1"
    ]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except FileNotFoundError:
        print(f"[ERROR] blastp binary not found in active conda environment.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] BLASTp search failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)

    hits={}
    for line in res.stdout.strip().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        qseqid = parts[0]
        hits[qseqid] = {
            "sseqid": parts[1],
            "pident": float(parts[2]),
            "length": int(parts[3]),
            "evalue": float(parts[4]),
            "bitscore": float(parts[5])
        }
    return hits

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("-----------------------------------------")
    print("    NOVELTY & DEREPLICATION SCREENING    ")
    print("-----------------------------------------")

    # 1. Ensure database is built and ready
    ensure_blast_database()

    # 2. Ingest safe candidates
    if not INPUT_FASTA.exists():
        print(f"[ERROR] Missing safe leads input: {INPUT_FASTA}", file=sys.stderr)
        sys.exit(1)

    candidates = parse_fasta(INPUT_FASTA)
    print(f"Loaded {len(candidates)} safe peptides from: {INPUT_FASTA.name}")

    # 3. Query local DRAMP database
    print("Executing BLASTp search against DRAMP database...")
    top_hits = run_blastp(INPUT_FASTA, DB_PREFIX)

    report_lines = []
    novel_peptides = []

    for header, seq in candidates:
        clean_id = header.split()[0]
        hit = top_hits.get(clean_id)

        if hit is None:
            pident = 0.0
            evalue = 1.0
            top_hit_id = "None"
            status = "NOVEL"
        else:
            pident = hit["pident"]
            evalue = hit["evalue"]
            top_hit_id = hit["sseqid"]
            if pident >= HOMOLOG_ID_THRESHOLD:
                status = "KNOWN_HOMOLOG"
            elif pident >= NOVEL_ID_THRESHOLD:
                status = "DISTANT_HOMOLOG"
            else:
                status = "NOVEL"

        report_lines.append(f"{clean_id}\t{len(seq)}\t{top_hit_id}\t{pident:.2f}\t{evalue:.2e}\t{status}\t{seq}")

        if status == "NOVEL":
            novel_peptides.append((header, seq))

    # 1. Write dereplication TSV report
    with open(OUTPUT_REPORT, "w", encoding="utf-8") as f:
        f.write("Candidate_ID\tLength\tTop_DRAMP_hit\tPercentage_Identity\tE_Value\tNovelty_Status\tSequence\n")
        f.write("\n".join(report_lines) + "\n")

    # 2. Write novel peptides FASTA
    with open(OUTPUT_NOVEL_FASTA, "w", encoding="utf-8") as f:
        for h, s in novel_peptides:
            f.write(f">{h}\n{s}\n")

    print("-----------------------------------------")
    print("          DEREPLICATION SUMMARY          ")
    print("-----------------------------------------")
    print(
        f"Total peptides analyzed: {len(candidates)}\n"
        f"Known homologs (>=80%) : {sum('KNOWN_HOMOLOG' in line for line in report_lines)}\n"
        f"Distant Homologs       : {sum('DISTANT_HOMOLOG' in line for line in report_lines)}\n"
        f"High-Confidence Novel  : {len(novel_peptides)} ({(len(novel_peptides)/len(candidates))*100:.1f}%)\n"
        f"Dereplication Report   : {OUTPUT_REPORT}\n"
        f"Novel peptides FASTA   : {OUTPUT_NOVEL_FASTA}\n"
    )

if __name__ == "__main__":
    main()