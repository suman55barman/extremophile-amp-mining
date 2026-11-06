import math
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple
import requests

BASE_DIR = Path("/home/suman_barman/projects/extremophile_amp_mining")
INPUT_FASTA = BASE_DIR / "data" / "04_novelty" / "novel_safe_amps.faa"

OUTPUT_DIR = BASE_DIR / "data" / "05_final_leads"
PDB_DIR = OUTPUT_DIR / "pdbs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PDB_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_REPORT = OUTPUT_DIR / "prioritized_leads_report.tsv"
OUTPUT_TOP_FASTA = OUTPUT_DIR / "top_candidate_leads.faa"

TOP_N_EXPORT = 5
ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

EISENBERG = {
    'I': 0.73, 'F': 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37,
    'M': 0.26, 'A': 0.25, 'G': 0.16, 'C': 0.04, 'Y': 0.02,
    'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40, 'E': -0.62,
    'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.76
}

BOMAN = {
    'L': -4.92, 'I': -4.92, 'V': -4.04, 'F': -2.98, 'M': -2.35,
    'W': -2.33, 'A': -1.81, 'C': -1.28, 'G': -0.94, 'Y': -0.14,
    'T': 2.57, 'S': 3.40, 'H': 4.66, 'Q': 5.54, 'K': 5.55,
    'N': 6.64, 'E': 6.81, 'D': 8.72, 'R': 9.38, 'P': 0.00
}

PKA_TABLE = {
    'N_TERM': 9.69, 'C_TERM': 2.34,
    'K': 10.53, 'R': 12.48, 'H':  6.00,
    'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10.07
}

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

def compute_net_charge(seq: str, ph: float = 7.4,) -> float:
    """Calculate formal net charge at pH 7.4 via Henderon-Hasselbalch."""
    charge = 1.0 / (1.0 + 10.0 ** (ph - PKA_TABLE['N_TERM']))
    charge -= 1.0 / (1.0 + 10.0 ** (PKA_TABLE['C_TERM'] - ph))

    for aa in seq:
        if aa in ('K', 'R', 'H'):
            pka = PKA_TABLE[aa]
            charge += 1.0 / (1.0 + 10.0 ** (ph - pka))
        elif aa in ('D', 'E', 'C', 'Y'):
            pka = PKA_TABLE[aa]
            charge -= 1.0 / (1.0 + 10.0 ** (pka - ph))
    return charge

def extract_consensus_prob(header: str) -> float:
    """Extracts upstream ML probability if present in header, else defaults to 0.5"""
    match = re.search(r"Prob=([0-9\.]+)", header)
    return float(match.group(1)) if match else 0.5

def compute_biophysics(header: str, seq: str) -> dict:
    """Compute hydrophobicity, Eisenberg hydrophobic moment, and Boman index."""
    seq_len = len(seq)
    mean_h = sum(EISENBERG.get(aa, 0.0) for aa in seq) / seq_len
    boman = sum(BOMAN.get(aa, 0.0) for aa in seq) / seq_len

    # Eisenberg hydrophobic moment for alpha-helix (delta = 100 deg)
    delta = 100.0 * (math.pi / 180.0)
    sin_sum = sum(EISENBERG.get(aa, 0.0) * math.sin(i * delta) for i, aa in enumerate(seq))
    cos_sum = sum(EISENBERG.get(aa, 0.0) * math.cos(i * delta) for i, aa in enumerate(seq))
    u_h = math.sqrt(sin_sum**2 + cos_sum**2) / seq_len

    net_charge = compute_net_charge(seq)
    prob_score = extract_consensus_prob(header)

    # Multi-parametric Lead Priority score (0 - 100)
    score = 0.0 
    if 2.0 <= net_charge <= 6.0:
        score += 35.0
    elif 1.0 <= net_charge < 2.0 or 6.0 < net_charge <= 8.0:
        score += 20.0
    elif net_charge > 0:
        score += 10.0

    if u_h >= 0.35:
        score += 35.0
    elif u_h >= 0.20:
        score += 20.0

    if -1.0 <= boman <= 2.5:
        score += 15.0

    if 15 <= seq_len <= 35:
        score += 15.0
    elif 10 <= seq_len < 15 or 35 < seq_len <= 50:
        score += 8.0

    return {
        "Length": seq_len,
        "net_charge": net_charge,
        "mean_hydrophobicity": mean_h,
        "hydrophobic_moment": u_h,
        "boman_index": boman,
        "consensus_prob": prob_score,
        "priority_score": score
    }

def query_esmfold_api(sequence: str, pdb_out_path: Path) -> float:
    """Sends sequence to Meta ESMFold API, writes PDB, and Calculatrs mean pLDDT."""
    header = {"Content-Type": "text/plain", "User-Agent": "Mozilla/5.0 Bioinformatics-Pipeline"}
    try:
        response = requests.post(ESMFOLD_API_URL, headers=header, data=sequence, timeout=90)
        if response.status_code != 200:
            print(f"   [WARN] ESMFold API returned status code {response.status_code}")
            return 0.0

        pdb_content = response.text
        with open(pdb_out_path, "w", encoding="utf-8") as f:
            f.write(pdb_content)

        # In ESMFold PDBs, B-factor column contains the pLDDT confidence (0-100)
        plddt_vals = []
        for line in pdb_content.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    b_factor = float(line[60:66].strip())
                    plddt_vals.append(b_factor)
                except ValueError:
                    continue

        return sum(plddt_vals) / len(plddt_vals) if plddt_vals else 0.0

    except Exception as e:
        print(f"   [WARN] ESMFold folding request failed: {e}")
        return 0.0

def main() -> None:
    print("------------------------------------------")
    print(" Biophysical, Ranking & ESMFold profiling ")
    print("------------------------------------------")

    if not INPUT_FASTA.exists():
        print(f"[ERROR] Missing Input: {INPUT_FASTA}", file=sys.stderr)
        sys.exit(1)

    candidates = parse_fasta(INPUT_FASTA)
    print(f"Ingested {len(candidates)} novel peptides from: {INPUT_FASTA.name}")

    evaluated = []
    for header, seq in candidates:
        metrics = compute_biophysics(header, seq)
        evaluated.append((header, seq, metrics))

    # Strict Multi-Tier Sorting
    evaluated.sort(
        key=lambda x: (
            x[2]["priority_score"],
            x[2]["hydrophobic_moment"],
            -abs(x[2]["net_charge"] - 3.5),
            x[2]["consensus_prob"]
        ),
        reverse=True
    )

    print("Candidates successfully sorted across multi-parameter biophysical hierarchy.")

    # Clean old PDBs to prevent mixed ranks
    for old_pdb in PDB_DIR.glob("rank*.pdb"):
        old_pdb.unlink(missing_ok=True)

    # 1. Fold Top N candidates using Meta ESMFold API
    print(f"\n Quering Meta ESMFold API for TOP {TOP_N_EXPORT} shortlisted leads...")
    top_leads = evaluated[:TOP_N_EXPORT]

    for rank, (header, seq, m) in enumerate(top_leads, start=1):
        clean_id = header.split()[0].replace("|","_")
        pdb_path = PDB_DIR / f"rank{rank}_{clean_id}.pdb"
        print(f" -> Folding [{rank}/{TOP_N_EXPORT}] ({clean_id[:30]}..., {len(seq)} AA...)")
        mean_plddt = query_esmfold_api(seq, pdb_path)
        m["mean_plddt"] = mean_plddt
        if mean_plddt > 0:
            print(f"   Saved: {pdb_path.name} | Mean pLDDT: {mean_plddt:.2f}")
        time.sleep(1.0) # Courtesy pause for public API endpoint

    # 2. Write Master Prioritized Report
    with open(OUTPUT_REPORT, "w", encoding="utf-8") as f:
        f.write("Rank\tCandidate_ID\tLength\tNet_Charge_pH7.4\tMean_Hydrophobicity\tHydrophobic_Moment\tBoman_Index\tConsensus_Prob\tPriority_Score\tMean_pLDDT\tSequence\n")
        for rank, (header, seq, m) in enumerate(evaluated, start=1):
            clean_id = header.split()[0]
            plddt_str = f"{m.get('mean_plddt', 0.0):.2f}" if "mean_plddt" in m else "N/A"
            f.write(
                f"{rank}\t{clean_id}\t{m['Length']}\t{m['net_charge']:.2f}\t"
                f"{m['mean_hydrophobicity']:.3f}\t{m['hydrophobic_moment']:.3f}\t"
                f"{m['boman_index']:.2f}\t{m['consensus_prob']:.4f}\t{m['priority_score']:.1f}\t"
                f"{plddt_str}\t{seq}\n"
            )

     # 3. Export Top N FASTA
    with open(OUTPUT_TOP_FASTA, "w", encoding="utf-8") as f:
        for rank, (header, seq, m) in enumerate(top_leads, start = 1):
            plddt_val = m.get('mean_plddt', 0.0)
            annotated = (
                f">{header}|Rank={rank}|Score={m['priority_score']:.1f}"
                f"|Charge={m['net_charge']:.2f}|uH={m['hydrophobic_moment']:.2f}|pLDDT={plddt_val:.1f}"
            )
            f.write(f"{annotated}\n{seq}\n")

    print("-------------------------------------------")
    print("           TOP SHORTLISTED LEADS           ")
    print("-------------------------------------------")
    for rank, (header, seq, m) in enumerate(top_leads, start=1):
        clean_id = header.split()[0][:35]
        plddt_disp = f"{m.get('mean_plddt', 0.0):.2f}" if m.get('mean_plddt', 0.0) > 0 else "Pending"
        print(f"[{rank}] {clean_id}...")
        print(f"   Length: {m['Length']} AA | Charge: +{m['net_charge']:.2f} | <uH>: {m['hydrophobic_moment']:.3f}")
        print(f"   Priority score: {m['priority_score']:.1f}/100 | Mean pLDDT: {plddt_disp}")
        print(f"   Sequence: {seq}\n")

    print(f"Prioritizatin Report : {OUTPUT_REPORT}")
    print(f"Top Leads FASTA      : {OUTPUT_TOP_FASTA}")
    print(f"PDB Structure folder : {PDB_DIR}\n")

if __name__ == "__main__":
    main()