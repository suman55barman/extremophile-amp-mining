import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import joblib
import numpy as np

BASE_DIR = Path("/home/suman_barman/projects/extremophile_amp_mining")
DATA_DIR = BASE_DIR / "data"
MODELS_DIR = BASE_DIR / "models"

HEMO_MODEL_PATH = MODELS_DIR / "hemolysis_svm_model.joblib"
TOXIC_MODEL_PATH = MODELS_DIR / "toxicity_svm_model.joblib"

# Input candidates: prefers prefers consensus output if available, falls back to validated smORFs
CONSENSUS_INPUT = DATA_DIR / "02_classifications" / "consensus" / "consensus_leads.faa"
FALLBACK_INPUT = DATA_DIR / "01_smorfs" / "clean_dedup_smorfs.faa"

OUTPUT_FASTA = DATA_DIR / "03_safety" / "safety_filtered_amps.faa"
OUTPUT_REPORT = DATA_DIR / "03_safety" / "safety_screening_report.tsv"

# Safety thresholds: Probabilities strictly <= threshold, are deemed "SAFE"
HEMO_THRESHOLD = 0.50
TOXIC_THRESHOLD = 0.50

STD_AMINO_ACID = "ACDEFGHIKLMNPQRSTVWY"
DIPEPTIDES = [a + b for a in STD_AMINO_ACID for b in STD_AMINO_ACID]

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

def parse_fasta(file_path: Path) -> List[Tuple[str, str]]:
    valid_aa = set(STD_AMINO_ACID)
    records = []

    if not file_path.exists():
        print(f"[ERROR] Misiing File: {file_path}", file=sys.stderr)
        return records
    
    with open(file_path, "r", encoding="utf-8") as f:
        curr_header = None
        seq_parts = []
        for line in f:
            clean = line.strip().replace("\r", "")
            if not clean:
                continue
            if clean.startswith(">"):
                if curr_header and seq_parts:
                    s = "".join(seq_parts).upper()
                    if 10 <= len(s) <= 100 and set(s).issubset(valid_aa):
                        records.append((curr_header, s))
                curr_header = clean[1:]
                seq_parts = []
            else:
                seq_parts.append(clean)

        if curr_header and seq_parts:
            s = "".join(seq_parts).upper()
            if 10 <= len(s) <= 100 and set(s).issubset(valid_aa):
                records.append((curr_header, s))

    return records

def extract_augmented_featured(seq: str) -> np.ndarray:
    seq_len = len(seq)
    mono = [seq.count(aa) / seq_len for aa in STD_AMINO_ACID]

    di_counts = {dp: 0 for dp in DIPEPTIDES}
    for i in range(seq_len - 1):
        pair = seq[i:i+2]
        if pair in di_counts:
            di_counts[pair] += 1
    di = [di_counts[dp] / (seq_len - 1) for dp in DIPEPTIDES]

    mean_h = sum(EISENBERG.get(aa, 0.0) for aa in seq) / seq_len
    boman = sum(BOMAN.get(aa, 0.0) for aa in seq) / seq_len

    delta = 100.0 * (math.pi / 180.0)
    sin_sum = sum(EISENBERG.get(aa, 0.0) * math.sin(i * delta) for i, aa in enumerate(seq))
    cos_sum = sum(EISENBERG.get(aa, 0.0) * math.cos(i * delta) for i, aa in enumerate(seq))
    u_h = math.sqrt(sin_sum**2 + cos_sum**2) / seq_len

    pos_density = (seq.count("K") + seq.count("R")) / seq_len
    descriptors = [mean_h, u_h, boman, pos_density]

    return np.array(mono + di + descriptors, dtype=np.float64)

def main() -> None:
    print("-----------------------------------------")
    print("  Hemolysis & Toxicity safety screening  ")
    print("-----------------------------------------")

    # Locate Inputs
    input_file = CONSENSUS_INPUT if CONSENSUS_INPUT.exists() else FALLBACK_INPUT
    if not input_file.exists():
        print(f"[ERROR] No candidate input found. Missing {CONSENSUS_INPUT} and {FALLBACK_INPUT}.", file=sys.stderr)
        sys.exit(1)

    if not HEMO_MODEL_PATH.exists() or not TOXIC_MODEL_PATH.exists():
        print(f"[ERROR] Trained model artifacts missing in {MODELS_DIR}.", file=sys.stderr)
        sys.exit(1)

    print(f"Ingesting candidate smORFs from: {input_file.name}")
    candidates = parse_fasta(input_file)
    print(f"Loaded {len(candidates)} candidate peptides for safety evaluation.")

    if not candidates:
        print("No sequences to screen. Exiting", file=sys.stderr)
        sys.exit(1)

    # Load serialized artifacts
    print("Loading Model C (Hemolysis SVM)...")
    hemo_pack = joblib.load(HEMO_MODEL_PATH)
    hemo_model = hemo_pack["model"]
    hemo_scaler = hemo_pack["scaler"]

    print("Loading Model D (Toxicity SVM)...")
    toxic_pack = joblib.load(TOXIC_MODEL_PATH)
    toxic_model = toxic_pack["model"]
    toxic_scaler = toxic_pack["scaler"]

    # Compute features
    print("Extracting 424-D representations across candidate smORFs...")
    X_raw = np.array([extract_augmented_featured(s) for _, s in candidates])

    # Model C predictions
    X_hemo_scaled = hemo_scaler.transform(X_raw)
    p_hemo = hemo_model.predict_proba(X_hemo_scaled)[:, 1]

    # Model C predictions
    X_toxic_scaled = toxic_scaler.transform(X_raw)
    p_toxic = toxic_model.predict_proba(X_toxic_scaled)[:, 1]

    # Filter and compile report
    safe_candidates = []
    report_rows = []

    for idx, (header, seq) in enumerate(candidates):
        h_prob = float(p_hemo[idx])
        t_prob = float(p_toxic[idx])

        is_safe = (h_prob <= HEMO_THRESHOLD) and (t_prob <= TOXIC_THRESHOLD)
        status = "SAFE" if is_safe else "FLAGGED"

        report_rows.append(
            f"{header}\t{len(seq)}\t{h_prob:.4f}\t{t_prob:.4f}\t{status}\t{seq}"
        )

        if is_safe:
            annotated_header = f">{header}|p_hemo={h_prob:.3f}|p_toxic={t_prob:.3f}"
            safe_candidates.append((annotated_header, seq))

    # Write TSV report
    with open(OUTPUT_REPORT, "w", encoding="utf-8") as f:
        f.write("Sequence_ID\tLength\tP_Hemolysis\tP_Toxicity\tSafety_status\tSequence\n")
        f.write("\n".join(report_rows) + "\n")

    # Write Safe FASTA
    with open(OUTPUT_FASTA, "w", encoding="utf-8") as f:
        for h, s in safe_candidates:
            f.write(f"{h}\n{s}\n")

    print("-----------------------------------------")
    print("        SAFETY SCREENING SUMMARY         ")
    print("-----------------------------------------")
    print(
        f"Total candiates screened  : {len(candidates)}\n"
        f"Passed both safety filters: {len(safe_candidates)} ({(len(safe_candidates) / len(candidates)) * 100:.1f}%)\n"
        f"Flagged for hemolysis     : {sum(p > HEMO_THRESHOLD for p in p_hemo)}\n"
        f"Output FASTA saved to     : {OUTPUT_FASTA}\n"
        f"TSV report saved to       : {OUTPUT_REPORT}\n"
    )

if __name__ == "__main__":
    main()