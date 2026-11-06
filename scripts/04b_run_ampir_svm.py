import sys
from pathlib import Path
from typing import Generator, List, Tuple
import numpy as np
import joblib

#Threshold cutoff: only peptides with probability >= 0.50 will be retained
PROBABILITY_CUTOFF = 0.50
STD_AMINO_ACID = "ACDEFGHIKLMNPQRSTVWY"
DIPEPTIDES = [a + b for a in STD_AMINO_ACID for b in STD_AMINO_ACID]

def stream_fasta(file_path: Path) -> Generator[Tuple[str, str], None, None]:
    """Yields (header, sequence) pairs from FASTA file"""
    with open(file_path, "r", buffering=64 * 1024) as f:
        header = None
        seq_chunks: List[str] = []
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith(">"):
                if header:
                    yield header, "".join(seq_chunks).replace("\r", "").replace("\n", "").upper()
                header = line_str[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line_str)
        if header:
            yield header, "".join(seq_chunks).replace("\r", "").replace("\n", "").upper()

def extract_features(sequence: str) -> np.ndarray:
    """EXtracts amino acid composition (20-dim) and dipeptide frequency (400-dim)."""
    seq = sequence.upper().strip()
    seq_len = len(seq)

    # 1. Monopeptide vector (20)
    if seq_len == 0:
        mono = [0.0 * len(STD_AMINO_ACID)]
    else:
        mono = [seq.count(aa) / seq_len for aa in STD_AMINO_ACID]

    # 2. Dipepide composition
    di_counts = {dp: 0 for dp in DIPEPTIDES}
    if seq_len > 1:
        for i in range(seq_len - 1):
            pair = seq[i:i+2]
            if pair in di_counts:
               di_counts[pair] += 1
        di = [di_counts[dp] / (seq_len -1) for dp in DIPEPTIDES]
    else:
        di = [0.0] * len(DIPEPTIDES)

    return np.array(mono + di, dtype=np.float64 )

def run_ampir(input_fasta: Path, model_file: Path, output_file: Path, cutoff: float) -> None:
    if not input_fasta.exists():
        print(f"[ERROR] Input file '{input_fasta}' not found.", file=sys.stderr)
        sys.exit(1)

    if not model_file.exists():
        print(f"[ERROR] Model file '{model_file}' not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading Model B weights from: {model_file.name}")
    artifact = joblib.load(model_file)
    scaler = artifact["scaler"]
    clf = artifact["model"]
    cv_auc = artifact.get("cv_roc_aur", 0.0)

    print(f"Model loaded successfully (Benchmark ROC-AUC: {cv_auc:.4f})")
    print(f"Scoring smORFs from: {input_fasta}")
    print(f"Applying probability filter threshold: >= {cutoff}")

    passed_records: List[Tuple[str, int, float, str]] = []
    total = 0

    for seq_id, seq in stream_fasta(input_fasta):
        if not seq:
            continue
        total += 1
        feat = extract_features(seq)
        if feat is None:
            continue
        feat_scaled = scaler.transform(feat.reshape(1,-1))
        prob = clf.predict_proba(feat_scaled)[0][1]
        if prob >= cutoff:
            passed_records.append((seq_id, len(seq), round(float(prob), 4), seq))

    # Rank candidates by descending SVM probability
    passed_records.sort(key=lambda item: item[2], reverse=True)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as out:
        out.write("Accession\tLength\tModelB_SVM_Prob\tSequence\n")
        for item in passed_records:
            clean_seq = item[3].strip()
            out.write(f"{item[0]}\t{item[1]}\t{item[2]:.4f}\t{clean_seq}\n")

    print(f"Processed & scored {total} total smORF sequences.")
    print(f"Filtered leads retrianed: {len(passed_records)} (Pass rate: {(len(passed_records) / max(total, 1)) * 100:.2f}%)")
    print(f"Results saved to: {output_file}")

    # Display Top 5 predications
    print("\n--- Top 5 Model B Candidates ---")
    for i, lead in enumerate(passed_records[:5], 1):
        print(f"{i}. {lead[0]} | Prob: {lead[2]:.4f} | Length: {lead[1]} AA")

def main() -> None:
    base_dir = Path(__file__).resolve().parent.parent
    input_fasta = base_dir / "data" / "01_smorfs" / "clean_dedup_smorfs.faa"
    model_file = base_dir / "models" / "amp_svm_model.joblib"
    output_tsv = base_dir / "data" / "02_classifications" / "ampir" / "ampir_predictions.tsv"

    run_ampir(input_fasta, model_file, output_tsv, cutoff=PROBABILITY_CUTOFF)

if __name__ == "__main__":
    main()