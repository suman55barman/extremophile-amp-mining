import math
import sys
from pathlib import Path
from typing import List
import joblib
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

BASE_DIR = Path("/home/suman_barman/projects/extremophile_amp_mining")
DB_DIR = BASE_DIR / "databases" / "safety_benchmarks"
MODEL_PATH = BASE_DIR / "models" / "hemolysis_svm_model.joblib"

POS_FASTA = DB_DIR / "swissprot_hemolytic.faa"
NEG_FASTA = DB_DIR / "swissprot_nonhemo_controls.faa"

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

def parse_fasta(file_path: Path) -> List[str]:
    valid_aa = set(STD_AMINO_ACID)
    sequences = []

    if not file_path.exists():
        print(f"[ERROR] Misiing File: {file_path}", file=sys.stderr)
        return sequences
    
    with open(file_path, "r", encoding="utf-8") as f:
        seq_parts = []
        for line in f:
            clean = line.strip().replace("\r", "")
            if not clean:
                continue
            if clean.startswith(">"):
                if seq_parts:
                    s = "".join(seq_parts).upper()
                    if 10 <= len(s) <= 100 and set(s).issubset(valid_aa):
                        sequences.append(s)
                seq_parts = []
            else:
                seq_parts.append(clean)

        if seq_parts:
            s = "".join(seq_parts).upper()
            if 10 <= len(s) <= 100 and set(s).issubset(valid_aa):
                sequences.append(s)

    return sequences

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
    MODEL_PATH.parent.mkdir(parents=True, exist_ok=True)

    print("-----------------------------------------------")
    print("  TRAINING: MODEL C (HEMOLYSIS / CYTOLYTIC SVM)  ")
    print("-----------------------------------------------")

    print(f"Loading positive benchamrks: {POS_FASTA.name}")
    pos_seqs = list(set(parse_fasta(POS_FASTA)))
    print(f"Loaded {len(pos_seqs)} unique positive sequences.")

    print(f"Loading negative controls: {NEG_FASTA.name}")
    neg_seqs = list(set(parse_fasta(NEG_FASTA)))
    print(f"Loaded {len(neg_seqs)} unique negative sequences.")

    if not pos_seqs or not neg_seqs:
        print("[ERROR] Training sequences could not be loaded. Check benchmark files", file=sys.stderr)
        sys.exit(1)

    # 1:1 CLass Balancing
    n_samples = min(len(pos_seqs), len(neg_seqs))
    np.random.seed(42)
    selected_pos = list(np.random.choice(pos_seqs, n_samples, replace=False))
    selected_neg = list(np.random.choice(neg_seqs, n_samples, replace=False))

    print(f"Balanced Training matrix: {n_samples * 2} samples ({n_samples} pos / {n_samples} neg)")

    print("Computing 424-D feature representations...")
    X = np.array([extract_augmented_featured(s) for s in selected_pos + selected_neg])
    y = np.array([1] * n_samples + [0] * n_samples)

    # Feature Scaling
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # RBF-SVM with Platt scaling probability calibration
    svm = SVC(kernel="rbf", C=2.0, gamma="scale", probability=True, random_state=42)

    print("Running 5-Fold Stratified Cross-validation on Hemolysis...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    scores = cross_val_score(svm, X_scaled, y, cv=cv, scoring="roc_auc", n_jobs=-1)
    mean_auc = float(np.mean(scores))
    std_auc = float(np.std(scores))
    print(f"5-Fold Cross-Validation ROC-AUC: {mean_auc:.4f} (+/- {std_auc:.4f})")

    print("Fitting final Hemolysis model...")
    svm.fit(X_scaled, y)

    artifact = {
        "model": svm,
        "scaler": scaler,
        "cv_roc_auc": mean_auc,
        "feature_dim": 424
    }
    joblib.dump(artifact, MODEL_PATH, compress=3)
    print(f"Model saved to: {MODEL_PATH}")

if __name__ == "__main__":
    main()