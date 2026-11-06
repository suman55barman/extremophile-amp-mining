import sys
import requests
from pathlib import Path
import numpy as np
import joblib
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score

STD_AMINO_ACID = "ACDEFGHIKLMNPQRSTVWY"
DIPEPTIDES = [a + b for a in STD_AMINO_ACID for b in STD_AMINO_ACID]
STANDARD_SET = set(STD_AMINO_ACID)

def extract_features(sequence: str) -> np.ndarray:
    """Extracts 20-D monopeptide and 400-D dipeptide composition vectors."""
    seq = sequence.upper()
    seq_len = len(seq)
    if seq_len == 0:
        return np.zeros(len(STD_AMINO_ACID) + len(DIPEPTIDES))

    # 20-dim monopeptide
    mono = [seq.count(aa) / seq_len for aa in STD_AMINO_ACID]

    # 400-dim dipeptide (20^2)
    di_counts = {dp: 0 for dp in DIPEPTIDES}
    if seq_len > 1:
        for i in range(seq_len - 1):
            pair = seq[i:i+2]
            if pair in di_counts:
                di_counts[pair] += 1
        di = [di_counts[dp] / (seq_len - 1) for dp in DIPEPTIDES]
    else: 
        di = [0.0] * len(DIPEPTIDES)
    return np.array(mono + di)

def load_clean_fasta(faa_path: Path):
    """Load FASTA sequences, filtering to canonical amino acids and lengths 10-100 AA."""
    seqs = []
    with open(faa_path, "r", buffering=64 * 1024) as f:
        cur = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    s = "".join(cur).upper()
                    if 10 <= len(s) <= 100 and set(s).issubset(STANDARD_SET):
                        seqs.append(s)
                cur = []
            else:
                cur.append(line)
        if cur:
            s = "".join(cur).upper()
            if 10 <= len(s) <= 100 and set(s).issubset(STANDARD_SET):
                seqs.append(s)
    return list(set(seqs))

def download_uniprot_file(url: str, dest_path: Path):
    """Stream FASTA stream directly using requests."""
    custom_headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) BioPipeline/1.0 (Python research tool)"}
    with requests.get(url, headers=custom_headers, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=64 * 1024):
                if chunk:
                    f.write(chunk)

def fetch_or_verify_benchmarks(benchmark_dir: Path):
    """
    Checks if benchmarks datasets exit locally.
    Downloads Swiss-Prot verified sets only if missing.
    """
    benchmark_dir.mkdir(parents=True, exist_ok=True)
    pos_file = benchmark_dir / "curated_amps.faa"
    neg_file = benchmark_dir / "curated_non_amps.faa"

    # Exact encoded stream URLs verified by curl
    amp_url = (
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&"
        "query=%28keyword%3AKW-0929%29+AND+%28reviewed%3Atrue%29+AND+%28length%3A%5B10+TO+100%5D%29"
    )
    non_amp_url = (
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&"
        "query=%28reviewed%3Atrue%29+AND+%28length%3A%5B10+TO+100%5D%29+NOT+%28keyword%3AKW-0929%29+NOT+%28keyword%3AKW-0800%29"
    )

    # Positive AMP dataset check / dowload
    if pos_file.exists() and pos_file.stat().st_size > 1024:
        print(f"Found existing local positive AMP database: {pos_file.name} (skipping download)")
    else:
        print(f"Downloading reviewed Swiss-Prot AMP benchmarks...")
        try:
            download_uniprot_file(amp_url, pos_file)
            print(f"Downloaded positive AMP dataset ({pos_file.stat().st_size // 1024} KB).")
        except Exception as e:
            pos_file.unlink(missing_ok=True)
            print(f"[ERROR] Failed to download positive AMPs: {e}", file=sys.stderr)
            sys.exit(1)

    # Negative non-AMP dataset check
    if neg_file.exists() and neg_file.stat().st_size > 1024:
        print(f"Found existing local negative non-AMP database: {neg_file.name} (skipping download)")
    else:
        print(f"Querying UniPort for reviewed Swiss-Prot Non-AMP controls...")
        try:
            download_uniprot_file(non_amp_url, neg_file)
            print(f"Downloaded negative Non-AMP dataset ({neg_file.stat().st_size // 1024} KB).")
        except Exception as e:
            neg_file.unlink(missing_ok=True)
            print(f"[ERROR] Failed to download negative Non-AMPs: {e}", file=sys.stderr)
            sys.exit(1)

    pos_seqs = load_clean_fasta(pos_file)
    neg_seqs = load_clean_fasta(neg_file)
    return pos_seqs, neg_seqs

def main():
    base_dir = Path(__file__).resolve().parent.parent
    benchmark_dir = base_dir / "databases" / "benchmarks"
    model_dir = base_dir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    model_file = model_dir / "amp_svm_model.joblib"

    print(f"Checking benchmark databases...")
    pos_seqs, neg_seqs = fetch_or_verify_benchmarks(benchmark_dir)

    print(
        f"\nLoaded positive AMP sequences: {len(pos_seqs)}\n"
        f"Loaded negative non-AMP sequences: {len(neg_seqs)}\n"
    )

    # Balance datasets to 1:1 ratio
    min_count = min(len(pos_seqs), len(neg_seqs))
    pos_seqs = pos_seqs[:min_count]
    neg_seqs = neg_seqs[:min_count]

    print(f"Extracting 420-D feature vectors ({min_count * 2} balanced samples)...")
    X = [extract_features(s) for s in pos_seqs] + [extract_features(s) for s in neg_seqs]
    y = [1] * min_count + [0] * min_count

    X = np.array(X)
    y = np.array(y)

    print("Normalizing features with StandardScaler...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    print("Evaluating 5-Fold Stratified Cross-Validation...")
    clf = SVC(kernel="rbf", C=2.0, gamma="scale", probability=True, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(clf, X_scaled, y, cv=cv, scoring="roc_auc", n_jobs=1)

    mean_auc = float(np.mean(cv_scores))
    std_auc = float(np.std(cv_scores))
    print(f"5-Fold Cross-Validation ROC-AUC: {mean_auc:.4f} (+/- {std_auc:.4f})")

    print("Fitting final SVM model on complete benchmark set...")
    clf.fit(X_scaled, y)

    # Serialize model artifact
    payload = {
        "scaler": scaler,
        "model": clf,
        "feature_dim": 420,
        "cv_roc_aur": mean_auc,
        "cv_std": std_auc
    }
    joblib.dump(payload, model_file, compress=3)
    print(f"Production model artifact saved to: {model_file}")

if __name__ == "__main__":
    main()