import subprocess
import sys
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = BASE_DIR / "scripts"

PIPELINE_STAGES = [
    ("Stage 1: smORFs Extraction", SCRIPTS_DIR / "01_extract_smorfs.py"),
    ("Stage 2: Sequence Deduplication & Filtering", SCRIPTS_DIR / "02_validate_dedup.py"),
    ("Stage 3: Macrel Random Forest Classification", SCRIPTS_DIR / "03_run_macrel.py"),
    ("Stage 4: Custom SVM Feature prediction", SCRIPTS_DIR / "04b_run_ampir_svm.py"),
    ("Stage 5: High-Confidence Consensus Merge", SCRIPTS_DIR / "05_consensus_merge.py"),
    ("Stage 6: Hemolysis & Cytotoxicity Safety filter", SCRIPTS_DIR / "06c_safety_filter.py"),
    ("Stage 7: DRAMP BLASTp Novelty Dereplication", SCRIPTS_DIR / "07_novelty_dereplication.py"),
    ("Stage 8: Biophysical Ranking & ESMFold 3D Modeling", SCRIPTS_DIR / "08_fold_and_biophysics.py"),
]

def run_stage(title: str, script_path: Path) -> None:
    print(f"\n{'-'*70}")
    print(f"STARTING: {title}")
    print(f"Script  : {script_path.name}")
    print(f"{'-'*70}")

    if not script_path.exists():
        print(f"[ERROR] Script not found: {script_path}", file=sys.stderr)
        sys.exit(1)

    start_time = time.time()
    result = subprocess.run([sys.executable, str(script_path)], cwd=str(BASE_DIR))
    elapsed = time.time() - start_time

    if result.returncode != 0:
        print(f"\n[FATAL] {title} failed with exit code {result.returncode}.", file=sys.stderr)
        sys.exit(result.returncode)

    print(f"COMPLETED: {title} ({elapsed:.1f}s)")

def main() -> None:
    print("-" * 70)
    print("   EXTREMOPHILE COMPUTATIONAL AMP BIOPROSPECTING ENGINE   ")
    print("-" * 70)
    total_start = time.time()

    for title, script in PIPELINE_STAGES:
        run_stage(title, script)

    total_time = time.time() - total_start
    print("\n" + "-" * 70)
    print(f"[SUCCESS] End-to-End pipeline finished in {total_time:.1f}s.")
    print("Final prioritized leads an 3D PDBs: data/05_final_leads/")
    print("-" * 70 + "\n")

if __name__ == "__main__":
    main()