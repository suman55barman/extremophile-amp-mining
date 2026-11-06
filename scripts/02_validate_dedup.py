import sys
import shutil
import tempfile
import subprocess
from pathlib import Path

STD_AMINO_ACID = set("ACDEFGHIKLMNPQRSTVWY")

def stream_fasta(file_path: Path):
    """Yiels (header, sequence) pairs from FASTA file"""
    with open(file_path, "r", buffering=64 * 1024) as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_chunks).upper()
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            yield header, "".join(seq_chunks).upper()

def run_validated_cdhit(
        input_path: Path,
        output_path: Path,
        identity: float = 1.0,
        min_len: int = 10,
        threads: int = 8
):
    if shutil.which("cd-hit") is None:
        print("[ERROR] 'cd-hit' executable not found in PATH.", file=sys.stderr)
        sys.exit(1)

    total = 0
    passed = 0
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        for header, seq in stream_fasta(input_path):
            total += 1
            if len(seq) < min_len:
                continue
            if not set(seq).issubset(STD_AMINO_ACID):
                continue
            tmp.write(f">{header}\n{seq}\n")
            passed += 1

    print(f"Total input sequences: {total}")
    print(f"Passed canonical peptides (>= {min_len} AA): {passed}")

    if passed == 0:
        print("No sequences passed validation. Exiting.")
        tmp_path.unlink(missing_ok=True)
        return
        
    cmd = [
        "cd-hit",
        "-i", str(tmp_path),
        "-o", str(output_path),
        "-c", str(identity),
        "-n", "4" if identity >= 0.7 else "2",
        "-l", str(min_len),
        "-s", "1.0",
        "-T", str(threads),
        "-M", "0",
        "-d", "0"
    ]
    
    print(f"Initializing CD-HIT stream (-c {identity}, -l {min_len}, -n 4, threads={threads})...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    tmp_path.unlink(missing_ok=True)

    if result.returncode != 0:
        print(f"[ERROR] CD-HIT failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    print(f"Output written to: {output_path}")

def main():
    input_file = Path("data/01_smorfs/extracted_smorfs.faa")
    output_file = Path("data/01_smorfs/clean_dedup_smorfs.faa")

    if not input_file.exists():
        print(f"[ERROR] Input file '{input_file}' not found in current directory!", file=sys.stderr)
        sys.exit(1)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    run_validated_cdhit(input_file, output_file, identity=1.0, min_len=10, threads=8)

if __name__ == "__main__":
    main()