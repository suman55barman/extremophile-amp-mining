import sys
import gzip
import shutil
import subprocess
from pathlib import Path

def run_macrel(input_fasta: Path, output_dir: Path, threads: int = 8):
    if shutil.which("macrel") is None:
        print("[ERROR] 'macrel' command not found in your active conda env.", file=sys.stderr)
        print("Please verify your conda env. has macrel installed.", file=sys.stderr)
        sys.exit(1)

    if output_dir.exists():
        shutil.rmtree(output_dir)

    cmd = [
        "macrel",
        "peptides",
        "--fasta", str(input_fasta.resolve()),
        "--output", str(output_dir.resolve()),
        "--tag", "macrel.out",
        "--threads", str(threads)
    ]

    print(f"Running Macrel classification on: {input_fasta}")
    print(f"Target directory: {output_dir} (Threads={threads})...")

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"[ERROR] Macrel execution failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    print(f"Macrel Classification completed successfully.")

def unpack_and_inspect_macrel_result(output_dir: Path):

    candidates = [
        f for f in output_dir.iterdir()
        if f.is_file() and "prediction" in f.name.lower() and not f.name.endswith(".md")
    ]

    if not candidates:
        print(f"Warning! Could not locate prediction file in {output_dir}")
        return

    gz_file = [f for f in candidates if f.name.endswith(".gz")]
    target_gz = gz_file[0] if gz_file else candidates[0]
    unpacked_tsv = output_dir / "macrel.out.prediction.tsv"

    if target_gz.name.endswith(".gz"):
        with gzip.open(target_gz, "rt") as src, open(unpacked_tsv, "w") as dst:
            dst.write(src.read())
        print(f"Unpacked plain TSV: {unpacked_tsv}")

    read_target = unpacked_tsv if unpacked_tsv.exists() else target_gz
    amp_count = 0
    non_hemo_amp = 0

    with open(read_target, "r") as f:
        header_line = f.readline()
        while header_line and header_line.startswith("#"):
            header_line = f.readline()

        if not header_line:
            print("Prediction file is empty.")
            return

        header = [h.strip() for h in header_line.split("\t")]

        hemo_idx = None
        for i,h in enumerate(header):
            if "hemo" in h.lower() and "prob" not in h.lower():
                hemo_idx = i
                break

        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue

            values = [v.strip() for v in line_str.split("\t")]

            amp_count += 1

            if hemo_idx is not None and len(values) > hemo_idx:
                val = values[hemo_idx].lower()
                if "non" in val:
                    non_hemo_amp += 1

    print(
        f"Total AMP candidates found: {amp_count}\n"
        f"Non-Hemo AMP candidates: {non_hemo_amp}\n"
        f"Full results saved at: {read_target}"
    )

def main():
    input_file = Path("data/01_smorfs/clean_dedup_smorfs.faa")
    output_directory = Path("data/02_classifications/macrel")

    if not input_file.exists():
        print(f"[ERROR] Input file '{input_file}' not found in current directory!", file=sys.stderr)
        sys.exit(1)

    run_macrel(input_file, output_directory, threads=8)
    unpack_and_inspect_macrel_result(output_directory)

if __name__ == "__main__":
    main()