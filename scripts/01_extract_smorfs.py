import glob
from collections import Counter
import math
import os
import re
import pyrodigal
from Bio import SeqIO

INPUT_DIR = "data/00_raw_genomes"
OUTPUT_FILE = "data/01_smorfs/extracted_smorfs.faa"
MIN_TRAINING_LENGTH = 1000000

CANONICAL_STARTS = {"ATG", "GTG", "TTG"}

def calculate_shannon_entropy(seq_str):
    if not seq_str:
        return 0.0
    length = len(seq_str)
    counts = Counter(seq_str)
    return -sum(
        (c / length) * math.log2(c / length) for c in counts.values() if c > 0
    )

total_extracted = 0
genome_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.fna")))

if not genome_files:
    print(f"Error! No .fna files found in directory '{INPUT_DIR}'.")
    exit(1)

print(f"Found {len(genome_files)} genome assemblies. Starting Extraction...\n")

with open(OUTPUT_FILE, "w") as out_faa:
    for fna_path in genome_files:
        genome_id = os.path.basename(fna_path).replace(".fna", "")
        records = list(SeqIO.parse(fna_path, "fasta"))

        if not records:
            print(f"Skipping empty file: {fna_path}")
            continue

        combined_seq = b"".join(bytes(r.seq) for r in records)
        if len(combined_seq) < MIN_TRAINING_LENGTH:
            print(f"Skipping {genome_id}: assembly size ({len(combined_seq)} bp) too small for training.")
            continue

        finder = pyrodigal.GeneFinder(min_gene=30, max_overlap=15)
        finder.train(combined_seq)

        genome_smorfs = 0

        for record in records:
            dna_bytes = bytes(record.seq)
            if len(dna_bytes) < 100:
                continue

            genes = finder.find_genes(dna_bytes)

            gene_coords = [(g.begin, g.end, g.strand, g.score) for g in genes]

            for i, gene in enumerate(genes):
                pep = gene.translate().rstrip("*")

                # 1. Microprotein length constraint
                if not (10 <= len(pep) <= 100 and "*" not in pep):
                    continue

                # 2.Shannon Entropy check: discard low-complexity tandem repeats (e.g. poly-L, poly-A)
                entropy = calculate_shannon_entropy(pep)
                min_entropy_threshold = 2.30 if len(pep) <= 30 else 2.7
                if entropy < min_entropy_threshold:
                    continue

                # 3. Initiation codon sanity check
                if gene.strand == 1:
                    start_codon = str(record.seq[gene.begin - 1 : gene.begin + 2]).upper()
                else:
                    start_codon = str(record.seq[gene.end - 3 : gene.end].reverse_complement()).upper()
                if start_codon not in CANONICAL_STARTS:
                    continue

                # 4. Thermodynamic RBS / Initiation Gate
                has_rbs = gene.rbs_motif is not None and str(gene.rbs_motif).strip() not in ("None", "")
                spacer_str = str(gene.rbs_spacer) if gene.rbs_spacer is not None else ""
                spacer_nums = [int(n) for n in re.findall(r"\d+", spacer_str)] if spacer_str else []
                valid_spacer = bool(spacer_nums and min(spacer_nums) >= 3 and max(spacer_nums) <= 15)
                if not ((has_rbs and valid_spacer) or gene.score >= 8.0):
                    continue

                # 5. Shadow ORF filter: Check for heavy antisense overlaps with higher-scoring genes
                is_shadow = False
                for other_begin, other_end, other_strand, other_score in gene_coords:
                    if other_strand != gene.strand and other_score > gene.score:
                        overlap = max(0, min(gene.end, other_end) - max(gene.begin, other_begin))
                        gene_len = gene.end - gene.begin
                        if gene_len > 0 and (overlap / gene_len) > 0.60:
                            is_shadow = True
                            break
                        if is_shadow:
                            continue

                rbs_str = gene.rbs_motif if has_rbs else "None"
                final_spacer_str = spacer_str if valid_spacer else "None"

                header = (
                    f">{genome_id}_{record.id}_{i+1} "
                    f"rbs={rbs_str} "
                    f"spacer={final_spacer_str} "
                    f"score={gene.score:.2f} "
                    f"start_codon={start_codon} "
                    f"entropy={entropy:.2f} "
                    f"start={gene.begin} "
                    f"end={gene.end} "
                    f"strand={gene.strand} "
                    f"len={len(pep)}"
                )

                out_faa.write(f"{header}\n{pep}\n")
                genome_smorfs += 1
                total_extracted += 1

        print(f"[{genome_id}] Retained {genome_smorfs} Ultra-refined smORFs")

print(f"\nStage 1 Complete! Extracted a total of {total_extracted} Ultra-refined smORFs into '{OUTPUT_FILE}'.")