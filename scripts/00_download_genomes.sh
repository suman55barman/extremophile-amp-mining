#!/bin/bash

# Create target directories
mkdir -p raw_data/zips raw_genomes

echo "[*] Launching automated retrieval for 32 extremophile genomes..."

while read -r acc || [ -n "$acc" ]; do
    # Trim whitespace, skip empty lines or comments
    acc=$(echo "$acc" | tr -d '\r\n ')
    [[ -z "$acc" || "$acc" =~ ^# ]] && continue

    out_zip="raw_data/zips/${acc}.zip"

    # Skip if FASTA already present
    if compgen -G "raw_genomes/${acc}*" > /dev/null; then
        echo "[=] Already present: $acc, skipping."
        continue
    fi

    echo "----------------------------------------------------"
    echo "[>] Downloading accession: $acc"
    
    datasets download genome accession "$acc" --include genome --filename "$out_zip"

    if [ -f "$out_zip" ]; then
        unzip -q -o "$out_zip" -d "raw_data/tmp_${acc}" 2>/dev/null
        find "raw_data/tmp_${acc}" -name "*_genomic.fna" -exec cp {} raw_genomes/ \; 2>/dev/null
        rm -rf "raw_data/tmp_${acc}" "$out_zip"
        echo "[✓] Extracted: $acc"
    else
        echo "[!] Warning: Download failed for $acc"
    fi

    # Small throttle to avoid NCBI rate limits
    sleep 1

done < scripts/accessions.txt

echo "===================================================="
TOTAL_FILES=$(ls -1 raw_genomes/*.fna 2>/dev/null | wc -l)
echo "[✓] Complete! FASTA assemblies ready in raw_genomes/: ${TOTAL_FILES} / 32"
