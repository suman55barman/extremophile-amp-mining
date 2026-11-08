set -euo pipefail

# 1. Resolve canonical repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "${SCRIPT_DIR}/../data" ]]; then
    ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [[ -d "${SCRIPT_DIR}/data" ]]; then
    ROOT_DIR="${SCRIPT_DIR}"
else
    ROOT_DIR="$(pwd)"
fi

OUTPUT_DIR="${ROOT_DIR}/data/00_raw_genomes"
TEMP_DIR="${ROOT_DIR}/.tmp_ncbi_downloads"

if [[ -f "${ROOT_DIR}/contigs/accessions.txt" ]]; then
    ACCESSION_FILE="${ROOT_DIR}/contigs/accessions.txt"
elif [[ -f "${ROOT_DIR}/accessions.txt" ]]; then
    ACCESSION_FILE="${ROOT_DIR}/accessions.txt"
else
    echo "[-] Error: accessions.txt not found in ${ROOT_DIR}/contigs/ or ${ROOT_DIR}/"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}" "${TEMP_DIR}"

# 2. Read accessions
mapfile -t ACCESSIONS < <(grep -E '^[A-Za-z0-9_]+' "${ACCESSION_FILE}" | tr -d '\r')
TOTAL=${#ACCESSIONS[@]}

echo "-------------------------------------------------------------------------------"
echo " Starting NCBI Download: ${TOTAL} Extremophile Genomes"
echo " Accession File:  ${ACCESSION_FILE}"
echo " Destination:     ${OUTPUT_DIR}"
echo "-------------------------------------------------------------------------------"

USE_NCBI_CLI=false
if command -v datasets &> /dev/null; then
    USE_NCBI_CLI=true
    echo "[*] Detected 'datasets' CLI. Using native NCBI toolchain."
else
    echo "[*] 'datasets' CLI not found. Using resilient Python REST streamer."
fi

COUNT=0
SKIPPED=0
DOWNLOADED=0
FAILED=0

for acc in "${ACCESSIONS[@]}"; do
    COUNT=$((COUNT + 1))
    FINAL_FNA="${OUTPUT_DIR}/${acc}.fna"

    # Match ANY file starting with the accession prefix (e.g. GCF_..._genomic.fna or GCF_....fna)
    MATCHES=($(compgen -G "${OUTPUT_DIR}/${acc}*.fna" || true))
    if [[ ${#MATCHES[@]} -gt 0 && -s "${MATCHES[0]}" ]]; then
        echo "[${COUNT}/${TOTAL}] [SKIPPED] ${acc} already present: $(basename "${MATCHES[0]}")"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo -n "[${COUNT}/${TOTAL}] [FETCHING] ${acc} ... "

    SUCCESS=false

    # Strategy A: NCBI datasets CLI with retry
    if [[ "${USE_NCBI_CLI}" == "true" ]]; then
        ZIP_FILE="${TEMP_DIR}/${acc}.zip"
        EXTRACT_DIR="${TEMP_DIR}/${acc}_extracted"
        for attempt in 1 2 3; do
            if datasets download genome accession "${acc}" --include genome --filename "${ZIP_FILE}" &> /dev/null; then
                unzip -q -o "${ZIP_FILE}" -d "${EXTRACT_DIR}"
                find "${EXTRACT_DIR}" -type f -name "*.fna" -exec mv {} "${FINAL_FNA}" \;
                rm -rf "${ZIP_FILE}" "${EXTRACT_DIR}"
                if [[ -s "${FINAL_FNA}" ]]; then
                    SUCCESS=true
                    break
                fi
            fi
            sleep $((attempt * 2))
        done
    fi

    # Strategy B: Python REST API fallback with retry & backoff
    if [[ "${SUCCESS}" == "false" ]]; then
        if python3 - "${acc}" "${FINAL_FNA}" << 'EOF'
import sys
import time
import urllib.request
import zipfile
import io

acc = sys.argv[1]
out_fna = sys.argv[2]
api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"
req = urllib.request.Request(api_url, headers={"User-Agent": "ExtremophileAMPMiner/1.0"})

success = False
for attempt in range(1, 4):
    try:
        with urllib.request.urlopen(req, timeout=90) as resp:
            zip_bytes = resp.read()
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as z:
            fna_names = [n for n in z.namelist() if n.endswith(".fna")]
            if not fna_names:
                sys.exit(1)
            target_fna = fna_names[0]
            with z.open(target_fna) as src, open(out_fna, "wb") as dst:
                dst.write(src.read())
        success = True
        break
    except Exception:
        time.sleep(attempt * 2)

sys.exit(0 if success else 1)
EOF
        then
            if [[ -s "${FINAL_FNA}" ]]; then
                SUCCESS=true
            fi
        fi
    fi

    if [[ "${SUCCESS}" == "true" ]]; then
        SIZE=$(du -h "${FINAL_FNA}" | cut -f1)
        echo "DONE (${SIZE})"
        DOWNLOADED=$((DOWNLOADED + 1))
    else
        echo "FAILED"
        rm -f "${FINAL_FNA}"
        FAILED=$((FAILED + 1))
    fi

    sleep 0.5
done

rm -rf "${TEMP_DIR}"

echo "-------------------------------------------------------------------------------"
echo " Download Run Complete!"
echo " Total: ${TOTAL} | Downloaded: ${DOWNLOADED} | Skipped: ${SKIPPED} | Failed: ${FAILED}"
echo " Target: ${OUTPUT_DIR}"
echo "-------------------------------------------------------------------------------"