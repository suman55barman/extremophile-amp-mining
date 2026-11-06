from pathlib import Path
import requests

DEST_DIR = Path("/home/suman_barman/projects/extremophile_amp_mining/databases/safety_benchmarks")
DEST_DIR.mkdir(parents=True, exist_ok=True)

HEADERS = {"User-Agent": "Mozilla/5.0 (Bioinformatics Pipeline; Linux x86_64)"}

def download_via_search(query: str, dest_file: Path, max_pages: int = 3) -> None:
    """Fetches up to max_pages * 500 sequences using the operational /search endpoint."""
    print(f"Querying /search for: {dest_file.name}...")
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "fasta",
        "size": 500
    }

    collected = 0
    with open(dest_file, "w", encoding="utf-8") as out_f:
        current_url = base_url
        current_params = params

        for page in range(1, max_pages + 1):
            r = requests.get(current_url, params=current_params, headers=HEADERS, timeout=30)
            if r.status_code != 200:
                print(f"[!] Warning: Page {page} returned status {r.status_code}")
                break

            text = r.text
            records = text.count(">")
            if records == 0:
                break

            out_f.write(text)
            collected += records
            print(f"    -> Page {page}: fetched {records} records (total: {collected})")

            # Check for next page in Link header
            if "next" in r.links:
                current_url = r.links["next"]["url"]
                current_params = None  # Next URL already includes pagination cursor
            else:
                break

    print(f"Successfully wrote {collected} records to {dest_file.name} ({dest_file.stat().st_size // 1024} KB)\n")


def main():
    # 1. Model C Positive: Curated Hemolytic / Cytolytic peptides across all taxa (15-90 AA)
    download_via_search(
        query="reviewed:true AND (keyword:KW-0205 OR keyword:KW-0354) AND length:[15 TO 90]",
        dest_file=DEST_DIR / "swissprot_hemolytic.faa",
        max_pages=3
    )

    # 2. Model D Positive: Curated Toxins & Venom peptides (KW-0800, 15-90 AA)
    download_via_search(
        query="reviewed:true AND keyword:KW-0800 AND length:[15 TO 90]",
        dest_file=DEST_DIR / "swissprot_toxins.faa",
        max_pages=3
    )

    # 3. Model C & D Negative: Multi-taxa neutral controls
    # Intracellular / cytoplasmic proteins with no toxin, cytolysis, or antimicrobial annotations
    download_via_search(
        query="reviewed:true NOT keyword:KW-0800 NOT keyword:KW-0205 NOT keyword:KW-0354 NOT keyword:KW-0929 AND length:[15 TO 80]",
        dest_file=DEST_DIR / "swissprot_nonhemo_controls.faa",
        max_pages=3
    )


if __name__ == "__main__":
    main()