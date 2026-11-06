import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

def normalize_seq(seq: str) -> str:
    """Strips whitespace, converts to uppercase, and standardizes leading M offset."""
    s = seq.strip().upper()
    return s[1:] if s.startswith("M") else s

def load_macrel_results(macrel_file: Path) -> Dict[str, dict]:
    """
    Loads macrel prediction output.
    Extracts accession, Macrel AMP classification, and cleaved sequences.
    """
    macrel_leads = {}
    if not macrel_file.exists():
        print(f"[ERROR] Macrel output not found at: {macrel_file}", file=sys.stderr)
        sys.exit(1)

    with open(macrel_file, "r", encoding="utf-8") as f:
        header_line = None
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue

            # Header row Check
            if line_str.startswith("Access") or line_str.startswith('# Access'):
                line_clean = line_str.lstrip('#').strip()
                header_line = line_clean.split("\t")
                if len(header_line) == 1: 
                    header_line = line_clean.split()
                header_line = [c.lower() for c in header_line]
                continue

            cols = line_str.split("\t")
            if len(cols) == 1:
                cols = line_str.split()

            if header_line:
                row = dict(zip(header_line, cols))
                acc = row.get("access") or row.get("accession") or cols[0]
                seq = row.get("sequence") or cols[1]
                amp_family = row.get("amp_family") or (cols[2] if len(cols) > 2 else "AMP")
                hemo = row.get("hemolytic") or (cols[4] if len(cols) > 4 else "N/A")
            else:
                acc = cols[0]
                seq = cols[1]
                amp_family = cols[2] if len(cols) > 2 else "AMP"
                hemo = cols[4] if len(cols) > 4 else "N/A"

            acc = acc.split()[0]
            seq = seq.strip().upper()

            macrel_leads[acc] = {
                "macrel_class": amp_family,
                "hemolytic": hemo,
                "mature_seq": seq,
                "norm_seq": normalize_seq(seq)
            } 

    return macrel_leads

def load_modelb_results(svm_file: Path) -> Dict[str, dict]:
    """Load Model B SVM predictions with full precursor sequences."""
    svm_leads = {}
    if not svm_file.exists():
        print(f"[ERROR] Model B output not found at: {svm_file}", file=sys.stderr)
        sys.exit(1)

    with open(svm_file, "r", encoding="utf-8") as f:
        header_line = None
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            cols = line_str.split("\t")
            if header_line is None:
                header_line = cols
                continue

            acc = cols[0].split()[0]
            length = int(cols[1])
            prob = float(cols[2])
            raw_seq = cols[3].strip().upper()

            svm_leads[acc] = {
                "length": length,
                "svm_prob": prob,
                "precursor_seq": raw_seq,
                "norm_seq": normalize_seq(raw_seq)
            }
    return svm_leads

def merge_and_tier_leads(
        macrel_leads: Dict[str, dict],
        svm_leads: Dict[str, dict],
        output_tsv: Path,
        output_faa: Path
) -> None:
    tier_1 = []
    tier_2 = []
    tier_3 = []

    # Cross-reference SVM leads against Macrel
    for acc, s_data in svm_leads.items():
        prob = s_data["svm_prob"]
        precursor = s_data["precursor_seq"]

        is_macrel = False
        mac_class = "N/A"
        mature_seq = precursor
        cleaved = False

        if acc in macrel_leads:
            is_macrel = True
            mac_class = macrel_leads[acc]["macrel_class"]
            mature_seq = macrel_leads[acc]["mature_seq"]
            cleaved = (mature_seq != precursor)
        else:
            # Check normalized match
            for m_acc, m_data in macrel_leads.items():
                if s_data["norm_seq"] == m_data["norm_seq"]:
                    is_macrel = True
                    mac_class = m_data["macrel_class"]
                    mature_seq = m_data["macrel_seq"]
                    cleaved = (mature_seq != precursor)
                    break

        if is_macrel:
            tier_1.append({
                "tier": "Tier_1_Consensus",
                "acc": acc,
                "len": len(mature_seq),
                "prob": prob,
                "macrel_class": mac_class,
                "cleaved": cleaved,
                "mature_seq": mature_seq,
                "precursor_seq": precursor
            })
        elif prob >= 0.85:
            tier_2.append({
                "tier": "Tier_2_HighSVM",
                "acc": acc,
                "len": len(precursor),
                "prob": prob,
                "macrel_class": "Non-Macrel",
                "cleaved": False,
                "mature_seq": precursor,
                "precursor_seq": precursor
            })

    # Capture any Macrel-only positives that fell below SVM threshold
    for acc, m_data in macrel_leads.items():
        if acc not in svm_leads and not any (m_data["norm_seq"] == s["norm_seq"] for s in svm_leads.values()):
            seq = m_data["mature_seq"]
            tier_3.append({
                "tier": "Tier_3_MacrelOnly",
                "acc": acc,
                "len": len(seq),
                "prob": 0.0,
                "macrel_class": m_data["macrel_class"],
                "cleaved": True,
                "mature_seq": seq,
                "precursor_seq": seq
                })

    tier_1.sort(key=lambda x: x["prob"], reverse=True)
    tier_2.sort(key=lambda x: x["prob"], reverse=True)
    tier_3.sort(key=lambda x: x["len"])

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w", encoding="utf-8") as tsv_out:
        tsv_out.write("Tier\tAccession\tLength\tModelB_SVM_Prob\tMacrel_Class\tN_Cleaved\tMature_Sequence\tFull_Precursor\n")
        for item in tier_1 + tier_2 + tier_3:
            tsv_out.write(
                f"{item['tier']}\t{item['acc']}\t{item['len']}\t{item['prob']:.4f}\t"
                f"{item['macrel_class']}\t{item['cleaved']}\t{item['mature_seq']}\t{item['precursor_seq']}\n"
            )

    # Write prioritized mature FASTA for dowstream biophysics and ESMFold
    with open(output_faa, "w", encoding="utf-8") as faa_out:
        for item in tier_1:
            faa_out.write(f">{item['acc']}|Tier1|Prob={item['prob']:.4f}|Cleaved={item['cleaved']}\n{item['mature_seq']}\n")
        for item in tier_2:
            faa_out.write(f">{item['acc']}|Tier2|Prob={item['prob']:.4f}\n{item['mature_seq']}\n")

    print(
        f"\n-----------------------------------------\n"
        f"          CONSENSUS MERGE SUMMARY          "
        f"\n-----------------------------------------\n"
        f"Tier 1 (Dual Consensus Hits): {len(tier_1)}\n"
        f"Tier 2 (High-Confidence SVM): {len(tier_2)}\n"
        f"Tier 3 (Macrel Only)        : {len(tier_3)}\n"
        f"Total Screened Candidates   : {len(tier_1)} + {len(tier_2)} + {len(tier_3)}"
        f"\n-----------------------------------------\n"
        f"Consensur TSV: {output_tsv}\n"
        f"Target FASTA : {output_faa}"
        f"\n-----------------------------------------\n"
    )

    if tier_1:
        print("--- Top tier 1 Consensus Leads ---")
        for i, lead in enumerate(tier_1[:5], 1):
            cleaved_tag = "[N-cleaved]" if lead["cleaved"] else "[Intact]"
            print(f" {i}. {lead['acc']} | SVM Prob: {lead['prob']:.4f} | {cleaved_tag} | Length: {lead['len']} AA")

def main() -> None:
    base_dir = Path(__file__).resolve().parent.parent
    macrel_file = base_dir / "data" / "02_classifications" / "macrel" / "macrel.out.prediction.tsv"
    svm_file = base_dir / "data" / "02_classifications" / "ampir" / "ampir_predictions.tsv"
    output_tsv = base_dir / "data" / "02_classifications" / "consensus" / "consensus_leads.tsv"
    output_faa = base_dir / "data" / "02_classifications" / "consensus" / "consensus_leads.faa"

    print("Ingesting Model B (SVM) predictions...")
    svm_leads = load_modelb_results(svm_file)
    print(f"Loaded {len(svm_leads)} SVM predictions.")

    print("Ingesting Macrel predictions...")
    macrel_leads = load_macrel_results(macrel_file)
    print(f"Loaded {len(macrel_leads)} Macrel positive predictions.")

    print("Reconciling N-terminal offsets and ranking tiers...")
    merge_and_tier_leads(macrel_leads, svm_leads, output_tsv, output_faa)

if __name__ == "__main__":
    main()
        