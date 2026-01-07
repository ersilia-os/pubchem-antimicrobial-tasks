import json
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import requests
import time

# --- CONFIG ---
MAX_WORKERS = 6
RETRIES = 3
DELAY = 2  # seconds between retries

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DISPLAY_DIR = DATA_RAW / "pubchem_bioassays" / "Display"
DISPLAY_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = DATA_PROCESSED / "04_display_summary.csv"

# --- Load existing summary if exists ---
if SUMMARY_CSV.exists():
    df_existing = pd.read_csv(SUMMARY_CSV)
    parsed_aids = set(df_existing["PubChem_AID"].astype(int))
else:
    df_existing = pd.DataFrame()
    parsed_aids = set()

# --- Define full set of AIDs (1 to 1.4M) ---
all_aids = set(range(1, 1_400_001))
already_downloaded = {int(f.stem) for f in DISPLAY_DIR.glob("*.json")}
pending_aids = sorted(already_downloaded - parsed_aids)

print(f"🧠 {len(already_downloaded):,} JSONs available.")
print(f"✅ {len(parsed_aids):,} already in 04_display_summary.csv")
print(f"📥 {len(pending_aids):,} remaining to parse and append")

# --- Helpers ---
def fetch_and_parse_json(aid, retries=RETRIES, delay=DELAY):
    try:
        data = json.loads((DISPLAY_DIR / f"{aid}.json").read_text(encoding="utf-8"))
        return parse_json(data)
    except Exception:
        return None

def extract_first_string(info):
    try:
        return info["Value"]["StringWithMarkup"][0]["String"]
    except:
        return None

def extract_taxid(info):
    url = info.get("URL", "")
    if "taxonomy" in url:
        try:
            return url.strip("/").split("/")[-1]
        except:
            return None
    return None

def parse_json(data):
    record = data.get("Record", {})
    aid = record.get("RecordNumber", None)

    result = {
        "PubChem_AID": aid,
        "ChEMBL_ID": None,
        "Compounds_Tested": None,
        "Compounds_Active": None,
        "Compounds_Inactive": None,
        "Target": None,
        "Assay_Type": None,
        "Assay_Format": None,
        "Assay_Organism": None,
        "Organism_TaxID": None,
        "Assay_Strain": None,
        "Organism_Target": None,
        "Protein_Target": None,
        "Source": None,
    }

    def walk(sections):
        for section in sections:
            heading = section.get("TOCHeading", "")

            if "Section" in section:
                walk(section["Section"])

            for info in section.get("Information", []):
                name = info.get("Name", "")
                val = extract_first_string(info)

                if heading in ["Source Information", "External ID"] and val and "CHEMBL" in val:
                    result["ChEMBL_ID"] = val.split("::")[-1]

                elif heading == "Tested Compounds":
                    value = info.get("Value", {}).get("Number", [None])[0]
                    if name == "All Compounds":
                        result["Compounds_Tested"] = value
                    elif name == "Active Compounds":
                        result["Compounds_Active"] = value
                    elif name == "Inactive Compounds":
                        result["Compounds_Inactive"] = value

                elif heading == "BioAssay Annotations":
                    if name == "Assay Type":
                        result["Assay_Type"] = val
                    elif name == "Assay Format":
                        result["Assay_Format"] = val
                    elif name == "Assay Organism":
                        result["Assay_Organism"] = val
                        result["Organism_TaxID"] = extract_taxid(info)
                    elif name == "Assay Strain":
                        result["Assay_Strain"] = val

                elif heading == "Protein Target":
                    if val:
                        result["Protein_Target"] = val

                elif heading == "Organism Target":
                    if val:
                        result["Organism_Target"] = val

                elif heading == "Source":
                    if val:
                        result["Source"] = val

                elif heading == "Target":
                    for subsection in section.get("Section", []):
                        for subinfo in subsection.get("Information", []):
                            val = extract_first_string(subinfo)
                            if val:
                                result["Target"] = val
                                break

    walk(record.get("Section", []))
    return result

# --- Main run ---
def main():
    if not pending_aids:
        print("✅ All JSONs already parsed.")
        return

    new_rows = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(fetch_and_parse_json, aid): aid for aid in pending_aids}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Parsing JSONs"):
            result = future.result()
            if result:
                new_rows.append(result)

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        updated_df = pd.concat([df_existing, new_df], ignore_index=True)
        updated_df.to_csv(SUMMARY_CSV, index=False)
        print(f"✅ Appended {len(new_df)} new entries → {SUMMARY_CSV}")
    else:
        print("⚠️ No new valid JSONs parsed.")

if __name__ == "__main__":
    main()