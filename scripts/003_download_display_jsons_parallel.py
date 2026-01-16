#!/usr/bin/env python3

import json
import time
import requests
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# ============================
# CONFIGURATION
# ============================

MAX_WORKERS = 6
RETRIES = 3
DELAY = 2  # seconds between retries

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

AIDS_CSV = DATA_PROCESSED / "05_filtered_aids.csv"
DISPLAY_JSON_DIR = DATA_RAW / "filtered_assays_v3" / "Display"
DISPLAY_JSON_DIR.mkdir(parents=True, exist_ok=True)

OUT_CSV = DATA_PROCESSED / "06_display_info.csv"

# ============================
# LOAD AIDs
# ============================

all_aids = pd.read_csv(AIDS_CSV)["AID"].dropna().astype(int).unique().tolist()

if OUT_CSV.exists():
    already_parsed = pd.read_csv(OUT_CSV)["PubChem_AID"].dropna().astype(int).unique().tolist()
else:
    already_parsed = []

pending_aids = sorted(set(all_aids) - set(already_parsed))
print(f"🔍 Total filtered AIDs: {len(all_aids)}")
print(f"✅ Already parsed: {len(already_parsed)}")
print(f"⏳ Remaining to process: {len(pending_aids)}\n")

# ============================
# FETCH + PARSE FUNCTIONS
# ============================

def fetch_and_parse_display_json(aid, retries=RETRIES, delay=DELAY):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/bioassay/{aid}/JSON/"
    json_path = DISPLAY_JSON_DIR / f"AID_{aid}_display.json"

    # Use cached version if available
    if json_path.exists():
        try:
            data = json.loads(json_path.read_text(encoding="utf-8"))
            return parse_json(data)
        except Exception:
            return None

    # Download from PubChem
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=15)
            r.raise_for_status()
            json_path.write_text(r.text, encoding="utf-8")
            return parse_json(r.json())
        except requests.exceptions.RequestException:
            time.sleep(delay)
        except Exception:
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
        "Tested_Substances": None,  # ✅ NEW
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

    def extract_first_string(info_obj):
        try:
            return info_obj["Value"]["StringWithMarkup"][0]["String"]
        except:
            return None

    def extract_taxid(info_obj):
        url = info_obj.get("URL", "")
        if "taxonomy" in url:
            return url.strip("/").split("/")[-1]
        return None

    def walk_sections(sections):
        for section in sections:
            heading = section.get("TOCHeading", "")

            # ✅ Extract Tested Substances
            if heading == "Tested Substances":
                for info in section.get("Information", []):
                    if info.get("Name") == "All Substances":
                        result["Tested_Substances"] = info.get("Value", {}).get("Number", [None])[0]

            # ✅ Target
            elif heading == "Target":
                for subsection in section.get("Section", []):
                    for info in subsection.get("Information", []):
                        val = extract_first_string(info)
                        if val:
                            result["Target"] = val
                            break

            # Other info
            for info in section.get("Information", []):
                name = info.get("Name", "")
                val = extract_first_string(info)

                if heading in ["Source Information", "External ID"] and val and "CHEMBL" in val:
                    result["ChEMBL_ID"] = val.split("::")[-1]

                elif heading == "Tested Compounds":
                    number = info.get("Value", {}).get("Number", [None])[0]
                    if name == "All Compounds":
                        result["Compounds_Tested"] = number
                    elif name == "Active Compounds":
                        result["Compounds_Active"] = number
                    elif name == "Inactive Compounds":
                        result["Compounds_Inactive"] = number

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

                elif heading == "Protein Target" and val:
                    result["Protein_Target"] = val

                elif heading == "Organism Target" and val:
                    result["Organism_Target"] = val

                elif heading == "Source" and val:
                    result["Source"] = val

            # Recurse into nested subsections
            if "Section" in section:
                walk_sections(section["Section"])

    walk_sections(record.get("Section", []))
    return result

# ============================
# MAIN EXECUTION
# ============================

def main():
    if not pending_aids:
        print("✅ All AIDs already processed.")
        return

    parsed = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(fetch_and_parse_display_json, aid): aid
            for aid in pending_aids
        }

        for fut in tqdm(as_completed(futures), total=len(futures), desc="Fetching display JSONs"):
            res = fut.result()
            if res:
                parsed.append(res)

                # Incremental write
                df_partial = pd.DataFrame([res])
                df_partial.to_csv(OUT_CSV, mode="a", header=not OUT_CSV.exists(), index=False)

    print(f"\n✅ Completed. New AIDs parsed: {len(parsed)}")
    print(f"📄 Output saved to: {OUT_CSV}")


if __name__ == "__main__":
    main()