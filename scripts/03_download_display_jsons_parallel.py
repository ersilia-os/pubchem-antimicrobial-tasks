import pandas as pd
import requests
import time
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# ---- CONFIGURATION ----
MAX_WORKERS = 6  # Number of parallel threads (adjust based on your system)
RETRIES = 3       # Number of retries on failure
DELAY = 2         # Delay between retries in seconds

PROJECT_ROOT = Path("/Users/maria/Documents/Ersilia/PubChem/pubchem-antimicrobial-tasks")
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_RAW = PROJECT_ROOT / "data" / "raw"
FILTERED_DISPLAY_DIR = DATA_RAW / "filtered_assays_v2" / "Display"
FILTERED_DISPLAY_DIR.mkdir(parents=True, exist_ok=True)

# ---- LOAD FILTERED AIDs ----
summary_description = pd.read_csv(DATA_PROCESSED / "filtered_description_with_organisms_v2_REBUILT.csv")
aids = summary_description["AID"].dropna().astype(int).unique().tolist()

# ---- FETCH + PARSE FUNCTION ----
def fetch_and_parse_display_json(aid, retries=RETRIES, delay=DELAY):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/bioassay/{aid}/JSON/"
    json_file = FILTERED_DISPLAY_DIR / f"AID_{aid}_display.json"

    if json_file.exists():
        try:
            data = json.loads(json_file.read_text(encoding="utf-8"))
            return parse_json(data)
        except Exception:
            return {"PubChem_AID": aid, "Status": "corrupted"}

    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            json_file.write_text(response.text, encoding="utf-8")
            data = response.json()
            return parse_json(data)
        except requests.exceptions.RequestException:
            time.sleep(delay)
        except Exception:
            return {"PubChem_AID": aid, "Status": "invalid_json"}

    return {"PubChem_AID": aid, "Status": "failed"}

# ---- PARSE JSON FUNCTION ----
def parse_json(data):
    record = data.get("Record", {})
    aid = record.get("RecordNumber", None)

    chembl_id = None
    compounds_all = None
    compounds_active = None
    compounds_inactive = None

    for section in record.get("Section", []):
        heading = section.get("TOCHeading", "")

        if heading == "Source Information":
            for info in section.get("Information", []):
                if info.get("Name") == "SourceID" and info.get("Value", {}).get("StringWithMarkup", []):
                    source_id = info["Value"]["StringWithMarkup"][0].get("String")
                    if source_id and "CHEMBL" in source_id:
                        chembl_id = source_id.split("::")[-1]  # Extract CHEMBLxxxxxxx only

        elif heading == "Tested Compounds":
            for info in section.get("Information", []):
                name = info.get("Name", "")
                value = info.get("Value", {}).get("Number", [None])[0]
                if name == "All Compounds":
                    compounds_all = value
                elif name == "Active Compounds":
                    compounds_active = value
                elif name == "Inactive Compounds":
                    compounds_inactive = value

    return {
        "PubChem_AID": aid,
        "ChEMBL_ID": chembl_id,
        "Compounds_All": compounds_all,
        "Compounds_Active": compounds_active,
        "Compounds_Inactive": compounds_inactive,
        "Status": "ok"
    }

# ---- PARALLEL DOWNLOAD + EXTRACTION ----
def main():
    results = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(fetch_and_parse_display_json, aid): aid for aid in aids}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading + Parsing"):
            result = future.result()
            results.append(result)

    df = pd.DataFrame(results)
    df = df[df["Status"] == "ok"].drop(columns="Status")
    df.to_csv(DATA_PROCESSED / "summary_display.csv", index=False)
    print(f"\n✅ Saved summary CSV to: {DATA_PROCESSED}")

if __name__ == "__main__":
    main()