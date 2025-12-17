import json
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# ---- Paths ----
PROJECT_ROOT = Path("/Users/maria/Documents/Ersilia/PubChem/pubchem-antimicrobial-tasks")
FILTERED_DISPLAY_DIR = PROJECT_ROOT / "data" / "raw" / "filtered_assays_v2" / "Display"
SUMMARY_CSV = PROJECT_ROOT / "data" / "processed" / "summary_display.csv"

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

    def extract_first_string(info_obj):
        try:
            return info_obj["Value"]["StringWithMarkup"][0]["String"]
        except:
            return None

    def extract_taxid(info_obj):
        url = info_obj.get("URL", "")
        if "taxonomy" in url:
            try:
                return url.strip("/").split("/")[-1]
            except:
                return None
        return None

    def walk_sections(sections):
        for section in sections:
            heading = section.get("TOCHeading", "")

            # ✅ Correct handling of the Target section (before looping through 'Information')
            if heading == "Target":
                for subsection in section.get("Section", []):
                    for info in subsection.get("Information", []):
                        val = extract_first_string(info)
                        if val:
                            result["Target"] = val
                            break  # stop after first match

            # ⬇️ Now process section's own info fields
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

            # Recurse into nested sections
            if "Section" in section:
                walk_sections(section["Section"])

    # 🔁 Start walking from top-level sections
    walk_sections(record.get("Section", []))

    return result

# ---- Rebuild CSV from all JSONs ----
def main():
    print(f"📂 Rebuilding summary from JSON files in {FILTERED_DISPLAY_DIR}")
    json_files = list(FILTERED_DISPLAY_DIR.glob("*.json"))
    print(f"🔍 Found {len(json_files):,} JSON files.")

    rows = []
    for file in tqdm(json_files, desc="Parsing JSONs"):
        try:
            data = json.loads(file.read_text(encoding="utf-8"))
            row = parse_json(data)
            if row.get("PubChem_AID") is not None:
                rows.append(row)
        except Exception as e:
            continue  # silently skip bad JSONs

    df = pd.DataFrame(rows)
    df.to_csv(SUMMARY_CSV, index=False)
    print(f"\n✅ Saved rebuilt summary to: {SUMMARY_CSV}")

if __name__ == "__main__":
    main()