import zipfile
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import gzip
import shutil

# --- Paths ---
PROJECT_ROOT = Path("/Users/maria/Documents/Ersilia/PubChem/pubchem-antimicrobial-tasks")
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_DIR = DATA_RAW / "pubchem_bioassays" / "data"
FILTERED_DATA_DIR = DATA_RAW / "filtered_assays_v2" / "Data"
FILTERED_DATA_DIR.mkdir(parents=True, exist_ok=True)
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

# --- Load filtered AIDs ---
df = pd.read_csv(DATA_PROCESSED / "filtered_description_with_organisms_v2_REBUILT.csv")
filtered_aids = set(df["AID"].dropna().astype(str))  # Make AIDs strings for filename matching

# --- Track what's found ---
found_aids = set()
not_found_aids = set(filtered_aids)  # Start assuming all are missing

# --- Process all ZIPs ---
zip_files = sorted(DATA_DIR.glob("*.zip"))
print(f"🔍 Scanning {len(zip_files)} ZIP files for matching AIDs...")

for zip_path in tqdm(zip_files, desc="Processing ZIPs"):
    with zipfile.ZipFile(zip_path, "r") as zf:
        for file_in_zip in zf.namelist():
            file_name = Path(file_in_zip).name
            aid = file_name.replace(".csv.gz", "")
            if aid in filtered_aids:
                try:
                    # Extract and decompress
                    with zf.open(file_in_zip) as compressed_file:
                        with gzip.open(compressed_file, 'rt') as gz_file:
                            output_csv_path = FILTERED_DATA_DIR / f"{aid}.csv"
                            with open(output_csv_path, 'w') as out_csv:
                                shutil.copyfileobj(gz_file, out_csv)

                    found_aids.add(aid)
                    not_found_aids.discard(aid)  # Remove from missing list

                except Exception as e:
                    print(f"⚠️ Error extracting AID {aid}: {e}")

# --- Save missing AIDs ---
if not_found_aids:
    with open(PROJECT_ROOT / "missing_filtered_aids.log", "w") as f:
        for aid in sorted(not_found_aids):
            f.write(f"{aid}\n")
    print(f"\n⚠️ Missing AIDs: {len(not_found_aids)}")
    print("📄 See missing_filtered_aids.log")

print("\n✅ Done!")
print(f"🗂️ Extracted CSVs saved to: {FILTERED_DATA_DIR}")
print(f"🧪 AIDs found: {len(found_aids)} / {len(filtered_aids)}")