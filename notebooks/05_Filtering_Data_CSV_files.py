# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: pamt
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 05. Filtering Data CSV files
#
# This notebook extracts the `*.csv.gz` data files for the **filtered AIDs** and stores them in a dedicated `Data/` folder.
#
# These CSV files contain the assay results, including compound identifiers and activity values.

# %% [markdown]
# ## 00. Setup

# %%
from pathlib import Path
import zipfile
import shutil
from tqdm import tqdm
import pandas as pd

# %%
# Paths
PROJECT_ROOT = Path("/Users/maria/Documents/Ersilia/PubChem/pubchem-antimicrobial-tasks")
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_DIR = DATA_RAW / "pubchem_bioassays" / "data"
FILTERED_DATA_DIR = DATA_RAW / "filtered_assays_v2" / "Data"
FILTERED_DATA_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## 01. Load list of filtered AIDs

# %%
# Load filtered assay metadata
filtered_descr = pd.read_csv(DATA_PROCESSED / "filtered_description_with_organisms_v2_REBUILT.csv")

# Extract list of AIDs (integers)
filtered_aids = sorted(filtered_descr["AID"].dropna().astype(int).unique())
print(f"Total filtered AIDs: {len(filtered_aids):,}")

# %% [markdown]
# ## 02. Map AIDs to ZIP files
#
# Each ZIP file is named as a range: `0000001_0001000.zip`
# For each AID, we can figure out which ZIP chunk it belongs to.

# %%
def get_zip_filename(aid):
    chunk_size = 1000
    start = (aid - 1) // chunk_size * chunk_size + 1
    end = start + chunk_size - 1
    return f"{start:07d}_{end:07d}.zip"

# Map AIDs to their ZIPs
aid_to_zip = {aid: get_zip_filename(aid) for aid in filtered_aids}

# Create ZIP → list of AIDs mapping
from collections import defaultdict
zip_to_aids = defaultdict(list)
for aid, zipfile_name in aid_to_zip.items():
    zip_to_aids[zipfile_name].append(aid)

print(f"Unique ZIP files to process: {len(zip_to_aids)}")

# %% [markdown]
# ## 03. Extract `.csv.gz` files for filtered AIDs

# %% [markdown]
# The script `02_run_parallel_extract_csv.py` has been used to unzip and save the data files from the filtered assays.

# %% [markdown]
# ## 04. Obtain counts of compounds and substances

# %%
# Prepare list of CSV files
csv_files = list(FILTERED_DATA_DIR.glob("*.csv"))

# Collect summary rows
summary_rows = []

print(f"Processing {len(csv_files)} files...")

for csv_file in tqdm(csv_files):
    try:
        df = pd.read_csv(csv_file, low_memory=False)
        
        aid = csv_file.stem  # file name without .csv → AID
        substances_count = df["PUBCHEM_SID"].nunique() if "PUBCHEM_SID" in df.columns else 0
        compound_count = df["PUBCHEM_CID"].nunique() if "PUBCHEM_CID" in df.columns else 0
        
        summary_rows.append({
            "AID": int(aid),
            "substances_count": substances_count,
            "compound_count": compound_count
        })

    except Exception as e:
        print(f"⚠️ Error processing {csv_file.name}: {e}")

# Create summary DataFrame
summary_data = pd.DataFrame(summary_rows).sort_values("AID")

# Save summary
summary_data.to_csv(DATA_PROCESSED / "summary_data.csv", index=False)
print(f"\n✔️ Summary saved to: {DATA_PROCESSED}")

summary_data.head(10)

# %%
