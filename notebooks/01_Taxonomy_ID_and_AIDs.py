# ---
# jupyter:
#   jupytext:
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
# # PubChem Taxonomy → BioAssay Master Table
#
# This notebook builds a unified table mapping:
#
# **Pathogen → Taxonomy ID (TaxID) → BioAssay ID (AID)**  
#
# This notebook:
# 1. Loads each pathogen’s exported file  
# 2. Extracts and expands pipe-separated AIDs  
# 3. Builds one unified TaxID–AID table  
# 4. Computes summary statistics  
# 5. Saves processed outputs for downstream bioactivity retrieval

# %% [markdown]
# ## 0. Setup

# %%
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import json
import requests

# %%
NOTEBOOK_DIR = Path().resolve()
DATA_RAW = NOTEBOOK_DIR.parent / "data" / "raw"
DATA_PROCESSED = NOTEBOOK_DIR.parent / "data" / "processed"
DATA_PROCESSED.mkdir(exist_ok=True)

# %%
pathogens = [
    "Acinetobacter baumannii", "Candida albicans", "Campylobacter",
    "Escherichia coli", "Enterococcus faecium", "Enterobacter",
    "Helicobacter pylori", "Klebsiella pneumoniae",
    "Mycobacterium tuberculosis", "Neisseria gonorrhoeae",
    "Pseudomonas aeruginosa", "Plasmodium falciparum",
    "Staphylococcus aureus", "Schistosoma mansoni",
    "Streptococcus pneumoniae"
]


# %% [markdown]
# ## 1. Function to Load and Expand One Pathogen File
#
# PubChem does not currently expose a stable API for retrieving
# organism-linked BioAssays directly from organism names.
# Therefore, the initial files used here were exported manually from:
#
# **PubChem → Search → "Organism name" → Taxonomy →  
# Actions → BioAssays → Download: *Summary (Search Results)***
#
# Each exported CSV contains:
# - `Taxonomy_ID`
# - `Taxonomy_Name`
# - `Linked_BioAssays` (pipe-separated AIDs)
#
# Each CSV contains a pipe-separated list of linked BioAssays.
# We expand these to one row per AID and annotate with the pathogen name.

# %%
def load_taxid_aid(pathogen: str) -> pd.DataFrame:
    """Load and expand the PubChem Taxonomy→BioAssay mapping for a pathogen,
    keeping ALL Taxonomy_ID rows even if no Linked_BioAssays exist.
    """
    filename = f"PubChem_taxonomy_text_{pathogen}.csv"
    filepath = DATA_RAW / filename

    df = pd.read_csv(filepath)

    # DO NOT DROP rows with missing Linked_BioAssays
    df["Linked_BioAssays"] = df["Linked_BioAssays"].fillna("")

    # Expand AIDs
    df["AID"] = df["Linked_BioAssays"].astype(str).str.split("|")
    df = df.explode("AID")

    # Clean up AID strings
    df["AID"] = df["AID"].astype(str).str.strip()
    df["AID"] = df["AID"].replace(["", "nan"], pd.NA)

    # Convert AID to numeric where possible (NaN stays)
    df["AID"] = pd.to_numeric(df["AID"], errors="coerce")

    df["Pathogen"] = pathogen

    return df[["Pathogen", "Taxonomy_ID", "Taxonomy_Name", "AID"]].drop_duplicates()

# %% [markdown]
# ## 2. Process All Pathogens

# %%
pathogen_tables = []

for pathogen in pathogens:
    print(f"Processing {pathogen}...")
    df_p = load_taxid_aid(pathogen)
    pathogen_tables.append(df_p)

df_all = pd.concat(pathogen_tables, ignore_index=True)
df_all = df_all[["Pathogen", "Taxonomy_ID", "Taxonomy_Name", "AID"]]
df_all.head()

# %% [markdown]
# ## 3. Summary Statistics

# %%
# Count number of AIDs per pathogen
summary_by_pathogen = (
    df_all.groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="Taxonomy_AIDs")
    .sort_values("Taxonomy_AIDs", ascending=False)
    .reset_index(drop=True)
)

summary_by_pathogen

# %%
x = np.arange(len(summary_by_pathogen))

values = summary_by_pathogen["Taxonomy_AIDs"].values
labels = summary_by_pathogen["Pathogen"].values

# Plot
plt.figure(figsize=(10, 5))

plt.bar(
    x, values,
    color="#AA96FA",
    edgecolor="k",
    zorder=2,
    label="Taxonomy_AIDS"
)

plt.yscale("log")
plt.ylim(1, values.max() * 1.5)   # space above bars

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Save Outputs

# %%
# Save final csv file with all Taxonomy_ID and AID pairs
df_all.to_csv(DATA_PROCESSED / "all_taxid_aid.csv", index=False)

print("Saved: all_taxid_aid.csv")

# %% [markdown]
# ## 5. Build Taxonomy Table

# %% [markdown]
# Each pathogen corresponds to multiple
# PubChem Taxonomy IDs (species, strains, variants, substrains).
#
# Here we extract all TaxID entries present in the manually exported CSV files
# and build:
#
# 1. A clean **taxonomy table** (`Pathogen`, `Taxonomy_ID`, `Taxonomy_Name`)
# 2. A **Python dictionary** mapping each pathogen to all its TaxIDs
#
# This dictionary will be used later to retrieve the full set of BioAssays
# linked to each organism from the full PubChem BioAssay dataset.

# %%
# Extract the taxonomy IDs and names per pathogen
df_taxonomy = (
    df_all[["Pathogen", "Taxonomy_ID", "Taxonomy_Name"]]
    .drop_duplicates()
    .sort_values(["Pathogen", "Taxonomy_ID"])
    .reset_index(drop=True)
)

df_taxonomy.head(20)

# %%
# Convert taxonomy table into dictionary
dict_taxonomy = (
    df_taxonomy.groupby("Pathogen")["Taxonomy_ID"]
    .apply(list)
    .to_dict()
)

dict_taxonomy

# %%
# Save taxonomy table and dictionary

df_taxonomy.to_csv(DATA_PROCESSED / "taxonomy_table.csv", index=False)

with open(DATA_PROCESSED / "dict_taxonomy.json", "w") as f:
    json.dump(dict_taxonomy, f, indent=2)

print("Saved: taxonomy_table.csv and dict_taxonomy.json")


# %% [markdown]
# ## 6. SDQ Query for Bioassays aid & taxids

# %%
def get_aids_for_taxid(taxid: int):
    """Return curated AIDs for a given TaxID from the PubChem taxonomy module."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/taxonomy/taxid/{taxid}/aids/JSON"
    r = requests.get(url)

    if r.status_code != 200:
        return []

    data = r.json()

    try:
        info = data["InformationList"]["Information"][0]
        return info.get("AID", [])
    except:
        return []


# %%
# For each pathogen, fetch curated AIDs (ETA: 25min)
pugrest_aids_counts = {}

for pathogen, taxids in dict_taxonomy.items():
    aids = set()
    for tid in taxids:
        aids.update(get_aids_for_taxid(tid))
    pugrest_aids_counts[pathogen] = len(aids)

pugrest_aids_counts

# %%
# Comparison with previous strategy
summary_by_pathogen["PUGREST_AIDs"] = summary_by_pathogen["Pathogen"].map(pugrest_aids_counts)
summary_by_pathogen = summary_by_pathogen.sort_values("Taxonomy_AIDs", ascending=False).reset_index(drop=True)

summary_by_pathogen

# %% [markdown]
# Let's also add the manually search bioassays count:

# %%
# Add UI (user interface) counts (PubChem website search)
ui_counts = {
    "Acinetobacter baumannii": 15778,
    "Candida albicans": 23814,
    "Campylobacter": 622,
    "Escherichia coli": 63263,
    "Enterococcus faecium": 3864,
    "Enterobacter": 4023,
    "Helicobacter pylori": 1670,
    "Klebsiella pneumoniae": 11883,
    "Mycobacterium tuberculosis": 25323,
    "Neisseria gonorrhoeae": 1019,
    "Pseudomonas aeruginosa": 26093,
    "Plasmodium falciparum": 24519,
    "Staphylococcus aureus": 59672,
    "Schistosoma mansoni": 1276,
    "Streptococcus pneumoniae": 9474
}

summary_by_pathogen["PubChem_Web_AIDs"] = summary_by_pathogen["Pathogen"].map(ui_counts)
summary_by_pathogen = summary_by_pathogen.sort_values("PubChem_Web_AIDs", ascending=False).reset_index(drop=True)

summary_by_pathogen

# %% [markdown]
# The two PubChem methods (CSV Taxonomy export vs PUG REST taxonomy→aids) return different and only partially overlapping assay sets, and the direction of the discrepancy varies by pathogen.

# %% [markdown]
# Let's also add the manually search bioassays count:

# %%
# Extract arrays
labels          = summary_by_pathogen["Pathogen"].values
PubChem_Web_AIDs = summary_by_pathogen["PubChem_Web_AIDs"].values
Taxonomy_AIDs    = summary_by_pathogen["Taxonomy_AIDs"].values
PUGREST_AIDs     = summary_by_pathogen["PUGREST_AIDs"].values

N = len(labels)
x = np.arange(N)

bar_width = 0.25

plt.figure(figsize=(14, 6))

# LEFT BAR → PubChem Web UI (gray)
plt.bar(
    x - bar_width,
    PubChem_Web_AIDs,
    width=bar_width,
    color="#D2D2D2",
    ec="black",
    zorder=2,
    label="AIDs (PubChem Website Search)"
)

# MIDDLE BAR → Taxonomy CSV (purple)
plt.bar(
    x,
    Taxonomy_AIDs,
    width=bar_width,
    color="#AA96FA",
    ec="black",
    zorder=2,
    label="AIDs (Taxonomy CSV Export)"
)

# RIGHT BAR → PUG REST (yellow)
plt.bar(
    x + bar_width,
    PUGREST_AIDs,
    width=bar_width,
    color="#FAD782",
    ec="black",
    zorder=2,
    label="AIDs (PUG REST taxonomy→aids)"
)

# Aesthetics
plt.yscale("log")
plt.ylim([1, PubChem_Web_AIDs.max() * 1.5])

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()

# %% [markdown]
# ## 7. PubChem’s Aid2Taxid.tsv.gz file

# %% [markdown]
# To obtain a complete list of all BioAssays linked to our taxonomy IDs, we now switch to PubChem’s official FTP file Aid2Taxid.tsv.gz, which contains the full mapping of every PubChem BioAssay to all associated NCBI Taxonomy IDs.
#
# This file was downloaded manually from:
#
# https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2Taxid.gz

# %%
# Load file

df_aid2tax = pd.read_csv(DATA_RAW / "Aid2Taxid.tsv", sep="\t")
df_aid2tax.head(10)

# %%
# 1. Build pathogen–TaxID lookup
rows = []
for pathogen, taxids in dict_taxonomy.items():
    for tax in taxids:
        rows.append({"Pathogen": pathogen, "TaxID": tax})

df_tax_lookup = pd.DataFrame(rows)

# 2. Merge and extract all matching AIDs
df_merged = df_aid2tax.merge(df_tax_lookup, on="TaxID", how="inner")
df_merged.head()

# %%
# 3. Count AIDs per pathogen
summary_aid2tax = (
    df_merged.groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="Aid2Taxid_AIDs")
    .sort_values("Aid2Taxid_AIDs", ascending=False)
    .reset_index(drop=True)
)

summary_aid2tax

# %%
aid2tax_dict = dict(zip(summary_aid2tax["Pathogen"], summary_aid2tax["Aid2Taxid_AIDs"]))
summary_by_pathogen["Aid2Taxid_AIDs"] = summary_by_pathogen["Pathogen"].map(aid2tax_dict)
summary_by_pathogen = summary_by_pathogen.sort_values("PubChem_Web_AIDs", ascending=False).reset_index(drop=True)

summary_by_pathogen

# %%
# Save summary_by_pathogen
filepath = DATA_PROCESSED / "summary_by_pathogen.csv"
summary_by_pathogen.to_csv(filepath, index=False)

# %%
# ---------------------------------------
# Resume support: track which ZIP chunks are done
# ---------------------------------------
zip_log = DATA_PROCESSED / "processed_zip_chunks.txt"

if zip_log.exists():
    processed_chunks = set(zip_log.read_text().splitlines())
else:
    processed_chunks = set()

# ---------------------------------------
# Loop to check ALL Description folders
# ---------------------------------------

all_records = []      # list of dicts: {"AID", "TaxIDs_detected", "ZipFolder"}
already_saved = set() # to avoid writing same XML twice
total_matched_global = 0

for zip_file in tqdm(zip_files, desc="Processing ZIP chunks"):
    
    zip_chunk = zip_file.stem

    # --------------- RESUME SUPPORT ---------------
    if zip_chunk in processed_chunks:
        print(f"⏩ Skipping {zip_chunk} (already processed)")
        continue
    # ----------------------------------------------

    print(f"\n=== Processing ZIP: {zip_chunk} ===")

    # Temporary extraction folder
    temp_extract = DESC_DIR / f"{zip_chunk}_tmp"
    temp_extract.mkdir(exist_ok=True)

    # 1) Extract contents
    with zipfile.ZipFile(zip_file, "r") as zf:
        zf.extractall(temp_extract)

    # 2) Locate XML files
    xml_files = list(temp_extract.rglob("*.xml.gz")) + list(temp_extract.rglob("*.xml"))
    print(f"Found XML files: {len(xml_files)}")

    matched_for_zip = []

    # 3) Scan XML files
    for xml_path in tqdm(xml_files, desc=f"Scanning XMLs in {zip_chunk}", leave=False):
        taxids = extract_taxids_from_descr(xml_path)

        if taxids & target_taxid_set:
            aid = int(Path(xml_path).stem.split(".")[0])

            matched_for_zip.append(aid)
            total_matched_global += 1

            all_records.append({
                "AID": aid,
                "TaxIDs_detected": sorted(list(taxids)),
                "ZipFolder": zip_chunk
            })

            # Save uncompressed XML only once
            if aid not in already_saved:
                dest_path = KEEP_DESC / f"{aid}.xml"

                if str(xml_path).endswith(".gz"):
                    with gzip.open(xml_path, "rb") as f_in:
                        with open(dest_path, "wb") as f_out:
                            f_out.write(f_in.read())
                else:
                    shutil.copy2(xml_path, dest_path)

                already_saved.add(aid)

    # 4) Clean up
    shutil.rmtree(temp_extract)
    print(f"✓ Finished ZIP {zip_chunk}, removed temp folder")

    # --------------- MARK ZIP AS DONE ---------------
    with open(zip_log, "a") as f:

        f.write(zip_chunk + "\n")
    processed_chunks.add(zip_chunk)
    # ------------------------------------------------

# ---------------------------------------
# Convert results to DataFrame
# ---------------------------------------
df_results = pd.DataFrame(all_records)
df_results.to_csv(DATA_PROCESSED / "filtered_assays_description_results.csv", index=False)

print("\n\n=====================================")
print("✔ ALL ZIP FILES PROCESSED SUCCESSFULLY")
print("✔ Filtered XML files saved in:", KEEP_DESC)
print("✔ Full results saved as: filtered_assays_description_results.csv")
print("Total matched AIDs:", df_results["AID"].nunique())
print("=====================================")

# %%
# Extract arrays
labels           = summary_by_pathogen["Pathogen"].values
PubChem_Web_AIDs = summary_by_pathogen["PubChem_Web_AIDs"].values
Taxonomy_AIDs    = summary_by_pathogen["Taxonomy_AIDs"].values
PUGREST_AIDs     = summary_by_pathogen["PUGREST_AIDs"].values
Aid2Taxid_AIDs   = summary_by_pathogen["Aid2Taxid_AIDs"].values

N = len(labels)
x = np.arange(N)

bar_width = 0.18   # narrower to fit 4 bars

plt.figure(figsize=(15, 6))

# 1. LEFT BAR → PubChem Web UI (gray)
plt.bar(
    x - 1.5*bar_width,
    PubChem_Web_AIDs,
    width=bar_width,
    color="#D2D2D2",
    ec="black",
    zorder=2,
    label="AIDs (PubChem Website Search)"
)

# 2. SECOND BAR → Taxonomy CSV (purple)
plt.bar(
    x - 0.5*bar_width,
    Taxonomy_AIDs,
    width=bar_width,
    color="#AA96FA",
    ec="black",
    zorder=2,
    label="AIDs (Taxonomy CSV Export)"
)

# 3. THIRD BAR → PUG REST (yellow)
plt.bar(
    x + 0.5*bar_width,
    PUGREST_AIDs,
    width=bar_width,
    color="#FAD782",
    ec="black",
    zorder=2,
    label="AIDs (PUG REST taxonomy→aids)"
)

# 4. RIGHT BAR → Aid2Taxid (salmon)
plt.bar(
    x + 1.5*bar_width,
    Aid2Taxid_AIDs,
    width=bar_width,
    color="#FAA08B",
    ec="black",
    zorder=2,
    label="AIDs (Aid2Taxid.tsv)"
)

# Aesthetics
plt.yscale("log")
plt.ylim([1, max(PubChem_Web_AIDs) * 1.5])

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()

# %% [markdown]
# The comparison shows that the **manual taxonomy-CSV** export captures the fewest assays, while both **PUG REST (taxonomy→aids)** and **Aid2Taxid.tsv** recover many more—almost identical to each other—because they rely on the same structured taxonomy annotations deposited in PubChem.
#
# However, both still fall below the counts shown on the **PubChem website**, meaning many assays mentioning a pathogen are not formally annotated with a `TaxID` and only appear in the website search because they contain the organism name somewhere in their free-text descriptions.

# %% [markdown]
# To capture all these additional assays, we now need to process the **full PubChem BioAssay dataset** and search for pathogen names directly **inside text fields**:
#
# - `BioAssay Name`
# - `BioAssay Types`
# - `Project Category`
# - `Source Name`
# - `Source ID`
# - `BioAssay Group`
#
# This free-text search will allow us to recover the assays that PubChem displays on the website but does not structurally tag with taxonomy IDs.

# %% [markdown]
# ## 7. Download PubChem Full BioAssay Dataset

# %% [markdown]
# To obtain a complete and unbiased mapping between BioAssays (AIDs) and taxonomy identifiers (TaxIDs), we now switch to PubChem’s official FTP dataset Aid2Taxid.gz.
#
# The file was downloaded manually from:
#
# https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/bioassays.tsv.gz 

# %% [markdown]
# This dataset contains one row per BioAssay (AID) together with a rich set of curated descriptors, including:
#
# Columns included:
# - `AID` – PubChem BioAssay identifier
# - `BioAssay Name` – title provided by the depositor
# - `Deposit Date` / `Modify Date`
# - `Source Name` – ChEMBL, Broad Institute, NIAID, etc.
# - `Source ID` – depositor-specific ID
# - `Substance Type` – small molecule / RNAi reagent
# - `Outcome Type` – Screening / Confirmatory / Summary / Literature
# - `Project Category`
# - `BioAssay Group`
# - `BioAssay Types` – Biochemical, Cell-based, Organism-based…
# - `Protein Accessions`
# - `UniProt IDs`
# - `Gene IDs`
# - `Target TaxIDs` – curated organism targets
# - `Taxonomy IDs` – additional taxonomy annotations
# - `Number of Tested SIDs / Active SIDs`
# - `Number of Tested CIDs / Active CIDs`

# %% [markdown]
#

# %%
bio_file = DATA_RAW / "bioassays.tsv"
df_bio = pd.read_csv(bio_file, sep="\t")


# %%
# 1. Clean taxonomy fields for matching

def to_list(field):
    """Convert pipe-separated taxid field to a clean list of strings."""
    if pd.isna(field):
        return []
    return [x.strip() for x in str(field).split("|") if x.strip() != ""]

df_bio["Target_TaxIDs_list"] = df_bio["Target TaxIDs"].apply(to_list)
df_bio["Taxonomy_IDs_list"] = df_bio["Taxonomy IDs"].apply(to_list)

# Load our pathogen taxonomy table
tax_table = pd.read_csv(DATA_PROCESSED / "taxonomy_table.csv")

# Map pathogen → list of TaxIDs
pathogen_taxids = (
    tax_table.groupby("Pathogen")["Taxonomy_ID"]
    .apply(lambda x: set(map(str, x)))
    .to_dict()
)

# 2. Function to extract AIDs for a pathogen

def extract_aids_for_pathogen(df, pathogen, taxid_set):
    """Return all AIDs where either taxonomy column matches the pathogen."""
    mask_target = df["Target_TaxIDs_list"].apply(
        lambda lst: any(t in lst for t in taxid_set)
    )
    mask_taxonomy = df["Taxonomy_IDs_list"].apply(
        lambda lst: any(t in lst for t in taxid_set)
    )
    subset = df[mask_target | mask_taxonomy]
    return set(subset["AID"].astype(int))


# 3. Extract AIDs for all pathogens

bioassay_aids = {}

for pathogen, taxids in pathogen_taxids.items():
    aids = extract_aids_for_pathogen(df_bio, pathogen, taxids)
    bioassay_aids[pathogen] = aids
    print(f"{pathogen:30s} → {len(aids)} AIDs")

bioassay_aids

# %%
# 4. Save output as CSV

rows = []
for pathogen, aids in bioassay_aids.items():
    for aid in aids:
        rows.append({"Pathogen": pathogen, "AID": aid})

df_out = pd.DataFrame(rows)
df_out.to_csv(DATA_PROCESSED / "bioassays_taxonomy_matched.csv", index=False)

print("Saved: bioassays_taxonomy_matched.csv")
