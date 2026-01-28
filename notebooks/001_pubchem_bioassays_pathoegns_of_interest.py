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
# # PubChem Pathogen-Linked BioAssay Mapping and Filtering
#
# This notebook builds a unified and reproducible workflow to identify and filter **PubChem BioAssays (AIDs)** relevant to pathogens of interest using multiple strategies.
#
# It:
# - Constructs a **Pathogen → Taxonomy ID** dictionary
# - Compares multiple methods to retrieve **AIDs per pathogen**, including:
#   - PubChem Website (UI) search
#   - Manual taxonomy exports
#   - PUG REST API queries
#   - PubChem's official `Aid2Taxid.tsv.gz` mapping
#   - Downloading and filtering the entire PubChem BioAssay dataset locally (~250k assays)
# - Extracts structured **Display metadata** (e.g., targets, organisms, compound counts)
# - Compares PubChem and ChEMBL bioassay coverage per pathogen
#
# The goal is to derive a high-confidence set of bioassays per pathogen that balances recall (UI alignment) and reproducibility (structured filtering), enabling downstream analysis and benchmarking.

# %% [markdown]
# ## 0. Setup

# %%
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import json
import requests
import re
from bs4 import BeautifulSoup
import seaborn as sns

# %%
# Project paths
NOTEBOOK_DIR = Path().resolve()

DATA_RAW = NOTEBOOK_DIR.parent / "data" / "raw"
DATA_PROCESSED = NOTEBOOK_DIR.parent / "data" / "processed"

PUBCHEM_DIR = DATA_RAW / "pubchem_bioassays"
DESC_DIR = PUBCHEM_DIR / "Description"
DATA_DIR = PUBCHEM_DIR / "Data"

KEEP_DIR = DATA_RAW / "filtered_assays"
KEEP_DESC = KEEP_DIR / "Description"
KEEP_DATA = KEEP_DIR / "Data"

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
# ## 01. Build Pathogen Taxonomy Table

# %% [markdown]
# ### 1.1. Manually download pathogen summary
# PubChem does not currently expose a stable API for retrieving organism-linked BioAssays directly from organism names. Therefore, the initial files used here were exported manually from:
#
# **PubChem → Search → "Organism name" → Taxonomy →  Download: *Summary (Search Results)**
#
# Each exported CSV (such as `PubChem_taxonomy_text_Acinetobacter baumannii.csv`) has been manually saved under `data\raw` and it contains:
# - `Taxonomy_ID`
# - `Taxonomy_Name`
#

# %%
abaumanii_summary = pd.read_csv(DATA_RAW / f"PubChem_taxonomy_text_Acinetobacter baumannii.csv")
abaumanii_summary.head(1)


# %% [markdown]
# ### 1.2. Keep only pathogen name and taxonomy
# Create a single Pathogen with associated taxonomies table

# %%
def pathogen_taxid(pathogen: str) -> pd.DataFrame:
    filename = f"PubChem_taxonomy_text_{pathogen}.csv"
    filepath = DATA_RAW / filename

    df = pd.read_csv(filepath)
    df["Pathogen"] = pathogen

    return df[["Pathogen", "Taxonomy_ID", "Taxonomy_Name"]].drop_duplicates()


# %%
pathogens_taxid_table = []

for pathogen in pathogens:
    print(f"Processing {pathogen}...")
    single_pathogen_taxid = pathogen_taxid(pathogen)
    pathogens_taxid_table.append(single_pathogen_taxid)

pathogens_taxid = pd.concat(pathogens_taxid_table, ignore_index=True)
pathogens_taxid = pathogens_taxid[["Pathogen", "Taxonomy_ID", "Taxonomy_Name"]]
pathogens_taxid

# %%
pathogens_taxid.to_csv(DATA_PROCESSED / "00_pathogens_taxid.csv", index=False)

# %% [markdown]
# ### 1.3. Manually clean taxonomy table
# Some taxonomy captures with PubChem query do not match the expected pathogens and have to be manually eliminated:

# %%
wrong_taxonomies = {
    "Acinetobacter baumannii": {
        "Acinetobacter calcoaceticus/baumannii complex",
    },
    "Candida albicans": {
        "Candida tropicalis",
    },
    "Campylobacter": {
        "Helicobacter pylori",
        "Helicobacter mustelae",
        "Aliarcobacter cryaerophilus",
        "Helicobacter cinaedi",
        "Helicobacter fennelliae",
        "Aliarcobacter butzleri",
        "Arcobacter nitrofigilis",
        "Helicobacter sp. CLO-3",
        "Firehammervirus CP220",
        "Firehammervirus CPt10",
        "Fletchervirus NCTC12673",
        "Fletchervirus CP81",
        "Fletchervirus CPX",
        "Firehammervirus CP21",
        "Fletchervirus CP30A",
        "Fletchervirus Los1",
    },
    "Escherichia coli": {
        "Tequintavirus AKFV33",
        "Enterobacteria phage CUS-3",
    },
    "Enterococcus faecium": {
        "Enterococcus casseliflavus",
    },
    "Enterobacter": {
        "Hafnia alvei",
        "Kosakonia radicincitans DSM 16656",
        "Kluyvera intermedia",
        "Cronobacter sakazakii",
        "Pluralibacter gergoviae",
        "Klebsiella aerogenes EA1509E",
        "Klebsiella aerogenes KCTC 2190",
        "Klebsiella aerogenes",
        "Pantoea agglomerans",
        "Lelliottia amnigena",
        "Lelliottia nimipressuralis",
        "Pluralibacter pyrinus",
        "Kosakonia radicincitans",
        "Siccibacter turicensis",
        "Franconibacter helveticus",
        "Franconibacter pulveris",
        "Kosakonia oryzae",
        "Kosakonia arachidis",
        "Kosakonia sacchari",
        "Kosakonia sacchari SP1",
        "Escherichia phage IME11",
        "Pluralibacter gergoviae ATCC 33028 = NBRC 105706",
        "Phytobacter massiliensis",
        "Atlantibacter hermannii",
        "Rahnella aquatilis",
        "Kosakonia cowanii",
        "Cronobacter sakazakii ATCC BAA-894",
        "Kosakonia oryzendophytica",
        "Kosakonia oryziphila",
        "Webervirus F20",
        "Franconibacter pulveris DSM 19144",
        "Hafnia phage Enc34",
        "Pseudenterobacter timonensis",
        "Karamvirus pg7",
        "Karamvirus cc31",
        "Slopekvirus eap3",
        "Eclunavirus EcL1",
        "Eapunavirus Eap1",
    },
    "Helicobacter pylori": {
        "Helicobacter mustelae",
    },
    "Klebsiella pneumoniae": {
        "Klebsiella michiganensis KCTC 1686",
        "Klebsiella variicola subsp. tropica",
    },
    "Mycobacterium tuberculosis": {
        "Mycobacterium avium",
        "Corynebacterium pseudotuberculosis",
    },
    "Pseudomonas aeruginosa": {
        "Pseudomonas virus Yua",
    },
    "Staphylococcus aureus": {
        "Dubowvirus dv11",
    },
}

# %%
# Remove rows with unwanted pathogens-taxonomy pairs
pathogens_taxid = pd.read_csv(DATA_PROCESSED / "00_pathogens_taxid.csv")

df = pathogens_taxid.copy()
keep = pd.Series(True, index=df.index)

for pathogen, taxonomy in wrong_taxonomies.items():
    # rows that match THIS pathogen AND have a name in THIS pathogen's bad set
    to_drop = (df["Pathogen"] == pathogen) & (df["Taxonomy_Name"].isin(taxonomy))
    
    # set those rows to False (drop them)
    keep &= ~to_drop

removed_pairs_df = df[keep].reset_index(drop=True)

# %%
# Remove rows with "phage" or "virus"
is_phage_or_virus = removed_pairs_df["Taxonomy_Name"].str.contains(
    r"phage|virus"
)

clean_df = removed_pairs_df[~is_phage_or_virus].reset_index(drop=True)

# %%
clean_df.to_csv(DATA_PROCESSED / "01_pathogens_taxid_cleaned.csv", index=False)

# %%
# Convert taxonomy table into dictionary
pathogens_taxid_cleaned_dict = (
    clean_df.groupby("Pathogen")["Taxonomy_ID"]
    .apply(list)
    .to_dict()
)

with open(DATA_PROCESSED / "02_pathogens_taxid_cleaned_dict.json", "w") as f:
    json.dump(pathogens_taxid_cleaned_dict, f, indent=2)

# %% [markdown]
# ## 02. Obtaining all Bioassays form PubcChem (AIDs) per pathogen

# %% [markdown]
# ### 2.1. AIDs per pathogen using PubChem user interface (UI)

# %%
# Add UI (user interface) counts (manually annotated 17.11.25)
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

AIDs_pathogen = (pd.DataFrame(ui_counts.items(), columns=["Pathogen", "UI_AIDs"])
    .sort_values("UI_AIDs", ascending=False)
    .reset_index(drop=True)
)

AIDs_pathogen

# %%
AIDs_pathogen.to_csv(DATA_PROCESSED / "03_aid_counts_per_pathogen.csv", index=False)


# %% [markdown]
# ### 2.2. AIDs per pathogen using PubChem Taxonomy search
#
# PubChem does not currently expose a stable API for retrieving organism-linked BioAssays directly from organism names. Therefore, the initial files used here were exported manually from:
#
# **PubChem → Search → "Organism name" → Taxonomy →  Download: *Summary (Search Results)**
#
# Each exported CSV (such as `PubChem_taxonomy_text_Acinetobacter baumannii.csv`) has been manually saved under `data\raw`.
#
# Each CSV contains a pipe-separated list of linked BioAssays. 
#
# We expand these to one row per AID and annotate with the pathogen name.
#
# Now, we will NOT filter out the wrong taxonomies we have manually filtered, in order to compare it with the UI_AIDs.

# %%
def taxid_aid(pathogen: str) -> pd.DataFrame:
    """Expand the PubChem Taxonomy→BioAssay mapping for a pathogen
    """

    df = pd.read_csv(DATA_RAW / f"PubChem_taxonomy_text_{pathogen}.csv")

    # Fill empty rows with empty strings "" (if no linked_bioassay, no rows for that pathogen to be counted)
    df["Linked_BioAssays"] = df["Linked_BioAssays"].fillna("")

    # Split into python list (from "1234|5678|9012" to ["1234", "5678", "9012"])
    df["AID"] = df["Linked_BioAssays"].astype(str).str.split("|")
    
    # Expand AIDs (one row per AID)
    df = df.explode("AID")

    df["Pathogen"] = pathogen

    return df[["Pathogen", "Taxonomy_ID", "Taxonomy_Name", "AID"]].drop_duplicates()


# %%
# Process all pathogens and create single table

pathogen_tables = []

for pathogen in pathogens:
    pathogen_table = taxid_aid(pathogen)
    pathogen_tables.append(pathogen_table)

all_pathogens = pd.concat(pathogen_tables, ignore_index=True)
all_pathogens = all_pathogens[["Pathogen", "Taxonomy_ID", "Taxonomy_Name", "AID"]]
all_pathogens.head()

# %%
# Count number of AIDs per pathogen
taxid_counts = (
    all_pathogens.groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="Taxonomy_AIDs")
    .sort_values("Taxonomy_AIDs", ascending=False)
    .reset_index(drop=True)
)

taxid_counts

# %%
# Merge with previous AID counts
AIDs_pathogen = AIDs_pathogen.merge(
    taxid_counts,
    on="Pathogen",
    how="left"   # keep all pathogens from AIDs_pathogen
)

AIDs_pathogen

# %%
AIDs_pathogen.to_csv(DATA_PROCESSED / "03_aid_counts_per_pathogen.csv", index=False)

# %%
# Extract arrays
labels = AIDs_pathogen["Pathogen"].values
UI_AIDs = AIDs_pathogen["UI_AIDs"].values
Taxonomy_AIDs = AIDs_pathogen["Taxonomy_AIDs"].values

N = len(labels)
x = np.arange(N)

bar_width = 0.25

plt.figure(figsize=(10, 4))

plt.bar(
    x - bar_width,
    UI_AIDs,
    width=bar_width,
    color="#D2D2D2",
    ec="black",
    zorder=2,
    label="PubChem Website Search"
)

plt.bar(
    x,
    Taxonomy_AIDs,
    width=bar_width,
    color="#AA96FA",
    ec="black",
    zorder=2,
    label="Manual Taxonomy"
)

# Aesthetics
plt.yscale("log")
plt.ylim([1, UI_AIDs.max() * 1.5])

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()


# %% [markdown]
# ### 2.3. SDQ Query for Bioassays aid & taxids
# PubChem provides a programmatic interface (PUG-REST) for retrieving BioAssay identifiers (AIDs) linked to NCBI Taxonomy IDs, but not directly from organism names.
#
# In this step, we will:
#
# 1.	Use the curated set of NCBI Taxonomy IDs associated with each pathogen (`pathogens_taxid_cleaned_dict.json`).
#
# 2.	Query PubChem’s taxonomy → taxid → AIDs endpoint via PUG-REST for each TaxID
#
# 3.	Retrieve the list of BioAssay IDs associated with each TaxID.
#
# 4.	Aggregate and deduplicate AIDs across all TaxIDs belonging to the same pathogen.
#
# This approach avoids manual UI exports and enables fully reproducible retrieval of organism-linked BioAssays directly from PubChem.
#
# Unlike the previous section, here we rely exclusively on curated Taxonomy IDs and programmatic access.

# %%
def get_aids_for_taxid(taxid: int):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/taxonomy/taxid/{taxid}/aids/JSON"

    try:
        r = requests.get(url, timeout=15)
    except requests.exceptions.RequestException:
        return []

    if r.status_code != 200:
        return []

    try:
        data = r.json()
        info = data["InformationList"]["Information"][0]
        return info.get("AID", [])
    except:
        return []


# %%
# For each pathogen, fetch curated AIDs (ETA: 25min)
pugrest_counts = {}

with open(DATA_PROCESSED / "02_pathogens_taxid_cleaned_dict.json", "r") as f:
    dict_taxonomy = json.load(f)

for pathogen, taxids in dict_taxonomy.items():
    aids = set()
    for tid in taxids:
        aids.update(get_aids_for_taxid(tid))
    pugrest_counts[pathogen] = len(aids)

pugrest_counts

# %%
# Merge with previous AID counts
pugrest_counts_df = (pd.DataFrame(pugrest_counts.items(), columns=["Pathogen", "PUGREST_AIDs"]))

AIDs_pathogen = AIDs_pathogen.merge(
    pugrest_counts_df,
    on="Pathogen",
    how="left"   # keep all pathogens from AIDs_pathogen
)

AIDs_pathogen

# %%
AIDs_pathogen.to_csv(DATA_PROCESSED / "03_aid_counts_per_pathogen.csv", index=False)

# %%
# Extract arrays
labels = AIDs_pathogen["Pathogen"].values
UI_AIDs = AIDs_pathogen["UI_AIDs"].values
Taxonomy_AIDs = AIDs_pathogen["Taxonomy_AIDs"].values
PUGREST_AIDs = AIDs_pathogen["PUGREST_AIDs"].values

N = len(labels)
x = np.arange(N)

bar_width = 0.25

plt.figure(figsize=(10, 4))

plt.bar(
    x - bar_width,
    UI_AIDs,
    width=bar_width,
    color="#D2D2D2",
    ec="black",
    zorder=2,
    label="PubChem Website Search"
)

plt.bar(
    x,
    Taxonomy_AIDs,
    width=bar_width,
    color="#AA96FA",
    ec="black",
    zorder=2,
    label="Manual Taxonomy"
)

plt.bar(
    x + bar_width,
    PUGREST_AIDs,
    width=bar_width,
    color="#FAD782",
    ec="black",
    zorder=2,
    label="PUGREST Taxonomy"
)

# Aesthetics
plt.yscale("log")
plt.ylim([1, UI_AIDs.max() * 1.5])

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()

# %% [markdown]
# ### 2.4. PubChem’s Aid2Taxid.tsv.gz file
#
# To obtain a complete list of all BioAssays linked to our taxonomy IDs, we now switch to PubChem’s official FTP file `Aid2Taxid.tsv.gz`, which contains the full mapping of every PubChem BioAssay to all associated NCBI Taxonomy IDs.
#
# This file was downloaded manually from:
#
# https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/Aid2Taxid.gz

# %%
# Load file
df_aid2tax = pd.read_csv(DATA_RAW / "Aid2Taxid.tsv", sep="\t")
df_aid2tax = df_aid2tax.rename(columns={"TaxID": "Taxonomy_ID"})
df_aid2tax.head(10)

# %%
# Build pathogen–TaxID lookup
pathogens_taxid = pd.read_csv(DATA_PROCESSED / "00_pathogens_taxid.csv")

df_merged = df_aid2tax.merge(
    pathogens_taxid[["Pathogen", "Taxonomy_ID"]],
    on="Taxonomy_ID",
    how="inner"
)

df_merged

# %%
# Count AIDs per pathogen
aid2tax_counts = (
    df_merged.groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="Aid2Taxid_AIDs")
    .sort_values("Aid2Taxid_AIDs", ascending=False)
    .reset_index(drop=True)
)

aid2tax_counts

# %%
# Merge with previous AID counts
AIDs_pathogen = AIDs_pathogen.merge(
    aid2tax_counts,
    on="Pathogen",
    how="left"   # keep all pathogens from AIDs_pathogen
)

AIDs_pathogen

# %%
AIDs_pathogen.to_csv(DATA_PROCESSED / "03_aid_counts_per_pathogen.csv", index=False)

# %%
# Extract arrays
labels = AIDs_pathogen["Pathogen"].values
UI_AIDs = AIDs_pathogen["UI_AIDs"].values
Taxonomy_AIDs = AIDs_pathogen["Taxonomy_AIDs"].values
PUGREST_AIDs = AIDs_pathogen["PUGREST_AIDs"].values
PUGREST_AIDs = AIDs_pathogen["PUGREST_AIDs"].values
Aid2Taxid_AIDs = AIDs_pathogen["Aid2Taxid_AIDs"].values

N = len(labels)
x = np.arange(N)

bar_width = 0.2

plt.figure(figsize=(12, 6))

plt.bar(
    x - bar_width,
    UI_AIDs,
    width=bar_width,
    color="#D2D2D2",
    ec="black",
    zorder=2,
    label="PubChem Website Search"
)

plt.bar(
    x,
    Taxonomy_AIDs,
    width=bar_width,
    color="#AA96FA",
    ec="black",
    zorder=2,
    label="Manual Taxonomy"
)

plt.bar(
    x + bar_width,
    PUGREST_AIDs,
    width=bar_width,
    color="#FAD782",
    ec="black",
    zorder=2,
    label="PUGREST Taxonomy"
)

plt.bar(
    x + 2*bar_width,
    Aid2Taxid_AIDs,
    width=bar_width,
    color="#FAA08B",
    ec="black",
    zorder=2,
    label="Aid2Taxid"
)

# Aesthetics
plt.yscale("log")
plt.ylim([1, UI_AIDs.max() * 1.5])

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()


# %% [markdown]
# The comparison shows that the **manual taxonomy** export captures the fewest assays, while both **PUG REST** and **Aid2Taxid.tsv** recover many more—almost identical to each other—because they rely on the same structured taxonomy annotations deposited in PubChem.
#
# However, both still fall below the counts shown on the **PubChem website**, meaning many assays mentioning a pathogen are not formally annotated with a `TaxID` and only appear in the website search because they contain the organism name somewhere in their free-text descriptions.

# %% [markdown]
# ### 2.5. Downloading ALL PubChem locally

# %% [markdown]
# #### 2.5.1. Understanding PubChem BioAssay FTP Folder Structure

# %% [markdown]
# PubChem stores all BioAssay files here:
#
# a) **Descriptions** (XML): https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/
#
# b) **Data** (CSV assay results): https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/
#
# Inside each folder files look like:

# %% [raw]
# CSV/
#     Description/
#         0000001_0001000.zip
#         0001001_0002000.zip
#         ...
#     Data/
#         0000001_0001000.zip
#         0001001_0002000.zip
#         ...

# %% [markdown]
# Each zip file corresponds to a range of AIDs. Example:
#
# `0000001_0001000.zip`  → contains assays with AID 1 to 1000
#
# In `Description/`, there are files such as:

# %% [markdown]
# 1.descr.xml.gz
#
# 2.descr.xml.gz
#
# 3.descr.xml.gz
#
# ...

# %% [markdown]
# Each `*.descr.xml` file follows the official PubChem BioAssay XML schema and may contain detailed metadata describing the assay. Depending on the assay, the XML can include:
#
# Always present (schema-required)
# - **Assay Name**: `<PC-AssayDescription_name>`
# - **Assay Description / Protocol text**: `<PC-AssayDescription_description>` (may include protocol steps, summary, conditions, etc.)
# - **Depositor / Source Information**: `<PC-AssayDescription_aid-source>`
#
#
# Present when deposited by submitter (schema-optional)
# - **Targets** (proteins, genes, taxonomy IDs):` <PC-AssayDescription_target>`
# - **Comments / Additional notes**: `<PC-AssayDescription_comment>`
# - **References** (PMIDs, DOIs, citation links): `<PC-AssayDescription_xref>`
# - **Relations to other assays**: `<PC-AssayDescription_relations>`

# %% [markdown]
# In `Data/`, there are files like:

# %% [markdown]
# 1.csv.gz
#
# 2.csv.gz
#
# 3.csv.gz
#
# ...

# %% [markdown]
# Each .csv contains the assay results:
#
# Columns 1–7:
# - `PUBCHEM_RESULT_TAG`
# - `PUBCHEM_SID`
# - `PUBCHEM_CID`
# - `PUBCHEM_ACTIVITY_OUTCOME`
# - `PUBCHEM_ACTIVITY_SCORE`
# - `PUBCHEM_ACTIVITY_URL`
# - `PUBCHEM_ASSAYDATA_COMMENT`
#
# Columns 8+:
# - depositor-defined results (IC50, % inhibition, etc.)

# %% [markdown]
# #### 2.5.2. Preparing download

# %%
# Calculate aprox zip file sizes to be downloaded

def get_total_ftp_size(url):
    print(f"Checking: {url}")
    r = requests.get(url).text

    # Matches lines like: "0000001_0001000.zip     120M"
    sizes = re.findall(r'\s+(\d+(?:\.\d+)?)([KMG])', r)

    total_bytes = 0
    for num, unit in sizes:
        num = float(num)
        if unit == "K": num *= 1e3
        if unit == "M": num *= 1e6
        if unit == "G": num *= 1e9
        total_bytes += num

    return total_bytes

def to_gb(bytes_val):
    return round(bytes_val / 1e9, 3)

desc_url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/"
data_url = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/"

desc_total = get_total_ftp_size(desc_url)
data_total = get_total_ftp_size(data_url)

print("\n=== FTP Total Size Summary ===")
print(f"Description/ total: {to_gb(desc_total)} GB")
print(f"Data/ total:         {to_gb(data_total)} GB")
print(f"Combined total:      {to_gb(desc_total + data_total)} GB")


# %%
# Count ZIPs in an FTP (File Transfer Protocol) PubChem directory
def list_zip_files(base_url):
    """Scrapes a PubChem FTP directory and returns list of .zip file names."""
    html = requests.get(base_url).text
    soup = BeautifulSoup(html, "html.parser")
    return [a.text for a in soup.find_all("a") if a.text.endswith(".zip")]

URL_DESC = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/"
URL_DATA = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/"

desc_zip_files = list_zip_files(URL_DESC)
data_zip_files = list_zip_files(URL_DATA)

len(desc_zip_files), len(data_zip_files)

# %% [markdown]
# #### 2.5.3. Downloading all PubChem Bioassays Description and Data

# %%
# Create a folder where everything will go

PUBCHEM_DIR = DATA_RAW / "pubchem_bioassays"
DESC_DIR = PUBCHEM_DIR / "Description"
DATA_DIR = PUBCHEM_DIR / "Data"

DESC_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# Due to the size and number of PubChem BioAssay CSV archives, downloading and managing these files directly inside a notebook is slow, fragile, and not reproducible.
#
# Therefore, all PubChem BioAssay CSV archives (Description and Data) are downloaded outside the notebook using a dedicated script `scripts/001_download_pubchem_bioassay_csv.py`.
#
# This script:
#
# - retrieves the full list of BioAssay CSV ZIP files from the PubChem FTP service,
# - supports resumable downloads via HTTP range requests
# - downloads files in parallel (default: 6 workers)
# - avoids re-downloading already completed files
# - stores all archives under:

# %% [raw]
# data/raw/pubchem_bioassays/
# ├── Description/
# └── Data/

# %%
# python scripts/001_download_pubchem_bioassay_csv.py --out data/raw/pubchem_bioassays --workers 6

# %% [markdown]
# #### 2.5.4. Filtering only bioassay files for pathogens of interest (TaxID and Organism)

# %% [markdown]
# We now filter the downloaded XML BioAssay descriptions to retain only those assays that match our pathogens of interest, using two complementary strategies:
# - Structured annotations (Taxonomy IDs) from:
#     - <PC-AssayTarget_tax-id>
#     - <PC-XRefData_taxonomy>
#     - <PC-AssayResultType_tax-id>
# - Free-text mentions from:
#     - <PC-AnnotatedXRef_comment> (often used for assay_organism)
#
# This dual filtering greatly increases recall of relevant AIDs for each pathogen.
#
# A dedicated script `scripts/002_filter_bioassay_descriptions.py` has been created to perform this filtering outside the notebook.
#
# It:
# - parses all uncompressed XML files in `data/raw/filtered_assays/Description/`
# - detects hits using both TaxIDs and free-text organism mentions
# - extracts additional metadata like ChEMBL IDs
# - outputs:
# 	- filtered_assays_description_taxid_organism.csv (AID-level metadata)
# 	- summary_xml_taxid_organism.csv (pathogen → AID counts)

# %% [markdown]
# #### 2.5.5. Counting downloaded AIDs per pathogen

# %%
# Load new filtered AIDs from script 002
Downloaded_AIDs = pd.read_csv(DATA_PROCESSED / "004b_filtered_aid_summary.csv")

Downloaded_AIDs

# %%
# Merge with previous AID counts
AIDs_pathogen = AIDs_pathogen.merge(
    Downloaded_AIDs,
    on="Pathogen",
    how="left"   # keep all pathogens from AIDs_pathogen
)

AIDs_pathogen

# %%
AIDs_pathogen.to_csv(DATA_PROCESSED / "03_aid_counts_per_pathogen.csv", index=False)

# %%
N = len(labels)
bar_width = 0.15
x = np.arange(N)

plt.figure(figsize=(14, 6))

# Center the 5 bars around x
offsets = [-2, -1, 0, 1, 2]  # for 5 bars
colors = ["#D2D2D2", "#AA96FA", "#FAD782", "#FAA08B", "#DC9FDC"]
datasets = [UI_AIDs, Taxonomy_AIDs, PUGREST_AIDs, Aid2Taxid_AIDs, Downloaded_AIDs]
labels_list = [
    "PubChem Website Search",
    "Manual Taxonomy",
    "PUGREST Taxonomy",
    "Aid2Taxid",
    "Downloaded"
]

for i, (data, color, label) in enumerate(zip(datasets, colors, labels_list)):
    plt.bar(
        x + offsets[i] * bar_width,
        data,
        width=bar_width,
        color=color,
        ec="black",
        zorder=2,
        label=label
    )

# Aesthetics
plt.yscale("log")
plt.ylim([1, UI_AIDs.max() * 1.5])

plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of AIDs (log scale)")

plt.grid(linestyle="--", zorder=1)
plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 9})

plt.tight_layout()
plt.show()

# %% [markdown]
# We select the **Downloaded strategy** for downstream analysis, as it most closely aligns with the PubChem UI results across pathogens while relying on a reproducible filtering approach based on both taxonomy IDs and organism mentions. With this refined set of AIDs, we now proceed to explore the **metadata** of each bioassay in detail.

# %% [markdown]
# ## 03. Descriptors of interest (Display files)

# %% [markdown]
# The information of each bioassay displayed using the UI is stored as a `Display` file and can be also download in JSON format. In this file we cand find descriptors of interest such as:
# - PubChem_AID	
# - ChEMBL_ID	
# - Compounds_Tested	
# - Compounds_Active	
# - Compounds_Inactive	
# - Target	
# - Assay_Type	
# - Assay_Format	
# - Assay_Organism	
# - Organism_TaxID	
# - Assay_Strain	
# - Organism_Target	
# - Protein_Target	
# - Source

# %% [markdown]
# `Display` JSON files provide detailed metadata for each bioassay. To enrich the filtered AIDs from `scripts/002_filter_bioassay_descriptions.py`, we use `scripts/003_download_display_jsons_parallel.py`, which:
# - downloads Display JSON files from PubChem for each filtered AID
# - supports incremental resume (skips already downloaded/parsed files)
# - uses parallel processing (default: 6 workers) for speed
# - extracts key metadata, including:
# 	- PubChem_AID, ChEMBL_ID
# 	- Compounds_Tested, Active, Inactive
# 	- Assay_Type, Format, Organism, Strain, TaxID
# 	- Target name/type, gene symbol, source
# - saves parsed metadata to `data/processed/06_display_info.csv`
# - stores raw JSON files in `data/raw/filtered_assays_v3/Display/`

# %%
# python scripts/003_download_display_jsons_parallel.py

# %%
display_info = pd.read_csv(DATA_PROCESSED / "06_display_info.csv")
display_info.head(10)

# %% [markdown]
# Are there any inconsistances between Target and Assay_Organism from 06_display_info.csv with the Pathogen from 05_filtered_aids.csv for each AID?

# %%
# Paths to the processed CSVs
description_info = pd.read_csv(DATA_PROCESSED / "05_filtered_aids.csv")
display_info = pd.read_csv(DATA_PROCESSED / "06_display_info.csv")

# Merge on AID (PubChem_AID in display_info)
df_merged = description_info.merge(display_info, left_on="AID", right_on="PubChem_AID", how="inner")

print(f"Merged shape: {df_merged.shape}")
df_merged.head()

# %%
df_merged.to_csv(DATA_PROCESSED / "07_filtered_aids_metadata.csv", index=False)

# %% [markdown]
#  Check if “Pathogen” matches either “Target” or “Assay_Organism”

# %%
# Add two columns to flag matching

# Normalize strings to lowercase for comparison
df_merged["pathogen_in_target"] = df_merged.apply(
    lambda row: str(row["Pathogen"]).lower() in str(row["Target"]).lower() if pd.notnull(row["Target"]) else False,
    axis=1
)

df_merged["pathogen_in_organism"] = df_merged.apply(
    lambda row: str(row["Pathogen"]).lower() in str(row["Assay_Organism"]).lower() if pd.notnull(row["Assay_Organism"]) else False,
    axis=1
)

# Optional: overall match
df_merged["pathogen_match_any"] = df_merged[["pathogen_in_target", "pathogen_in_organism"]].any(axis=1)

# Preview
df_merged[["AID", "Pathogen", "Target", "Assay_Organism", "pathogen_in_target", "pathogen_in_organism", "pathogen_match_any"]].head()

# %%
# Count matches
summary = df_merged["pathogen_match_any"].value_counts(normalize=True) * 100
print("% of AIDs where pathogen matches either target or organism:\n")
print(summary)

# %%
non_matching = df_merged[df_merged["pathogen_match_any"] == False]

# Display relevant columns
non_matching[["AID", "Pathogen", "Target", "Assay_Organism"]]

# %%
from IPython.display import display
import pandas as pd

# Show only the relevant columns
non_matching_display = non_matching[["AID", "Pathogen", "Target", "Assay_Organism"]]

# Scrollable display
with pd.option_context('display.max_rows', None, 'display.max_colwidth', None):
    display(non_matching_display)

# %% [markdown]
# **What should we do? Remove the non-matching ones? Correct them manually?**

# %% [markdown]
# ## 04. PubChem-Chembl comparison

# %%
pubchem_assays = pd.read_csv(DATA_PROCESSED / "07_filtered_aids_metadata.csv")
pubchem_assays

# %%
chembl_assays = pd.read_csv(DATA_PROCESSED / "ChEMBL_bioassays/assays_ChEMBL_all.csv") # This file was produced by Arnau
chembl_assays

# %%
chembl_pathogens = sorted(chembl_assays["pathogen"].dropna().unique())
chembl_pathogens

# %%
# Define mapping from ChEMBL codes to canonical names
chembl_to_full = {
    "abaumannii": "Acinetobacter baumannii",
    "calbicans": "Candida albicans",
    "ecoli": "Escherichia coli",
    "efaecium": "Enterococcus faecium",
    "enterobacter": "Enterobacter",
    "hpylori": "Helicobacter pylori",
    "kpneumoniae": "Klebsiella pneumoniae",
    "mtuberculosis": "Mycobacterium tuberculosis",
    "ngonorrhoeae": "Neisseria gonorrhoeae",
    "paeruginosa": "Pseudomonas aeruginosa",
    "pfalciparum": "Plasmodium falciparum",
    "saureus": "Staphylococcus aureus",
    "smansoni": "Schistosoma mansoni",
    "spneumoniae": "Streptococcus pneumoniae",
}

# %%
# Map ChEMBL shorthand to canonical names
chembl_assays["Pathogen"] = chembl_assays["pathogen"].map(pathogen_map)

# Drop rows where the mapping failed (optional)
chembl_assays = chembl_assays.dropna(subset=["Pathogen"])

# %% [markdown]
# ### 4.1. Bioassay counts

# %%
# PubChem
pubchem_counts = (
    pubchem_assays
    .groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="PubChem_Assays")
)

# ChEMBL (now with Pathogen column)
chembl_counts = (
    chembl_assays
    .groupby("Pathogen")["assay_id"]
    .nunique()
    .reset_index(name="ChEMBL_Assays")
)

# Merge + format
summary = (
    pubchem_counts
    .merge(chembl_counts, on="Pathogen", how="outer")
    .fillna(0)
)

summary[["PubChem_Assays", "ChEMBL_Assays"]] = summary[["PubChem_Assays", "ChEMBL_Assays"]].astype(int)
summary["Total"] = summary["PubChem_Assays"] + summary["ChEMBL_Assays"]
summary = summary.sort_values("Total", ascending=False).reset_index(drop=True)

# %%
# Data
labels = summary["Pathogen"].values
pubchem_vals = summary["PubChem_Assays"].values
chembl_vals = summary["ChEMBL_Assays"].values

N = len(labels)
x = np.arange(N)

# ⏩ Bar styling
bar_width = 0.35
bar_spacing = 0.2

color_pubchem = "#AA96FA"  # purple 
color_chembl  = "#BEE6B4"    # Green

# Plot
plt.figure(figsize=(14, 6))

plt.bar(x - bar_spacing, pubchem_vals, width=bar_width, color=color_pubchem,
        ec="black", zorder=2, label="PubChem Assays")

plt.bar(x + bar_spacing, chembl_vals, width=bar_width, color=color_chembl,
        ec="black", zorder=2, label="ChEMBL Assays")

# Add labels
for i in range(N):
    plt.text(x[i] - bar_spacing, pubchem_vals[i] + 1, str(pubchem_vals[i]),
             ha='center', va='bottom', fontsize=8)
    plt.text(x[i] + bar_spacing, chembl_vals[i] + 1, str(chembl_vals[i]),
             ha='center', va='bottom', fontsize=8)

# Aesthetics
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of Assays")
plt.title("Number of Assays per Pathogen (PubChem vs ChEMBL)")
plt.grid(axis="y", linestyle="--", zorder=1)
plt.ylim(0, max(pubchem_vals.max(), chembl_vals.max()) * 1.2)

plt.legend(loc="upper right", framealpha=1, edgecolor="k", prop={"size": 10})
plt.tight_layout()
plt.show()

# %% [markdown]
# ### 4.1. % Chembl_ID present in PubChem bioassays

# %%
# 1. Filter rows with ChEMBL IDs
aids_with_chembl = pubchem_assays[pubchem_assays["ChEMBL_ID"].notna()]

# 2. Group total AIDs and AIDs with ChEMBL by pathogen
total_aids = pubchem_assays.groupby("Pathogen")["AID"].nunique().reset_index(name="Total_AIDs")
with_chembl = aids_with_chembl.groupby("Pathogen")["AID"].nunique().reset_index(name="AIDs_with_ChEMBL")

# 3. Merge and compute percentage
summary = total_aids.merge(with_chembl, on="Pathogen", how="left")
summary["AIDs_with_ChEMBL"] = summary["AIDs_with_ChEMBL"].fillna(0).astype(int)
summary["Percent_ChEMBL_IDs"] = 100 * summary["AIDs_with_ChEMBL"] / summary["Total_AIDs"]
summary = summary.sort_values("Percent_ChEMBL_IDs", ascending=False).reset_index(drop=True)

# %%
# Data
labels = summary["Pathogen"].values
percent_chembl = summary["Percent_ChEMBL_IDs"].values

N = len(labels)
x = np.arange(N)

# Plot settings
bar_width = 0.6
color_chembl = "#AA96FA"

plt.figure(figsize=(14, 6))

bars = plt.bar(x, percent_chembl, width=bar_width,
               color=color_chembl, ec="black", zorder=2,
               label="% AIDs with ChEMBL ID")

# Add labels on top
for i in range(N):
    plt.text(x[i], percent_chembl[i] + 2, f"{percent_chembl[i]:.1f}%",
             ha='center', va='bottom', fontsize=9, zorder=3)

# Aesthetics
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Percentage (%)")
plt.ylim(0, 110)
plt.title("Percentage of AIDs with ChEMBL ID per Pathogen")
plt.grid(axis="y", linestyle="--", zorder=1)
plt.tight_layout()
plt.show()

# %% [markdown]
# ### 4.2. Number of compounds per assay

# %% [markdown]
# In PubChem terminology, a **substance** is a chemical sample description provided by a single source and a **compound** is a normalized chemical structure representation found in one or more contributed substances.
#
# PubChem calls the community-provided sample descriptions *substances*.  Each record found in the PubChem **Substance database** (https://www.ncbi.nlm.nih.gov/pcsubstance) contains information provided by an individual contributor about a particular chemical substance.  Substance records are independent of each other.  Two different Substance records (from the same or different providers) could provide different information about the same chemical structure.  
#
# For example, one substance record may give information about the biological role of aspirin, while another may give information about a research grade sample of aspirin.  The Substance database maintains the provenance of chemical substance information in PubChem.  It helps users see who provided what.  
#
# As a result, there may be many substance records about a given molecule, presenting a problem for users who are interested in an aggregated view of information on the molecule.  This is where the PubChem Compound database (https://www.ncbi.nlm.nih.gov/pccompound) comes into play.
#
# The **Compound database** is derived from the chemical structure contents found in the Substance database.  Each chemical is computationally examined with a series of validation and normalization steps.  This process results in a normalized representation of the chemical structure for a substance record.  
#
# Chemical substances in the Substance database that are not completely described or that fail normalization procedures are not included in the Compound database.  Those substances in the Substance database that pass chemical structure normalization procedures are linked to a “compound” record in the Compound database.  If two substances refer to the same chemical structure, they point to the same compound.  This allows data from different Substance data providers to be aggregated through a common Compound record.  However, also having separate substance records is still valuable to users, who, for example, might be interested in the provenance of a substance or a particular state of the chemical (e.g., a different tautomeric form).  
#
# In essence, a primary purpose of the PubChem Compound database is to provide a “non-redundant” view of the depositor-contributed chemical structure contents stored in the PubChem Substance database.

# %%
pubchem_assays = pd.read_csv(DATA_PROCESSED / "07_filtered_aids_metadata.csv")
chembl_assays = pd.read_csv(DATA_PROCESSED / "ChEMBL_bioassays/assays_ChEMBL_all.csv") 

# %%
# From PubChem
pubchem_sub = pubchem_assays[
    [
        "Pathogen",
        "PubChem_AID",
        "ChEMBL_ID",
        "Compounds_Tested",
        "Compounds_Active",
        "Compounds_Inactive",
        "Tested_Substances"
    ]
]

# Keep only assays with a ChEMBL_ID in PubChem
pubchem_sub = pubchem_sub[
    pubchem_sub["ChEMBL_ID"].notna()
]

# Keep only assays with a compound counts available
pubchem_sub = pubchem_sub[
    pubchem_sub["Compounds_Tested"].notna()
]

# From ChEMBL
chembl_sub = chembl_assays[
    ["assay_id", "cpds"]
]

# Rename assay_id → ChEMBL_ID for merging
chembl_sub = chembl_sub.rename(columns={"assay_id": "ChEMBL_ID"})

# %%
# Merge on ChEMBL_ID (We exclude assays only present either in PubChem or in Chembl selection)

compounds = (
    pubchem_sub
    .merge(
        chembl_sub,
        on="ChEMBL_ID",
        how="inner"
    )
)

# %%
# Ensure numeric types

num_cols = ["Compounds_Tested", "Tested_Substances", "cpds"]

for col in num_cols:
    compounds[col] = pd.to_numeric(compounds[col], errors="coerce")

# Compute compound difference
compounds["Compounds_diff"] = (
    compounds["Compounds_Tested"] - compounds["cpds"]
)

compounds["Substances_diff"] = (
    compounds["Tested_Substances"] - compounds["cpds"]
)

compounds

# %%
# How many assays have different compound counts
assays_compunds_diff = compounds[compounds["Compounds_diff"] != 0]
assays_substances_diff = compounds[compounds["Substances_diff"] != 0]

print(len(assays_compunds_diff))
print(len(assays_substances_diff))

# %%
assays_compunds_diff.groupby("Pathogen").size().sort_values(ascending=False)

# %%
assays_substances_diff.groupby("Pathogen").size().sort_values(ascending=False)


# %%
# Classification of bioassays depending on whether cpds matches either Compounds_Tested or Tested_Substances.

def classify_assay(row):
    cpds = row["cpds"]
    comp = row["Compounds_Tested"]
    subs = row["Tested_Substances"]

    if pd.isna(cpds):
        return "no_cpds"

    matches_comp = not pd.isna(comp) and cpds == comp
    matches_subs = not pd.isna(subs) and cpds == subs

    if matches_comp and matches_subs:
        return "cpds_match_both"
    elif matches_comp:
        return "cpds_match_only_compounds"
    elif matches_subs:
        return "cpds_match_only_substances"
    else:
        return "cpds_mismatch"

compounds["cpds_match_class"] = compounds.apply(classify_assay, axis=1)

# %%
compounds["cpds_match_class"].value_counts()

# %%
# Count assays per pathogen and class
counts = compounds.groupby(["Pathogen", "cpds_match_class"]).size().unstack(fill_value=0)

# Define order and colors
match_order = [
    "cpds_match_both",
    "cpds_match_only_compounds",
    "cpds_match_only_substances",
    "cpds_mismatch",
]
colors = {
    "cpds_match_both": "#BEE6B4",             # Green
    "cpds_match_only_compounds": "#AA96FA",   # Purple
    "cpds_match_only_substances": "#FAD782",  # Yellow
    "cpds_mismatch": "#FAA08B",               # Orange
}

# Reorder columns and normalize to percentage
counts = counts.reindex(columns=match_order, fill_value=0)
percentages = counts.div(counts.sum(axis=1), axis=0) * 100

# Plot
labels = percentages.index.tolist()
N = len(labels)
x = np.arange(N)
bar_width = 0.8

plt.figure(figsize=(14, 6))

bottom = np.zeros(N)

for match_class in match_order:
    plt.bar(
        x,
        percentages[match_class],
        bottom=bottom,
        width=bar_width,
        color=colors[match_class],
        ec="black",
        label=match_class.replace("_", " ").title(),
        zorder=2
    )
    bottom += percentages[match_class].values

# Aesthetics
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Percentage of Assays (%)")
plt.title("Pubchem (compound/substance) vs Chembl (cpds) counts match")
plt.ylim(0, 100)
plt.grid(axis="y", linestyle="--", zorder=1)
plt.legend(loc="lower left", framealpha=1, edgecolor="k", prop={"size": 9})
plt.tight_layout()
plt.show()

# %%
# Filter for mismatch cases
mismatch = compounds[
    (compounds["cpds_match_class"] == "cpds_mismatch") &
    compounds["Compounds_diff"].notna()
]

# Plot
plt.figure(figsize=(14, 6))

sns.stripplot(
    data=mismatch,
    x="Pathogen",
    y="Compounds_diff",
    jitter=True,
    alpha=0.7,
    size=6,
    linewidth=0.5,
    edgecolor="black",
    color="#FAA08B"  # Match the orange used for mismatch
)

plt.axhline(0, linestyle="--", color="gray")
plt.title("Compound Count Difference per Assay (Only Mismatches)", fontsize=14)
plt.ylabel("Compound Count Difference (PubChem − ChEMBL)", fontsize=12)
plt.xlabel("Pathogen", fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()

# %%
# Filter and sort
mismatch_sorted = compounds[
    (compounds["cpds_match_class"] == "cpds_mismatch") &
    (compounds["Compounds_diff"].notna())
].sort_values("Compounds_diff", ascending=False)

# Show the full DataFrame with all columns
mismatch_sorted.reset_index(drop=True, inplace=True)
mismatch_sorted


# %% [markdown]
# ### 4.3. Active vs Inactive Compounds

# %%
# Classify each assay
def classify_activity(row):
    active = pd.notna(row["Compounds_Active"])
    inactive = pd.notna(row["Compounds_Inactive"])

    if active and inactive:
        return "both"
    elif active:
        return "only_active"
    elif inactive:
        return "only_inactive"
    else:
        return "none"

# Apply classification
compounds["activity_info_class"] = compounds.apply(classify_activity, axis=1)

# Count per pathogen and class
counts = compounds.groupby(["Pathogen", "activity_info_class"]).size().unstack(fill_value=0)

# Define class order and colors
class_order = ["both", "only_active", "only_inactive", "none"]
colors = {
    "both": "#BEE6B4",         # Green
    "only_active": "#AA96FA",  # Purple
    "only_inactive": "#FAD782",# Yellow
    "none": "#FAA08B",         # Red-ish
}

# Reindex and convert to percentages
counts = counts.reindex(columns=class_order, fill_value=0)
percentages = counts.div(counts.sum(axis=1), axis=0) * 100

# Plotting
labels = percentages.index.tolist()
N = len(labels)
x = np.arange(N)
bar_width = 0.8

plt.figure(figsize=(14, 6))
bottom = np.zeros(N)

for cls in class_order:
    plt.bar(
        x,
        percentages[cls],
        bottom=bottom,
        width=bar_width,
        color=colors[cls],
        edgecolor="black",
        label=cls.replace("_", " ").title(),
        zorder=2
    )
    bottom += percentages[cls].values

# Aesthetics
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Percentage of Assays (%)")
plt.title("Presence of Active/Inactive Compound Info in Pubchem")
plt.ylim(0, 100)
plt.grid(axis="y", linestyle="--", zorder=1)
plt.legend(loc="upper left", framealpha=1, edgecolor="k", prop={"size": 9})
plt.tight_layout()
plt.show()

# %%
compounds.to_csv(DATA_PROCESSED / "08_pubchem_vs_chembl_assays.csv", index=False)
