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
# # Downloading ALL PubChem locally

# %% [markdown]
# ## 00. Setup

# %%
import requests
import re
from bs4 import BeautifulSoup
from tqdm import tqdm
from pathlib import Path
import os

# %%
# Project paths
NOTEBOOK_DIR = Path().resolve()
DATA_RAW = NOTEBOOK_DIR.parent / "data" / "raw"
DATA_PROCESSED = NOTEBOOK_DIR.parent / "data" / "processed"


# %% [markdown]
# ### Understand the PubChem BioAssay FTP Folder Structure

# %% [markdown]
# PubChem stores all BioAssay files here:
#
# **Descriptions** (XML): https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/
#
# **Data** (CSV assay results): https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/
#
# Inside each folder files look like:

# %% [raw] vscode={"languageId": "raw"}
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

# %% [raw]
# 1.descr.xml.gz
# 2.descr.xml.gz
# 3.descr.xml.gz
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

# %% [raw]
# 1.csv.gz
# 2.csv.gz
# 3.csv.gz
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
#

# %% [markdown]
# ### Calculate aprox file size

# %%
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

# %% [markdown]
# ## 1. Download all BioAssays

# %%
# Create a folder where everything will go

PUBCHEM_DIR = DATA_RAW / "pubchem_bioassays"
DESC_DIR = PUBCHEM_DIR / "Description"
DATA_DIR = PUBCHEM_DIR / "Data"

DESC_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

PUBCHEM_DIR, DESC_DIR, DATA_DIR


# %%
# List ZIPs in an FTP (File Transfer Protocol) PubChem directory
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


# %%
# Function to download all .zip files

def download_zip(url, output_path):
    """
    Download a ZIP file with *full resume support* and a progress bar.

    Behaviors:
    - If the file does NOT exist → start fresh.
    - If the file exists but is incomplete → resume from last byte.
    - If the file is already complete → skip.
    - Download in 1 MB chunks to protect RAM.
    - Shows a live tqdm progress bar.
    """

    local_path = Path(output_path)   # Transform filesystem path to Path object → gives .exists(), .stat(), .name, etc.

    # ---------------------------------------------------------
    # 1. DEFAULT VALUES (fresh download unless resume detected)
    # ---------------------------------------------------------
    mode = "wb"        # "wb" = write binary (start new file)
    headers = {}       # HTTP headers (empty until resume needed)

    # ---------------------------------------------------------
    # 2. CHECK IF FILE ALREADY EXISTS → POSSIBLE RESUME
    # ---------------------------------------------------------
    if local_path.exists():
        existing_size = local_path.stat().st_size
        print(f"→ Found existing file ({existing_size} bytes): {local_path.name}")

        # HEAD request = ask server for metadata WITHOUT downloading the file (file size, last modified date, etc.)
        head = requests.head(url)

        total_size = (
            int(head.headers["Content-Length"])
            if "Content-Length" in head.headers
            else None                                   # If not provided, resume may not work
        )

        # If file is already complete → skip download
        if total_size is not None and existing_size >= total_size:
            print(f"✓ Already fully downloaded: {local_path.name}")
            return

        # Otherwise → resume from last byte downloaded
        headers["Range"] = f"bytes={existing_size}-"    # Request remaining bytes
        mode = "ab"                                     # Append mode (do not overwrite)

    # ---------------------------------------------------------
    # 3. FRESH DOWNLOAD CASE (file not previously downloaded or not able to resume)
    # ---------------------------------------------------------
    else:
        print(f"↓ Starting fresh download: {local_path.name}")

    # ---------------------------------------------------------
    # 4. STREAM DOWNLOAD WITH PROGRESS BAR
    # ---------------------------------------------------------
    with requests.get(url, headers=headers, stream=True) as r: # stream=True means download in chunks, not all at once (prevent RAM overload)
        r.raise_for_status()   # Stop if URL not accessible

        # Determine total size for progress bar
        if "Content-Length" in r.headers:
            total_size = int(r.headers["Content-Length"])
        else:
            total_size = None

        # tqdm settings
        chunk_size = 1024 * 1024   # 1 MB chunks
        progress = tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            initial=(existing_size if local_path.exists() else 0),
            desc=f"Downloading {local_path.name}",
        )

        # Download loop
        with open(local_path, mode) as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    progress.update(len(chunk))   # Update progress bar

        progress.close()

    print(f"✔ Done: {local_path.name}")


# %%
# Download ALL Description ZIP files (ETA: 69 min)

print("\n==============================")
print("Downloading DESCRIPTION files")
print("==============================")

for zip_name in tqdm(desc_zip_files):
    url = URL_DESC + zip_name
    out = DESC_DIR / zip_name
    download_zip(url, out)

print("\n🎉 ALL PubChem BioAssay Description ZIP files downloaded!")

# %%
# Download ALL Data ZIP files (ETA: 135 min)

print("\n==============================")
print("Downloading DATA files")
print("==============================")

for zip_name in tqdm(data_zip_files):
    url = URL_DATA + zip_name
    out = DATA_DIR / zip_name
    download_zip(url, out)

print("\n🎉 ALL PubChem BioAssay Data ZIP files downloaded!")

# %%
