
# 🦠 Antimicrobial binary ML tasks from PubChem 💊

Get antimicrobial tasks from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) framed as binary classifications.

This repository is currently **WORK IN PROGRESS**. ⚠️🚧

---

## Background

PubChem is the largest public repository of bioassay data, including **hundreds of thousands of biological assays (AIDs)** that test compounds against pathogens, proteins, and organisms. However, extracting **pathogen-specific antimicrobial assays** from PubChem is non-trivial due to:

- inconsistent taxonomy annotations,
- partial or missing NCBI Taxonomy IDs,
- heavy reliance on free-text fields,
- multiple overlapping access methods (UI search, taxonomy search, APIs, FTP dumps).

This project builds a **reproducible, high-recall pipeline** to:
1. identify PubChem bioassays relevant to antimicrobial pathogens,
2. extract standardized assay metadata,
3. prepare the ground for **binary ML task construction**.

---

## Tracking details

The project is tracked by Git (mainly for code) and DVC (mainly for data):

* Tracked by Git and linked to a Github repository: only src, scripts and notebooks.
* Tracked by DVC and linked to a Google Drive folder inside "Projects/<<Repository name>>".

---

## Repository structure

```
pubchem-antimicrobial-tasks/
│
├── LICENSE
├── README.md
├── .gitignore
├── install.sh
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   ├── PubChem_taxonomy_text_*.csv
│   │   ├── Aid2Taxid.tsv
│   │   └── pubchem_bioassays/
│   │       ├── Description/
│   │       └── Data/
│   │
│   └── processed/
│       ├── 00_pathogens_taxid.csv
│       ├── 01_pathogens_taxid_cleaned.csv
│       ├── 02_pathogens_taxid_cleaned_dict.json
│       ├── 03_aid_counts_per_pathogen.csv
│       ├── 004b_filtered_aid_summary.csv
│       ├── 05_filtered_aids.csv
│       ├── 06_display_info.csv
│       ├── 07_filtered_aids_metadata.csv
│       └── 08_pubchem_vs_chembl_assays.csv
│
├── scripts/
│   ├── 001_build_pathogen_taxonomy.py
│   ├── 002_download_pubchem_bioassay_csv.py
│   ├── 002_filter_bioassay_descriptions.py
│   └── 003_download_display_jsons_parallel.py
│
├── notebooks/
│   └── 001_pubchem_bioassays_pathogens_of_interest.ipynb
│
├── assets/
├── output/
│   ├── results/
│   └── plots/
│
├── src/
├── tools/
├── docs/
├── tmp/
│
└── .git/
```

📌 Empty folders are preserved with `.gitkeep`.

---

## Project motivation and goal

The ultimate goal is to build **binary antimicrobial ML tasks** (active vs inactive) from PubChem assays, aligned with ChEMBL where possible.

To do this reliably, we must first answer:

- *Which PubChem assays truly target a pathogen?*
- *How many compounds were tested per assay?*
- *How consistent are PubChem and ChEMBL compound counts?*
- *Which assays contain enough metadata for ML?*

This repository focuses on building the **data foundation** needed to answer these questions.

---

## 🚀 Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/ersilia-os/pubchem-antimicrobial-tasks.git
cd pubchem-antimicrobial-tasks
```

### 2. Install dependencies
```bash
bash install.sh
```
or
```bash
pip install -r requirements.txt
```

---

## Pipeline overview

The **final pipeline** consists of **three main stages**, implemented via scripts and documented in `notebooks/001_pubchem_bioassays_pathoegns_of_interest.ipynb`.

---

### 1️⃣ Build Pathogen Taxonomy Table

### Purpose
The first step is to create a high-confidence mapping between each pathogen and its associated **NCBI Taxonomy IDs**, which is required to query PubChem BioAssays in a structured way.

PubChem doesn’t expose a stable API for this, so we use a hybrid manual + automated approach.

### Inputs

The first step is to build a **curated mapping** between each **pathogen name** and the associated **NCBI Taxonomy IDs**, which are essential for programmatically retrieving PubChem BioAssays.

**Manual Step: Download Pathogen Taxonomy Summaries**

From the PubChem website:

PubChem → Search → Organism Name → “Taxonomy” → Export → Summary (CSV)

You should repeat this for each of the 15 pathogens listed below:

````
Acinetobacter baumannii
Candida albicans
Campylobacter
Escherichia coli
Enterococcus faecium
Enterobacter
Helicobacter pylori
Klebsiella pneumoniae
Mycobacterium tuberculosis
Neisseria gonorrhoeae
Pseudomonas aeruginosa
Plasmodium falciparum
Staphylococcus aureus
Schistosoma mansoni
Streptococcus pneumoniae
````

Each file is named like:

```bash
data/raw/PubChem_taxonomy_text_<Pathogen>.csv
```

Example:

```bash
PubChem_taxonomy_text_Acinetobacter baumannii.csv
```

Manually downloaded .csv files must be saved to:

```
data/raw/
```

### Scripts
Once all .csv files are in data/raw/, run:

```bash
python scripts/001_build_pathogen_taxonomy.py
```

### Outputs

- `00_pathogens_taxid.csv` → Raw merged table from all .csv manual exports
- `01_pathogens_taxid_cleaned.csv` → After filtering out irrelevant taxonomies (e.g. viruses, phages, other organisms)
- `02_pathogens_taxid_cleaned_dict.json` → Python dictionary format: {Pathogen: [TaxID1, TaxID2, ...]}

These are used downstream to retrieve AIDs via multiple strategies.

---

### 2️⃣ Downloading ALL PubChem BioAssays locally

### Purpose

We evaluated several strategies to retrieve pathogen-linked bioassays:

- **PubChem UI** search captures the most assays by searching free text across all fields, including references and comments.

- Other methods like **Taxonomy exports**, **PUG REST API**, and **Aid2Taxid** use structured links via NCBI Taxonomy IDs — but miss many relevant assays.

Ultimately, we selected the **Download All** strategy:

- It retrieves ~50% of the UI assay count
- But it’s fully automated, reproducible, and uniquely allows:
  - Matching both Taxonomy IDs
  - Searching structured metadata fields (e.g., *Target, Assay Organism, Strain*)

This comparison is detailed in `notebooks/001_pubchem_bioassays_pathogens_of_interest.ipynb`.

### Scripts

To download all .zip files from the PubChem BioAssay FTP (Description + Data) run:

```bash
scripts/002_download_pubchem_bioassay_csv.py
```

Then, to parse the downloaded files and filter for pathogen-linked assays using both Taxonomy IDs and organism mentions run:

```bash
scripts/003_filter_bioassay_descriptions.py
```

### Outputs
- `*.xml` → Filtered BioAssay XML files matching pathogens (one per AID)
- `04_filtered_aid_summary.csv` → Summary table of all matched AIDs with `AID`, `Pathogen`, `ChEMBLid` & `ZipFolder`
- `05_filtered_aids.csv` → Long-format table with one row per (AID, Pathogen) pair
- `processed_zips.txt` → List of completed ZIP chunks, used to resume work without reprocessing

---

### 3️⃣ Descriptors of interest (Display files)

### Script
- `003_download_display_jsons_parallel.py`

### Output
- `06_display_info.csv`
- `07_filtered_aids_metadata.csv`

---

## PubChem ↔ ChEMBL comparison

Compound and assay-level comparison between PubChem and ChEMBL, including mismatch analysis.

### Output
- `08_pubchem_vs_chembl_assays.csv`

---

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South.  
Ersilia’s main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia).

![Ersilia Logo](assets/Ersilia_Brand.png)
