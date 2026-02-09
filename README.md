
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

The project is tracked by Git (for code) and EOSVC (for data). The data folder does not have a version control, so be careful when downloading and uploading it. Read the eosvc [documentation](https://github.com/ersilia-os/eosvc) for more information.

---

## Repository structure

```
pubchem-antimicrobial-tasks/
│
├── LICENSE
├── README.md
├── .gitignore
│
├── data/
|   ├── config/ _data only available upon manual download or originating from another source, do not delete or edit_
│   ├── raw/ _data automatically downloaded via scripts from the internet_
│   └── processed/ _data from curation processes_
├── scripts/ _numbered scripts to run the pipeline_
│
├── notebooks/ _easy to look at checks and tests_
├── assets/
├── output/
│   ├── results/ _data cleaned for modelling_
│   └── plots/ _analysis results_
│
└── .git/
```
---

## Project motivation and goal

The ultimate goal is to build **binary antimicrobial ML tasks** (active vs inactive) from PubChem assays, without duplicating the efforts done in ChEMBL data curation. Therefore, we only keep assays in PubChem that:

1. Contain more than 100 datapoints
2. Contain more than 10% of data difference (number of compounds) between PubChem and ChEMBL (with more data in PubChem) or
3. Are not present in ChEMBL at all

### Pathogens of interest

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

---

## 🚀 Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/ersilia-os/pubchem-antimicrobial-tasks.git
cd pubchem-antimicrobial-tasks
eosvc download --path data
```
---

## Pipeline overview

### 1. Obtain config data
The data necessary to run the pipeline comes from either manual curation or from other pipelines at Ersilia. You can get it by simply doing an eosvc download, and it will be located in 'data/config'
In short:
* taxnonomy_raw: downloaded Taxonomy for each pathogen name from the PubChem website: PubChem → Search → Organism Name → “Taxonomy” → Export → Summary (CSV)
* bioassays_summary: likewise, downloaded Bioassays for each pathogen name from the PubChem website: PubChem → Search → Organism Name → “Bioassays” → Export → Summary (CSV)
* bioassays_selected_manually: manual curation of assays that would be discarded based on automatic curation but we do want to include (curation from processed/bioassays_summary/bioassays_{pathogen_name}_manual_check.csv)
* chembl_mappings: assays in chembl (assays.csv) and individual assay lists per each pathogen of interest, with the total number of compounds per assay. Comes from the [chembl-antimicrobial-tasks](https://github.com/ersilia-os/chembl-antimicrobial-tasks)

### 2. Run scripts sequentially
The scripts are numbered and when run sequentially, will generate the necessary files:
* '00_preprocess_bioassays.py': curates the manually downloaded lists of taxonomy names associated with a pathogen and keeps the right ones. Uses the taxonomies to curate the bioassays for each pathogen associated to the right taxid or targettaxid (if both fields empty, add to a manual_check list).
* '01_downoad_pubchem_bioassays.py': automatically downloads the Bioassays (Data and Description) and the entire bioassays.csv file in the data/raw/bioassays folder.
* '02_extract_bioassays.py': compares pubchem and chembl and keeps a list of pubchem AIDs to consider for individual dataset extraction.


## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South.  
Ersilia’s main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia).

![Ersilia Logo](assets/Ersilia_Brand.png)
