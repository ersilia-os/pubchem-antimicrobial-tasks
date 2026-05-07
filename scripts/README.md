# Scripts

Pipeline for extracting PubChem antimicrobial bioassays and filtering out data already present in ChEMBL, yielding ML-ready binary classification tasks (active vs inactive) for 15 pathogens.

## Pathogens

Acinetobacter baumannii, Campylobacter, Candida albicans, Enterobacter, Enterococcus faecium, Escherichia coli, Helicobacter pylori, Klebsiella pneumoniae, Mycobacterium tuberculosis, Neisseria gonorrhoeae, Plasmodium falciparum, Pseudomonas aeruginosa, Schistosoma mansoni, Staphylococcus aureus, Streptococcus pneumoniae.

Pathogen list and shared constants (`MIN_COMPOUNDS`, `CHEMBL_MISMATCH_THRESHOLD`) are defined in `src/default.py`.

## Pipeline Overview

```
Manual setup
  ├── Download taxonomy & bioassay CSVs from PubChem → data/config/taxonomy_raw/, data/config/bioassays_summary/
  ├── Obtain ChEMBL assay mappings → data/config/chembl_mappings/
  └── Download ChEMBL target dictionary → data/config/chembl_mappings/target_dictionary.csv

00_preprocess_bioassays.py
  └── Filter bioassays by pathogen taxonomy, keep assays with ≥MIN_COMPOUNDS compounds
       → data/processed/00_bioassays_summary/bioassays_{pathogen}.csv
       → output/00_preprocess_bioassays/00_total_aids.png

01_download_bioassays_tsv.py
  └── Download bioassays.tsv metadata file from NCBI FTP
       → data/raw/01_bioassays/bioassays.tsv

02_bioassays_not_in_chembl.py
  └── Compare PubChem assays against ChEMBL; keep assays absent from or underrepresented in ChEMBL
       → data/processed/02_bioassays_to_keep/aids_{pathogen}.csv

03_download_data_zips.py
  └── Download only the Data ZIP blocks that contain the selected AIDs
       → data/raw/03_data_zips/

04_extract_pubchem_bioassays.py
  └── Extract and clean the selected assay CSVs from ZIP blocks; convert to InChIKey
       → data/processed/04_extracted_bioassays/{pathogen}/{aid}.csv
       → data/processed/04_unique_cids/

05_curate_bioassays.py
  └── Copy final selected assay files to output directory; generate per-pathogen report figures and cross-pathogen summary table
       → output/05_curate_bioassays/{pathogen}/{aid}.csv
       → output/05_curate_bioassays/{pathogen}/{aid}_meta.csv
       → output/05_curate_bioassays/{pathogen}/05_{pathogen}_report.png
       → output/05_curate_bioassays/05_summary_table.csv

06_annotate_selected_bioassays.py
  └── Annotate each selected AID with target type, ChEMBL ID, and compound counts; generate summary plots and per-assay task files
       → output/06_annotate_selected_bioassays/06_summary.csv
       → output/06_annotate_selected_bioassays/06_activity_outcomes.png
       → output/06_annotate_selected_bioassays/06_target_types.png
       → output/06_selected_bioassays/{code}/{aid}.csv

08_annotate_bioassays.py
  └── Annotate each AID with target type (ChEMBL + PubChem) and activity type (PubChem readout columns)
       → data/processed/08_annotated_assays/summaries.csv
       → data/processed/08_annotated_assays/summaries_organism.csv

09_filter_bioassays_to_model.py
  └── Apply quality thresholds, label assays A/B, flag non-antimicrobial counter-screens
       → data/processed/09_bioassays_to_model/bioassays_to_model.csv
```

`manual_download_individual_aid.py` is a utility for downloading a single assay via the PubChem REST API; it is not part of the main batch pipeline.

## Scripts

### 00_preprocess_bioassays.py

Curates taxonomy data and filters raw bioassay summaries to retain assays relevant to each pathogen.

**Inputs**
- `data/config/taxonomy_raw/PubChem_taxonomy_{pathogen}.csv` — raw taxonomy files (manual download from PubChem interface)
- `data/config/bioassays_summary/PubChem_bioassay_{pathogen}.csv` — raw bioassay summaries (manual download from PubChem interface)
- `data/config/bioassays_selected_manually/{pathogen}.csv` — optional manually curated inclusions (datasets that are relevant but not selected from the interface)

**Logic**
1. Remove phage/virus entries from taxonomy; retain bacteria, fungi, and parasites.
2. Match pathogen species names against taxonomy IDs (`targettaxid`, `taxids`).
3. Retain assays with ≥`MIN_COMPOUNDS` compounds; flag assays with no taxonomy info for manual review.
4. Optionally merge manually selected assays.

**Outputs**
- `data/processed/00_taxonomy_processed/taxonomy_{pathogen}.csv`
- `data/processed/00_bioassays_summary/bioassays_{pathogen}.csv`
- `data/processed/00_bioassays_summary/summary.csv`
- `output/00_preprocess_bioassays/00_total_aids.png` — bar charts of AID counts and compound counts, with and without the ≥`MIN_COMPOUNDS` filter

---

### 01_download_bioassays_tsv.py

Downloads `bioassays.tsv.gz` from the NCBI FTP server and unpacks it, delivering the assay metadata for all PubChem bioassays (AID, Source ID, compound counts). Skips the download if the file already exists.

---

### 02_bioassays_not_in_chembl.py

Compares PubChem assays against ChEMBL to identify new or significantly larger assays, avoiding duplication with the ChEMBL pipeline.

**Inputs**
- `data/processed/00_bioassays_summary/bioassays_{pathogen}.csv` — from script 00
- `data/raw/01_bioassays/bioassays.tsv` — from script 01
- `data/config/chembl_mappings/assays.csv` — ChEMBL–PubChem AID mappings, obtained from the [chembl-antimicrobial-tasks](https://github.com/ersilia-os/chembl-antimicrobial-tasks) repository
- `data/config/chembl_mappings/compounds_per_assay_{pathogen}.csv` — ChEMBL compound counts per assay, obtained from the [chembl-antimicrobial-tasks](https://github.com/ersilia-os/chembl-antimicrobial-tasks) repository

**Outputs**
- `data/processed/02_bioassays_to_keep/aids_{pathogen}.csv` — final AID selection
- `data/processed/02_bioassays_to_keep/aids_not_in_chembl_{pathogen}.csv` — AIDs with no ChEMBL counterpart
- `data/processed/02_bioassays_to_keep/summary.csv` — per-pathogen statistics
- `data/processed/02_bioassays_to_keep/chembl_assays_in_pubchem_{pathogen_code}.csv` — ChEMBL assays with PubChem counterpart
- `data/processed/02_bioassays_to_keep/chembl_cpds_mismatch_{pathogen}.csv` — Assays with mismatched data between ChEMBL and PubChem
- `output/02_bioassays_not_in_chembl/02_aids_to_keep.png` — bar chart of selected AIDs per pathogen
- `output/02_bioassays_not_in_chembl/02_chembl_pubchem_overview.png` — scatter/bar overview comparing PubChem vs ChEMBL counts

**Selection criteria** — an assay is kept if it meets **either** of:
- Not present in ChEMBL and has ≥`MIN_COMPOUNDS` compounds.
- Present in ChEMBL but PubChem contains >`CHEMBL_MISMATCH_THRESHOLD` (10%) more compounds AND PubChem count ≥`MIN_COMPOUNDS`.

---

### 03_download_data_zips.py

Downloads only the PubChem Data ZIP blocks that contain the AIDs selected by script 02. Each block covers 1000 AIDs (e.g. `0743001_0744000.zip` contains AIDs 743001–744000). Designed to save download time and disk space

---

### 04_extract_pubchem_bioassays.py

Extracts and cleans the selected assay CSVs from the downloaded ZIP blocks, then aggregates compound statistics across pathogens.

**Inputs**
- `data/raw/03_data_zips/` — ZIP blocks from script 03
- `data/processed/02_bioassays_to_keep/aids_{pathogen}.csv` — final AID list from script 02

**Outputs**
- `data/raw/04_unzipped/{pathogen}/` — raw extracted CSVs
- `data/processed/04_extracted_bioassays/{pathogen}/{aid}.csv` — cleaned assay data
- `data/processed/04_extracted_bioassays/{pathogen}/{aid}_meta.csv` — assay metadata
- `data/processed/04_extracted_bioassays/{pathogen}/summary.csv` — per-AID compound counts
- `data/processed/04_unique_cids/unique_cids_{pathogen}.csv` — list of unique CIDs (compound identifier) selected per pathogen
- `data/processed/04_unique_cids/unique_cids_all_pathogens.csv` — list of unique CIDs (compound identifier) selected across all pathogens
- `data/processed/04_unique_cids/summary.csv` — summary of number of compounds and assays selected per pathogen

**Cleaned CSV columns:** `sid`, `cid`, `smiles`, `inchikey`, `activity`

Activity encoding: `1` active, `0` inactive, `-1` inconclusive, `2` unspecified, `3` chemical probe.

---

### 05_curate_bioassays.py

Copies the selected assay files to the output directory, then generates per-pathogen report figures and a cross-pathogen summary table comparing PubChem and ChEMBL datasets. Needs to have the [chembl-antimicrobial-tasks](https://github.com/ersilia-os/chembl-antimicrobial-tasks) cloned in the same root as this repository.

**Inputs**
- `data/processed/02_bioassays_to_keep/aids_{pathogen}.csv` — from script 02
- `data/processed/02_bioassays_to_keep/summary.csv` — from script 02
- `data/processed/04_extracted_bioassays/{pathogen}/` — from script 04
- `data/processed/04_unique_cids/unique_cids_{pathogen}.csv` — from script 04
- `../../chembl-antimicrobial-tasks/output/{code}/20_all_smiles.csv` — ChEMBL SMILES (companion repo)

**Outputs**
- `output/05_curate_bioassays/{pathogen}/{aid}.csv`
- `output/05_curate_bioassays/{pathogen}/{aid}_meta.csv`
- `output/05_curate_bioassays/{pathogen}/05_{pathogen}_report.png` — per-pathogen figure (Panel A: compound counts; Panel B: % selected AIDs)
- `output/05_curate_bioassays/05_summary_table.csv` — cross-pathogen summary

**Figure panels**
- Panel A: bar chart of unique compound counts for ChEMBL, all PubChem (non-unique across assays), and novel PubChem (selected assays only).
- Panel B: percentage of total pathogen AIDs that are not in ChEMBL vs. mismatched (present in ChEMBL but with >10% more compounds in PubChem).

**Summary table columns**

| Column | Description |
|--------|-------------|
| `pathogen` | PubChem pathogen name |
| `chembl_code` | ChEMBL output directory code for this pathogen |
| `n_chembl_compounds` | Number of unique compounds in ChEMBL for this pathogen (by InChIKey) |
| `n_pubchem_all` | Total compound count across all PubChem bioassays for this pathogen (non-unique across assays) |
| `n_pubchem_novel` | Number of unique compounds in the selected PubChem assays (by InChIKey) |
| `pct_pubchem_novel` | `n_pubchem_novel` as a percentage of `n_pubchem_all` |
| `n_bioassay_aids` | Total number of PubChem AIDs associated with this pathogen after script 00 filtering |
| `n_not_in_chembl` | Total number of AIDs with no ChEMBL link, regardless of compound count |
| `n_not_in_chembl_min` | AIDs with no ChEMBL link that also meet the ≥`MIN_COMPOUNDS` threshold (i.e. selected) |
| `n_mismatched` | Number of selected AIDs present in ChEMBL but with >10% more compounds in PubChem |
| `pct_not_in_chembl` | `n_not_in_chembl` as a percentage of `n_bioassay_aids` |
| `pct_not_in_chembl_min` | `n_not_in_chembl_min` as a percentage of `n_bioassay_aids` |
| `pct_mismatched` | `n_mismatched` as a percentage of `n_bioassay_aids` |

---

### 06_annotate_selected_bioassays.py

Annotates every selected AID with target type, protein ID, ChEMBL ID, compound counts, and a quality label.

**Inputs**
- `data/processed/02_bioassays_to_keep/aids_{pathogen}.csv` — from script 02
- `data/processed/02_bioassays_to_keep/chembl_assays_in_pubchem_{code}.csv` — from script 02
- `data/processed/04_extracted_bioassays/{pathogen}/summary.csv` — from script 04
- `data/config/bioassays_summary/PubChem_bioassay_{pathogen}.csv` — manual download
- `data/config/bioassays_summary/bioassays_manually_excluded.csv` — manual curation; AIDs with `discard=1` are excluded before annotation
- `../../chembl-antimicrobial-tasks/output/{code}/18_assays_master.csv` — companion repo

**Outputs**
- `output/06_annotate_selected_bioassays/06_summary.csv`
- `output/06_selected_bioassays/{code}/{aid}.csv` — active + inactive compounds only (`smiles`, `bin`); one file per non-discarded assay

**Output columns**

| Column | Description |
|---|---|
| `pathogen_code` | Short pathogen code (e.g. `abaumannii`, `mtuberculosis`) |
| `aid` | PubChem Assay ID |
| `chembl_id` | Linked ChEMBL assay ID, if the assay had a compound mismatch with ChEMBL |
| `target_type` | `single_protein`, `organism`, `discarded`, or `other` |
| `protein_id` | UniProt accession(s) from `18_assays_master.csv` (ChEMBL assays) or `protacxns` from PubChem (novel assays) |
| `cids` | Number of unique compounds |
| `actives` | Active compound count |
| `inactives` | Inactive compound count |
| `inconclusive` | Inconclusive compound count |
| `unspecified` | Unspecified compound count |
| `chemical_probe` | Chemical probe compound count |
| `label` | `A` ; `B` ; `discarded` if denominator is zero |

Label mimicks the ChEMBL decision tree:
| Label | Condition |
|---|---|
| `A` | Large class-balanced dataset: `cids ≥ 1000`, `actives ≥ 50`, `ratio < 0.5` |
| `B` | Active-enriched dataset: `actives ≥ 100`, `ratio ≥ 0.5` |
| discarded | Neither condition met |

**target_type logic**
- **ChEMBL assays** (present in `chembl_assays_in_pubchem_{code}.csv`): read `target_type_curated_extra` from `18_assays_master.csv` → `SINGLE PROTEIN` → `single_protein`; `ORGANISM` → `organism`; `DISCARDED` → `discarded`; anything else → `other`
- **Novel PubChem assays**: `protacxns` or `geneids` populated → `single_protein`; both empty → `organism`

---

### manual_download_individual_aid.py

Utility for downloading and processing a single PubChem assay via the REST API (max 10 000 SIDs per request). Not part of the main batch pipeline; intended for manual curation of specific assays.

---

## Output Format

Each final assay CSV contains one compound per row:

| Column | Description |
|--------|-------------|
| `sid` | PubChem Substance ID |
| `cid` | PubChem Compound ID |
| `smiles` | Canonical SMILES |
| `inchikey` | Standard InChIKey |
| `activity` | Binary label: `1` active, `0` inactive |

---

External data: PubChem FTP/REST API, ChEMBL assay mappings (from the companion `chembl-antimicrobial-tasks` repository).
