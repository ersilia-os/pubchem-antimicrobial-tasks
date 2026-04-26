# Scripts

Pipeline for extracting PubChem antimicrobial bioassays and filtering out data already present in ChEMBL, yielding ML-ready binary classification tasks (active vs inactive) for 15 pathogens.

## Pathogens

Acinetobacter baumannii, Campylobacter, Candida albicans, Enterobacter, Enterococcus faecium, Escherichia coli, Helicobacter pylori, Klebsiella pneumoniae, Mycobacterium tuberculosis, Neisseria gonorrhoeae, Plasmodium falciparum, Pseudomonas aeruginosa, Schistosoma mansoni, Staphylococcus aureus, Streptococcus pneumoniae.

Pathogen list and shared constants (`MIN_COMPOUNDS`, `CHEMBL_MISMATCH_THRESHOLD`) are defined in `src/default.py`.

## Pipeline Overview

```
Manual setup
  ├── Download taxonomy & bioassay CSVs from PubChem → data/config/taxonomy_raw/, data/config/bioassays_summary/
  └── Obtain ChEMBL assay mappings → data/config/chembl_mappings/

00_preprocess_bioassays.py
  └── Filter bioassays by pathogen taxonomy, keep assays with ≥MIN_COMPOUNDS compounds
       → data/processed/bioassays_summary/bioassays_{pathogen}.csv

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
  └── Copy final selected assay files to output directory
       → output/results/{pathogen}/{aid}.csv

06_plot_results.py
  └── Per-pathogen report figures (compound counts, % selected AIDs) and cross-pathogen summary table
       → output/plots/06_{pathogen}_report.png
       → output/06_summary_table.csv

07_summary_activity.py
  └── Stacked bar chart of activity outcome counts (active, inactive, inconclusive, unspecified, chemical probe) per pathogen
       → output/plots/07_summary_activity.png
```

`manual_download_individual_aid.py` is a utility for downloading a single assay via the PubChem REST API; it is not part of the main batch pipeline.

## Scripts

### 00_preprocess_bioassays.py

Curates taxonomy data and filters raw bioassay summaries to retain assays relevant to each pathogen.

**Inputs**
- `data/config/taxonomy_raw/PubChem_taxonomy_{pathogen}.csv` — raw taxonomy files (manual download)
- `data/config/bioassays_summary/PubChem_bioassay_{pathogen}.csv` — raw bioassay summaries (manual download)
- `data/config/bioassays_selected_manually/{pathogen}.csv` — optional manually curated inclusions

**Outputs**
- `data/processed/00_taxonomy_processed/taxonomy_{pathogen}.csv`
- `data/processed/00_bioassays_summary/bioassays_{pathogen}.csv` — filtered assay lists
- `data/processed/00_bioassays_summary/summary.csv`

**Logic**
1. Remove phage/virus entries from taxonomy; retain bacteria, fungi, and parasites.
2. Match pathogen species names against taxonomy IDs (`targettaxid`, `taxids`).
3. Retain assays with ≥`MIN_COMPOUNDS` compounds; flag assays with no taxonomy info for manual review.
4. Optionally merge manually selected assays.

---

### 01_download_bioassays_tsv.py

Downloads `bioassays.tsv.gz` from the NCBI FTP server and unpacks it. Skips the download if the file already exists.

**Usage**
```bash
python scripts/01_download_bioassays_tsv.py --out data/raw/bioassays
```

**Output**
- `data/raw/01_bioassays/bioassays.tsv` — assay metadata for all PubChem bioassays (AID, Source ID, compound counts)

---

### 02_bioassays_not_in_chembl.py

Compares PubChem assays against ChEMBL to identify new or significantly larger assays, avoiding duplication with the ChEMBL pipeline. `bioassays.tsv` is loaded once before the pathogen loop.

**Inputs**
- `data/processed/00_bioassays_summary/bioassays_{pathogen}.csv` — from script 00
- `data/raw/01_bioassays/bioassays.tsv` — from script 01
- `data/config/chembl_mappings/assays.csv` — ChEMBL–PubChem AID mappings
- `data/config/chembl_mappings/compounds_per_assay_{pathogen}.csv` — ChEMBL compound counts per assay

**Outputs**
- `data/processed/02_bioassays_to_keep/aids_{pathogen}.csv` — final AID selection
- `data/processed/02_bioassays_to_keep/summary.csv` — per-pathogen statistics
- `data/processed/02_bioassays_to_keep/chembl_assays_in_pubchem_{code}.csv` — named using ChEMBL code (e.g. `mtuberculosis`)
- `data/processed/02_bioassays_to_keep/chembl_cpds_mismatch_{pathogen}.csv`

**Selection criteria** — an assay is kept if it meets **either** of:
- Not present in ChEMBL and has ≥`MIN_COMPOUNDS` compounds.
- Present in ChEMBL but PubChem contains >`CHEMBL_MISMATCH_THRESHOLD` (10%) more compounds AND PubChem count ≥`MIN_COMPOUNDS`.

---

### 03_download_data_zips.py

Downloads only the PubChem Data ZIP blocks that contain the AIDs selected by script 02. Each block covers 1000 AIDs (e.g. `0743001_0744000.zip` contains AIDs 743001–744000). Already-complete files are skipped; partial files are resumed.

**Usage**
```bash
python scripts/03_download_data_zips.py \
    --aids-dir data/processed/02_bioassays_to_keep \
    --out data/raw/03_data_zips
```

**Inputs**
- `data/processed/02_bioassays_to_keep/aids_*.csv` — from script 02

**Output**
- `data/raw/03_data_zips/` — ZIP blocks covering only the selected AIDs

**Features:** targeted block selection, resumable downloads, parallel workers, exponential backoff retries.

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
- `data/processed/04_unique_cids/unique_cids_{pathogen}.csv`
- `data/processed/04_unique_cids/unique_cids_all_pathogens.csv`
- `data/processed/04_unique_cids/summary.csv`

**Cleaned CSV columns:** `sid`, `cid`, `smiles`, `inchikey`, `activity`

Activity encoding: `1` active, `0` inactive, `-1` inconclusive, `2` unspecified, `3` chemical probe.

SMILES are converted to InChIKeys via RDKit for normalisation.

---

### 05_curate_bioassays.py

Final step: copies the selected assay files to the output directory, ready for ML task creation.

**Inputs**
- `data/processed/02_bioassays_to_keep/aids_{pathogen}.csv` — from script 02
- `data/processed/04_extracted_bioassays/{pathogen}/` — from script 04

**Outputs**
- `output/05_results/{pathogen}/{aid}.csv`
- `output/05_results/{pathogen}/{aid}_meta.csv`

---

### 06_plot_results.py

Generates per-pathogen report figures and a cross-pathogen summary table comparing PubChem and ChEMBL datasets. ChEMBL InChIKeys are computed from SMILES on the first run and cached to avoid recomputation.

**Inputs**
- `data/processed/02_bioassays_to_keep/summary.csv` — from script 02
- `data/processed/04_unique_cids/unique_cids_{pathogen}.csv` — from script 04
- `../../chembl-antimicrobial-tasks/output/{code}/20_all_smiles.csv` — ChEMBL SMILES (companion repo)

**Outputs**
- `output/plots/06_{pathogen}_report.png` — per-pathogen figure (Panel A: compound counts; Panel B: % selected AIDs)
- `output/06_summary_table.csv` — cross-pathogen summary

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

### 07_summary_activity.py

Aggregates activity outcome counts from the per-AID summary CSVs produced by script 04 and plots a stacked bar chart showing the breakdown per pathogen.

**Inputs**
- `data/processed/04_extracted_bioassays/{pathogen}/summary.csv` — per-AID activity counts from script 04

**Outputs**
- `output/plots/07_summary_activity.png` — stacked bar chart with one bar per pathogen, broken down by activity outcome

**Activity categories (stacked order)**

| Category | Encoded value | Colour |
|---|---|---|
| Active | `1` | mint |
| Inactive | `0` | purple |
| Inconclusive | `-1` | orange |
| Unspecified | `2` | gray |
| Chemical probe | `3` | yellow |

---

### manual_download_individual_aid.py

Utility for downloading and processing a single PubChem assay via the REST API (max 10 000 SIDs per request). Not part of the main batch pipeline; intended for manual curation of specific assays.

**Hardcoded defaults:** AID 743156, pathogen `campylobacter`.

**Outputs**
- `output/results/{pathogen}/{aid}.csv`
- `output/results/{pathogen}/{aid}_meta.csv`

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

## Dependencies

```
pandas
rdkit
requests
beautifulsoup4
tqdm
matplotlib
stylia
```

External data: PubChem FTP/REST API, ChEMBL assay mappings (from the companion `chembl-antimicrobial-tasks` repository).
