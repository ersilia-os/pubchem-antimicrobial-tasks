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
# # PubChem BioAssays for *Mycobacterium tuberculosis*
#
# This notebook processes PubChem BioAssay data associated with *Mycobacterium tuberculosis* and its related taxa.
#
# 1. Loads the exported TaxID–AID table  
# 2. Cleans the dataset and removes taxa with no BioAssays  
# 3. Expands the pipe-separated AID lists into individual assay records  
# 4. Produces a tidy TaxID–AID table  
# 5. Summarizes:
#    - number of unique taxonomy entries  
#    - number of unique BioAssays  
#    - number of assays per taxon  
#
# This serves as the **first step** for building a complete organism-level
# BioAssay retrieval pipeline for antimicrobial pathogens.

# %% [markdown]
# ## 0. Setup

# %%
import pandas as pd
from pathlib import Path

NOTEBOOK_DIR = Path().resolve()
DATA_RAW = NOTEBOOK_DIR.parent / "data" / "raw"
DATA_PROCESSED = NOTEBOOK_DIR.parent / "data" / "processed"
DATA_PROCESSED.mkdir(exist_ok=True)

# %% [markdown]
# ## 1. Load the TaxID–AID Mapping
#
# PubChem currently does not provide a stable programmatic endpoint to query
# taxonomy-linked BioAssays directly from an organism name. Therefore, the initial dataset used here was downloaded manually from the
# PubChem interface:
#
# Search → “Mycobacterium tuberculosis” → Taxonomy →  Actions → BioAssays → Download: Summary (Search Results).
#
# The downloaded file contains:
# - PubChem Taxonomy IDs (TaxID)
# - Scientific names
# - Pipe-separated BioAssay IDs (AID) linked to each taxon

# %%
df_raw = pd.read_csv(DATA_DIR / "PubChem_taxonomy_text_Mycobacterium tuberculosis.csv")
df_raw.head()

# %% [markdown]
# ## 2. Select Relevant Columns & Remove Empty Entries

# %%
df = (
    df_raw[["Taxonomy_ID", "Taxonomy_Name", "Linked_BioAssays"]]
    .dropna(subset=["Linked_BioAssays"])
)

df.head()

# %% [markdown]
# ## 3. Expand BioAssays into Individual AIDs

# %%
df["AID"] = df["Linked_BioAssays"].str.split("|")
df_expanded = df.explode("AID").copy()
df_expanded["AID"] = df_expanded["AID"].astype(int)

df_final = df_expanded[["Taxonomy_ID", "Taxonomy_Name", "AID"]].drop_duplicates()
df_final.head()

# %% [markdown]
# ## 4. Count Assays per Taxonomy ID

# %%
assays_per_tax = (
    df_final.groupby(["Taxonomy_ID", "Taxonomy_Name"])["AID"]
    .nunique()
    .reset_index(name="N_Assays")
    .sort_values("N_Assays", ascending=False)
    .reset_index(drop=True)
)

assays_per_tax

# %% [markdown]
# ## 5. Summary

# %%
summary = pd.DataFrame({
    "Metric": ["Unique Taxonomy IDs", "Unique Assay IDs"],
    "Value": [
        df_final["Taxonomy_ID"].nunique(),
        df_final["AID"].nunique()
    ]
})

summary

# %% [markdown]
# ## 6. Comparing donwloaded with website

# %%
# 1. Curated AIDs from your exported CSV (only TaxID 1773)
df_curated_1773 = df_final[df_final["Taxonomy_ID"] == 1773].copy()
curated_aids_1773 = set(df_curated_1773["AID"].astype(int))

print("Curated AIDs (TaxID 1773):", len(curated_aids_1773))

# 2. AIDs from the UI-downloaded assay table
df_ui = pd.read_csv(DATA_RAW / "pubchem_taxid_1773_bioassay.csv")
ui_aids = set(df_ui["aid"].astype(int))

print("UI AIDs:", len(ui_aids))

# 3. Compute differences
missing_from_curated = ui_aids - curated_aids_1773
extra_in_curated = curated_aids_1773 - ui_aids

print("AIDs missing from curated:", len(missing_from_curated))
print("Curated AIDs missing from UI:", len(extra_in_curated))

# %% [markdown]
# We are missing 9769 AIDs in the original curated list.

# %%
# 4. Extract metadata for missing ones
df_missing = df_ui[df_ui["aid"].isin(missing_from_curated)]
df_missing.head()


# %% [markdown]
# All missing AIDs are ChEMBL-based BioAssays???

# %%
# 5. Extract metadata for curated ones
df_found = df_ui[df_ui["aid"].isin(curated_aids_1773)].copy()

# Show sources
print("\nSources of missing AIDs:")
print(df_missing["aidsrcname"].value_counts())

print("\nSources of curated (UI-matched) AIDs:")
print(df_found["aidsrcname"].value_counts())

# %% [markdown]
# No, so the bast majority of assays are shared with ChEMBL, this is not the issue.
