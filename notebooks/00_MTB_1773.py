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
# Hello

# %% [markdown]
# # Mycobacterium tuberculosis (TaxID 1773) – BioAssay Comparison
#
# Objective:
# Compare BioAssays obtained from:
# 1. **PubChem Web UI export** for TaxID 1773
# 2. **Aid2Taxid.tsv explicit taxonomy mapping**
#
# We identify assays that appear in the Web UI download but *not* in Aid2Taxid,
# and inspect metadata to understand the differences.

# %% [markdown]
# ## 0. Setup

# %%
import pandas as pd
from pathlib import Path

NOTEBOOK_DIR = Path().resolve()
DATA_RAW = NOTEBOOK_DIR.parent / "data" / "raw"
DATA_PROCESSED = NOTEBOOK_DIR.parent / "data" / "processed"

NOTEBOOK_DIR, DATA_RAW

# %% [markdown]
# ## 1. Load UI-exported BioAssays for MTB
#
# File downloaded manually from:
#
# **PubChem → Taxonomy → Mycobacterium tuberculosis (TaxID 1773) → Actions → BioAssays → Summary**

# %%
df_ui = pd.read_csv(DATA_RAW / "pubchem_taxid_1773_bioassay.csv")

df_ui.head()

# %% [markdown]
# ## 2. Load Aid2Taxid.tsv mapping file

# %%
df_map = pd.read_csv(DATA_RAW / "Aid2Taxid.tsv", sep="\t")

df_map.head()

# %% [markdown]
# ## 3. Define the MTB TaxID
#
# We only examine **TaxID = 1773**.

# %%
mtb_taxid = "1773"

# %% [markdown]
# ## 4. Extract AIDs explicitly linked to MTB in Aid2Taxid

# %%
df_map_mtb = df_map[df_map["TaxID"].astype(str) == mtb_taxid]

aid2tax_aids = set(df_map_mtb["AID"].astype(int))
len(aid2tax_aids)

# %% [markdown]
# ## 5. Extract UI-linked AIDs

# %%
ui_aids = set(df_ui["aid"].astype(int))
len(ui_aids)

# %% [markdown]
# ## 6. Compare the two sets

# %%
missing_from_aid2tax = ui_aids - aid2tax_aids
extra_in_aid2tax = aid2tax_aids - ui_aids

summary = pd.DataFrame({
    "Metric": ["UI AIDs", "Aid2Taxid AIDs", "Missing (UI→Aid2Taxid)", "Unexpected (Aid2Taxid→UI)"],
    "Count": [len(ui_aids), len(aid2tax_aids), len(missing_from_aid2tax), len(extra_in_aid2tax)]
})

summary

# %% [markdown]
# ## 7. Inspect Missing AIDs (in UI but *not* in Aid2Taxid)

# %%
df_missing = df_ui[df_ui["aid"].isin(missing_from_aid2tax)].copy()

df_missing.head()

# %% [markdown]
# ### 7.1 Where do the missing assays come from?

# %%
df_missing["aidsrcname"].value_counts()

# %%
df_missing["aidtype"].value_counts()

# %% [markdown]
# ### 7.2 Evaluate TaxIDs inside the `taxids` column

# %%
def extract_taxids(taxid_field):
    """Return list of TaxIDs (strings) from the pipe-separated taxid field."""
    if pd.isna(taxid_field):
        return []
    return str(taxid_field).split("|")


# %%
df_missing["parsed_taxids"] = df_missing["taxids"].apply(extract_taxids)

df_missing["has_mtb"] = df_missing["parsed_taxids"].apply(
    lambda lst: mtb_taxid in lst
)

df_missing["has_mtb"].value_counts()

# %% [markdown]
# ### 7.3 Inspect AIDs that DO contain the MTB TaxID

# %%
df_missing_with_mtb = df_missing[df_missing["has_mtb"]]

df_missing_with_mtb[["aid", "aidsrcname", "aidname", "taxids", "parsed_taxids"]].head(20)

# %% [markdown]
# ### 7.4 Inspect AIDs that do NOT contain the MTB TaxID

# %%
df_missing_no_mtb = df_missing[~df_missing["has_mtb"]]

df_missing_no_mtb[["aid", "aidsrcname", "aidname", "parsed_taxids"]].head(20)

# %% [markdown]
# ### 7.5 Are missing assays mapped to ANY of our 15 pathogens?

# %%
# Load pathogen taxonomy table
tax_table = pd.read_csv(DATA_PROCESSED / "taxonomy_table.csv")

pathogen_tax_map = (
    tax_table.groupby("Pathogen")["Taxonomy_ID"]
    .apply(lambda x: set(map(str, x)))
    .to_dict()
)

def match_pathogens(tid_list):
    matches = []
    for pathogen, ids in pathogen_tax_map.items():
        if any(t in tid_list for t in ids):
            matches.append(pathogen)
    return matches

df_missing_no_mtb["other_pathogen_hits"] = df_missing_no_mtb["parsed_taxids"].apply(match_pathogens)

df_missing_no_mtb[["aid", "parsed_taxids", "other_pathogen_hits"]].head(20)

# %% [markdown]
# ### 7.6 Missing values in TaxIDs

# %%
df_missing["taxids"].isna().sum()

# %% [markdown]
# ## 8. Summary of missing categories

# %%
summary_missing = pd.DataFrame({
    "Category": [
        "Missing AIDs (UI→Aid2Taxid)",
        "Contain MTB TaxID",
        "Match Other Pathogens",
        "Match No Known Pathogens",
        "No TaxID Present"
    ],
    "Count": [
        len(df_missing),
        df_missing["has_mtb"].sum(),
        sum(df_missing_no_mtb["other_pathogen_hits"].apply(len) > 0),
        sum(df_missing_no_mtb["other_pathogen_hits"].apply(len) == 0),
        sum(df_missing["taxids"].isna())
    ]
})

summary_missing

# %% [markdown]
# ## 9 Comparing with downloaded PubChem taxid and taxonomy search

# %%
# Load your XML-derived filtered assays
df_xml = pd.read_csv(DATA_PROCESSED / "filtered_assays_description_results.csv")

df_xml.head()


# %%
# Parse the TaxIDs inside the XML-derived file

def parse_list(s):
    if pd.isna(s): return []
    # stored as a Python-like string: "['1773', '1234']"
    s = s.strip("[]")
    items = [x.strip(" '\"") for x in s.split(",") if x.strip()]
    return items

df_xml["Parsed_TaxIDs"] = df_xml["TaxIDs_detected"].apply(parse_list)

df_xml.head()

# %%
# Extract AIDs explicitly linked to MTB

xml_aids = set(
    df_xml[df_xml["Parsed_TaxIDs"].apply(lambda lst: mtb_taxid in lst)]["AID"]
)

len(xml_aids)

# We already have the UI-linked AIDs as 'ui_aids'

# %%
# Compare the 2 sets

missing_from_xml = ui_aids - xml_aids
extra_in_xml = xml_aids - ui_aids

summary = pd.DataFrame({
    "Metric": ["UI AIDs", "XML AIDs", "Missing (UI→XML)", "Unexpected (AXML→UI)"],
    "Count": [len(ui_aids), len(xml_aids), len(missing_from_xml), len(extra_in_xml)]
})

summary

# %%
# Inspect Missing AIDs (in UI but not in xml)

df_missing_xml = df_ui[df_ui["aid"].isin(missing_from_xml)].copy()

df_missing_xml.head()

# %%
df_missing_xml["aidsrcname"].value_counts()

# %%
df_missing_xml["aidtype"].value_counts()

# %%
# Evaluate TaxIDs inside the `taxids` column
df_missing_xml["parsed_taxids"] = df_missing_xml["taxids"].apply(extract_taxids)

df_missing_xml["has_mtb"] = df_missing_xml["parsed_taxids"].apply(
    lambda lst: mtb_taxid in lst
)

df_missing_xml["has_mtb"].value_counts()

# %% [markdown]
# Compared to Aid2Taxid, we are getting ALL assays with the taxid correct! :)

# %%
# Inspect AIDs that DO contain the MTB TaxID
df_missing_xml_with_mtb = df_missing_xml[df_missing_xml["has_mtb"]]

df_missing_xml_with_mtb[["aid", "aidsrcname", "aidname", "taxids", "parsed_taxids"]].head(20)

# %%
