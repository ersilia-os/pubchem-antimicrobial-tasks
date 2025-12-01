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
# # Neisseria gonorrhoeae – BioAssay Comparison
#
# Objective:
# Compare BioAssays obtained from:
# 1. **PubChem Web UI export**
# 2. **Aid2Taxid.tsv explicit taxonomy mapping**
#
# We identify assays that appear in the UI but *not* in Aid2Taxid,
# and inspect their metadata to understand why.

# %% [markdown]
# ## 0. Setup

# %%
import pandas as pd
from pathlib import Path

NOTEBOOK_DIR = Path().resolve()
DATA_RAW = NOTEBOOK_DIR.parent / "data" / "raw"
DATA_PROCESSED = NOTEBOOK_DIR.parent / "data" / "processed"

DATA_RAW

# %% [markdown]
# ## 1. Load UI-exported BioAssays
#
# This file was manually downloaded from:
# PubChem → Search “Neisseria gonorrhoeae” → Taxonomy → Actions → BioAssays → Summary (Search Results)

# %%
ui_file = DATA_RAW / "PubChem_bioassay_text_Neisseria gonorrhoeae.csv"
df_ui = pd.read_csv(ui_file)

df_ui.head()

# %% [markdown]
# ## 2. Load Aid2Taxid.tsv mapping file

# %%
aid2taxid_file = DATA_RAW / "Aid2Taxid.tsv"
df_map = pd.read_csv(aid2taxid_file, sep="\t")

df_map.head()

# %% [markdown]
# ## 3. Identify all Neisseria TaxIDs
#
# Using the taxonomy table created earlier (taxonomy_table.csv)

# %%
tax_table = pd.read_csv(DATA_PROCESSED / "taxonomy_table.csv")

neisseria_taxids = tax_table.loc[
    tax_table["Pathogen"] == "Neisseria gonorrhoeae",
    "Taxonomy_ID"
].tolist()

neisseria_taxids

# %% [markdown]
# ## 4. Extract AIDs explicitly linked to Neisseria in Aid2Taxid

# %%
df_map_neis = df_map[df_map["TaxID"].isin(neisseria_taxids)]

aid2tax_aids = set(df_map_neis["AID"].astype(int))
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
df_missing = df_ui[df_ui["aid"].isin(missing_from_aid2tax)]

df_missing.head(5)

# %% [markdown]
# ### 7.1 Where do they come from?

# %%
df_missing["aidsrcname"].value_counts()

# %%
df_missing["aidtype"].value_counts()


# %% [markdown]
# ### 7.2 Do they show Neisseria TaxIDs inside the `taxids` column?

# %%
def extract_taxids(taxid_field):
    """Return a list of TaxIDs from the pipe-separated taxid field."""
    if pd.isna(taxid_field):
        return []
    return str(taxid_field).split("|")

def contains_neisseria(taxid_field):
    """Return True if any Neisseria taxid is in the field."""
    tid_list = extract_taxids(taxid_field)
    return any(t in tid_list for t in map(str, neisseria_taxids))


# %%
df_missing = df_missing.copy()
df_missing["parsed_taxids"] = df_missing["taxids"].apply(extract_taxids)
df_missing["contains_neisseria"] = df_missing["parsed_taxids"].apply(
    lambda lst: any(t in lst for t in map(str, neisseria_taxids))
)

df_missing["contains_neisseria"].value_counts()

# %% [markdown]
# ### 7.3 Inspect the Neisseria-containing AIDs

# %%
df_missing_neis = df_missing[df_missing["contains_neisseria"]]

df_missing_neis[["aid", "aidsrcname", "aidname", "taxids","parsed_taxids"]].head(20)

# %% [markdown]
# ### 7.4 Inspect Non-Neisseria AIDs from the UI

# %%
df_missing_non = df_missing[~df_missing["contains_neisseria"]]

df_missing_non[["aid", "aidsrcname", "aidname", "parsed_taxids"]].head(20)

# %% [markdown]
# ### 7.5 Check Whether These TaxIDs Belong to ANY of Your 15 Pathogens

# %%
# Load your pathogen taxonomy table
pathogen_tax_map = (
    tax_table.groupby("Pathogen")["Taxonomy_ID"]
    .apply(lambda x: set(map(str, x)))
    .to_dict()
)

def match_other_pathogens(tid_list):
    """Return pathogen names whose taxids appear in the list."""
    matches = []
    for pathogen, tids in pathogen_tax_map.items():
        if any(t in tid_list for t in tids):
            matches.append(pathogen)
    return matches

df_missing_non["other_pathogen_hits"] = df_missing_non["parsed_taxids"].apply(match_other_pathogens)

df_missing_non[["aid", "taxids", "parsed_taxids", "other_pathogen_hits"]].head(20)

# %% [markdown]
# ### 7.6. Check missing values in taxids

# %%
df_missing["taxids"].isna().sum()

# %% [markdown]
# ### 7.6 Summary Statistics

# %%
summary_missing = pd.DataFrame({
    "Category": [
        "Missing AIDs (UI→Aid2Taxid)",
        "Contain Neisseria TaxID",
        "Match Other Pathogens",
        "Match No Known Pathogens",
        "Do Not Contain Any Taxid"
    ],
    "Count": [
        len(df_missing),
        df_missing["contains_neisseria"].sum(),
        sum(df_missing_non["other_pathogen_hits"].apply(len) > 0),
        sum(df_missing_non["other_pathogen_hits"].apply(len) == 0),
        sum(df_missing["taxids"].isna())
    ]
})

summary_missing

# %% [markdown]
# ## 8. Comparing with XML (downloaded PubChem taxid and taxonomy search)

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

df_xml["parsed_taxids"] = df_xml["TaxIDs_detected"].apply(parse_list)

df_xml.head()

# %%
# Extract AIDs explicitly linked to Neisseria

neisseria_taxids_set = set(map(str, neisseria_taxids))

xml_aids = set(
    df_xml[
        df_xml["parsed_taxids"].apply(
            lambda lst: bool(neisseria_taxids_set.intersection(lst))
        )
    ]["AID"]
)

len(xml_aids)

# We already have the UI-linked AIDs as 'ui_aids'

# %%
# Compare the 2 datasets
missing_from_xml = ui_aids - xml_aids
extra_in_xml = xml_aids - ui_aids

summary = pd.DataFrame({
    "Metric": ["UI AIDs", "XML AIDs", "Missing (UI→XML)", "Unexpected (XML→UI)"],
    "Count": [len(ui_aids), len(xml_aids), len(missing_from_xml), len(extra_in_xml)]
})

summary

# %%
# Inspect Missing AIDs (in UI but not in XML)

df_missing_xml = df_ui[df_ui["aid"].isin(missing_from_xml)]

df_missing_xml.head(5)

# %%
df_missing_xml["aidsrcname"].value_counts()

# %%
df_missing_xml["aidtype"].value_counts()

# %%
# Do they show Neisseria TaxIDs inside the `taxids` column?

df_missing_xml.loc[:, "parsed_taxids"] = df_missing_xml["taxids"].apply(extract_taxids)

df_missing_xml.loc[:, "contains_neisseria"] = df_missing_xml["parsed_taxids"].apply(
    lambda lst: any(t in lst for t in map(str, neisseria_taxids))
)

df_missing_xml["contains_neisseria"].value_counts()

# %% [markdown]
# Compared to Aid2Taxid, we are getting ALL assays with the taxid correct! :)
