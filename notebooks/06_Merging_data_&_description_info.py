# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
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
# # 06. Merge assay metadata with compound counts
#
# This script merges two sources:
# - `summary_data.csv`: compound and substance counts per AID
# - `filtered_description_with_organisms_v2_REBUILT.csv`: metadata including pathogen and ChEMBL ID
#
# Output: `summary_bioassays.csv`

# %% [markdown]
# ## 00. Setup

# %%
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# Define paths
PROJECT_ROOT = Path("/Users/maria/Documents/Ersilia/PubChem/pubchem-antimicrobial-tasks")
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

# %% [markdown]
# ## 01. Load input files

# %%
summary_data = pd.read_csv(DATA_PROCESSED / "summary_data.csv")
summary_data = summary_data.sort_values('AID', ascending = True ).reset_index(drop = True)
summary_data.head(5)

# %%
# Check for an specific filtered AID
summary_data[summary_data["AID"]== 1626]

# %%
summary_data[summary_data["AID"]== 1332]

# %%
summary_description = pd.read_csv(DATA_PROCESSED / "filtered_description_with_organisms_v2_REBUILT.csv")
summary_description = summary_description.sort_values('AID', ascending = True ).reset_index(drop = True)
summary_description

# %%
# Check for an specific filtered AID
summary_description[summary_description["AID"]== 1626]
summary_description[summary_description["AID"]== 1332]

# %% [markdown]
# ## 02. Merge on AID

# %%
summary_bioassays = summary_data.merge(
    summary_description[["AID", "Pathogen", "ChEMBLid"]],
    on="AID",
    how="left"
)

summary_bioassays = summary_bioassays[["Pathogen", "AID", "compound_count", "substances_count", "ChEMBLid"]]

summary_bioassays = summary_bioassays.rename(columns={
    "AID": "PubChem_AID",
    "substances_count": "PubChem_substances",
    "compound_count": "PubChem_compounds",
    "ChEMBLid": "ChEMBL_ID"
})

summary_bioassays

# %% [markdown]
# ## 03. Save output

# %%
PUBCHEM = DATA_PROCESSED / "PubChem_bioassays"
PUBCHEM.mkdir(parents=True, exist_ok=True)

summary_bioassays.to_csv(PUBCHEM / "summary_bioassays.csv", index=False)

# %% [markdown]
# ## 04. Merge with ChEMBL dataframes

# %% [markdown]
# Some assays have more than one pathogen detected, because we have searched for taxid and assay organism.

# %%
summary_bioassays["Pathogen"].unique()

# %%
# Convert Pathogen column to string, just in case
summary_bioassays["Pathogen"] = summary_bioassays["Pathogen"].astype(str)

# Convert to list (split by comma + space)
summary_bioassays["Pathogen_list"] = summary_bioassays["Pathogen"].str.split(", ")

# Find rows with >1 pathogen
multi_pathogen_rows = summary_bioassays[summary_bioassays["Pathogen_list"].apply(len) > 1].copy()

# Create a clean string version for display or export (optional)
multi_pathogen_rows["Pathogen_str"] = multi_pathogen_rows["Pathogen_list"].apply(lambda x: ", ".join(x))

# Show count and preview
print(f"Multi-pathogen assays: {len(multi_pathogen_rows)}")
multi_pathogen_rows[["PubChem_AID", "Pathogen_str"]]

# %%
# Get list of AIDs with multiple pathogens
aids_to_exclude = set(multi_pathogen_rows["PubChem_AID"])

# Filter df to exclude those AIDs
summary_bioassays = summary_bioassays[~summary_bioassays["PubChem_AID"].isin(aids_to_exclude)]
summary_bioassays = summary_bioassays.drop(columns = ['Pathogen_list'])
summary_bioassays

# %%
# Prepare list of pathogen files (from Arnau's ChEMBL project)
pathogen_dir = DATA_PROCESSED / "ChEMBL_bioassays"
pathogen_files = list(pathogen_dir.glob("*_ChEMBL_data.csv"))

# Initialize empty list to collect all [assay_id, cpds] rows
cpds_rows = []

# Process each pathogen file
for file in pathogen_files:
    df = pd.read_csv(file, low_memory=False)
    if "assay_chembl_id" in df.columns and "compound_chembl_id" in df.columns:
        tmp = df.groupby("assay_chembl_id")["compound_chembl_id"].nunique().reset_index()
        tmp = tmp.rename(columns={"assay_chembl_id": "ChEMBL_ID", "compound_chembl_id": "ChEMBL_compounds"})
        cpds_rows.append(tmp)

# Concatenate all assay-level compound counts
df_cpds = pd.concat(cpds_rows).drop_duplicates(subset=["ChEMBL_ID"])

# Merge into summary
summary_bioassays = summary_bioassays.merge(df_cpds, on="ChEMBL_ID", how="left")

# Fill missing counts with 0 and convert to integer
summary_bioassays["ChEMBL_compounds"] = summary_bioassays["ChEMBL_compounds"].fillna(0).astype(int)

# Save updated summary
summary_bioassays.to_csv(DATA_PROCESSED / "summary_bioassays_with_chembl_compounds.csv", index=False)

summary_bioassays.head()

# %% [markdown]
# ## 06. Check number of compounds in PubChem vs ChEMBL

# %%
# Checking if compounds numbers are equal in PubChem vs ChEMBL

summary_bioassays["same_compound_count"] = (
    summary_bioassays["PubChem_compounds"] == summary_bioassays["ChEMBL_compounds"]
)

summary_bioassays["compound_diff"] = (
    summary_bioassays["PubChem_compounds"] - summary_bioassays["ChEMBL_compounds"]
)

summary_bioassays.to_csv(DATA_PROCESSED / "summary_bioassays_with_chembl_compounds.csv", index=False)

# %%
summary_bioassays

# %%
# Group by pathogen
summary_pubchem_chembl = summary_bioassays.groupby("Pathogen").agg(
    PubChem_AIDs=("PubChem_AID", "nunique"),
    ChEMBL_IDs=("ChEMBL_ID", lambda x: x.notna().sum()),
    Matching_Compounds=("same_compound_count", "sum"),
    Total_with_ChEMBL=("same_compound_count", "count")  # Only rows with ChEMBL compound counts
).reset_index()


# Compute percentage of ChEMBL IDs
summary_pubchem_chembl["Percent_ChEMBL_IDs"] = (
    100 * summary_pubchem_chembl["ChEMBL_IDs"] / summary_pubchem_chembl["PubChem_AIDs"]
).round(1)

# Compute percentage of matched compounds (among only those with ChEMBL IDs)
summary_pubchem_chembl["Percent_Matching_Compounds"] = (
    100 * summary_pubchem_chembl["Matching_Compounds"] / summary_pubchem_chembl["Total_with_ChEMBL"]
).round(1)

# Drop intermediate column
summary_pubchem_chembl.drop(columns="Total_with_ChEMBL", inplace=True)

# Reorder columns
summary = summary_pubchem_chembl[
    ["Pathogen", "PubChem_AIDs", "ChEMBL_IDs", "Percent_ChEMBL_IDs", 
     "Matching_Compounds", "Percent_Matching_Compounds"]
]

summary_pubchem_chembl

# %%
# Data
summary = summary_pubchem_chembl.copy()
labels = summary["Pathogen"].values
percent_chembl = summary["Percent_ChEMBL_IDs"].values
percent_matching = summary["Percent_Matching_Compounds"].values

N = len(labels)
x = np.arange(N)

# ⏩ Adjust spacing: smaller bar width + larger offset
bar_width = 0.25
bar_spacing = 0.15  # controls the separation inside each group

# Colors
color_chembl = "#AA96FA"    # Purple
color_match  = "#BEE6B4"    # Green

# Plot
plt.figure(figsize=(16, 6))

# 1. Bar: Percent_ChEMBL_IDs
bars1 = plt.bar(x - bar_spacing, percent_chembl, width=bar_width,
                color=color_chembl, ec="black", zorder=2, label="% ChEMBL ID coverage")

# 2. Bar: Percent_Matching_Compounds
bars2 = plt.bar(x + bar_spacing, percent_matching, width=bar_width,
                color=color_match, ec="black", zorder=2, label="% Matching compound counts")

# Add labels on top of bars
for i in range(N):
    plt.text(x[i] - bar_spacing, percent_chembl[i] + 2, f"{percent_chembl[i]:.1f}%",
             ha='center', va='bottom', fontsize=9, zorder=3)
    plt.text(x[i] + bar_spacing, percent_matching[i] + 2, f"{percent_matching[i]:.1f}%",
             ha='center', va='bottom', fontsize=9, zorder=3)

# Aesthetics
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Percentage (%)")
plt.ylim(0, 110)
plt.title("Percent of ChEMBL-Mapped Assays and Matching Compounds by Pathogen")
plt.grid(linestyle="--", zorder=1, axis="y")

# Legend at bottom
plt.legend(
    loc="lower center",
    bbox_to_anchor=(0.18, 0.05),
    ncol=2,
    framealpha=1,
    edgecolor="k",
    prop={"size": 10}
)

plt.tight_layout()
plt.show()

# %%
# Get only rows with ChEMBL_id

ChEMBL_id = summary_bioassays[
    summary_bioassays["ChEMBL_ID"].notna()
].dropna(
    subset=["compound_diff", "Pathogen", "PubChem_compounds", "ChEMBL_compounds"]
).copy()

# Get top 10 assays with largest absolute difference (more compunds in top_PubChem)
top_PubChem = ChEMBL_id.loc[
    ChEMBL_id["compound_diff"].abs().nlargest(10).index,
    ["Pathogen", "PubChem_AID", "ChEMBL_ID", "PubChem_compounds", "ChEMBL_compounds", "compound_diff"]
]

# Clean index for display
top_PubChem = top_PubChem.reset_index(drop=True)

top_PubChem

# %%
# Reuse ChEMBL_id-filtered DataFrame
ChEMBL_id = summary_bioassays[
    summary_bioassays["ChEMBL_ID"].notna()
].dropna(
    subset=["compound_diff", "Pathogen", "PubChem_compounds", "ChEMBL_compounds"]
).copy()

# Plot: Compound Count Difference per Assay (only those with ChEMBL ID)
plt.figure(figsize=(16, 6))
sns.stripplot(
    data=ChEMBL_id,
    x="Pathogen",
    y="compound_diff",
    jitter=True,
    alpha=0.7,
    size=6,
    linewidth=0.5,
    edgecolor="black",
    color="#AA96FA"  # match with previous purple bar
)

plt.axhline(0, linestyle="--", color="gray")
plt.title("Compound Count Difference per Assay (PubChem − ChEMBL)", fontsize=14)
plt.ylabel("Compound Count Difference", fontsize=12)
plt.xlabel("Pathogen", fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()

# %%
# Bottom 10: ChEMBL >> PubChem
bottom_ChEMBL = ChEMBL_id.loc[
    ChEMBL_id["compound_diff"].nsmallest(10).index,
    ["Pathogen", "PubChem_AID", "ChEMBL_ID", "PubChem_compounds", "ChEMBL_compounds", "compound_diff"]
].reset_index(drop=True)

bottom_ChEMBL

# %% [markdown]
# ### New file to compare number of compunds

# %%
# Lod new file
assays_ChEMBL = pd.read_csv(DATA_PROCESSED / "ChEMBL_bioassays" / "assays_ChEMBL.csv")
assays_ChEMBL.head()

# %%
# Merge data from PubChem and from ChEMBL
import pandas as pd

# --- Load original dataframes ---
df1 = summary_bioassays.copy()       # Pathogen, PubChem_AID, PubChem_compounds, ChEMBL_ID
df2 = assays_ChEMBL.copy()           # assay_id, cpds

# --- Standardize ChEMBL IDs ---

# Remove "CHEMBL" prefix in assays_ChEMBL
df2["ChEMBL_ID_clean"] = df2["assay_id"].astype(str).str.replace("CHEMBL", "", regex=False)

# Remove any accidental prefix in summary_bioassays (safe to run)
df1["ChEMBL_ID_clean"] = df1["ChEMBL_ID"].astype(str).str.replace("CHEMBL", "", regex=False)

# --- Select required columns ---
df1_small = df1[[
    "Pathogen", "PubChem_AID", "PubChem_compounds", "ChEMBL_ID", "ChEMBL_ID_clean"
]]

df2_small = df2[[
    "assay_id", "cpds", "ChEMBL_ID_clean"
]]

# --- Merge ---
summary_bioassays_2 = df1_small.merge(
    df2_small,
    on="ChEMBL_ID_clean",
    how="left"
)

# --- Optional clean-up: drop helper column ---
summary_bioassays_2.drop(columns=["ChEMBL_ID_clean"], inplace=True)

# --- Preview ---
summary_bioassays_2.head()


# %%
# Checking if compounds numbers are equal in PubChem vs ChEMBL

summary_bioassays_2["same_compound_count"] = (
    summary_bioassays_2["PubChem_compounds"] == summary_bioassays_2["cpds"]
)

summary_bioassays_2["compound_diff"] = (
    summary_bioassays_2["PubChem_compounds"] - summary_bioassays_2["cpds"]
)

summary_bioassays_2.to_csv(DATA_PROCESSED / "summary_bioassays_with_chembl_compounds_2.csv", index=False)

# %%
summary_bioassays_2.head()

# %%
# Group by pathogen
summary_pubchem_chembl_2 = summary_bioassays_2.groupby("Pathogen").agg(
    PubChem_AIDs=("PubChem_AID", "nunique"),
    ChEMBL_IDs=("ChEMBL_ID", lambda x: x.notna().sum()),
    Matching_Compounds=("same_compound_count", "sum"),
    Total_with_ChEMBL=("same_compound_count", "count")  # Only rows with ChEMBL compound counts
).reset_index()


# Compute percentage of ChEMBL IDs
summary_pubchem_chembl_2["Percent_ChEMBL_IDs"] = (
    100 * summary_pubchem_chembl_2["ChEMBL_IDs"] / summary_pubchem_chembl_2["PubChem_AIDs"]
).round(1)

# Compute percentage of matched compounds (among only those with ChEMBL IDs)
summary_pubchem_chembl_2["Percent_Matching_Compounds"] = (
    100 * summary_pubchem_chembl_2["Matching_Compounds"] / summary_pubchem_chembl_2["Total_with_ChEMBL"]
).round(1)

# Drop intermediate column
summary_pubchem_chembl_2.drop(columns="Total_with_ChEMBL", inplace=True)

# Reorder columns
summary = summary_pubchem_chembl_2[
    ["Pathogen", "PubChem_AIDs", "ChEMBL_IDs", "Percent_ChEMBL_IDs", 
     "Matching_Compounds", "Percent_Matching_Compounds"]
]

summary_pubchem_chembl_2

# %%
# Step 1: Filter to rows with ChEMBL_ID
df_valid = summary_bioassays_2[summary_bioassays_2["ChEMBL_ID"].notna()].copy()

# Step 2: Check for duplicate ChEMBL_IDs
duplicated_ids = df_valid["ChEMBL_ID"].duplicated(keep=False)
num_duplicates = duplicated_ids.sum()

if num_duplicates > 0:
    print(f"⚠️ Found {num_duplicates} duplicated ChEMBL_ID entries. Dropping duplicates to ensure uniqueness.")
    # Optional: print them
    # print(df_valid[duplicated_ids].sort_values("ChEMBL_ID"))
    
    # Drop duplicates keeping the first occurrence
    df_valid = df_valid.drop_duplicates(subset="ChEMBL_ID", keep="first")

# %%
# Find duplicated ChEMBL_IDs (show all duplicates, not just the second+ occurrence)
duplicates = summary_bioassays_2[summary_bioassays_2["ChEMBL_ID"].duplicated(keep=False)]

# Sort for easier inspection
duplicates_sorted = duplicates.sort_values("ChEMBL_ID")

# Display the first few rows
duplicates_sorted.head(20)

# %%
# Data
summary = summary_pubchem_chembl_2.copy()
labels = summary["Pathogen"].values
percent_chembl = summary["Percent_ChEMBL_IDs"].values
percent_matching = summary["Percent_Matching_Compounds"].values

N = len(labels)
x = np.arange(N)

# ⏩ Adjust spacing: smaller bar width + larger offset
bar_width = 0.25
bar_spacing = 0.15  # controls the separation inside each group

# Colors
color_chembl = "#AA96FA"    # Purple
color_match  = "#BEE6B4"    # Green

# Plot
plt.figure(figsize=(16, 6))

# 1. Bar: Percent_ChEMBL_IDs
bars1 = plt.bar(x - bar_spacing, percent_chembl, width=bar_width,
                color=color_chembl, ec="black", zorder=2, label="% ChEMBL ID coverage")

# 2. Bar: Percent_Matching_Compounds
bars2 = plt.bar(x + bar_spacing, percent_matching, width=bar_width,
                color=color_match, ec="black", zorder=2, label="% Matching compound counts")

# Add labels on top of bars
for i in range(N):
    plt.text(x[i] - bar_spacing, percent_chembl[i] + 2, f"{percent_chembl[i]:.1f}%",
             ha='center', va='bottom', fontsize=9, zorder=3)
    plt.text(x[i] + bar_spacing, percent_matching[i] + 2, f"{percent_matching[i]:.1f}%",
             ha='center', va='bottom', fontsize=9, zorder=3)

# Aesthetics
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Percentage (%)")
plt.ylim(0, 110)
plt.title("Percent of ChEMBL-Mapped Assays and Matching Compounds by Pathogen")
plt.grid(linestyle="--", zorder=1, axis="y")

# Legend at bottom
plt.legend(
    loc="lower center",
    bbox_to_anchor=(0.18, 0.05),
    ncol=2,
    framealpha=1,
    edgecolor="k",
    prop={"size": 10}
)

plt.tight_layout()
plt.show()

# %%
# Get only rows with ChEMBL_id

ChEMBL_id = summary_bioassays_2[
    summary_bioassays_2["ChEMBL_ID"].notna()
].dropna(
    subset=["compound_diff", "Pathogen", "PubChem_compounds", "cpds"]
).copy()

# Get top 10 assays with largest absolute difference (more compunds in top_PubChem)
top_PubChem = ChEMBL_id.loc[
    ChEMBL_id["compound_diff"].abs().nlargest(10).index,
    ["Pathogen", "PubChem_AID", "ChEMBL_ID", "PubChem_compounds", "cpds", "compound_diff"]
]

# Clean index for display
top_PubChem = top_PubChem.reset_index(drop=True)

top_PubChem

# %%
# Reuse ChEMBL_id-filtered DataFrame
ChEMBL_id = summary_bioassays_2[
    summary_bioassays_2["ChEMBL_ID"].notna()
].dropna(
    subset=["compound_diff", "Pathogen", "PubChem_compounds", "cpds"]
).copy()

# Plot: Compound Count Difference per Assay (only those with ChEMBL ID)
plt.figure(figsize=(16, 6))
sns.stripplot(
    data=ChEMBL_id,
    x="Pathogen",
    y="compound_diff",
    jitter=True,
    alpha=0.7,
    size=6,
    linewidth=0.5,
    edgecolor="black",
    color="#AA96FA"  # match with previous purple bar
)

plt.axhline(0, linestyle="--", color="gray")
plt.title("Compound Count Difference per Assay (PubChem − ChEMBL)", fontsize=14)
plt.ylabel("Compound Count Difference", fontsize=12)
plt.xlabel("Pathogen", fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.grid(axis="y", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 07. ChEMBL_id found in PubChem

# %%
# Paths
pathogen_dir = DATA_PROCESSED / "ChEMBL_bioassays"

# Prepare set for faster lookup
arnau_chembl_ids = set(summary_bioassays["ChEMBL_ID"].dropna().astype(str))

# Mapping from short names in filenames → full pathogen names in summary_bioassays
pathogen_name_map = {
    "ecoli": "Escherichia coli",
    "saureus": "Staphylococcus aureus",
    "mtuberculosis": "Mycobacterium tuberculosis",
    "calbicans": "Candida albicans",
    "paeruginosa": "Pseudomonas aeruginosa",
    "kpneumoniae": "Klebsiella pneumoniae",
    "efaecium": "Enterococcus faecium",
    "pfalciparum": "Plasmodium falciparum",
    "abaumannii": "Acinetobacter baumannii",
    "enterobacter": "Enterobacter",
    "campylobacter": "Campylobacter",
    "smansoni": "Schistosoma mansoni",
    "hpylori": "Helicobacter pylori",
    "ngonorrhoeae": "Neisseria gonorrhoeae",
    "spneumoniae": "Streptococcus pneumoniae"
}

# Loop through pathogen files
rows = []
for file in pathogen_dir.glob("*_ChEMBL_data.csv"):
    df = pd.read_csv(file)
    short_name = file.stem.split("_ChEMBL_data")[0]
    pathogen_name = pathogen_name_map.get(short_name, short_name)  # fallback if not in dict

    file_chembl_ids = df["assay_chembl_id"].dropna().astype(str).unique()
    total = len(file_chembl_ids)
    matched = sum(aid in arnau_chembl_ids for aid in file_chembl_ids)
    percent = round(100 * matched / total, 1) if total > 0 else 0.0

    rows.append({
        "Pathogen": pathogen_name,
        "ChEMBL_assays_in_file": total,
        "Matched_in_summary": matched,
        "Percent_matched": percent
    })

chembl_match_summary = pd.DataFrame(rows).sort_values("Percent_matched", ascending=False)
chembl_match_summary

# %%
import pandas as pd
from pathlib import Path

# Paths
PROJECT_ROOT = Path("/Users/maria/Documents/Ersilia/PubChem/pubchem-antimicrobial-tasks")
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
summary_file = DATA_PROCESSED / "summary_bioassays.csv"
pathogen_dir = DATA_PROCESSED / "ChEMBL_bioassays"

# Clean up ChEMBL IDs
summary_bioassays["ChEMBL_ID"] = summary_bioassays["ChEMBL_ID"].astype(str)
summary_bioassays.loc[summary_bioassays["ChEMBL_ID"] == "nan", "ChEMBL_ID"] = None

# Mapping from filename to full pathogen names
pathogen_name_map = {
    "ecoli": "Escherichia coli",
    "saureus": "Staphylococcus aureus",
    "mtuberculosis": "Mycobacterium tuberculosis",
    "calbicans": "Candida albicans",
    "paeruginosa": "Pseudomonas aeruginosa",
    "kpneumoniae": "Klebsiella pneumoniae",
    "efaecium": "Enterococcus faecium",
    "pfalciparum": "Plasmodium falciparum",
    "abaumannii": "Acinetobacter baumannii",
    "enterobacter": "Enterobacter",
    "campylobacter": "Campylobacter",
    "smansoni": "Schistosoma mansoni",
    "hpylori": "Helicobacter pylori",
    "ngonorrhoeae": "Neisseria gonorrhoeae",
    "spneumoniae": "Streptococcus pneumoniae"
}

rows = []

for file in pathogen_dir.glob("*_ChEMBL_data.csv"):
    df = pd.read_csv(file, low_memory=False)
    short_name = file.stem.split("_ChEMBL_data")[0]
    pathogen = pathogen_name_map.get(short_name, short_name)

    # ChEMBL IDs from the file
    file_chembl_ids = df["assay_chembl_id"].dropna().astype(str).unique()
    total_chembl = len(file_chembl_ids)

    # ChEMBL IDs in summary for this pathogen
    summary_chembl_ids = summary_bioassays[
        (summary_bioassays["Pathogen"] == pathogen) &
        (summary_bioassays["ChEMBL_ID"].notna())
    ]["ChEMBL_ID"].astype(str).unique()

    # Matched = those in both
    matched = len(set(file_chembl_ids) & set(summary_chembl_ids))
    only_chembl = total_chembl - matched

    # PubChem assays for this pathogen with no ChEMBL_ID
    only_pubchem = summary_bioassays[
        (summary_bioassays["Pathogen"] == pathogen) &
        (summary_bioassays["ChEMBL_ID"].isna())
    ]["PubChem_AID"].nunique()

    percent_matched = round(100 * matched / total_chembl, 1) if total_chembl > 0 else 0.0

    rows.append({
        "Pathogen": pathogen,
        "ChEMBL_assays_in_file": total_chembl,
        "Matched_in_summary": matched,
        "Only_ChEMBL": only_chembl,
        "Only_PubChem": only_pubchem,
        "Both": matched,
        "Percent_matched": percent_matched
    })

# Final dataframe
id_overlap_summary = pd.DataFrame(rows).sort_values("Pathogen").reset_index(drop=True)

# Show it
id_overlap_summary

# %%
import matplotlib.pyplot as plt
import numpy as np

# Reuse your existing dataframe
df = id_overlap_summary.copy()

# Values
labels         = df["Pathogen"].values
matched        = df["Matched_in_summary"].values
only_chembl    = df["Only_ChEMBL"].values
only_pubchem   = df["Only_PubChem"].values

N = len(labels)
x = np.arange(N)

bar_width = 0.25
bar_spacing = 0  # controls the separation inside each group

# Colors
color_matched      = "#AA96FA"  # purple
color_only_chembl  = "#BEE6B4"  # green
color_only_pubchem = "#8DC7FA"  # light blue

# Plot
plt.figure(figsize=(16, 6))

plt.bar(x - spacing, matched, width=bar_width,
        color=color_matched, edgecolor="black", zorder=2, label="Matched (both IDs)")

plt.bar(x, only_chembl, width=bar_width,
        color=color_only_chembl, edgecolor="black", zorder=2, label="Only in ChEMBL")

plt.bar(x + spacing, only_pubchem, width=bar_width,
        color=color_only_pubchem, edgecolor="black", zorder=2, label="Only in PubChem")

# Labels and layout
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Number of Assays")
plt.title("Overlap of Assays by Pathogen (ChEMBL vs PubChem)", fontsize=14)
plt.grid(axis="y", linestyle="--", alpha=0.6, zorder=1)

# Legend at bottom
plt.legend(
    loc="lower center",
    bbox_to_anchor=(0.2, 0.85),
    ncol=3,
    framealpha=1,
    edgecolor="k",
    prop={"size": 10}
)

plt.tight_layout()
plt.show()

# %%
