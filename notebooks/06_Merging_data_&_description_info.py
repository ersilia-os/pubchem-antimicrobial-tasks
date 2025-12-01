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
summary_description = pd.read_csv(DATA_PROCESSED / "filtered_description_with_organisms_v2_REBUILT.csv")
summary_description = summary_description.sort_values('AID', ascending = True ).reset_index(drop = True)
summary_description

# %%
summary_description["Has_ChEMBL"] = summary_description["ChEMBLid"].notna()

summary_chembl = (
    summary_description[summary_description["Has_ChEMBL"]]
    .groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="AIDs_with_ChEMBL")
    .sort_values("AIDs_with_ChEMBL", ascending=False)
)

summary_chembl

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
#

# %% [markdown]
# ## 06. Check number of compounds in PubChem vs ChEMBL

# %%
# Checking if compounds numbers are equal in PubChem vs ChEMBL

summary_bioassays["same_compound_count"] = (
    summary_bioassays["PubChem_compounds"] == summary_bioassays["ChEMBL_compounds"]
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
# Filter rows where compound_diff is not zero and not null
df_check = summary_bioassays[
    (summary_bioassays["compound_diff"].notna()) &
    (summary_bioassays["compound_diff"] != 0)
]

# Check how many of those have ChEMBL_compounds == 0
count_zero_chembl = (df_check["ChEMBL_compounds"] == 0).sum()
total_diff = len(df_check)

print(f"Total rows with non-zero compound_diff: {total_diff}")
print(f"Rows where ChEMBL_compounds == 0: {count_zero_chembl}")
print(f"✅ All zero? {count_zero_chembl == total_diff}")

# %%
# Get top 10 assays with smallest absolute difference (more compunds in ChEMBL)
top_ChEMBL = ChEMBL_id.loc[
    ChEMBL_id["compound_diff"].abs().nsmallest(10).index,
    ["Pathogen", "PubChem_AID", "ChEMBL_ID", "PubChem_compounds", "ChEMBL_compounds", "compound_diff"]
]

# Clean index for display
top_ChEMBL = top_ChEMBL.reset_index(drop=True)

top_ChEMBL
