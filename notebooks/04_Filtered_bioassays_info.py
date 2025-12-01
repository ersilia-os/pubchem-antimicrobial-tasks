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
# # Assay Statistics Summary (Notebook 04)
#
# This notebook summarizes the filtered PubChem assays:
# - Total assays per pathogen
# - Number and % of assays with ChEMBL ID
# - Global % with ChEMBL ID

# %%
import pandas as pd
from pathlib import Path

# %%
# Set paths
DATA_PROCESSED = Path("../data/processed")

# %% [markdown]
# ## 1. Load filtered assay metadata (v2)

# %%
# Load dataframe
df = pd.read_csv(DATA_PROCESSED / "filtered_description_with_organisms_v2_REBUILT.csv")

# %% [markdown]
# Some assays might have more than one pathogen detected, because we have searched for taxid and assay organism.

# %%
# Convert Pathogen string to list
df["Pathogen_list"] = df["Pathogen"].str.split(", ")

# Find rows where the list has more than 1 pathogen
df_multi = df[df["Pathogen_list"].apply(len) > 1].copy()

# Convert the list to a comma-separated string to allow drop_duplicates
df_multi["Pathogen_str"] = df_multi["Pathogen_list"].apply(lambda x: ", ".join(x))

# Show unique (AID, Pathogen_str) pairs
multi_pathogen_assays = df_multi[["AID", "Pathogen_str"]].reset_index(drop=True)

# Display result
multi_pathogen_assays

# %% [markdown]
# For the time being, we will ignore these assays until manually checked for true pathogen.

# %%
# Get list of AIDs with multiple pathogens
aids_to_exclude = set(multi_pathogen_assays["AID"])

# Filter df to exclude those AIDs
df = df[~df["AID"].isin(aids_to_exclude)]

# %% [markdown]
# ## 2. Count total assays per pathogen

# %%
summary_assays = (
    df.groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="Total_Assays")
    .sort_values("Total_Assays", ascending=False)
    .reset_index(drop=True)
)

summary_assays

# %% [markdown]
# ## 3. Count assays with ChEMBL ID

# %%
df["Has_ChEMBL"] = df["ChEMBLid"].notna()

summary_chembl = (
    df[df["Has_ChEMBL"]]
    .groupby("Pathogen")["AID"]
    .nunique()
    .reset_index(name="AIDs_with_ChEMBL")
    .sort_values("AIDs_with_ChEMBL", ascending=False)
)

summary_chembl

# %% [markdown]
# ## 4. Merge and compute percentages

# %%
# Merge the two summaries
df_summary = summary_assays.merge(summary_chembl, on="Pathogen", how="left")

# Fill missing values with 0 (i.e., no ChEMBL ID found)
df_summary["AIDs_with_ChEMBL"] = df_summary["AIDs_with_ChEMBL"].fillna(0).astype(int)

# Calculate percentage of assays with ChEMBL ID
df_summary["Percent_with_ChEMBL"] = (
    100 * df_summary["AIDs_with_ChEMBL"] / df_summary["Total_Assays"]
).round(1)

df_summary

# %% [markdown]
# ## 5. Global statistics

# %%
total_assays = df_summary["Total_Assays"].sum()
total_with_chembl = df_summary["AIDs_with_ChEMBL"].sum()
global_percentage = 100 * total_with_chembl / total_assays

print(f"Global summary:")
print(f"- Total filtered assays: {total_assays:,}")
print(f"- With ChEMBL ID: {total_with_chembl:,} ({global_percentage:.1f}%)")

# %% [markdown]
# ## 6. Llist of ChEMBL ids

# %%
# Filter rows with a ChEMBL ID
chembl_pairs = df[df["ChEMBLid"].notna()][["AID", "ChEMBLid"]]
chembl_pairs.head()

# %%
# Save to CSV
output_path = DATA_PROCESSED / "aid_chembl_pairs.csv"
chembl_pairs.to_csv(output_path, index=False)
