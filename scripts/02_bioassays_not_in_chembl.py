import os
import sys
import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
import stylia as st

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import pathogens

datapath = os.path.join(root, "..", "data")
configpath = os.path.join(datapath, "config")

aids_dict = {}

for pathogen in pathogens:
    # read list of aids for this pathogen
    bioassay_aids = set(
        pd.read_csv(
            os.path.join(datapath, "processed", "bioassays_summary", f"bioassays_{pathogen}.csv"),
            low_memory=False,
        )["aid"].astype(int)
        .tolist()
    )
    print(f"Total number of AIDs for {pathogen}:", len(bioassay_aids))
    
    # Extract raw bioassay information
    df = pd.read_csv(os.path.join(datapath,"raw", "bioassays", "bioassays.tsv"), sep="\t", low_memory=False)
    df = df[df["AID"].isin(bioassay_aids)]
    df["AID"] = df["AID"].astype(int)
    df["Source ID"] = df["Source ID"].fillna("").astype(str)
    df["Number of Tested CIDs"] = pd.to_numeric(df["Number of Tested CIDs"], errors="coerce")
    print("AIDS without molecules:", len(df[df["Number of Tested CIDs"]==0]))
    df = df[df["Number of Tested CIDs"]>0]
    aids_no_chembl_source = set(df.loc[~df["Source ID"].str.startswith("CHEMBL"), "AID"].tolist())
    print("Number of AIDs without CHEMBL source:", len(aids_no_chembl_source))

    # Extract ChEMBL assays linked to PubChem AIDs from the ChEMBL assay file
    # obtain ChEMBL assay ids for this pathogen
    chembl = pd.read_csv(os.path.join(configpath, "chembl_mappings", f"compounds_per_assay_{pathogen.lower()}.csv"))
    chembl_ids = set(chembl["assay_chembl_id"].tolist())
    print("Number of ChEMBL assays linked to this pathogen:", len(chembl_ids))
    assays = pd.read_csv(os.path.join(configpath, "chembl_mappings","assays.csv"))
    assays = assays[assays["chembl_id"].isin(chembl_ids)]
    chembl_ids_with_pubchem_source = set(assays[assays["src_assay_id"].isin(bioassay_aids)]["chembl_id"].tolist())
    pubchem_aids_linked_in_chembl = set(assays[assays["src_assay_id"].isin(bioassay_aids)]["src_assay_id"].tolist())
    print("Number of ChEMBL IDS with Pubchem Source:", len(chembl_ids_with_pubchem_source), len(pubchem_aids_linked_in_chembl))
    aids_in_chembl=aids_no_chembl_source.intersection(pubchem_aids_linked_in_chembl)
    aids_not_in_chembl = aids_no_chembl_source.difference(pubchem_aids_linked_in_chembl)
    print("Number of AIDs with ChEMBL source not captured by PubChem:", len(aids_in_chembl))
    print("Number of PubChem AIDS never linked to ChEMBL", len(aids_not_in_chembl))

    # get mapping of PubChem AIDs linked in ChEMBL
    assays_dict = {}

    for aid in aids_in_chembl:
        chembl_id = assays.loc[assays["src_assay_id"] == aid, "chembl_id"].iloc[0]
        assays_dict[aid] = chembl_id
    
    # Obtain the AIDS that are significantly mismatched in terms of number of compounds tested between PubChem and ChEMBL
    df_chembl_0 = df[df["Source ID"].str.startswith("CHEMBL")]
    df_chembl_1 = df[df["AID"].isin(aids_in_chembl)].copy()
    df_chembl_1["Source ID"] = df_chembl_1["AID"].map(assays_dict).fillna(df_chembl_1["Source ID"])
    df_chembl = pd.concat([df_chembl_0, df_chembl_1], axis=0)
    df_chembl = df_chembl.drop_duplicates() # just in case
    df_chembl.rename(columns={"Number of Tested CIDs": "pubchem_cpds"}, inplace=True)

    merged = df_chembl.merge(
        chembl.rename(columns={"n_compounds": "chembl_cpds"}),
        left_on="Source ID",
        right_on="assay_chembl_id",
        how="left",
    )

    # AIDs where CHEMBL Source ID exists AND chembl mapping exists AND cpds mismatch
    cpds_mismatch = merged.loc[merged["pubchem_cpds"] != merged["chembl_cpds"],
        ["AID", "Source ID", "chembl_cpds", "pubchem_cpds"]
    ]

    cpds_mismatch["mismatch"] = (cpds_mismatch["pubchem_cpds"]-cpds_mismatch["chembl_cpds"])/cpds_mismatch["pubchem_cpds"]

    print("Number of AIDs with CHEMBL source and CPDs mismatch:", len(cpds_mismatch))
    cpds_mismatch.to_csv(
        os.path.join(datapath, "processed", "bioassays_to_keep", f"chembl_cpds_mismatch_{pathogen.lower()}.csv"),
        index=False,
    )

    # Keep rules:

    # 1. Compounds mismatch is over 0.1 or Chembl ID reports 0 compounds
    # 2. Compound number is over 100 in either the mismatched or the not in ChEMBL data

    cpds_mismatch = cpds_mismatch[(cpds_mismatch["mismatch"] > 0.1) | (cpds_mismatch["chembl_cpds"].isna())]
    cpds_mismatch  =  cpds_mismatch[cpds_mismatch["pubchem_cpds"]>100]
    aids_mismatched = set(cpds_mismatch["AID"].tolist())
    print("Number of mismatched assays kept:", len(aids_mismatched))

    df_ = df[df["AID"].isin(aids_not_in_chembl)]
    df_ = df_[df_["Number of Tested CIDs"]>100]
    aids_not_in_chembl_keep = set(df_["AID"].tolist())

    #final list of AIDS
    aids_to_use = set.union(aids_mismatched, aids_not_in_chembl_keep)
    print("Final AIDS to consider:", len(aids_to_use))
    aids_dict[pathogen] = [len(bioassay_aids), len(aids_no_chembl_source),len(aids_not_in_chembl), len(aids_in_chembl),len(aids_to_use), len(aids_mismatched), len(aids_not_in_chembl_keep)]
    filepath = os.path.join(datapath, "processed", "bioassays_to_keep", f"aids_{pathogen.lower()}.csv")
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aid"]) 
        for aid in aids_to_use:
            writer.writerow([aid])

df_plot = pd.DataFrame.from_dict(
    aids_dict,
    orient="index",
    columns=["bioassay_aids", "aids_no_chembl_source", "aids_not_in_chembl", "aids_in_chembl","aids_to_use", "aids_mismatched", "aids_not_in_chembl_keep"]
)

df_plot.index.name = "pathogen"
df_plot = df_plot.reset_index()
df_plot.to_csv(os.path.join(datapath, "processed", "bioassays_to_keep", "summary.csv"), index=False)