import os
import sys
import pandas as pd
import csv
import numpy as np
import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import pathogens, MIN_COMPOUNDS, CHEMBL_MISMATCH_THRESHOLD, PATHOGEN_TO_CODE

datapath = os.path.join(root, "..", "data")
configpath = os.path.join(datapath, "config")
outpath = os.path.join(datapath, "processed", "02_bioassays_to_keep")
plotpath = os.path.join(root, "..", "output", "02_bioassays_not_in_chembl")

os.makedirs(outpath, exist_ok=True)
os.makedirs(plotpath, exist_ok=True)

# Load once — file is large (~hundreds of MB)
df_bioassays_raw = pd.read_csv(
    os.path.join(datapath, "raw", "01_bioassays", "bioassays.tsv"),
    sep="\t",
    low_memory=False,
)
df_bioassays_raw["AID"] = df_bioassays_raw["AID"].astype(int)
df_bioassays_raw["Source ID"] = df_bioassays_raw["Source ID"].fillna("").astype(str)
df_bioassays_raw["Number of Tested CIDs"] = pd.to_numeric(
    df_bioassays_raw["Number of Tested CIDs"], errors="coerce"
)

aids_dict = {}

for pathogen in pathogens:
    # read list of aids for this pathogen
    df_bioassay = pd.read_csv(
        os.path.join(datapath, "processed", "00_bioassays_summary", f"bioassays_{pathogen.lower()}.csv"),
        low_memory=False,
    )
    bioassay_aids = set(df_bioassay["aid"].astype(int).tolist())
    print(f"Total number of AIDs for {pathogen}:", len(bioassay_aids))
    cpd_count_all_bioassays = pd.to_numeric(df_bioassay["cnt"], errors="coerce").sum()

    df = df_bioassays_raw[df_bioassays_raw["AID"].isin(bioassay_aids)].copy()
    print("AIDS without molecules:", len(df[df["Number of Tested CIDs"] == 0]))
    df = df[df["Number of Tested CIDs"] > 0]
    aids_no_chembl_source = set(df.loc[~df["Source ID"].str.startswith("CHEMBL"), "AID"].tolist())
    print("Number of AIDs without CHEMBL source:", len(aids_no_chembl_source))
    cpd_count_aids_no_chembl_source = (
        df.loc[~df["Source ID"].str.startswith("CHEMBL"), "Number of Tested CIDs"]
        .astype(float)
        .sum()
    )

    # Extract ChEMBL assays linked to PubChem AIDs from the ChEMBL assay file
    chembl = pd.read_csv(
        os.path.join(configpath, "chembl_mappings", f"compounds_per_assay_{pathogen.lower()}.csv")
    )
    chembl_ids = set(chembl["assay_chembl_id"].tolist())
    print("Number of ChEMBL assays linked to this pathogen:", len(chembl_ids))
    cpd_count_chembl = pd.to_numeric(chembl["n_compounds"], errors="coerce").sum()
    assays = pd.read_csv(os.path.join(configpath, "chembl_mappings", "assays.csv"), low_memory=False)
    assays["src_assay_id"] = pd.to_numeric(assays["src_assay_id"], errors="coerce")
    bioassay_aids_set = set(bioassay_aids)
    assays = assays[assays["chembl_id"].isin(chembl_ids)]
    chembl_ids_with_pubchem_source = set(
        assays[assays["src_assay_id"].isin(bioassay_aids_set)]["chembl_id"].tolist()
    )
    pubchem_aids_linked_in_chembl = set(
        assays[assays["src_assay_id"].isin(bioassay_aids_set)]["src_assay_id"].dropna().astype(int).tolist()
    )
    print(
        "Number of ChEMBL assays with PubChem source:", len(chembl_ids_with_pubchem_source),
        "| PubChem AIDs linked in ChEMBL:", len(pubchem_aids_linked_in_chembl),
    )
    aids_in_chembl = aids_no_chembl_source.intersection(pubchem_aids_linked_in_chembl)
    aids_not_in_chembl = aids_no_chembl_source.difference(pubchem_aids_linked_in_chembl)
    print("AIDs not sourced from ChEMBL but present in ChEMBL:", len(aids_in_chembl))
    print("Number of PubChem AIDS never linked to ChEMBL", len(aids_not_in_chembl))

    cpd_count_aids_in_chembl = (
        df.loc[df["AID"].isin(aids_in_chembl), "Number of Tested CIDs"].astype(int).sum()
    )
    cpd_count_aids_not_in_chembl = (
        df.loc[df["AID"].isin(aids_not_in_chembl), "Number of Tested CIDs"].astype(int).sum()
    )

    # get mapping of PubChem AIDs linked in ChEMBL
    assays_dict = {}
    for aid in aids_in_chembl:
        chembl_id = assays.loc[assays["src_assay_id"] == aid, "chembl_id"].iloc[0]
        assays_dict[aid] = chembl_id

    # Obtain the AIDs that are significantly mismatched in number of compounds
    # between PubChem and ChEMBL
    df_chembl_0 = df[df["Source ID"].str.startswith("CHEMBL")]
    df_chembl_1 = df[df["AID"].isin(aids_in_chembl)].copy()
    df_chembl_1["Source ID"] = df_chembl_1["AID"].map(assays_dict).fillna(
        df_chembl_1["Source ID"]
    )
    df_chembl = pd.concat([df_chembl_0, df_chembl_1], axis=0).drop_duplicates()
    df_chembl.rename(columns={"Number of Tested CIDs": "pubchem_cpds"}, inplace=True)

    merged = df_chembl.merge(
        chembl.rename(columns={"n_compounds": "chembl_cpds"}),
        left_on="Source ID",
        right_on="assay_chembl_id",
        how="left",
    )

    cpds_mismatch = merged.loc[
        merged["pubchem_cpds"] != merged["chembl_cpds"],
        ["AID", "Source ID", "chembl_cpds", "pubchem_cpds"],
    ].copy()
    cpds_mismatch["mismatch"] = (
        (cpds_mismatch["pubchem_cpds"] - cpds_mismatch["chembl_cpds"])
        / cpds_mismatch["pubchem_cpds"]
    )

    print("Number of AIDs with CHEMBL source and CPDs mismatch:", len(cpds_mismatch))
    cpds_mismatch.to_csv(
        os.path.join(
            outpath,
            f"chembl_cpds_mismatch_{pathogen.lower()}.csv",
        ),
        index=False,
    )

    # Keep rules:
    # 1. Compound mismatch is over CHEMBL_MISMATCH_THRESHOLD or ChEMBL reports 0 compounds
    # 2. PubChem compound count exceeds MIN_COMPOUNDS
    cpds_mismatch = cpds_mismatch[
        (cpds_mismatch["mismatch"] > CHEMBL_MISMATCH_THRESHOLD)
        | (cpds_mismatch["chembl_cpds"].isna())
    ]
    cpds_mismatch = cpds_mismatch[cpds_mismatch["pubchem_cpds"] > MIN_COMPOUNDS]
    cpds_mismatch.to_csv(
        os.path.join(
            outpath,
            f"chembl_assays_in_pubchem_{PATHOGEN_TO_CODE[pathogen]}.csv",
        ),
        index=False,
    )
    aids_mismatched = set(cpds_mismatch["AID"].tolist())
    print("Number of mismatched assays kept:", len(aids_mismatched))

    df_ = df[df["AID"].isin(aids_not_in_chembl)]
    df_ = df_[df_["Number of Tested CIDs"] > MIN_COMPOUNDS]
    aids_not_in_chembl_keep = set(df_["AID"].tolist())

    df_[["AID", "Number of Tested CIDs"]].rename(
        columns={"AID": "aid", "Number of Tested CIDs": "n_compounds"}
    ).to_csv(
        os.path.join(outpath, f"aids_not_in_chembl_{pathogen.lower()}.csv"),
        index=False,
    )

    # final list of AIDs
    aids_to_use = set.union(aids_mismatched, aids_not_in_chembl_keep)
    print("Final AIDS to consider:", len(aids_to_use))
    cpd_count_aids_to_use = (
        df.loc[df["AID"].isin(aids_to_use), "Number of Tested CIDs"].astype(int).sum()
    )
    cpd_count_aids_mismatched = (
        df.loc[df["AID"].isin(aids_mismatched), "Number of Tested CIDs"].astype(int).sum()
    )
    cpd_count_aids_not_in_chembl_keep = (
        df.loc[df["AID"].isin(aids_not_in_chembl_keep), "Number of Tested CIDs"].astype(int).sum()
    )

    aids_dict[pathogen] = [
        len(bioassay_aids), cpd_count_all_bioassays,
        len(chembl_ids), cpd_count_chembl,
        len(aids_no_chembl_source), cpd_count_aids_no_chembl_source,
        len(aids_not_in_chembl), cpd_count_aids_not_in_chembl,
        len(aids_in_chembl), cpd_count_aids_in_chembl,
        len(aids_to_use), cpd_count_aids_to_use,
        len(aids_mismatched), cpd_count_aids_mismatched,
        len(aids_not_in_chembl_keep), cpd_count_aids_not_in_chembl_keep,
    ]

    filepath = os.path.join(
        outpath, f"aids_{pathogen.lower()}.csv"
    )
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aid"])
        for aid in aids_to_use:
            writer.writerow([aid])

df_plot = pd.DataFrame.from_dict(
    aids_dict,
    orient="index",
    columns=[
        "bioassay_aids", "cpd_count_bioassay_aids",
        "chembl_ids", "cpd_count_chembl_ids",
        "aids_no_chembl_source", "cpd_count_no_chembl_source",
        "aids_not_in_chembl", "cpd_count_not_in_chembl",
        "aids_in_chembl", "cpd_count_in_chembl",
        "aids_to_use", "cpd_count_to_use",
        "aids_mismatched", "cpd_count_mismatched",
        "aids_not_in_chembl_keep", "cpd_count_not_in_chembl_keep",
    ],
)

df_plot.index.name = "pathogen"
df_plot = df_plot.reset_index()
df_plot.to_csv(
    os.path.join(outpath, "summary.csv"), index=False
)


# --- Plots ---
# Format: slide | Style: ersilia — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("ersilia")


def plot_aids_to_keep(df, outdir):
    nc = stylia.NamedColors()
    x = np.arange(len(df["pathogen"]))
    width = 0.25

    fig, axs = stylia.create_figure(1, 1)
    ax = axs.next()
    ax.bar(x - width, df["aids_to_use"], width, label="AIDs to use", color=nc.plum)
    ax.bar(x, df["aids_mismatched"], width, label="Mismatched", color=nc.mint)
    ax.bar(x + width, df["aids_not_in_chembl_keep"], width, label="Not in ChEMBL", color=nc.blue)
    ax.set_xticks(x)
    ax.set_xticklabels(df["pathogen"], rotation=45, ha="right")
    ax.legend()
    stylia.label(ax, title="AIDs kept per pathogen", xlabel="", ylabel="Number of AIDs")

    stylia.save_figure(os.path.join(outdir, "02_aids_to_keep.png"))


def plot_chembl_pubchem_overview(df, outdir):
    nc = stylia.NamedColors()
    y = np.arange(len(df))

    df = df.copy()
    df["pct_aids_not_in_chembl"] = df["aids_not_in_chembl"] / df["bioassay_aids"] * 100
    df["pct_aids_mismatched"] = df["aids_mismatched"] / df["bioassay_aids"] * 100
    df["pct_cpds_not_in_chembl"] = df["cpd_count_not_in_chembl"] / df["cpd_count_bioassay_aids"] * 100
    df["pct_cpds_mismatched"] = df["cpd_count_mismatched"] / df["cpd_count_bioassay_aids"] * 100
    df["pct_cpds_to_use"] = df["cpd_count_to_use"] / df["cpd_count_bioassay_aids"] * 100

    fig, axs = stylia.create_figure(1, 3)

    ax = axs.next()
    ax.barh(y - 0.2, df["pct_aids_not_in_chembl"], height=0.4, label="Not linked to ChEMBL", color=nc.plum)
    ax.barh(y + 0.2, df["pct_aids_mismatched"], height=0.4, label="Mismatched", color=nc.mint)
    ax.set_yticks(y)
    ax.set_yticklabels(df["pathogen"])
    ax.legend()
    stylia.label(ax, title="% of assay IDs (vs total)", xlabel="", ylabel="")

    ax = axs.next()
    ax.barh(y - 0.2, df["pct_cpds_not_in_chembl"], height=0.4, label="Not linked to ChEMBL", color=nc.plum)
    ax.barh(y + 0.2, df["pct_cpds_mismatched"], height=0.4, label="Mismatched", color=nc.mint)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.legend()
    stylia.label(ax, title="% of compounds in assays (vs total)", xlabel="", ylabel="")

    ax = axs.next()
    ax.barh(y, df["pct_cpds_to_use"], color=nc.blue)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    stylia.label(ax, title="% of compounds used (vs total)", xlabel="", ylabel="")

    stylia.save_figure(os.path.join(outdir, "02_chembl_pubchem_overview.png"))


plot_aids_to_keep(df_plot, plotpath)
plot_chembl_pubchem_overview(df_plot, plotpath)
