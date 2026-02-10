import os
import sys
import pandas as pd

import zipfile
import gzip
import shutil

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(root, "..", "src"))
from default import pathogens

datapath = os.path.join(root, "..", "data")
outpath = os.path.join(root, "..", "output")


def extract_aid_csv(datadir,aid,outdir):
    """
    Given a PubChem BioAssay AID and a local PubChem /Data directory containing
    0000001_0001000.zip blocks, extract the AID's .csv.gz and gunzip it to .csv.

    Returns the path to the produced .csv file.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    start = ((aid - 1) // 1000) * 1000 + 1
    end = start + 999
    block = f"{start:07d}_{end:07d}"
    zip_path = os.path.join(datadir,f"{block}.zip")

    if not os.path.exists(zip_path):
        raise FileNotFoundError(f"ZIP block not found: {zip_path}")

    csv_path = os.path.join(outdir,f"{aid}.csv")
    member = f"{block}/{aid}.csv.gz"
    with zipfile.ZipFile(zip_path, "r") as zf:
        if member not in zf.namelist():
            raise FileNotFoundError(f"{member} not found inside {zip_path}")
        with zf.open(member) as zipped_gz, gzip.open(zipped_gz, "rb") as fin, open(csv_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    return csv_path

def clean_bioassay_csv(inputpath):
    df = pd.read_csv(inputpath)
    print(df.head())
    is_meta = df["PUBCHEM_RESULT_TAG"].fillna("").astype(str).str.startswith("RESULT")
    meta_df = df.loc[is_meta].copy()
    results_df = df.loc[~is_meta].copy()
    results_df = results_df[["PUBCHEM_SID","PUBCHEM_CID","PUBCHEM_EXT_DATASOURCE_SMILES","PUBCHEM_ACTIVITY_OUTCOME"]]
    results_df.rename(columns={"PUBCHEM_SID":"sid",
                               "PUBCHEM_CID":"cid",
                               "PUBCHEM_EXT_DATASOURCE_SMILES":"smiles",
                               "PUBCHEM_ACTIVITY_OUTCOME":"activity"}, inplace=True)
    results_df["sid"] = pd.to_numeric(results_df["sid"], errors="coerce").astype("Int64")
    results_df["cid"] = pd.to_numeric(results_df["cid"], errors="coerce").astype("Int64")
    results_df = results_df[~results_df["smiles"].isna()]
    results_df["activity"] = results_df["activity"].str.strip().str.capitalize()
    results_df = results_df[results_df["activity"].isin(["Active", "Inactive"])]
    results_df["activity"] = results_df["activity"].map({"Inactive": 0, "Active": 1})
    return results_df, meta_df
    

datadir = os.path.join(datapath, "raw", "bioassays", "Data")

for p in pathogens:
    aids = pd.read_csv(os.path.join(datapath,"processed","bioassays_to_keep", f"aids_{p.lower()}.csv"))["aid"].tolist()
    outdir_zips = os.path.join(datapath, "raw", "unzipped", f"{p}")
    outdir_clean = os.path.join(outpath, "results", f"{p}")
    if not os.path.exists(outdir_zips):
        os.makedirs(outdir_zips)
    if not os.path.exists(outdir_clean):
        os.makedirs(outdir_clean)
    for a in aids:
        print(a)
        csv_path = extract_aid_csv(datadir,a,outdir_zips)
        results_df, meta_df = clean_bioassay_csv(csv_path)
        results_df.to_csv(os.path.join(outdir_clean, f"{a}.csv"), index=False)
        meta_df.to_csv(os.path.join(outdir_clean, f"{a}_meta.csv"), index=False)


