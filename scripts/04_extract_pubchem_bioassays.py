import os
import sys
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import zipfile
import gzip
import shutil
from rdkit import Chem

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(root, "..", "src"))
from default import pathogens

datapath = os.path.join(root, "..", "data")
outpath = os.path.join(root, "..", "output")


def _aid_to_block(aid):
    start = ((aid - 1) // 1000) * 1000 + 1
    end = start + 999
    return f"{start:07d}_{end:07d}"


def extract_block_aids(datadir, aids, outdir):
    """
    Extract multiple AIDs that share the same ZIP block in a single pass.
    Opens the ZIP once, extracts all requested members.

    Returns a dict mapping aid -> csv_path (or aid -> None if not found).
    """
    os.makedirs(outdir, exist_ok=True)
    if not aids:
        return {}

    block = _aid_to_block(aids[0])
    zip_path = os.path.join(datadir, f"{block}.zip")
    results = {}

    if not os.path.exists(zip_path):
        for aid in aids:
            results[aid] = None
        print(f"ZIP block not found: {zip_path}")
        return results

    with zipfile.ZipFile(zip_path, "r") as zf:
        names = set(zf.namelist())
        for aid in aids:
            member = f"{block}/{aid}.csv.gz"
            csv_path = os.path.join(outdir, f"{aid}.csv")
            if member not in names:
                print(f"{member} not found inside {zip_path}")
                results[aid] = None
            else:
                with zf.open(member) as zipped_gz, gzip.open(zipped_gz, "rb") as fin, open(csv_path, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
                results[aid] = csv_path

    return results


def _smiles_to_inchikey(s):
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None


def clean_bioassay_csv(inputpath):
    df = pd.read_csv(inputpath)
    is_meta = df["PUBCHEM_RESULT_TAG"].fillna("").astype(str).str.startswith("RESULT")
    meta_df = df.loc[is_meta].copy()
    results_df = df.loc[~is_meta].copy()
    results_df = results_df[
        ["PUBCHEM_SID", "PUBCHEM_CID", "PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_ACTIVITY_OUTCOME"]
    ]
    results_df.rename(
        columns={
            "PUBCHEM_SID": "sid",
            "PUBCHEM_CID": "cid",
            "PUBCHEM_EXT_DATASOURCE_SMILES": "smiles",
            "PUBCHEM_ACTIVITY_OUTCOME": "activity",
        },
        inplace=True,
    )
    results_df["sid"] = pd.to_numeric(results_df["sid"], errors="coerce").astype("Int64")
    results_df["cid"] = pd.to_numeric(results_df["cid"], errors="coerce").astype("Int64")
    sids = results_df["sid"].tolist()
    cids = results_df[~results_df["cid"].isna()]["cid"].tolist()
    results_df = results_df[~results_df["smiles"].isna()]
    results_df["inchikey"] = results_df["smiles"].apply(_smiles_to_inchikey)
    results_df = results_df[~results_df["inchikey"].isna()]
    results_df["activity"] = results_df["activity"].str.strip().str.lower()
    actives = len(results_df[results_df["activity"] == "active"])
    inactives = len(results_df[results_df["activity"] == "inactive"])
    inconclusive = len(results_df[results_df["activity"] == "inconclusive"])
    unspecified = len(results_df[results_df["activity"] == "unspecified"])
    chem_probe = len(results_df[results_df["activity"] == "chemical probe"])
    results_df["activity"] = results_df["activity"].map(
        {"inactive": 0, "active": 1, "inconclusive": -1, "unspecified": 2, "chemical probe": 3}
    )
    summary_dict = {
        "sids": len(sids),
        "cids": len(cids),
        "smiles": len(results_df["smiles"]),
        "actives": actives,
        "inactives": inactives,
        "inconclusive": inconclusive,
        "unspecified": unspecified,
        "chem_probe": chem_probe,
    }
    return results_df, meta_df, summary_dict


def unique_cids_smiles_counts(
    input_dir,
    aids,
    cid_col="cid",
    smiles_col="smiles",
    inchikey_col="inchikey",
):
    cid_to_smiles = {}
    cid_to_inchikey = {}
    cid_to_count = {}

    for aid in tqdm(aids, desc="Collecting CID+SMILES+InChIKey+counts"):
        path = os.path.join(input_dir, f"{aid}.csv")
        if not os.path.exists(path):
            continue

        df = pd.read_csv(
            path,
            usecols=[cid_col, smiles_col, inchikey_col],
            dtype={cid_col: "Int64", smiles_col: str, inchikey_col: str},
        ).dropna(subset=[cid_col, smiles_col, inchikey_col])

        for cid, smi, ik in zip(df[cid_col], df[smiles_col], df[inchikey_col]):
            cid = int(cid)
            if cid not in cid_to_smiles:
                cid_to_smiles[cid] = smi
                cid_to_inchikey[cid] = ik
            cid_to_count[cid] = cid_to_count.get(cid, 0) + 1

    out_df = pd.DataFrame(
        {
            "cid": list(cid_to_smiles.keys()),
            "smiles": [cid_to_smiles[c] for c in cid_to_smiles],
            "inchikey": [cid_to_inchikey[c] for c in cid_to_smiles],
            "n_occurrences": [cid_to_count[c] for c in cid_to_smiles],
        }
    )
    return out_df


datadir = os.path.join(datapath, "raw", "03_data_zips")
unique_cids_path = os.path.join(datapath, "processed", "04_unique_cids")
os.makedirs(unique_cids_path, exist_ok=True)

unique_cids_all = {}
all_dfs = []

for p in pathogens:
    print(p)
    # Read the filtered AID list produced by 02_bioassays_not_in_chembl.py
    aids = pd.read_csv(
        os.path.join(datapath, "processed", "02_bioassays_to_keep", f"aids_{p.lower()}.csv")
    )["aid"].tolist()

    outdir_zips = os.path.join(datapath, "raw", "04_unzipped", f"{p.lower()}")
    outdir_clean = os.path.join(datapath, "processed", "04_extracted_bioassays", f"{p.lower()}")
    os.makedirs(outdir_zips, exist_ok=True)
    os.makedirs(outdir_clean, exist_ok=True)

    pathogen_dict = {}
    aids_with_data = []

    # Separate already-done AIDs from those that still need extraction
    aids_to_extract = []
    for a in aids:
        clean_csv = os.path.join(outdir_clean, f"{a}.csv")
        if os.path.exists(clean_csv):
            aids_with_data.append(a)
            # Rebuild summary_dict from the existing file
            df_existing = pd.read_csv(clean_csv)
            activity_counts = df_existing["activity"].value_counts()
            pathogen_dict[a] = {
                "sids": len(df_existing),
                "cids": df_existing["cid"].notna().sum(),
                "smiles": df_existing["smiles"].notna().sum(),
                "actives": int(activity_counts.get(1, 0)),
                "inactives": int(activity_counts.get(0, 0)),
                "inconclusive": int(activity_counts.get(-1, 0)),
                "unspecified": int(activity_counts.get(2, 0)),
                "chem_probe": int(activity_counts.get(3, 0)),
            }
        else:
            aids_to_extract.append(a)

    if aids_to_extract:
        print(f"  {len(aids_with_data)} AIDs already processed, extracting {len(aids_to_extract)} remaining")

    # Group remaining AIDs by ZIP block and extract each block once
    by_block = defaultdict(list)
    for a in aids_to_extract:
        by_block[_aid_to_block(a)].append(a)

    for block, block_aids in tqdm(by_block.items(), desc="  ZIP blocks"):
        extracted = extract_block_aids(datadir, block_aids, outdir_zips)
        for a, csv_path in extracted.items():
            if csv_path:
                results_df, meta_df, summary_dict = clean_bioassay_csv(csv_path)
                results_df.to_csv(os.path.join(outdir_clean, f"{a}.csv"), index=False)
                meta_df.to_csv(os.path.join(outdir_clean, f"{a}_meta.csv"), index=False)
                os.remove(csv_path)
                aids_with_data.append(a)
            else:
                summary_dict = {
                    "sids": 0, "cids": 0, "smiles": 0,
                    "actives": 0, "inactives": 0,
                    "inconclusive": 0, "unspecified": 0, "chem_probe": 0,
                }
            pathogen_dict[a] = summary_dict

    df_summary = pd.DataFrame.from_dict(
        pathogen_dict,
        orient="index",
        columns=["sids", "cids", "smiles", "actives", "inactives", "inconclusive", "unspecified", "chem_probe"],
    )
    df_summary.index.name = "aid"
    df_summary = df_summary.reset_index()
    df_summary.to_csv(os.path.join(outdir_clean, "summary.csv"), index=False)

    cid2smi_df = unique_cids_smiles_counts(outdir_clean, aids_with_data)
    unique_cids_all[p.lower()] = [len(cid2smi_df), len(aids_with_data)]
    cid2smi_df.to_csv(
        os.path.join(unique_cids_path, f"unique_cids_{p.lower()}.csv"),
        index=False,
    )
    all_dfs.append(cid2smi_df)

df_cids = pd.DataFrame.from_dict(unique_cids_all, orient="index", columns=["n_cid", "n_aid"])
df_cids.index.name = "pathogen"
df_cids = df_cids.reset_index()
df_cids.to_csv(os.path.join(unique_cids_path, "summary.csv"), index=False)

df_all = pd.concat(all_dfs, ignore_index=True)
df_global = (
    df_all.groupby("cid", as_index=False)
    .agg({"smiles": "first", "inchikey": "first", "n_occurrences": "sum"})
)
df_global = df_global.sort_values("n_occurrences", ascending=False)
df_global.to_csv(
    os.path.join(unique_cids_path, "unique_cids_all_pathogens.csv"),
    index=False,
)
