import os
import sys
import pandas as pd
from tqdm import tqdm
from itertools import combinations
import zipfile
import gzip
import shutil
from rdkit import Chem

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
            print(f"{member} not found inside {zip_path}")
            return
        with zf.open(member) as zipped_gz, gzip.open(zipped_gz, "rb") as fin, open(csv_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    return csv_path

def _smiles_to_inchikey(s):
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol)
    except:
        return None

def clean_bioassay_csv(inputpath):
    df = pd.read_csv(inputpath)
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
    sids = results_df["sid"].tolist()
    cids = results_df[~results_df["cid"].isna()]["cid"].tolist()
    results_df = results_df[~results_df["smiles"].isna()]
    results_df["inchikey"] = results_df["smiles"].apply(_smiles_to_inchikey)
    results_df = results_df[~results_df["inchikey"].isna()]
    results_df["activity"] = results_df["activity"].str.strip().str.lower()
    actives = len(results_df[results_df["activity"]=="active"])
    inactives = len(results_df[results_df["activity"]=="inactive"])
    inconclusive = len(results_df[results_df["activity"]=="inconclusive"])
    unspecified = len(results_df[results_df["activity"]=="unspecified"])
    chem_probe = len(results_df[results_df["activity"]=="chemical probe"])
    results_df["activity"] = results_df["activity"].map({"inactive": 0, "active": 1, "inconclusive": -1, "unspecified": 2, "chemical probe":3})
    summary_dict = {"sids":len(sids), "cids":len(cids), "smiles":len(results_df["smiles"]),
                    "actives": actives,"inactives":inactives,
                    "inconclusive": inconclusive, "unspecified":unspecified,
                    "chem_probe": chem_probe}
    return results_df, meta_df, summary_dict

def _aid_to_cids(input_path):
    cid_col = "cid"
    df = pd.read_csv(input_path, usecols=[cid_col], dtype={cid_col: "Int64"})[cid_col]
    return set(df.dropna().astype(int).unique())

def compare_aid_cid_overlap(input_dir, aids, out_csv,chunk_size=200_000):
    aid_sets: dict[int, set[int]] = {}
    for aid in tqdm(aids, desc="Loading CID sets"):
        p = os.path.join(input_dir,f"{aid}.csv")
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing file: {aid}")
        cids = _aid_to_cids(p)
        aid_sets[aid] = cids
    items = list(aid_sets.items())
    n = len(items)
    total_pairs = n * (n - 1) // 2
    cols = ["aid_1", "aid_2", "n_cids_1", "n_cids_2", "intersection", "ratio_min", "jaccard"]
    buffer = []
    wrote_header = False

    for (aid1, c1), (aid2, c2) in tqdm(list(combinations(items, 2)),total=total_pairs, desc="Computing overlaps"):
        n1, n2 = len(c1), len(c2)
        inter = len(c1 & c2)
        union = len(c1 | c2)
        ratio_min = inter / min(len(c1), len(c2)) if min(len(c1), len(c2)) else 0.0
        jaccard = inter / union if union else 0.0
        buffer.append([aid1, aid2, n1, n2, inter, ratio_min, jaccard])
        if len(buffer) >= chunk_size:
            pd.DataFrame(buffer, columns=cols).to_csv(out_csv, mode="a", header=not wrote_header, index=False)
            wrote_header = True
            buffer.clear()
    if buffer:
        pd.DataFrame(buffer, columns=cols).to_csv(out_csv, mode="a", header=not wrote_header, index=False)

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

            # store first occurrence
            if cid not in cid_to_smiles:
                cid_to_smiles[cid] = smi
                cid_to_inchikey[cid] = ik

            # increment occurrence count
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

def count_unique_cids(input_dir, aids, cid_col="cid", smi_col ="smiles"):
    cid2smi = {}
    for aid in tqdm(aids, desc="Collecting CIDs+SMILES"):
        path = os.path.join(input_dir, f"{aid}.csv")
        if not os.path.exists(path):
            continue 
        df = pd.read_csv(path, usecols=[cid_col, smi_col], dtype={cid_col: "Int64", smi_col: str})
        df = df[[cid_col,smi_col]].dropna()
        for cid, smi in zip(df[cid_col], df[smi_col]):
            cid = int(cid)
            if cid not in cid2smi:
                cid2smi[cid] = smi
    return cid2smi

datadir = os.path.join(datapath, "raw", "bioassays", "Data")
unique_cids_all = {}
all_dfs = []
for p in pathogens:
    print(p)
    aids = pd.read_csv(os.path.join(datapath,"processed","bioassays_summary", f"bioassays_{p.lower()}.csv"))["aid"].tolist()
    outdir_zips = os.path.join(datapath, "raw", "unzipped", f"{p.lower()}")
    outdir_clean = os.path.join(datapath, "processed","extracted_bioassays", f"{p.lower()}")
    if not os.path.exists(outdir_zips):
        os.makedirs(outdir_zips)
    if not os.path.exists(outdir_clean):
        os.makedirs(outdir_clean)
    pathogen_dict = {}
    aids_with_data = []
    for a in aids:
        csv_path= extract_aid_csv(datadir,a,outdir_zips)
        if csv_path:
            aids_with_data += [a]
    """
    for a in aids:
        csv_path= extract_aid_csv(datadir,a,outdir_zips)
        if csv_path:
            results_df, meta_df, summary_dict = clean_bioassay_csv(csv_path)
            results_df.to_csv(os.path.join(outdir_clean, f"{a}.csv"), index=False)
            meta_df.to_csv(os.path.join(outdir_clean, f"{a}_meta.csv"), index=False)
            aids_with_data += [a]
        else:
            summary_dict = {"sids":0, "cids":0,"smiles":0,"actives":0,"inactives":0,"inconclusive":0, "unspecified":0,"chem_probe":0}
        pathogen_dict[a]=summary_dict
    df_summary = pd.DataFrame.from_dict(
        pathogen_dict,
        orient="index",
        columns=["sids", "cids", "smiles","actives","inactives","inconclusive", "unspecified","chem_probe"]
        )
    df_summary.index.name = "aid"
    df_summary = df_summary.reset_index()
    df_summary.to_csv(os.path.join(outdir_clean, "summary.csv"), index=False)
    #overlap_csv = os.path.join(outdir_clean, "overlap.csv")
    #compare_aid_cid_overlap(outdir_clean, aids_with_data, overlap_csv)
    """
    cid2smi_df=unique_cids_smiles_counts(outdir_clean,aids_with_data)
    unique_cids_all[p.lower()] = [len(cid2smi_df), len(aids_with_data)]
    cid2smi_df.to_csv(os.path.join(datapath, "processed", "unique_cids", f"unique_cids_{p.lower()}.csv"), index=False)
    all_dfs.append(cid2smi_df)

df_cids = pd.DataFrame.from_dict(unique_cids_all,orient="index", columns = ["n_cid", "n_aid"])
df_cids.index.name = "pathogen"
df_cids = df_cids.reset_index()
df_cids.to_csv(os.path.join(datapath, "processed", "unique_cids", "summary.csv"), index=False)

df_all = pd.concat(all_dfs, ignore_index=True)
df_global = (
    df_all.groupby("cid", as_index=False)
          .agg({
              "smiles": "first",
              "inchikey": "first",
              "n_occurrences": "sum"
          })
)
df_global = df_global.sort_values("n_occurrences", ascending=False)
df_global.to_csv(os.path.join(datapath, "processed", "unique_cids", "unique_cids_all_pathogens.csv"), index=False)