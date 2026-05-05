"""
Step 08 — Annotate selected bioassays with target and activity type information.

Reads the selected AIDs from previous steps and annotates each assay using two
complementary signals: ChEMBL curated metadata (for assays with a ChEMBL entry)
and PubChem bioassay summary data (for all assays).

Assay-type classification (single_protein / organism / other):
  - ChEMBL: derived from target_type in target_dictionary joined via assays.tid.
      SINGLE PROTEIN → single_protein
      ORGANISM / NON-MOLECULAR → organism
      everything else → other
      no tid / no match → None (fall back to PubChem signal)
  - PubChem: derived from presence of protacxns or geneids fields.

Activity-type classification (from _meta.csv column headers, organism assays only):
  ic50, ec50, ki, kd, mic, inhibition, fluorescence, luminescence, absorbance,
  doseresponse, other.

Inputs:
    data/processed/04_extracted_bioassays/{pathogen}/summary.csv
    data/processed/02_bioassays_to_keep/chembl_assays_in_pubchem_{code}.csv
    data/config/chembl_mappings/assays.csv
    data/config/target_dictionary.csv
    data/config/bioassays_summary/PubChem_bioassay_{pathogen}.csv
    output/results/{pathogen}/{aid}_meta.csv

Outputs:
    data/processed/08_annotated_assays/summaries.csv          (all AIDs, annotated)
    data/processed/08_annotated_assays/summaries_organism.csv (organism-level AIDs only)
"""

import csv
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.append(str(ROOT / "src"))
from default import PATHOGEN_TO_CODE, pathogens  # noqa: E402


_PUBCHEM_STANDARD_COLS = {
    "PUBCHEM_RESULT_TAG", "PUBCHEM_SID", "PUBCHEM_CID",
    "PUBCHEM_EXT_DATASOURCE_SMILES", "PUBCHEM_ACTIVITY_OUTCOME",
    "PUBCHEM_ACTIVITY_SCORE", "PUBCHEM_ACTIVITY_URL", "PUBCHEM_ASSAYDATA_COMMENT",
}

# Each entry: (category, [(pattern, use_word_boundary), ...])
# Word-boundary patterns are precompiled; plain strings are matched as substrings.
_ACTIVITY_KEYWORDS: list[tuple[str, list[tuple[re.Pattern | str, bool]]]] = [
    ("ic50",         [("ic50", False), ("ic90", False)]),
    ("ec50",         [("ec50", False), ("ac50", False), ("pac50", False), ("pic50", False)]),
    ("ki",           [(re.compile(r"\bki\b", re.IGNORECASE), True)]),
    ("kd",           [(re.compile(r"\bkd\b", re.IGNORECASE), True)]),
    ("mic",          [(re.compile(r"\bmic\b", re.IGNORECASE), True)]),
    ("inhibition",   [("inhibition", False), ("pctinhib", False), ("pct inh", False),
                      ("% inhib", False), ("% activ", False), ("growth inhib", False),
                      ("percent_response", False), ("percentresponse", False),
                      ("% survival", False), ("inhib_", False)]),
    ("fluorescence", [("fluorescence", False)]),
    ("luminescence", [("luminescence", False)]),
    ("absorbance",   [("absorbance", False)]),
    ("doseresponse", [("hill slope", False), ("hill_slope", False),
                      ("dose-response", False), ("dose_response", False),
                      ("absac", False), ("potency", False)]),
]

_TARGET_TYPE_MAP = {
    "SINGLE PROTEIN": "single_protein",
    "ORGANISM": "organism",
    "NON-MOLECULAR": "organism",
}


def _classify_target_type(target_type) -> str | None:
    if pd.isna(target_type):
        return None
    return _TARGET_TYPE_MAP.get(target_type, "other")


def _classify_pubchem(protacxns, geneids) -> str:
    if not pd.isna(protacxns) or not pd.isna(geneids):
        return "single_protein"
    return "organism"


def _classify_activity(col_names: list[str]) -> str:
    joined = " | ".join(col_names).lower()
    for category, keywords in _ACTIVITY_KEYWORDS:
        for kw, word_boundary in keywords:
            if word_boundary:
                if any(kw.search(col) for col in col_names):
                    return category
            else:
                if kw in joined:
                    return category
    return "other"


def build_summaries() -> pd.DataFrame:
    dfs = []
    for pathogen in pathogens:
        code = PATHOGEN_TO_CODE[pathogen]
        path = ROOT / "data" / "processed" / "04_extracted_bioassays" / pathogen.lower() / "summary.csv"
        if not path.exists():
            print(f"Skipping {pathogen}: summary.csv not found.")
            continue
        df = pd.read_csv(path)
        df.insert(0, "pathogen", pathogen)
        df.insert(1, "code", code)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def annotate_assay_types(df: pd.DataFrame) -> pd.DataFrame:
    # --- build (pathogen, aid) → chembl_id map ---
    aid_to_chembl: dict[tuple[str, int], str] = {}
    for pathogen in pathogens:
        code = PATHOGEN_TO_CODE[pathogen]
        local = ROOT / "data" / "processed" / "02_bioassays_to_keep" / f"chembl_assays_in_pubchem_{code}.csv"
        if not local.exists():
            continue
        mapping = pd.read_csv(local, usecols=["AID", "Source ID"])
        for _, row in mapping.iterrows():
            aid_to_chembl[(pathogen, int(row["AID"]))] = row["Source ID"]

    # --- load target_dictionary: tid → target_type ---
    tid_to_target_type: dict[int, str] = {}
    target_dict_path = ROOT / "data" / "config" / "target_dictionary.csv"
    if target_dict_path.exists():
        print("Loading target_dictionary.csv…")
        for _, row in pd.read_csv(target_dict_path, usecols=["tid", "target_type"]).iterrows():
            tid_to_target_type[int(row["tid"])] = row["target_type"]

    # --- load assays.csv filtered to relevant ChEMBL IDs ---
    relevant_chembl_ids = set(aid_to_chembl.values())
    chembl_info: dict[str, str | None] = {}
    assays_path = ROOT / "data" / "config" / "chembl_mappings" / "assays.csv"
    if relevant_chembl_ids and assays_path.exists():
        print(f"Loading assays.csv (filtering to {len(relevant_chembl_ids)} ChEMBL IDs)…")
        cols = ["chembl_id", "tid"]
        for chunk in pd.read_csv(assays_path, usecols=cols, chunksize=100_000):
            for _, row in chunk[chunk["chembl_id"].isin(relevant_chembl_ids)].iterrows():
                tid = int(row["tid"]) if not pd.isna(row["tid"]) else None
                target_type = tid_to_target_type.get(tid) if tid is not None else None
                chembl_info[row["chembl_id"]] = target_type

    # --- load PubChem bioassay summaries ---
    pubchem_target: dict[tuple[str, int], tuple] = {}
    for pathogen in pathogens:
        local = ROOT / "data" / "config" / "bioassays_summary" / f"PubChem_bioassay_{pathogen}.csv"
        if not local.exists():
            continue
        for _, row in pd.read_csv(local, usecols=["aid", "aidname", "aiddesc", "protacxns", "geneids"]).iterrows():
            pubchem_target[(pathogen, int(row["aid"]))] = (row["protacxns"], row["geneids"], row["aidname"], row["aiddesc"])

    # --- annotate ---
    chembl_id_col, chembl_target_type_col = [], []
    target_type_chembl, target_type_pubchem = [], []
    pubchem_name_col, pubchem_desc_col = [], []
    pubchem_protacxns_col, pubchem_geneids_col = [], []

    for _, row in df.iterrows():
        key = (row["pathogen"], int(row["aid"]))

        chembl_id = aid_to_chembl.get(key)
        if chembl_id and chembl_id in chembl_info:
            target_type = chembl_info[chembl_id]
            chembl_id_col.append(chembl_id)
            chembl_target_type_col.append(target_type)
            target_type_chembl.append(_classify_target_type(target_type))
        else:
            for lst in (chembl_id_col, chembl_target_type_col, target_type_chembl):
                lst.append(None)

        if key in pubchem_target:
            prot, gene, name, desc = pubchem_target[key]
            pubchem_protacxns_col.append(prot)
            pubchem_geneids_col.append(gene)
            pubchem_name_col.append(name)
            pubchem_desc_col.append(desc)
            target_type_pubchem.append(_classify_pubchem(prot, gene))
        else:
            pubchem_protacxns_col.append(None)
            pubchem_geneids_col.append(None)
            pubchem_name_col.append(None)
            pubchem_desc_col.append(None)
            target_type_pubchem.append(None)

    df = df.copy()
    df["target_type_chembl"] = target_type_chembl
    df["target_type_pubchem"] = target_type_pubchem
    df["chembl_id"] = chembl_id_col
    df["chembl_target_type"] = chembl_target_type_col
    df["pubchem_name"] = pubchem_name_col
    df["pubchem_description"] = pubchem_desc_col
    df["pubchem_protacxns"] = pubchem_protacxns_col
    df["pubchem_geneids"] = pubchem_geneids_col

    chembl_n = sum(1 for v in target_type_chembl if v is not None)
    print(f"  ChEMBL annotation: {chembl_n}/{len(df)} AIDs")
    print(f"  PubChem annotation: {sum(1 for v in target_type_pubchem if v is not None)}/{len(df)} AIDs")
    return df


def filter_organism_assays(df: pd.DataFrame) -> pd.DataFrame:
    mask = (df["target_type_pubchem"] == "organism") & (
        df["target_type_chembl"].isna() | (df["target_type_chembl"] == "organism")
    )
    organism_df = df[mask].copy()
    print(f"  Organism assays: {len(organism_df)}/{len(df)}")
    return organism_df


def annotate_activity_types(df: pd.DataFrame) -> pd.DataFrame:
    readout_cols_col, activity_type_pubchem_col = [], []

    for _, row in df.iterrows():
        pathogen = row["pathogen"]
        aid = int(row["aid"])
        meta_path = ROOT / "output" / "results" / pathogen / f"{aid}_meta.csv"

        if not meta_path.exists():
            readout_cols_col.append(None)
            activity_type_pubchem_col.append(None)
            continue

        with open(meta_path) as fh:
            header = next(csv.reader(fh))

        assay_cols = [c for c in header if c not in _PUBCHEM_STANDARD_COLS]
        readout_cols_col.append("|".join(assay_cols) if assay_cols else None)
        activity_type_pubchem_col.append(_classify_activity(assay_cols) if assay_cols else "other")

    df = df.copy()
    df["pubchem_readout_columns"] = readout_cols_col
    df["activity_type_pubchem"] = activity_type_pubchem_col

    classified = sum(1 for v in activity_type_pubchem_col if v is not None)
    print(f"  Activity type (PubChem): {classified}/{len(df)} AIDs classified")
    return df


out_dir = ROOT / "data" / "processed" / "08_annotated_assays"
out_dir.mkdir(parents=True, exist_ok=True)

print("Building summaries…")
summaries = build_summaries()
summaries_path = out_dir / "summaries.csv"

print("Annotating assay types…")
summaries = annotate_assay_types(summaries)
summaries.to_csv(summaries_path, index=False)
print(f"  Saved {len(summaries)} rows → {summaries_path}")

print("Filtering organism assays…")
organism = filter_organism_assays(summaries)

print("Annotating activity types…")
organism = annotate_activity_types(organism)
organism_path = out_dir / "summaries_organism.csv"
organism.to_csv(organism_path, index=False)
print(f"  Saved {len(organism)} rows → {organism_path}")
