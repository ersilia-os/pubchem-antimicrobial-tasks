"""
Step 09 — Filter annotated bioassays to a modelling-ready table.

Takes the organism-level annotated assays from step 08, drops metadata-only
columns not needed for modelling, and adds an active ratio column.

Inputs:
    data/processed/08_annotated_assays/summaries_organism.csv

Outputs:
    data/processed/09_bioassays_to_model/bioassays_to_model.csv
"""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent.parent
sys.path.append(str(ROOT / "src"))

MIN_CPDS_CONDITION_A = 1000
MIN_POSITIVES_CONDITION_A = 50
MIN_POSITIVES_CONDITION_B = 100

# AIDs manually identified as counter-screens or non-antimicrobial assays
_NON_ANTIMICROBIAL_AIDS = {
    2327,    # mammalian fibroblast toxicity counter-screen (calbicans)
    588517,  # compound fluorescence interference counter-screen (calbicans)
    588335,  # biochemical artifact counter-screen (mtuberculosis)
    527,     # quorum sensing / virulence pathway inhibition, not bacterial killing (saureus)
    1159583, # hypoxia-regulated fluorescent biosensor assay, not bacterial killing (mtuberculosis)
    488966,  # bacterial capsule biogenesis / virulence factor, not growth inhibition (ecoli)
    463173,  # teichoic acid synthesis / virulence factor, non-essential for in vitro viability (saureus)
}

_DROP_COLS = [
    "pathogen", "sids", "smiles", "inconclusive", "unspecified",
    "chem_probe", "chembl_id", "chembl_target_type",
    "pubchem_protacxns", "pubchem_geneids",
]

in_path = ROOT / "data" / "processed" / "08_annotated_assays" / "summaries_organism.csv"
out_dir = ROOT / "data" / "processed" / "09_bioassays_to_model"
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(in_path)
print(f"Loaded {len(df)} organism assays.")

df = df.drop(columns=[c for c in _DROP_COLS if c in df.columns])

total = df["actives"] + df["inactives"]
df["ratio"] = (df["actives"] / total.where(total > 0)).round(3)

cond_a = (df["cids"] >= MIN_CPDS_CONDITION_A) & (df["actives"] >= MIN_POSITIVES_CONDITION_A) & (df["ratio"] < 0.5)
cond_b = (df["actives"] >= MIN_POSITIVES_CONDITION_B) & (df["ratio"] >= 0.5)
df["label"] = None
df.loc[cond_a, "label"] = "A"
df.loc[cond_b, "label"] = "B"
print(f"  Label A: {cond_a.sum()} | Label B: {cond_b.sum()} | discarded: {(~cond_a & ~cond_b).sum()}")
df = df[df["label"].notna()].reset_index(drop=True)
print(f"  Retained {len(df)} assays.")

df["keep"] = ~df["aid"].isin(_NON_ANTIMICROBIAL_AIDS)
print(f"  keep=False: {(~df['keep']).sum()} non-antimicrobial assays flagged.")

_COL_ORDER = [
    "code", "aid", "cids",
    "actives", "inactives", "ratio", "label", "keep",
    "activity_type_pubchem",
    "target_type_chembl", "target_type_pubchem",
    "pubchem_name", "pubchem_description", "pubchem_readout_columns",
]
df = df[[c for c in _COL_ORDER if c in df.columns]]

out_path = out_dir / "bioassays_to_model.csv"
df.to_csv(out_path, index=False)
print(f"Saved {len(df)} rows → {out_path}")
