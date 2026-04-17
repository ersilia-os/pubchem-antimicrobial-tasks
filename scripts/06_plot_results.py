"""
Generate per-pathogen report figures and a summary table comparing
PubChem (novel) and ChEMBL compound datasets.

Requires:
  - data/processed/02_bioassays_to_keep/summary.csv  (script 02)
  - data/processed/04_unique_cids/unique_cids_{pathogen}.csv  (script 04)
  - ../../chembl-antimicrobial-tasks/output/{code}/20_all_smiles.csv
  - ../../chembl-antimicrobial-tasks/output/{code}/20_general_datasets.csv

Outputs:
  - output/plots/06_{pathogen}_report.png  — per-pathogen figure
  - output/06_summary_table.csv            — cross-pathogen summary

ChEMBL InChIKeys are cached in data/processed/06_chembl_inchikeys/ to
avoid recomputing them on subsequent runs.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia
from rdkit import Chem
from tqdm import tqdm

root = Path(__file__).resolve().parent
sys.path.append(str(root.parent / "src"))
from default import pathogens, PATHOGEN_TO_CODE

datapath = root.parent / "data"
outpath = root.parent / "output"
chembl_root = root.parent.parent / "chembl-antimicrobial-tasks"

plotpath = outpath / "plots"
cachepath = datapath / "processed" / "06_chembl_inchikeys"

plotpath.mkdir(parents=True, exist_ok=True)
cachepath.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def _smiles_to_inchikey(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None


def load_chembl_inchikeys(code: str) -> set:
    """Load (or compute and cache) InChIKeys for all ChEMBL SMILES."""
    cache_file = cachepath / f"{code}.txt"
    if cache_file.exists():
        with open(cache_file) as f:
            return {line.strip() for line in f if line.strip()}

    smiles_path = chembl_root / "output" / code / "20_all_smiles.csv"
    if not smiles_path.exists():
        return set()

    print(f"  Computing InChIKeys for {code} ChEMBL SMILES (cached after first run)...")
    df = pd.read_csv(smiles_path)
    iks = [
        ik for ik in tqdm(
            df["smiles"].apply(_smiles_to_inchikey), desc=f"  {code}", leave=False
        )
        if ik is not None
    ]
    with open(cache_file, "w") as f:
        f.write("\n".join(iks))
    return set(iks)


def load_pubchem_inchikeys(pathogen: str) -> set:
    path = datapath / "processed" / "04_unique_cids" / f"unique_cids_{pathogen.lower()}.csv"
    if not path.exists():
        return set()
    df = pd.read_csv(path, usecols=["inchikey"])
    return set(df["inchikey"].dropna())



def load_pubchem_summary() -> pd.DataFrame:
    path = datapath / "processed" / "02_bioassays_to_keep" / "summary.csv"
    return pd.read_csv(path).set_index("pathogen")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Format: print | Style: ersilia — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("ersilia")


def plot_compound_counts(ax, n_chembl, n_pubchem_all, n_pubchem_novel):
    nc = stylia.NamedColors()
    ax.bar(
        ["ChEMBL", "All\nPubChem", "Novel\nPubChem"],
        [n_chembl, n_pubchem_all, n_pubchem_novel],
        color=[nc.purple, nc.blue, nc.mint],
    )
    stylia.label(ax, title="Compound counts", xlabel="", ylabel="Count", abc="A")


def plot_selected_aids_pct(ax, pct_not_in_chembl, pct_mismatched):
    nc = stylia.NamedColors()
    ax.bar(
        ["Not in\nChEMBL", "Mismatch"],
        [pct_not_in_chembl, pct_mismatched],
        color=[nc.mint, nc.orange],
    )
    stylia.label(ax, title="Selected AIDs (% of total)", xlabel="", ylabel="%", abc="B")


def plot_pathogen_report(pathogen: str, code: str, row: pd.Series, outdir: Path):
    chembl_iks = load_chembl_inchikeys(code)
    pubchem_iks = load_pubchem_inchikeys(pathogen)

    n_chembl_compounds = len(chembl_iks)
    n_pubchem_compounds = len(pubchem_iks)
    n_pubchem_all = int(row.get("cpd_count_bioassay_aids", 0))
    n_bioassay_aids = int(row.get("bioassay_aids", 0))
    n_not_in_chembl = int(row.get("aids_not_in_chembl", 0))
    n_not_in_chembl_min = int(row.get("aids_not_in_chembl_keep", 0))
    n_mismatched = int(row.get("aids_mismatched", 0))

    pct_not_in_chembl = round(n_not_in_chembl / n_bioassay_aids * 100, 1) if n_bioassay_aids > 0 else 0.0
    pct_not_in_chembl_min = round(n_not_in_chembl_min / n_bioassay_aids * 100, 1) if n_bioassay_aids > 0 else 0.0
    pct_mismatched = round(n_mismatched / n_bioassay_aids * 100, 1) if n_bioassay_aids > 0 else 0.0

    _, axs = stylia.create_figure(1, 2, width=0.5)
    plot_compound_counts(axs.next(), n_chembl_compounds, n_pubchem_all, n_pubchem_compounds)
    plot_selected_aids_pct(axs.next(), pct_not_in_chembl_min, pct_mismatched)
    stylia.save_figure(str(outdir / f"06_{pathogen.lower()}_report.png"))

    return {
        "pathogen": pathogen,
        "chembl_code": code,
        "n_chembl_compounds": n_chembl_compounds,
        "n_pubchem_all": n_pubchem_all,
        "n_pubchem_novel": n_pubchem_compounds,
        "pct_pubchem_novel": (
            round(n_pubchem_compounds / n_pubchem_all * 100, 1)
            if n_pubchem_all > 0 else 0.0
        ),
        "n_bioassay_aids": n_bioassay_aids,
        "n_not_in_chembl": n_not_in_chembl,
        "n_not_in_chembl_min": n_not_in_chembl_min,
        "n_mismatched": n_mismatched,
        "pct_not_in_chembl": pct_not_in_chembl,
        "pct_not_in_chembl_min": pct_not_in_chembl_min,
        "pct_mismatched": pct_mismatched,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

pubchem_summary = load_pubchem_summary()
rows = []

for pathogen in pathogens:
    code = PATHOGEN_TO_CODE[pathogen]
    print(f"\n{pathogen} ({code})")
    row = pubchem_summary.loc[pathogen] if pathogen in pubchem_summary.index else pd.Series(dtype=float)
    result = plot_pathogen_report(pathogen, code, row, plotpath)
    rows.append(result)

summary_df = pd.DataFrame(rows)
summary_path = outpath / "06_summary_table.csv"
summary_df.to_csv(summary_path, index=False)
print(f"\nSummary table saved: {summary_path}")
