"""
Step 06 — Annotate selected bioassays with target type and compound counts.

For each selected AID produces a summary table with:
  - compound counts (cids, actives, inactives, …) from script 04 summaries
  - chembl_id if the assay was matched to ChEMBL with a compound mismatch
  - target_type:
      · ChEMBL assays → read target_type_curated_extra from the companion
        chembl-antimicrobial-tasks repo (output/{code}/18_assays_master.csv):
          SINGLE PROTEIN → single_protein
          ORGANISM       → organism
          DISCARDED      → discarded
      · Novel PubChem assays → geneids or protacxns populated → single_protein,
        otherwise → organism
  - protein_id: uniprot_accession from 18_assays_master.csv for ChEMBL assays,
    protacxns from PubChem bioassay summary for novel assays
  - label: A if actives / (actives + inactives) < 0.5, B if >= 0.5, None if
    denominator is zero

Inputs:
    data/processed/02_bioassays_to_keep/aids_{pathogen}.csv
    data/processed/02_bioassays_to_keep/chembl_assays_in_pubchem_{code}.csv
    data/processed/04_extracted_bioassays/{pathogen}/summary.csv
    data/config/bioassays_summary/PubChem_bioassay_{pathogen}.csv
    data/config/bioassays_summary/bioassays_manually_excluded.csv
    ../../chembl-antimicrobial-tasks/output/{code}/18_assays_master.csv

Outputs:
    output/06_annotate_selected_bioassays/06_summary.csv
    output/06_annotate_selected_bioassays/06_activity_outcomes.png
    output/06_annotate_selected_bioassays/06_target_types.png
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia

ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(ROOT / "src"))
from default import PATHOGEN_TO_CODE, pathogens  # noqa: E402

CHEMBL_ROOT = ROOT.parent / "chembl-antimicrobial-tasks"

_TARGET_TYPE_MAP = {
    "SINGLE PROTEIN": "single_protein",
    "ORGANISM":       "organism",
    "DISCARDED":      "discarded",
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_chembl_map() -> dict[tuple[str, int], str]:
    """Return {(pathogen, aid): chembl_id} for mismatched assays."""
    result: dict[tuple[str, int], str] = {}
    for pathogen in pathogens:
        code = PATHOGEN_TO_CODE[pathogen]
        path = ROOT / "data" / "processed" / "02_bioassays_to_keep" / f"chembl_assays_in_pubchem_{code}.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path, usecols=["AID", "Source ID"])
        for _, row in df.iterrows():
            result[(pathogen, int(row["AID"]))] = str(row["Source ID"])
    return result


def load_chembl_master(code: str) -> pd.DataFrame:
    """Load 18_assays_master.csv for a pathogen code, indexed by assay_id."""
    path = CHEMBL_ROOT / "output" / code / "18_assays_master.csv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, usecols=["assay_id", "target_type_curated_extra", "uniprot_accession"])
    df = df.set_index("assay_id")
    return df


def load_pubchem_info() -> dict[tuple[str, int], tuple[str | None, str | None]]:
    """Return {(pathogen, aid): (protacxns, geneids)} from PubChem bioassay summaries."""
    result: dict[tuple[str, int], tuple] = {}
    for pathogen in pathogens:
        path = ROOT / "data" / "config" / "bioassays_summary" / f"PubChem_bioassay_{pathogen}.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path, usecols=["aid", "protacxns", "geneids"], low_memory=False)
        for _, row in df.iterrows():
            prot = None if pd.isna(row["protacxns"]) else str(row["protacxns"])
            gene = None if pd.isna(row["geneids"]) else str(row["geneids"])
            result[(pathogen, int(row["aid"]))] = (prot, gene)
    return result


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def classify_novel(protacxns: str | None, geneids: str | None) -> str:
    if protacxns is not None or geneids is not None:
        return "single_protein"
    return "organism"


def compute_label(cids: int | None, actives: int | None, inactives: int | None) -> str:
    if cids is None or actives is None or inactives is None:
        return "discarded"
    denom = actives + inactives
    ratio = actives / denom if denom > 0 else None
    if cids >= 1000 and actives >= 50 and ratio is not None and ratio < 0.5:
        return "A"
    if actives >= 100 and ratio is not None and ratio >= 0.5:
        return "B"
    return "discarded"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

out_dir = ROOT / "output" / "06_annotate_selected_bioassays"
out_dir.mkdir(parents=True, exist_ok=True)

excluded_path = ROOT / "data" / "config" / "bioassays_summary" / "bioassays_manually_excluded.csv"
excluded: set[tuple[str, int]] = set()
if excluded_path.exists():
    excl_df = pd.read_csv(excluded_path)
    excl_df = excl_df[excl_df["discard"] == 1]
    for _, row in excl_df.iterrows():
        excluded.add((str(row["pathogen"]), int(row["aid"])))
    print(f"Loaded {len(excluded)} manually excluded AIDs.")

print("Loading ChEMBL ID map…")
chembl_map = load_chembl_map()

print("Loading PubChem bioassay info…")
pubchem_info = load_pubchem_info()

# Cache ChEMBL master tables per code
chembl_masters: dict[str, pd.DataFrame] = {}

rows = []

for pathogen in pathogens:
    code = PATHOGEN_TO_CODE[pathogen]

    aids_path = ROOT / "data" / "processed" / "02_bioassays_to_keep" / f"aids_{pathogen.lower()}.csv"
    if not aids_path.exists():
        print(f"  Skipping {pathogen}: aids file not found.")
        continue
    aids = pd.read_csv(aids_path)["aid"].tolist()

    summary_path = ROOT / "data" / "processed" / "04_extracted_bioassays" / pathogen.lower() / "summary.csv"
    summary_df = pd.read_csv(summary_path).set_index("aid") if summary_path.exists() else pd.DataFrame()

    if code not in chembl_masters:
        chembl_masters[code] = load_chembl_master(code)
    master = chembl_masters[code]

    print(f"  {pathogen}: {len(aids)} AIDs")

    for aid in aids:
        aid = int(aid)
        key = (pathogen, aid)
        if key in excluded:
            continue

        # Compound counts
        if aid in summary_df.index:
            s = summary_df.loc[aid]
            cids       = int(s.get("cids", 0))
            actives    = int(s.get("actives", 0))
            inactives  = int(s.get("inactives", 0))
            inconc     = int(s.get("inconclusive", 0))
            unspec     = int(s.get("unspecified", 0))
            chem_probe = int(s.get("chem_probe", 0))
        else:
            cids = actives = inactives = inconc = unspec = chem_probe = None

        # ChEMBL ID (only for mismatched assays)
        chembl_id = chembl_map.get(key)

        # target_type and protein_id
        if chembl_id and not master.empty and chembl_id in master.index:
            row_m = master.loc[chembl_id]
            raw_tt = row_m.get("target_type_curated_extra")
            target_type = _TARGET_TYPE_MAP.get(str(raw_tt), "other") if pd.notna(raw_tt) else "other"
            prot_raw = row_m.get("uniprot_accession")
            protein_id = None if pd.isna(prot_raw) or str(prot_raw).strip() == "" else str(prot_raw)
        else:
            prot, gene = pubchem_info.get(key, (None, None))
            target_type = classify_novel(prot, gene)
            protein_id = prot  # protacxns from PubChem

        label = compute_label(cids, actives, inactives)

        rows.append({
            "pathogen_code":  code,
            "aid":            aid,
            "chembl_id":      chembl_id,
            "target_type":    target_type,
            "protein_id":     protein_id,
            "cids":           cids,
            "actives":        actives,
            "inactives":      inactives,
            "inconclusive":   inconc,
            "unspecified":    unspec,
            "chemical_probe": chem_probe,
            "label":          label,
        })

out_df = pd.DataFrame(rows)
out_path = out_dir / "06_summary.csv"
out_df.to_csv(out_path, index=False)
print(f"\nSaved {len(out_df)} rows → {out_path}")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

stylia.set_format("print")
stylia.set_style("ersilia")

nc = stylia.NamedColors()
short_names = [p.replace("_", " ") for p in pathogens]

# --- Plot 1: activity outcomes stacked bar (from script 07) ---

ACT_CATS   = ["actives", "inactives", "inconclusive", "unspecified", "chem_probe"]
ACT_LABELS = ["Active", "Inactive", "Inconclusive", "Unspecified", "Chemical probe"]
ACT_COLORS = [nc.mint, nc.purple, nc.orange, nc.gray, nc.yellow]

act_rows = []
for p in pathogens:
    summary_path = ROOT / "data" / "processed" / "04_extracted_bioassays" / p.lower() / "summary.csv"
    if not summary_path.exists():
        act_rows.append({c: 0 for c in ACT_CATS})
        continue
    df = pd.read_csv(summary_path)
    act_rows.append({c: int(df[c].sum()) if c in df.columns else 0 for c in ACT_CATS})

act_agg = pd.DataFrame(act_rows)

fig, axs = stylia.create_figure(1, 1)
ax = axs.next()
x = np.arange(len(pathogens))
bottoms = np.zeros(len(pathogens))
for cat, label, color in zip(ACT_CATS, ACT_LABELS, ACT_COLORS):
    values = act_agg[cat].values.astype(float)
    ax.bar(x, values, bottom=bottoms, color=color, label=label)
    bottoms += values
ax.set_xticks(x)
ax.set_xticklabels(short_names, rotation=45, ha="right")
ax.legend()
stylia.label(ax, xlabel="", ylabel="Compound measurements", title="Activity outcome per pathogen")
stylia.save_figure(str(out_dir / "06_activity_outcomes.png"))
print(f"Saved: {out_dir / '06_activity_outcomes.png'}")

# --- Plot 2: single_protein vs organism per pathogen ---

tt_agg = (
    out_df[out_df["target_type"].isin(["single_protein", "organism"])]
    .groupby(["pathogen_code", "target_type"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=["organism", "single_protein"], fill_value=0)
)
# reorder to match pathogens list
code_order = [PATHOGEN_TO_CODE[p] for p in pathogens if PATHOGEN_TO_CODE[p] in tt_agg.index]
tt_agg = tt_agg.reindex(code_order)

fig, axs = stylia.create_figure(1, 1)
ax = axs.next()
x = np.arange(len(tt_agg))
width = 0.4
ax.bar(x - width / 2, tt_agg["organism"].values,      width, color=nc.blue,   label="Organism")
ax.bar(x + width / 2, tt_agg["single_protein"].values, width, color=nc.purple, label="Single protein")
ax.set_xticks(x)
ax.set_xticklabels(tt_agg.index, rotation=45, ha="right")
ax.legend()
stylia.label(ax, xlabel="", ylabel="Number of assays", title="Target type per pathogen")
stylia.save_figure(str(out_dir / "06_target_types.png"))
print(f"Saved: {out_dir / '06_target_types.png'}")
