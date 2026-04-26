"""
Plot per-pathogen activity breakdown across all selected bioassays.

Reads the summary.csv files produced by script 04 (one per pathogen) and
produces a stacked bar chart showing actives, inactives, inconclusive,
unspecified, and chemical probe counts.

Output: output/plots/07_summary_activity.png
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia

root = Path(__file__).resolve().parent
sys.path.append(str(root.parent / "src"))
from default import pathogens

datapath = root.parent / "data"
outpath = root.parent / "output"
plotpath = outpath / "plots"
plotpath.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load and aggregate activity counts per pathogen
# ---------------------------------------------------------------------------

CATEGORIES = ["actives", "inactives", "inconclusive", "unspecified", "chem_probe"]
LABELS = ["Active", "Inactive", "Inconclusive", "Unspecified", "Chemical probe"]

rows = []
for p in pathogens:
    summary_path = datapath / "processed" / "04_extracted_bioassays" / p.lower() / "summary.csv"
    if not summary_path.exists():
        rows.append({"pathogen": p, **{c: 0 for c in CATEGORIES}})
        continue
    df = pd.read_csv(summary_path)
    row = {"pathogen": p}
    for c in CATEGORIES:
        row[c] = int(df[c].sum()) if c in df.columns else 0
    rows.append(row)

agg = pd.DataFrame(rows).set_index("pathogen")

# Short labels for x-axis (replace underscores with spaces)
short_names = [p.replace("_", " ") for p in agg.index]

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

# Format: print | Style: ersilia — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("ersilia")

nc = stylia.NamedColors()
COLORS = [nc.mint, nc.purple, nc.orange, nc.gray, nc.yellow]


def plot_activity_stacked(ax, agg, short_names, categories, labels, colors):
    x = np.arange(len(agg))
    bottoms = np.zeros(len(agg))
    for cat, label, color in zip(categories, labels, colors):
        values = agg[cat].values.astype(float)
        ax.bar(x, values, bottom=bottoms, color=color, label=label)
        bottoms += values
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=45, ha="right")
    ax.legend()
    stylia.label(ax, xlabel="", ylabel="Compound measurements", title="Activity outcome per pathogen")


fig, axs = stylia.create_figure(1, 1)
plot_activity_stacked(axs.next(), agg, short_names, CATEGORIES, LABELS, COLORS)
stylia.save_figure(str(plotpath / "07_summary_activity.png"))
print(f"Saved: {plotpath / '07_summary_activity.png'}")
