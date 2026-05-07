import os
import sys
import pandas as pd
import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import pathogens, MIN_COMPOUNDS

datapath = os.path.join(root, "..", "data")
configpath = os.path.join(datapath, "config")
taxpath = os.path.join(datapath, "processed", "00_taxonomy_processed")
bioassaypath = os.path.join(datapath, "processed", "00_bioassays_summary")

plotpath = os.path.join(root, "..", "output", "00_preprocess_bioassays")

os.makedirs(taxpath, exist_ok=True)
os.makedirs(bioassaypath, exist_ok=True)
os.makedirs(plotpath, exist_ok=True)

pathogen_aids = {}

for p in pathogens:
    df = pd.read_csv(os.path.join(configpath, "taxonomy_raw", f"PubChem_taxonomy_{p}.csv"))
    print(p, len(df))
    keep = []
    for t in df["Taxonomy_Name"].to_list():
        if pd.isna(t):
            raise ValueError(f"NaN found in Taxonomy_Name for {p} — check the raw taxonomy file.")
        if "phage" in t.lower() or "virus" in t.lower():
            continue
        if "_" in p:
            species = p.split("_")[1]
            if species.lower() in t.lower():
                keep.append(t)
        else:
            if p.lower() in t.lower():
                keep.append(t)
    df = df[df["Taxonomy_Name"].isin(keep)]
    print(p, len(df))
    df.to_csv(os.path.join(taxpath, f"taxonomy_{p}.csv"), index=False)
print("Taxonomy preprocessing completed.")

for p in pathogens:
    tax = pd.read_csv(os.path.join(taxpath, f"taxonomy_{p}.csv"))
    bioassays = pd.read_csv(os.path.join(configpath, "bioassays_summary", f"PubChem_bioassay_{p}.csv"), low_memory=False)
    print(p, len(bioassays))
    valid_target_ids = set(tax["Taxonomy_ID"].astype(float))
    bioassays1 = bioassays[bioassays["targettaxid"].isin(valid_target_ids)]
    print("Bioassays1:", len(bioassays1))

    # For rows without targettaxid, check the pipe-separated `taxids` field.
    valid_tax_ids = set(tax["Taxonomy_ID"].astype(str).tolist())
    def taxids_match(cell):
        if pd.isna(cell):
            return False
        for part in str(cell).split("|"):
            if part.strip() in valid_tax_ids:
                return True
        return False
    has_no_target = bioassays["targettaxid"].isna()
    bioassays2 = bioassays[has_no_target & bioassays["taxids"].apply(taxids_match)]
    print("Bioassays2:", len(bioassays2))

    # for rows without targettax id nor taxid, keep them for manual checking
    no_tax = bioassays[has_no_target & bioassays["taxids"].isna()]
    no_tax = no_tax[no_tax["cnt"] >= MIN_COMPOUNDS]
    print(f"No taxid info with >100 mols: {len(no_tax)} assays need manual checking")
    if len(no_tax)>0:
        no_tax.to_csv(os.path.join(bioassaypath, f"bioassays_{p.lower()}_manual_check.csv"), index=False)
    if os.path.exists(os.path.join(configpath, "bioassays_selected_manually", f"{p}.csv")):
        manual = pd.read_csv(os.path.join(configpath, "bioassays_selected_manually", f"{p}.csv"))
        print(f"Adding {len(manual)} manually selected assays for {p}")
        bioassays3 = bioassays[bioassays["aid"].isin(manual["aid"].tolist())]
        print("Bioassays3:", len(bioassays3))
        bioassays_final = pd.concat([bioassays1, bioassays2, bioassays3], axis=0).drop_duplicates()
    else:
        bioassays_final = pd.concat([bioassays1, bioassays2], axis=0).drop_duplicates()
    print(p, len(bioassays_final))
    bioassays_final_100 = bioassays_final[bioassays_final["cnt"] >= MIN_COMPOUNDS]
    pathogen_aids[p]=  [len(bioassays_final), bioassays_final["cnt"].sum(),len(bioassays_final_100), bioassays_final_100["cnt"].sum()]
    bioassays_final.to_csv(os.path.join(bioassaypath, f"bioassays_{p.lower()}.csv"), index=False)


df_plot = pd.DataFrame.from_dict(
    pathogen_aids,
    orient="index",
    columns=["n_bioassays", "total_cnt", "n_bioassays_100", "total_cnt_100"],
)
df_plot.index.name = "pathogen"
df_plot = df_plot.reset_index()
df_plot.to_csv(os.path.join(bioassaypath, "summary.csv"), index=False)


# --- Plots ---
# Format: slide | Style: ersilia — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("ersilia")


def plot_bioassay_summary(df, outdir):
    fig, axs = stylia.create_figure(2, 2, width=0.5)
    nc = stylia.NamedColors()
    ax = axs.next()
    ax.barh(df["pathogen"], df["n_bioassays"], color=nc.purple)
    stylia.label(ax, title="AIDs number", xlabel="", ylabel="")

    ax = axs.next()
    ax.barh(df["pathogen"], df["total_cnt"], color=nc.pink)
    ax.tick_params(axis="y", labelleft=False)
    stylia.label(ax, title="Total compound count", xlabel="", ylabel="")

    ax = axs.next()
    ax.barh(df["pathogen"], df["n_bioassays_100"], color=nc.orange)
    stylia.label(ax, title=f"AIDs number (≥{MIN_COMPOUNDS} cpds)", xlabel="", ylabel="")

    ax = axs.next()
    ax.barh(df["pathogen"], df["total_cnt_100"], color=nc.yellow)
    ax.tick_params(axis="y", labelleft=False)
    stylia.label(ax, title=f"Total compound count (≥{MIN_COMPOUNDS} cpds)", xlabel="", ylabel="")

    stylia.save_figure(os.path.join(outdir, "00_total_aids.png"))


plot_bioassay_summary(df_plot, plotpath)