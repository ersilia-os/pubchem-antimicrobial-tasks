import os
import sys
import pandas as pd
import stylia as st
import matplotlib.pyplot as plt

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import pathogens

datapath = os.path.join(root, "..", "data")
configpath = os.path.join(datapath, "config")

pathogen_aids = {}

for p in pathogens:
    df = pd.read_csv(os.path.join(configpath, "taxonomy_raw", f"PubChem_taxonomy_{p}.csv"))
    print(p, len(df))
    keep = []
    for t in df["Taxonomy_Name"].to_list():
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
    df.to_csv(os.path.join(datapath, "processed", "taxonomy_processed", f"taxonomy_{p}.csv"), index=False)
print("Taxonomy preprocessing completed.")

for p in pathogens:
    tax = pd.read_csv(os.path.join(datapath, "processed", "taxonomy_processed", f"taxonomy_{p}.csv"))
    bioassays = pd.read_csv(os.path.join(configpath, "bioassays_summary", f"PubChem_bioassay_{p}.csv"), low_memory=False)
    print(p, len(bioassays))
    bioassays1 = bioassays[bioassays["targettaxid"].isin(tax["Taxonomy_ID"].tolist())]
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
    no_tax = no_tax[no_tax["cnt"]>100]
    print(f"No taxid info with >100 mols: {len(no_tax)} assays need manual checking")
    if len(no_tax)>0:
        no_tax.to_csv(os.path.join(datapath, "processed", "bioassays_summary", f"bioassays_{p}_manual_check.csv"), index=False)
    
    if os.path.exists(os.path.join(configpath, "bioassays_selected_manually", f"{p}.csv")):
        manual = pd.read_csv(os.path.join(configpath, "bioassays_selected_manually", f"{p}.csv"))
        print(f"Adding {len(manual)} manually selected assays for {p}")
        bioassays3 = bioassays[bioassays["aid"].isin(manual["aid"].tolist())]
        print("Bioassays3:", len(bioassays3))
        bioassays_final = pd.concat([bioassays1, bioassays2, bioassays3], axis=0).drop_duplicates()
    else:
        bioassays_final = pd.concat([bioassays1, bioassays2], axis=0).drop_duplicates()
    print(p, len(bioassays_final))
    bioassays_final_100 = bioassays_final[bioassays_final["cnt"]>100]
    pathogen_aids[p]=  [len(bioassays_final), bioassays_final["cnt"].sum(),len(bioassays_final_100), bioassays_final_100["cnt"].sum()]
    bioassays_final.to_csv(os.path.join(datapath, "processed", "bioassays_summary", f"bioassays_{p}.csv"), index=False)


df_plot = pd.DataFrame.from_dict(
    pathogen_aids,
    orient="index",
    columns=["n_bioassays", "total_cnt", "n_bioassays_100", "total_cnt_100"],
)
df_plot.index.name = "pathogen"
df_plot = df_plot.reset_index()
df_plot.to_csv(os.path.join(datapath, "processed", "bioassays_summary","summary.csv"), index=False)

fig, axs = st.create_figure(2,2)
ax = axs.next()
ax.barh(df_plot["pathogen"], df_plot["n_bioassays"])
st.label(ax, title = "AIDs Number", xlabel = "", ylabel = "")
ax = axs.next()
ax.barh(df_plot["pathogen"], df_plot["total_cnt"])
ax.tick_params(axis="y", labelleft=False)
st.label(ax, title = "Total compound count", xlabel = "", ylabel = "")
ax = axs.next()
ax.barh(df_plot["pathogen"], df_plot["n_bioassays_100"])
st.label(ax, title = "AIDs Number (100 max)", xlabel = "", ylabel = "")
ax = axs.next()
ax.barh(df_plot["pathogen"], df_plot["total_cnt_100"])
ax.tick_params(axis="y", labelleft=False)
st.label(ax, title = "Total compound count (100 max)", xlabel = "", ylabel = "")
plt.tight_layout()
st.save_figure("../output/plots/00_total_aids.png")