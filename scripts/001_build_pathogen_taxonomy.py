import pandas as pd
from pathlib import Path
import json

# Define project directories
ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"

# List of pathogens
pathogens = [
    "Acinetobacter baumannii", "Candida albicans", "Campylobacter",
    "Escherichia coli", "Enterococcus faecium", "Enterobacter",
    "Helicobacter pylori", "Klebsiella pneumoniae",
    "Mycobacterium tuberculosis", "Neisseria gonorrhoeae",
    "Pseudomonas aeruginosa", "Plasmodium falciparum",
    "Staphylococcus aureus", "Schistosoma mansoni",
    "Streptococcus pneumoniae"
]

# Step 1: Load and merge all raw pathogen taxonomy files
def load_pathogen_taxids(pathogen):
    file = DATA_RAW / f"PubChem_taxonomy_text_{pathogen}.csv"
    df = pd.read_csv(file)
    df["Pathogen"] = pathogen
    return df[["Pathogen", "Taxonomy_ID", "Taxonomy_Name"]]

dfs = [load_pathogen_taxids(p) for p in pathogens]
merged = pd.concat(dfs, ignore_index=True).drop_duplicates()

# Save raw merged output
DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
merged.to_csv(DATA_PROCESSED / "00_pathogens_taxid.csv", index=False)

# Step 2: Clean unwanted taxonomies
# Manual cleaning dictionary
wrong_taxonomies = {
    "Acinetobacter baumannii": {"Acinetobacter calcoaceticus/baumannii complex"},
    "Candida albicans": {"Candida tropicalis"},
    "Campylobacter": {
        "Helicobacter pylori", "Helicobacter mustelae",
        "Aliarcobacter cryaerophilus", "Helicobacter cinaedi",
        "Helicobacter fennelliae", "Aliarcobacter butzleri",
        "Arcobacter nitrofigilis", "Helicobacter sp. CLO-3",
        "Firehammervirus CP220", "Firehammervirus CPt10",
        "Fletchervirus NCTC12673", "Fletchervirus CP81",
        "Fletchervirus CPX", "Firehammervirus CP21",
        "Fletchervirus CP30A", "Fletchervirus Los1",
    },
    "Escherichia coli": {"Tequintavirus AKFV33", "Enterobacteria phage CUS-3"},
    "Enterococcus faecium": {"Enterococcus casseliflavus"},
    "Enterobacter": {
        "Hafnia alvei", "Kosakonia radicincitans DSM 16656", "Kluyvera intermedia",
        "Cronobacter sakazakii", "Pluralibacter gergoviae", "Klebsiella aerogenes",
        "Pantoea agglomerans", "Lelliottia amnigena", "Franconibacter helveticus",
        "Kosakonia oryzae", "Kosakonia arachidis", "Kosakonia sacchari", 
        "Kosakonia cowanii", "Kosakonia oryzendophytica", "Franconibacter pulveris",
        "Hafnia phage Enc34", "Karamvirus pg7", "Slopekvirus eap3",
        "Escherichia phage IME11", "Phytobacter massiliensis",
        "Kosakonia radicincitans", "Pluralibacter pyrinus",
        "Pseudenterobacter timonensis", "Rahnella aquatilis",
        "Webervirus F20", "Kosakonia sacchari SP1",
        "Pluralibacter gergoviae ATCC 33028 = NBRC 105706",
        "Klebsiella aerogenes EA1509E", "Klebsiella aerogenes KCTC 2190",
        "Kosakonia oryziphila", "Karamvirus cc31", "Eclunavirus EcL1", "Eapunavirus Eap1"
    },
    "Helicobacter pylori": {"Helicobacter mustelae"},
    "Klebsiella pneumoniae": {"Klebsiella michiganensis KCTC 1686", "Klebsiella variicola subsp. tropica"},
    "Mycobacterium tuberculosis": {"Mycobacterium avium", "Corynebacterium pseudotuberculosis"},
    "Pseudomonas aeruginosa": {"Pseudomonas virus Yua"},
    "Staphylococcus aureus": {"Dubowvirus dv11"},
}

# Remove unwanted taxonomy matches
mask = pd.Series(True, index=merged.index)
for pathogen, bad_names in wrong_taxonomies.items():
    mask &= ~((merged["Pathogen"] == pathogen) & merged["Taxonomy_Name"].isin(bad_names))

cleaned = merged[mask].copy()

# Remove general phage/virus matches
cleaned = cleaned[~cleaned["Taxonomy_Name"].str.contains("phage|virus", case=False, na=False)]

# Save cleaned table
cleaned.to_csv(DATA_PROCESSED / "01_pathogens_taxid_cleaned.csv", index=False)

# Step 3: Export as dictionary
taxid_dict = (
    cleaned.groupby("Pathogen")["Taxonomy_ID"]
    .apply(list)
    .to_dict()
)

with open(DATA_PROCESSED / "02_pathogens_taxid_cleaned_dict.json", "w") as f:
    json.dump(taxid_dict, f, indent=2)

print("✅ Finished building pathogen taxonomy files.")