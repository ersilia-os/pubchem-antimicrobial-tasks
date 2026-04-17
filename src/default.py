PATHOGEN_TO_CODE = {
    "Acinetobacter_baumannii": "abaumannii",
    "Candida_albicans": "calbicans",
    "Campylobacter": "campylobacter",
    "Escherichia_coli": "ecoli",
    "Enterococcus_faecium": "efaecium",
    "Enterobacter": "enterobacter",
    "Helicobacter_pylori": "hpylori",
    "Klebsiella_pneumoniae": "kpneumoniae",
    "Mycobacterium_tuberculosis": "mtuberculosis",
    "Neisseria_gonorrhoeae": "ngonorrhoeae",
    "Pseudomonas_aeruginosa": "paeruginosa",
    "Plasmodium_falciparum": "pfalciparum",
    "Staphylococcus_aureus": "saureus",
    "Schistosoma_mansoni": "smansoni",
    "Streptococcus_pneumoniae": "spneumoniae",
}

pathogens = [
    "Acinetobacter_baumannii", "Candida_albicans", "Campylobacter",
    "Escherichia_coli", "Enterococcus_faecium", "Enterobacter",
    "Helicobacter_pylori", "Klebsiella_pneumoniae",
    "Mycobacterium_tuberculosis", "Neisseria_gonorrhoeae",
    "Pseudomonas_aeruginosa", "Plasmodium_falciparum",
    "Staphylococcus_aureus", "Schistosoma_mansoni",
    "Streptococcus_pneumoniae"
]

# Minimum number of compounds an assay must contain to be included
MIN_COMPOUNDS = 100

# Minimum fractional surplus of PubChem compounds over ChEMBL to retain a
# matched assay (i.e. PubChem has at least 10% more compounds than ChEMBL)
CHEMBL_MISMATCH_THRESHOLD = 0.1