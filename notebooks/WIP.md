# Work in progress 🚧
## 1. How to download BioAssays from PubChem?

A. Check how PubChem works → [PubChem_info.md](pubchem-antimicrobial-tasks/notebooks/PubChem_info.md)

B. Look for taxonomy ID for the 15 selected pathogen. Make a Python script that:

- takes a list of pathogen species
- queries PubChem for all taxonomies
- retrieves all AIDs
- downloads summaries or data
- stores everything in CSV / JSON

Currently, the following list of pathogens is processed:
````
["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter", "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa", "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
```