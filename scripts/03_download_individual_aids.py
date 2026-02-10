import sys
import os
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(root, "..", "src"))
from pubchem_aid_parser import PubChemBioAssayRecord
from json2df_pubchem import PubChemBioAssayJsonConverter
from default import pathogens

datapath = os.path.join("..", "data")
outpath = os.path.join("..", "output")

tmp_dir = os.path.join("..", "tmp")
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

# This method only works with assays of less than 10000 SIDS: <Details>Assay record retrieval is limited to 10000 SIDs</Details>

p = "campylobacter"
aid = 743156

assert p in pathogens

if not os.path.exists(os.path.join(outpath,"results", p)):
    os.mkdir(os.path.join(outpath,"results", p))

print(aid)
c = PubChemBioAssayRecord(aid)
c.save_json(tmp_dir)
c = PubChemBioAssayJsonConverter(tmp_dir, f"PUBCHEM{aid}.json")
df = c.get_all_results()
c.save_df(df, os.path.join(outpath, "results", f"{p.lower()}"))
c.get_description(os.path.join(outpath, "results", f"{p.lower()}"))
