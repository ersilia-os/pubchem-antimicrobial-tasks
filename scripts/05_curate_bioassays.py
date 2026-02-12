import os
import sys
import pandas as pd

import zipfile
import gzip
import shutil

root = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(root, "..", "src"))
from default import pathogens

datapath = os.path.join(root, "..", "data")
outpath = os.path.join(root, "..", "output")

datadir = os.path.join(datapath, "processed", "extracted_bioassays")
destdir = os.path.join(outpath, "results")

for p in pathogens:
    aids = pd.read_csv(os.path.join(datapath,"processed","bioassays_to_keep", f"aids_{p.lower()}.csv"))["aid"].tolist()
    pathogen_dir = os.path.join(destdir, p.lower())
    os.makedirs(pathogen_dir, exist_ok=True)
    print(p)
    for a in aids:
        src_csv = os.path.join(datadir, p.lower(), f"{a}.csv")
        src_meta = os.path.join(datadir, p.lower(), f"{a}_meta.csv")
        if os.path.exists(src_csv):
            shutil.copy(src_csv, pathogen_dir)
        else:
            print(f"{a} assay does not exist")
        if os.path.exists(src_meta):
            shutil.copy(src_meta, pathogen_dir)
        else:
            print(f"{a} meta assay does not exist")



