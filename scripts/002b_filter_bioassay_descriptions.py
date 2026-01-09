#!/usr/bin/env python3

import gzip
import shutil
import zipfile
import time
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
from tqdm.auto import tqdm
from bs4 import BeautifulSoup
import json


# ===================================================================
# ------------------------ XML PARSING ------------------------------
# ===================================================================

def extract_info_from_descr(xml_path):
    if str(xml_path).endswith(".gz"):
        with gzip.open(xml_path, "rb") as f:
            xml_content = f.read()
    else:
        with open(xml_path, "rb") as f:
            xml_content = f.read()

    soup = BeautifulSoup(xml_content, "lxml-xml")

    taxids = set()
    for node in soup.find_all("PC-AssayTarget_tax-id"):
        txt = node.text.strip()
        if txt.isdigit():
            taxids.add(txt)
    for node in soup.find_all("PC-XRefData_taxonomy"):
        txt = node.text.strip()
        if txt.isdigit():
            taxids.add(txt)
    for node in soup.find_all("PC-AssayResultType_tax-id"):
        txt = node.text.strip()
        if txt.isdigit():
            taxids.add(txt)

    assay_organisms = []
    for node in soup.find_all("PC-AnnotatedXRef_comment"):
        txt = node.text.strip()
        if txt:
            assay_organisms.append(txt)
    assay_organisms = [x for x in assay_organisms if len(x) > 3]

    chembl_id = None
    for db in soup.find_all("PC-DBTracking"):
        name = db.find("PC-DBTracking_name")
        if name and name.text.strip() == "ChEMBL":
            obj = db.find("Object-id_str")
            if obj:
                chembl_id = obj.text.strip()
                break

    return {
        "TaxIDs": sorted(taxids),
        "AssayOrganism": assay_organisms,
        "ChEMBL": chembl_id,
    }


def detect_matching_pathogens(info, taxid_to_pathogen, pathogens):
    matched = set()

    for tid in info["TaxIDs"]:
        if tid in taxid_to_pathogen:
            matched.add(taxid_to_pathogen[tid])

    for org in info["AssayOrganism"]:
        org_low = org.lower()
        for pathogen in pathogens:
            if pathogen.lower() in org_low:
                matched.add(pathogen)

    return sorted(matched)


def worker_process_xml(args):
    xml_path_str, taxid_to_pathogen, pathogens = args
    xml_path = Path(xml_path_str)

    try:
        info = extract_info_from_descr(xml_path)
        pathogens_hit = detect_matching_pathogens(info, taxid_to_pathogen, pathogens)
        if not pathogens_hit:
            return None

        aid = int(xml_path.stem.split(".")[0])
        return {
            "xml_path": xml_path_str,
            "AID": aid,
            "Pathogens": pathogens_hit,
            "ChEMBL": info["ChEMBL"]
        }
    except Exception:
        return None


def process_one_zip(zip_file, DESC_DIR, KEEP_DESC_V3, N_WORKERS,
                    taxid_to_pathogen, pathogens):

    zip_chunk = zip_file.stem
    temp = DESC_DIR / f"{zip_chunk}_tmp"
    temp.mkdir(exist_ok=True)

    with zipfile.ZipFile(zip_file, "r") as zf:
        zf.extractall(temp)

    xml_files = [str(p) for p in temp.rglob("*.xml")] + \
                [str(p) for p in temp.rglob("*.xml.gz")]

    tasks = [(fp, taxid_to_pathogen, pathogens) for fp in xml_files]
    records = []

    with Pool(N_WORKERS) as pool:
        for result in pool.imap_unordered(worker_process_xml, tasks, chunksize=20):
            if result is None:
                continue

            aid = result["AID"]
            xml_path = Path(result["xml_path"])
            dst = KEEP_DESC_V3 / f"{aid}.xml"

            if xml_path.suffix == ".gz":
                with gzip.open(xml_path, "rb") as fi, open(dst, "wb") as fo:
                    fo.write(fi.read())
            else:
                shutil.copy2(xml_path, dst)

            records.append({
                "AID": aid,
                "Pathogen": ", ".join(result["Pathogens"]),
                "ChEMBLid": result["ChEMBL"],
                "ZipFolder": zip_chunk
            })

    shutil.rmtree(temp)
    return records


def main():
    PROJECT_ROOT = Path(__file__).resolve().parents[1]

    DATA_RAW = PROJECT_ROOT / "data" / "raw"
    DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

    PUBCHEM_DIR = DATA_RAW / "pubchem_bioassays"
    DESC_DIR = PUBCHEM_DIR / "Description"

    # NEW: Filtered output for v3
    KEEP_V3 = DATA_RAW / "filtered_assays_v3"
    KEEP_DESC_V3 = KEEP_V3 / "Description"
    KEEP_DESC_V3.mkdir(parents=True, exist_ok=True)

    OUT_CSV = DATA_PROCESSED / "004b_filtered_aid_summary.csv"
    ZIP_LOG = DATA_PROCESSED / "processed_zips_v3.txt"

    if ZIP_LOG.exists():
        processed_chunks = set(ZIP_LOG.read_text().splitlines())
    else:
        processed_chunks = set()

    # NEW: Load cleaned taxonomy dict
    with open(DATA_PROCESSED / "02_pathogens_taxid_cleaned_dict.json") as f:
        dict_taxonomy = json.load(f)

    taxid_to_pathogen = {
        str(tid): pathogen
        for pathogen, tids in dict_taxonomy.items()
        for tid in tids
    }
    pathogens = list(dict_taxonomy.keys())

    zip_files = sorted(DESC_DIR.glob("*.zip"))
    N_ZIPS = len(zip_files)
    N_WORKERS = min(6, cpu_count())

    all_records = []

    print(f"🔧 Using {N_WORKERS} worker processes")
    print(f"📦 Total ZIP files: {N_ZIPS}\n")

    start_time = time.time()
    ema_speed = None

    for i, zip_file in enumerate(zip_files, start=1):
        zip_chunk = zip_file.stem
        if zip_chunk in processed_chunks:
            print(f"⏩ Skipping {zip_chunk} (already processed)")
            continue

        zip_start = time.time()

        recs = process_one_zip(
            zip_file, DESC_DIR, KEEP_DESC_V3, N_WORKERS,
            taxid_to_pathogen, pathogens
        )
        all_records.extend(recs)

        with open(ZIP_LOG, "a") as f:
            f.write(zip_chunk + "\n")
        processed_chunks.add(zip_chunk)

        elapsed = time.time() - zip_start
        ema_speed = elapsed if ema_speed is None else 0.3 * elapsed + 0.7 * ema_speed
        zips_left = N_ZIPS - i
        eta_sec = zips_left * ema_speed

        print(f"[{i}/{N_ZIPS}] ZIP {zip_chunk} ✅ {len(recs)} matches | ETA {eta_sec/60:.1f} min")

    df_final = pd.DataFrame(all_records).drop_duplicates()
    df_final.to_csv(OUT_CSV, index=False)

    summary = (
        df_final.groupby("Pathogen")["AID"]
        .nunique()
        .reset_index(name="Downloaded_AIDs")
        .sort_values("Downloaded_AIDs", ascending=False)
    )
    summary.to_csv(OUT_CSV, index=False)

    total_minutes = (time.time() - start_time) / 60
    print("\n✅ Extraction complete!")
    print(f"📝 Summary saved to: {OUT_CSV}")
    print(f"📂 Filtered XMLs saved to: {KEEP_DESC_V3}")
    print(f"⏱ Total time: {total_minutes:.2f} min")


if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    main()