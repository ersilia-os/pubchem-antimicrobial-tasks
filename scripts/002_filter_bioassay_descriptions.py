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
    """
    Extract:
    - TaxIDs (3 sources)
    - Assay organism strings
    """
    if xml_path.suffix == ".gz":
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

    assay_organisms = [
        n.text.strip()
        for n in soup.find_all("PC-AnnotatedXRef_comment")
        if len(n.text.strip()) > 3
    ]

    return {
        "TaxIDs": sorted(taxids),
        "AssayOrganism": assay_organisms,
    }


# ===================================================================
# ------------------------ PATHOGEN MATCHING ------------------------
# ===================================================================

def detect_matching_pathogens(info, taxid_to_pathogen, pathogens):
    """
    Returns list of (TaxID | None, Pathogen)
    """
    matched = set()

    # Match by TaxID
    for tid in info["TaxIDs"]:
        if tid in taxid_to_pathogen:
            matched.add((tid, taxid_to_pathogen[tid]))

    # Match by organism name
    for org in info["AssayOrganism"]:
        org_low = org.lower()
        for p in pathogens:
            if p.lower() in org_low:
                matched.add((None, p))

    return sorted(matched)


# ===================================================================
# ------------------------ WORKER FUNCTION --------------------------
# ===================================================================

def worker(args):
    xml_path_str, taxid_to_pathogen, pathogens = args
    xml_path = Path(xml_path_str)

    try:
        info = extract_info_from_descr(xml_path)
        hits = detect_matching_pathogens(info, taxid_to_pathogen, pathogens)

        if not hits:
            return []

        aid = int(xml_path.stem.split(".")[0])

        return [
            (aid, taxid, pathogen, xml_path_str)
            for taxid, pathogen in hits
        ]

    except Exception:
        return []


# ===================================================================
# ------------------------------- MAIN ------------------------------
# ===================================================================

def main():
    PROJECT_ROOT = Path(__file__).resolve().parents[1]

    DATA_RAW = PROJECT_ROOT / "data" / "raw"
    DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

    DESC_DIR = DATA_RAW / "pubchem_bioassays" / "Description"
    KEEP_DESC = DATA_RAW / "filtered_assays" / "Description"
    KEEP_DESC.mkdir(parents=True, exist_ok=True)

    OUT_PATH = DATA_PROCESSED / "04_aid_xml_filtered.csv"

    # ---------------------------------------------------
    # Load curated taxonomy
    # ---------------------------------------------------
    with open(DATA_PROCESSED / "02_pathogens_taxid_cleaned_dict.json") as f:
        dict_taxonomy = json.load(f)

    taxid_to_pathogen = {
        str(tid): pathogen
        for pathogen, tids in dict_taxonomy.items()
        for tid in tids
    }
    pathogens = list(dict_taxonomy.keys())

    # ---------------------------------------------------
    # ZIP files 
    # ---------------------------------------------------
    zip_files = sorted(DESC_DIR.glob("*.zip"))
    # Example for testing:
    # zip_files = [DESC_DIR / "2062001_2063000.zip"]

    total_zips = len(zip_files)
    N_WORKERS = min(cpu_count(), 6)

    # ---------------------------------------------------
    # Resume support
    # ---------------------------------------------------
    if OUT_PATH.exists() and OUT_PATH.stat().st_size > 0:
        df_existing = pd.read_csv(OUT_PATH)
        completed = set(df_existing[["AID", "ZipFolder"]].apply(tuple, axis=1))
        all_records = df_existing.to_dict("records")
    else:
        completed = set()
        all_records = []

    print(f"🔧 Workers: {N_WORKERS}")
    print(f"📦 Total ZIPs: {total_zips}\n")

    start_time = time.time()
    ema_speed = None

    # ---------------------------------------------------
    # Main loop
    # ---------------------------------------------------
    for i, zip_file in enumerate(zip_files, 1):
        zip_chunk = zip_file.stem

        # Skip ZIP if already fully processed
        if any(z == zip_chunk for (_, z) in completed):
            print(f"⏩ Skipping {zip_chunk} (already processed)")
            continue

        temp_dir = DESC_DIR / f"{zip_chunk}_tmp"
        temp_dir.mkdir(exist_ok=True)

        try:
            with zipfile.ZipFile(zip_file, "r") as zf:
                zf.extractall(temp_dir)
        except zipfile.BadZipFile:
            print(f"❌ Corrupted ZIP skipped: {zip_chunk}")
            shutil.rmtree(temp_dir)
            continue

        xml_files = (
            list(temp_dir.rglob("*.xml")) +
            list(temp_dir.rglob("*.xml.gz"))
        )

        tasks = [(str(p), taxid_to_pathogen, pathogens) for p in xml_files]

        zip_start = time.time()
        new_rows = []

        with Pool(N_WORKERS) as pool:
            for results in pool.imap_unordered(worker, tasks, chunksize=20):
                for aid, taxid, pathogen, xml_path_str in results:

                    key = (aid, zip_chunk)
                    if key in completed:
                        continue

                    # Save XML (once per AID)
                    dst = KEEP_DESC / f"{aid}.xml"
                    src = Path(xml_path_str)

                    if not dst.exists():
                        if src.suffix == ".gz":
                            with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                                fo.write(fi.read())
                        else:
                            shutil.copy2(src, dst)

                    new_rows.append({
                        "AID": aid,
                        "Pathogen": pathogen,
                        "TaxID": taxid,
                        "ZipFolder": zip_chunk
                    })

                    completed.add(key)

        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            print(f"⚠️ Could not delete temp folder {temp_dir}: {e}")

        all_records.extend(new_rows)
        pd.DataFrame(all_records).drop_duplicates().to_csv(OUT_PATH, index=False)

        elapsed = time.time() - zip_start
        ema_speed = elapsed if ema_speed is None else 0.3 * elapsed + 0.7 * ema_speed
        zips_left = total_zips - i
        eta_sec = zips_left * ema_speed

        print(
            f"[{i}/{total_zips}] {zip_chunk} ✅ "
            f"{len(new_rows)} new | ETA {eta_sec/60:.1f} min"
        )

    total_minutes = (time.time() - start_time) / 60
    print("\n✅ Finished filtering bioassays.")
    print(f"📝 Output saved to: {OUT_PATH}")
    print(f"📂 XMLs saved to: {KEEP_DESC}")
    print(f"⏱ Total time: {total_minutes:.1f} min")


if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    main()