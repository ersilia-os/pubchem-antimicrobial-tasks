import gzip
import shutil
import zipfile
import time
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
from tqdm.auto import tqdm
from bs4 import BeautifulSoup


# ===================================================================
# ------------------------ CONFIG -----------------------------------
# ===================================================================

def setup_paths():
    PROJECT_ROOT = Path(__file__).resolve().parents[1]

    DATA_RAW = PROJECT_ROOT / "data" / "raw"
    DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

    PUBCHEM_DIR = DATA_RAW / "pubchem_bioassays"
    DESC_DIR = PUBCHEM_DIR / "Description"
    DATA_DIR = PUBCHEM_DIR / "Data"

    KEEP_DIR = DATA_RAW / "filtered_assays"
    KEEP_DESC = KEEP_DIR / "Description"
    KEEP_DATA = KEEP_DIR / "Data"

    KEEP_DESC.mkdir(parents=True, exist_ok=True)
    KEEP_DATA.mkdir(parents=True, exist_ok=True)

    return PROJECT_ROOT, DATA_RAW, DATA_PROCESSED, DESC_DIR, DATA_DIR, KEEP_DESC, KEEP_DATA


# ===================================================================
# ------------------------ XML PARSING ------------------------------
# ===================================================================

def extract_info_from_descr(xml_path: Path):
    if str(xml_path).endswith(".gz"):
        with gzip.open(xml_path, "rb") as f:
            xml_content = f.read()
    else:
        with open(xml_path, "rb") as f:
            xml_content = f.read()

    soup = BeautifulSoup(xml_content, "lxml-xml")

    # TaxIDs
    taxids = set()
    for node in soup.find_all("PC-XRefData_taxonomy"):
        txt = node.text.strip()
        if txt.isdigit():
            taxids.add(txt)
    for node in soup.find_all("PC-AssayTarget_tax-id"):
        txt = node.text.strip()
        if txt.isdigit():
            taxids.add(txt)

    # Assay organism
    assay_organisms = []
    for node in soup.find_all("PC-AnnotatedXRef_comment"):
        txt = node.text.strip()
        if txt:
            assay_organisms.append(txt)
    assay_organisms = [x for x in assay_organisms if len(x) > 3]

    # ChEMBL
    chembl_id = None
    for db in soup.find_all("PC-DBTracking"):
        name = db.find("PC-DBTracking_name")
        if name and name.text.strip() == "ChEMBL":
            obj = db.find("Object-id_str")
            if obj:
                chembl_id = obj.text.strip()
                break

    return {
        "TaxIDs": sorted(list(taxids)),
        "AssayOrganism": assay_organisms,
        "ChEMBL": chembl_id,
    }


# ===================================================================
# ------------------------ PATHOGEN MATCHING ------------------------
# ===================================================================

def detect_matching_pathogens(info, taxid_to_pathogen, pathogens):
    matched = set()

    # by TaxID
    for tid in info["TaxIDs"]:
        if tid in taxid_to_pathogen:
            matched.add(taxid_to_pathogen[tid])

    # by organism name
    for org in info["AssayOrganism"]:
        org_low = org.lower()
        for pathogen in pathogens:
            if pathogen.lower() in org_low:
                matched.add(pathogen)

    return sorted(matched)


# ===================================================================
# ------------------------ WORKER FUNCTION --------------------------
# ===================================================================

def worker_process_xml(args):
    xml_path_str, taxid_to_pathogen, pathogens = args
    xml_path = Path(xml_path_str)

    try:
        info = extract_info_from_descr(xml_path)
        pathogens_hit = detect_matching_pathogens(info, taxid_to_pathogen, pathogens)

        if not pathogens_hit:
            return None

        stem = xml_path.stem
        aid_str = stem.split(".")[0]
        aid = int(aid_str)

        return {
            "xml_path": xml_path_str,
            "AID": aid,
            "Pathogens": pathogens_hit,
            "ChEMBL": info["ChEMBL"],
        }

    except Exception:
        return None


# ===================================================================
# ------------------------ PROCESS ONE ZIP --------------------------
# ===================================================================

def process_one_zip(zip_file, already_saved, DESC_DIR, KEEP_DESC, N_WORKERS,
                    taxid_to_pathogen, pathogens, LOG_SPEED, SPEED_LOG):

    zip_chunk = zip_file.stem
    t0 = time.time()

    temp = DESC_DIR / f"{zip_chunk}_tmp"
    temp.mkdir(exist_ok=True)

    with zipfile.ZipFile(zip_file, "r") as zf:
        zf.extractall(temp)

    xml_files = [str(p) for p in temp.rglob("*.xml")] + \
                [str(p) for p in temp.rglob("*.xml.gz")]

    tasks = [(fp, taxid_to_pathogen, pathogens) for fp in xml_files]

    records = []

    with Pool(processes=N_WORKERS) as pool:
        for result in pool.imap_unordered(worker_process_xml, tasks, chunksize=20):
            if result is None:
                continue

            aid = result["AID"]
            xml_path = Path(result["xml_path"])

            if aid not in already_saved:
                dst = KEEP_DESC / f"{aid}.xml"
                if dst.exists():
                    dst.unlink()

                if xml_path.suffix == ".gz":
                    with gzip.open(xml_path, "rb") as fi, open(dst, "wb") as fo:
                        fo.write(fi.read())
                else:
                    shutil.copy2(xml_path, dst)

                already_saved.add(aid)

            records.append({
                "AID": aid,
                "Pathogen": ", ".join(result["Pathogens"]),
                "ChEMBLid": result["ChEMBL"],
                "ZipFolder": zip_chunk,
            })

    shutil.rmtree(temp)
    elapsed = time.time() - t0

    if LOG_SPEED:
        with open(SPEED_LOG, "a") as f:
            f.write(f"{zip_chunk}\t{len(xml_files)}\t{len(records)}\t{elapsed:.2f}\n")

    print(f"✓ ZIP {zip_chunk}: {len(xml_files)} XMLs, {len(records)} matches, {elapsed:.1f}s")
    return records, already_saved


# ===================================================================
# ------------------------------- MAIN ------------------------------
# ===================================================================

def main():

    # Paths
    PROJECT_ROOT, DATA_RAW, DATA_PROCESSED, DESC_DIR, DATA_DIR, KEEP_DESC, KEEP_DATA = setup_paths()

    # Config
    DEBUG = False
    LOG_SPEED = True
    N_WORKERS = min(6, cpu_count())
    SPEED_LOG = DATA_PROCESSED / "parallel_speed_log.txt"
    ZIP_LOG = DATA_PROCESSED / "processed_zip_chunks_with_organisms.txt"
    OUT_CSV = DATA_PROCESSED / "filtered_description_with_organisms.csv"

    # Load taxonomy & pathogens
    tax_df = pd.read_csv(DATA_PROCESSED / "taxonomy_table.csv")

    dict_taxonomy = (
        tax_df.groupby("Pathogen")["Taxonomy_ID"]
        .apply(lambda s: list(map(str, s)))
        .to_dict()
    )

    taxid_to_pathogen = {}
    for pathogen, taxids in dict_taxonomy.items():
        for tid in taxids:
            taxid_to_pathogen[str(tid)] = pathogen

    pathogens = list(dict_taxonomy.keys())

    # List zip files
    zip_files = sorted(DESC_DIR.glob("*.zip"))

    # Resume
    processed_chunks = set()
    if ZIP_LOG.exists():
        processed_chunks = set(ZIP_LOG.read_text().splitlines())

    # AIDs already saved
    already_saved = set(int(f.stem) for f in KEEP_DESC.glob("*.xml") if f.stem.isdigit())

    all_records = []
    start = time.time()

    for zip_file in tqdm(zip_files, desc="Processing ZIP chunks"):
        zip_chunk = zip_file.stem

        if zip_chunk in processed_chunks:
            continue

        recs, already_saved = process_one_zip(
            zip_file, already_saved,
            DESC_DIR, KEEP_DESC, N_WORKERS,
            taxid_to_pathogen, pathogens,
            LOG_SPEED, SPEED_LOG
        )

        all_records.extend(recs)

        with open(ZIP_LOG, "a") as zf:
            zf.write(zip_chunk + "\n")
        processed_chunks.add(zip_chunk)

    # Save results
    df_final = pd.DataFrame(all_records).drop_duplicates()
    df_final.to_csv(OUT_CSV, index=False)

    print("\n=====================================")
    print("✔ ALL ZIP FILES PROCESSED")
    print(f"✔ Total unique AIDs: {df_final['AID'].nunique()}")
    print(f"✔ Output CSV: {OUT_CSV}")
    print(f"✔ XML saved in: {KEEP_DESC}")
    print(f"⏱ Total time: {(time.time()-start)/3600:.2f} hours")
    print("=====================================")


# ===================================================================
# ----------------------- ENTRY POINT -------------------------------
# ===================================================================

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()