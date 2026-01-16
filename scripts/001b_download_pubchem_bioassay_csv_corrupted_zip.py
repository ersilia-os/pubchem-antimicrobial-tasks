#!/usr/bin/env python3
"""
Recovery script: re-download corrupted PubChem BioAssay Description ZIPs.

- Deletes corrupted files listed in data/processed/corrupted_description_zips.txt
- Downloads them again with resume and retry support

Example:
  python scripts/001b_download_pubchem_bioassay_csv_corrupted_zip.py --out data/raw/pubchem_bioassays --desc --workers 6
"""

from __future__ import annotations

import argparse
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import requests
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# -------------------------------------------------------------------
# URLs
# -------------------------------------------------------------------

URL_DESC = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/"

# -------------------------------------------------------------------
# HTTP session
# -------------------------------------------------------------------

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "pubchem-bioassay-downloader/1.0 (contact: you@example.com)",
    "Accept": "*/*",
})

# -------------------------------------------------------------------
# Data structures
# -------------------------------------------------------------------

@dataclass(frozen=True)
class DownloadJob:
    url: str
    out_path: Path

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

def _get_text_with_retries(url: str, timeout: float, retries: int) -> str:
    for attempt in range(retries + 1):
        try:
            r = SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            return r.text
        except requests.RequestException as e:
            if attempt == retries:
                raise
            time.sleep(1.0 * (2 ** attempt))


def list_zip_files(base_url: str, timeout: float = 120.0, retries: int = 5) -> List[str]:
    html = _get_text_with_retries(base_url, timeout=timeout, retries=retries)
    soup = BeautifulSoup(html, "html.parser")
    return [a.text.strip() for a in soup.find_all("a") if a.text.strip().endswith(".zip")]


def download_zip_resumable(
    url: str,
    out_path: Path,
    timeout: float = 120.0,
    max_retries: int = 5,
    chunk_size: int = 1024 * 1024,
    polite_sleep: float = 0.0,
) -> Tuple[str, str, int]:

    out_path.parent.mkdir(parents=True, exist_ok=True)
    filename = out_path.name

    headers = {}
    mode = "wb"

    for attempt in range(max_retries + 1):
        try:
            with SESSION.get(url, headers=headers, stream=True, timeout=timeout) as r:
                if r.status_code != 200:
                    return (filename, "failed", 0)

                bytes_written = 0
                with open(out_path, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            bytes_written += len(chunk)

            if polite_sleep:
                time.sleep(polite_sleep)

            return (filename, "downloaded", bytes_written)

        except (requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectTimeout,
                requests.exceptions.ConnectionError):
            if attempt == max_retries:
                break
            time.sleep(0.5 * (2 ** attempt))

    return (filename, "failed", 0)

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--timeout", type=float, default=120.0)
    ap.add_argument("--retries", type=int, default=5)
    ap.add_argument("--polite-sleep", type=float, default=0.0)
    args = ap.parse_args()

    out_root = args.out
    desc_dir = out_root / "Description"

    # Load corrupted ZIP list
    corrupted_list_path = Path("data/processed/corrupted_description_zips.txt")
    if not corrupted_list_path.exists():
        print("❌ No corrupted_description_zips.txt found. Aborting.")
        return 1

    with open(corrupted_list_path) as f:
        corrupted_zips = sorted({line.strip() for line in f if line.strip()})

    if not corrupted_zips:
        print("✅ List is empty — nothing to re-download.")
        return 0

    print(f"♻️ Recovery mode: {len(corrupted_zips)} corrupted ZIPs will be re-downloaded")
    print(f"🧹 Deleting old copies...")

    for name in corrupted_zips:
        local_path = desc_dir / name
        if local_path.exists():
            try:
                local_path.unlink()
                print(f"🗑️  Deleted: {name}")
            except Exception as e:
                print(f"⚠️  Failed to delete {name}: {e}")

    print(f"\n🌐 Fetching file list from PubChem...")
    available = list_zip_files(URL_DESC, timeout=args.timeout, retries=args.retries)
    available_set = set(available)

    jobs: List[DownloadJob] = []
    for name in corrupted_zips:
        if name not in available_set:
            print(f"⚠️  Skipping missing file on server: {name}")
            continue
        url = URL_DESC + name
        out_path = desc_dir / name
        jobs.append(DownloadJob(url, out_path))

    print(f"\n📥 Starting download of {len(jobs)} files...")

    results = {"downloaded": 0, "failed": 0}

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = [
            ex.submit(
                download_zip_resumable,
                job.url,
                job.out_path,
                args.timeout,
                args.retries,
                1024 * 1024,
                args.polite_sleep,
            )
            for job in jobs
        ]

        for fut in tqdm(as_completed(futures), total=len(futures), desc="Downloading ZIPs"):
            _, status, _ = fut.result()
            results[status] += 1

    print("\n=== Summary ===")
    for k, v in results.items():
        print(f"{k:>10}: {v}")

    print("\n✅ Recovery complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())