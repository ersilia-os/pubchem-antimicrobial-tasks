#!/usr/bin/env python3
"""
Download only the PubChem Bioassay Data ZIP blocks required by the selected AIDs.

Reads the aids_*.csv files produced by 02_bioassays_not_in_chembl.py, computes
which ZIP blocks contain those AIDs, and downloads only those blocks — avoiding
the full PubChem archive (hundreds of GB).

Each block covers 1000 AIDs (e.g. 0743001_0744000.zip contains AIDs 743001-744000).
Already-complete files are skipped; partial files are resumed.

Usage:
    python scripts/03_download_data_zips.py \\
        --aids-dir data/processed/bioassays_to_keep \\
        --out data/raw/bioassays/Data

"""

from __future__ import annotations

import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Set, Tuple

import requests
from bs4 import BeautifulSoup
from tqdm import tqdm

URL_DATA = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/"

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "pubchem-bioassay-downloader/1.0"})


# ---------------------------------------------------------------------------
# AID → block helpers
# ---------------------------------------------------------------------------

def aid_to_block(aid: int) -> str:
    start = ((aid - 1) // 1000) * 1000 + 1
    end = start + 999
    return f"{start:07d}_{end:07d}.zip"


def load_aids(aids_dir: Path) -> List[int]:
    aids: List[int] = []
    for csv_path in sorted(aids_dir.glob("aids_*.csv")):
        with open(csv_path, newline="") as f:
            for row in csv.DictReader(f):
                try:
                    aids.append(int(row["aid"]))
                except (KeyError, ValueError):
                    pass
    return aids


def required_blocks(aids: List[int]) -> Set[str]:
    return {aid_to_block(aid) for aid in aids}


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def _get_text(url: str, timeout: float, retries: int) -> str:
    last_err: Optional[Exception] = None
    for attempt in range(retries + 1):
        try:
            r = SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            return r.text
        except (requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectTimeout,
                requests.exceptions.ConnectionError) as e:
            last_err = e
            if attempt == retries:
                raise
            time.sleep(1.0 * (2 ** attempt))
    raise last_err  # pragma: no cover


def list_remote_zips(base_url: str, timeout: float, retries: int) -> Set[str]:
    html = _get_text(base_url, timeout=timeout, retries=retries)
    soup = BeautifulSoup(html, "html.parser")
    return {a.text.strip() for a in soup.find_all("a") if a.text.strip().endswith(".zip")}


@dataclass(frozen=True)
class DownloadJob:
    url: str
    out_path: Path


def download_zip_resumable(
    url: str,
    out_path: Path,
    timeout: float,
    max_retries: int,
    chunk_size: int = 1024 * 1024,
) -> Tuple[str, str]:
    """Download with resume support. Returns (filename, status)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    existing = out_path.stat().st_size if out_path.exists() else 0

    remote_size: Optional[int] = None
    try:
        head = SESSION.head(url, timeout=timeout, allow_redirects=True)
        if head.status_code == 200 and "Content-Length" in head.headers:
            remote_size = int(head.headers["Content-Length"])
    except requests.RequestException:
        pass

    if remote_size is not None and existing >= remote_size > 0:
        return (out_path.name, "skipped")

    headers = {"Range": f"bytes={existing}-"} if existing > 0 else {}
    mode = "ab" if existing > 0 else "wb"
    resumed = existing > 0

    for attempt in range(max_retries + 1):
        try:
            with SESSION.get(url, headers=headers, stream=True, timeout=timeout) as r:
                if r.status_code not in (200, 206):
                    return (out_path.name, "failed")
                with open(out_path, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
            return (out_path.name, "resumed" if resumed else "downloaded")
        except (requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectTimeout,
                requests.exceptions.ConnectionError):
            if attempt == max_retries:
                break
            time.sleep(0.5 * (2 ** attempt))
        except requests.RequestException:
            break

    return (out_path.name, "failed")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

root = Path(__file__).resolve().parent
aids_dir = root.parent / "data" / "processed" / "02_bioassays_to_keep"
out_dir = root.parent / "data" / "raw" / "03_data_zips"

WORKERS = 6
TIMEOUT = 120.0
RETRIES = 5


def main() -> int:
    out_dir.mkdir(parents=True, exist_ok=True)

    aids = load_aids(aids_dir)
    if not aids:
        print(f"No AIDs found in {aids_dir}. Did 02_bioassays_not_in_chembl.py run?")
        return 1

    needed = required_blocks(aids)
    print(f"{len(aids)} AIDs require {len(needed)} ZIP block(s)")

    print(f"Listing remote blocks: {URL_DATA}")
    available = list_remote_zips(URL_DATA, timeout=TIMEOUT, retries=RETRIES)
    missing = needed - available
    if missing:
        print(f"Warning: {len(missing)} block(s) not found on server: {sorted(missing)[:5]}")

    to_download = sorted(needed & available)
    jobs = [DownloadJob(URL_DATA + fn, out_dir / fn) for fn in to_download]

    if not jobs:
        print("Nothing to download.")
        return 0

    print(f"Downloading {len(jobs)} block(s) -> {out_dir.resolve()}")
    results = {"downloaded": 0, "resumed": 0, "skipped": 0, "failed": 0}

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futures = [
            ex.submit(download_zip_resumable, job.url, job.out_path, TIMEOUT, RETRIES)
            for job in jobs
        ]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="ZIPs"):
            _filename, status = fut.result()
            results[status] += 1

    print("\n=== Summary ===")
    for k, v in results.items():
        print(f"{k:>10}: {v}")

    if results["failed"] > 0:
        print("\nSome downloads failed. Re-run to resume partial files.")
        return 2

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
