#!/usr/bin/env python3
"""
Download PubChem bioassays.tsv.gz from NCBI FTP and unpack it.

The unpacked TSV is required by 02_bioassays_not_in_chembl.py to compare
assay metadata (Source ID, compound counts) against ChEMBL.

Skips the download if bioassays.tsv already exists.

Usage:
    python scripts/01_download_bioassays_tsv.py
"""

from __future__ import annotations

import gzip
import os
import shutil
import time
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError

URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/bioassays.tsv.gz"

root = Path(__file__).resolve().parent
out_dir = root.parent / "data" / "raw" / "01_bioassays"


def download(url: str, dst: Path, timeout: float = 120.0, retries: int = 5) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r, open(dst, "wb") as f:
                shutil.copyfileobj(r, f)
            return
        except (URLError, HTTPError, TimeoutError, OSError) as e:
            last_err = e
            if attempt < retries:
                time.sleep(min(2 ** (attempt - 1), 10))
    raise RuntimeError(f"Download failed after {retries} attempts: {url}") from last_err


def main() -> int:
    out_dir.mkdir(parents=True, exist_ok=True)
    gz_path = out_dir / "bioassays.tsv.gz"
    tsv_path = out_dir / "bioassays.tsv"

    if tsv_path.exists():
        print(f"Already exists, skipping: {tsv_path.resolve()}")
        return 0

    print(f"Downloading {URL} ...")
    download(URL, gz_path)
    print("Unpacking ...")
    with gzip.open(gz_path, "rb") as fin, open(tsv_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    gz_path.unlink()
    print(f"Ready: {tsv_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
