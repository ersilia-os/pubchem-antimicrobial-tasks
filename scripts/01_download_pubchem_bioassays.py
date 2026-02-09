#!/usr/bin/env python3
"""
01_download_pubchem_bioassay_csv.py

Download PubChem BioAssay CSV ZIP archives (Description + Data) from NCBI FTP over HTTPS with:
- directory listing (HTML) → discovers all .zip files
- resumable downloads via HTTP Range requests
- parallel downloads (default: 6 workers)
- retries + exponential backoff for flaky network/slow server
- polite User-Agent header
- progress bar over files processed

Example:
  python scripts/01_download_pubchem_bioassays.py --out ../data/raw/bioassays --workers 6
"""

from __future__ import annotations

import argparse
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
import gzip
import shutil
import urllib.request
from urllib.error import URLError, HTTPError
import requests
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


URL_DESC = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Description/"
URL_DATA = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/CSV/Data/"
URL_BIOASSAYS = "https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/bioassays.tsv.gz"


# Reuse TCP connections + set a polite UA
SESSION = requests.Session()
SESSION.headers.update(
    {
        "User-Agent": "pubchem-bioassay-downloader/1.0 (contact: you@example.com)",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    }
)


@dataclass(frozen=True)
class DownloadJob:
    url: str
    out_path: Path


def parse_human_size_listing(url: str, timeout: float = 120.0, retries: int = 5) -> Optional[float]:
    """
    Best-effort: scrape directory HTML and sum sizes like "120M", "2.3G".
    Returns total bytes or None if pattern not found.
    """
    html = _get_text_with_retries(url, timeout=timeout, retries=retries)
    sizes = re.findall(r"\s+(\d+(?:\.\d+)?)([KMG])", html)
    if not sizes:
        return None

    total_bytes = 0.0
    for num, unit in sizes:
        x = float(num)
        if unit == "K":
            x *= 1e3
        elif unit == "M":
            x *= 1e6
        elif unit == "G":
            x *= 1e9
        total_bytes += x
    return total_bytes


def _get_text_with_retries(url: str, timeout: float, retries: int) -> str:
    """GET a URL and return text, with retries/backoff for timeouts/connection issues."""
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


def list_zip_files(base_url: str, timeout: float = 120.0, retries: int = 5) -> List[str]:
    """Return list of .zip file names from an NCBI FTP HTTPS directory listing (HTML)."""
    html = _get_text_with_retries(base_url, timeout=timeout, retries=retries)
    soup = BeautifulSoup(html, "html.parser")
    return [a.text.strip() for a in soup.find_all("a") if a.text.strip().endswith(".zip")]


def download_zip_resumable(
    url: str,
    out_path: Path,
    timeout: float = 120.0,
    max_retries: int = 5,
    chunk_size: int = 1024 * 1024,  # 1 MB
    polite_sleep: float = 0.0,
) -> Tuple[str, str, int]:
    """
    Download a ZIP file with resume support.
    Returns: (filename, status, bytes_written)
      status in {"skipped", "downloaded", "resumed", "failed"}
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    filename = out_path.name

    existing_size = out_path.stat().st_size if out_path.exists() else 0

    # Determine remote size (best-effort)
    remote_size: Optional[int] = None
    try:
        head = SESSION.head(url, timeout=timeout, allow_redirects=True)
        if head.status_code == 200 and "Content-Length" in head.headers:
            remote_size = int(head.headers["Content-Length"])
    except requests.RequestException:
        remote_size = None

    # If already complete, skip
    if remote_size is not None and existing_size >= remote_size and existing_size > 0:
        return (filename, "skipped", 0)

    headers = {}
    mode = "wb"
    resumed = False
    if existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"
        mode = "ab"
        resumed = True

    last_err: Optional[Exception] = None
    for attempt in range(max_retries + 1):
        try:
            with SESSION.get(url, headers=headers, stream=True, timeout=timeout) as r:
                # 200 = fresh, 206 = partial content (resume)
                if r.status_code not in (200, 206):
                    return (filename, "failed", 0)

                bytes_written = 0
                with open(out_path, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            bytes_written += len(chunk)

            if polite_sleep:
                time.sleep(polite_sleep)

            return (filename, "resumed" if resumed else "downloaded", bytes_written)

        except (requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectTimeout,
                requests.exceptions.ConnectionError) as e:
            last_err = e
            if attempt == max_retries:
                break
            time.sleep(0.5 * (2 ** attempt))
            continue
        except requests.RequestException as e:
            # other requests errors
            last_err = e
            break

    return (filename, "failed", 0)

def _download_with_retries(url: str, dst: Path, timeout: float, retries: int) -> None:
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            dst.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(url, timeout=timeout) as r, open(dst, "wb") as f:
                shutil.copyfileobj(r, f)
            return
        except (URLError, HTTPError, TimeoutError, OSError) as e:
            last_err = e
            if attempt < retries:
                # simple backoff
                time.sleep(min(2 ** (attempt - 1), 10))
            else:
                raise RuntimeError(f"Failed to download after {retries} attempts: {url}") from last_err


def download_and_unpack_bioassays_tsv(
    out_dir: Path,
    *,
    url: str = URL_BIOASSAYS,
    overwrite: bool = False,
    timeout: float = 120.0,
    retries: int = 5,
) -> Path:
    """
    Download PubChem BioAssay 'bioassays.tsv.gz' and unpack to 'bioassays.tsv'.
    Returns the path to the unpacked TSV.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gz_path = out_dir / "bioassays.tsv.gz"
    tsv_path = out_dir / "bioassays.tsv"

    if tsv_path.exists() and not overwrite:
        return tsv_path

    if overwrite or (not gz_path.exists()):
        _download_with_retries(url, gz_path, timeout=timeout, retries=retries)

    # Unpack (streaming)
    with gzip.open(gz_path, "rb") as fin, open(tsv_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)

    return tsv_path

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output directory root (e.g., data/raw/pubchem_bioassays)",
    )
    ap.add_argument("--workers", type=int, default=6, help="Parallel download workers (default: 6)")
    ap.add_argument("--desc", action="store_true", help="Download Description zips")
    ap.add_argument("--data", action="store_true", help="Download Data zips")
    ap.add_argument("--timeout", type=float, default=120.0, help="HTTP timeout seconds (default: 120)")
    ap.add_argument("--retries", type=int, default=5, help="Retries for listing and downloads (default: 5)")
    ap.add_argument("--polite-sleep", type=float, default=0.0, help="Sleep seconds after each file (default: 0)")
    ap.add_argument("--bioassays", action="store_true", help="Download and unpack Extras/bioassays.tsv.gz")
    ap.add_argument("--bioassays-url", type=str, default=URL_BIOASSAYS, help="Override bioassays.tsv.gz URL")
    ap.add_argument("--bioassays-overwrite", action="store_true", help="Redownload/re-unpack bioassays even if present")
    args = ap.parse_args()

    # default to both if neither selected
    if not args.desc and not args.data:
        args.desc = True
        args.data = True

    out_root: Path = args.out
    desc_dir = out_root / "Description"
    data_dir = out_root / "Data"
    
    if args.bioassays:
        print(f"Downloading/unpacking bioassays.tsv.gz -> {out_root}")
        try:
            tsv_path = download_and_unpack_bioassays_tsv(
                out_root,
                url=args.bioassays_url,
                overwrite=args.bioassays_overwrite,
                timeout=args.timeout,
                retries=args.retries,
            )
            print(f"Bioassays TSV ready: {tsv_path.resolve()}")
        except Exception as e:
            print(f"Bioassays download/unpack failed: {e}")
            return 2

    jobs: List[DownloadJob] = []

    if args.desc:
        print(f"Listing: {URL_DESC}")
        desc_files = list_zip_files(URL_DESC, timeout=args.timeout, retries=args.retries)
        total = parse_human_size_listing(URL_DESC, timeout=args.timeout, retries=args.retries)
        if total is not None:
            print(f"Description estimated total: {total/1e9:.3f} GB")
        jobs += [DownloadJob(URL_DESC + fn, desc_dir / fn) for fn in desc_files]

    if args.data:
        print(f"Listing: {URL_DATA}")
        data_files = list_zip_files(URL_DATA, timeout=args.timeout, retries=args.retries)
        total = parse_human_size_listing(URL_DATA, timeout=args.timeout, retries=args.retries)
        if total is not None:
            print(f"Data estimated total:        {total/1e9:.3f} GB")
        jobs += [DownloadJob(URL_DATA + fn, data_dir / fn) for fn in data_files]

    if not jobs:
        print("No jobs to run.")
        return 0

    print(f"\nTotal ZIPs to process: {len(jobs)}")
    print(f"Output root: {out_root.resolve()}")
    print(f"Workers: {args.workers}")
    print(f"Timeout: {args.timeout}s | Retries: {args.retries}\n")

    results = {"downloaded": 0, "resumed": 0, "skipped": 0, "failed": 0}

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

        for fut in tqdm(as_completed(futures), total=len(futures), desc="ZIPs processed"):
            filename, status, _nbytes = fut.result()
            results[status] += 1

    print("\n=== Summary ===")
    for k, v in results.items():
        print(f"{k:>10}: {v}")

    if results["failed"] > 0:
        print("\nSome downloads failed. Re-run the script; resume support will pick up partial files.")
        return 2

    print("\nAll done ✅")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())