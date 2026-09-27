"""Streaming downloads with resume, caching and a SHA-256 manifest.

Files are streamed to ``<name>.part`` and renamed on completion so an interrupted
download can be resumed via an HTTP Range request. Files already present in
``data/raw`` are skipped. Sizes, SHA-256 hashes, the download date and the
STRING/UniProt release versions are logged to ``data/raw/MANIFEST.txt``.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import os
from typing import Dict, Iterable, Optional

import requests
from tqdm import tqdm

CHUNK = 1 << 20  # 1 MiB


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def _filename_for(url: str) -> str:
    return url.split("?")[0].rstrip("/").split("/")[-1]


def fetch_uniprot_release(cfg: dict) -> str:
    """Read the current UniProt KB release string from reldate.txt (if reachable)."""
    url = cfg.get("urls", {}).get("uniprot_reldate")
    if not url:
        return "unknown"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        first = r.text.strip().splitlines()[0]
        # e.g. "UniProt Knowledgebase Release 2026_03 consists of:"
        return first.replace(" consists of:", "").strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def download_one(url: str, dest_dir: str, filename: Optional[str] = None,
                 force: bool = False) -> str:
    """Download *url* into *dest_dir* with resume support. Returns the file path."""
    os.makedirs(dest_dir, exist_ok=True)
    filename = filename or _filename_for(url)
    dest = os.path.join(dest_dir, filename)
    part = dest + ".part"

    if os.path.exists(dest) and not force:
        print(f"[skip] {filename} already present ({os.path.getsize(dest):,} bytes)")
        return dest

    resume_from = os.path.getsize(part) if os.path.exists(part) else 0
    headers = {}
    if resume_from:
        headers["Range"] = f"bytes={resume_from}-"
        print(f"[resume] {filename} from byte {resume_from:,}")

    with requests.get(url, stream=True, headers=headers, timeout=120) as r:
        # If the server ignores the Range request, restart from scratch.
        if resume_from and r.status_code == 200:
            resume_from = 0
        r.raise_for_status()

        total = int(r.headers.get("Content-Length", 0))
        if total:
            total += resume_from

        mode = "ab" if resume_from else "wb"
        with open(part, mode) as fh, tqdm(
            total=total or None, initial=resume_from, unit="B",
            unit_scale=True, desc=filename,
        ) as bar:
            for chunk in r.iter_content(chunk_size=CHUNK):
                if chunk:
                    fh.write(chunk)
                    bar.update(len(chunk))

    os.replace(part, dest)
    print(f"[done] {filename} ({os.path.getsize(dest):,} bytes)")
    return dest


def write_manifest(dest_dir: str, files: Dict[str, str], versions: Dict[str, str]) -> str:
    """Write a MANIFEST.txt with sizes, hashes, date and release versions."""
    manifest = os.path.join(dest_dir, "MANIFEST.txt")
    now = _dt.datetime.now(_dt.timezone.utc).isoformat()
    with open(manifest, "w") as fh:
        fh.write("# STRING network pipeline download manifest\n")
        fh.write(f"# generated: {now}\n")
        for k, v in versions.items():
            fh.write(f"# {k}: {v}\n")
        fh.write("#\n")
        fh.write("# sha256  size_bytes  filename  source_url\n")
        for url, path in files.items():
            if not os.path.exists(path):
                continue
            digest = _sha256(path)
            size = os.path.getsize(path)
            fh.write(f"{digest}  {size}  {os.path.basename(path)}  {url}\n")
    print(f"[manifest] wrote {manifest}")
    return manifest


# Which config URL keys to fetch for a standard `make data` run.
DEFAULT_KEYS: Iterable[str] = (
    "links",
    "links_detailed",
    "links_physical",
    "links_physical_detailed",
    "protein_info",
    "uniprot_idmapping",
    "string_aliases",
    "reviewed_list",
)


def download_all(cfg: dict, keys: Optional[Iterable[str]] = None,
                 force: bool = False) -> Dict[str, str]:
    """Download all configured source files and refresh the manifest."""
    raw = cfg["paths"]["raw"]
    urls = cfg["urls"]
    keys = list(keys or DEFAULT_KEYS)

    fetched: Dict[str, str] = {}
    for key in keys:
        url = urls.get(key)
        if not url:
            print(f"[warn] no URL configured for '{key}', skipping")
            continue
        # The reviewed list has no natural filename; give it one.
        fname = None
        if key == "reviewed_list":
            fname = f"reviewed_{cfg['taxon']}.list"
        try:
            path = download_one(url, raw, filename=fname, force=force)
            fetched[url] = path
        except Exception as exc:  # noqa: BLE001
            print(f"[error] failed to download {key}: {exc}")

    # Build the manifest from *every* known source file present in raw, so a
    # partial download does not shrink the manifest.
    all_files: Dict[str, str] = {}
    for key, url in urls.items():
        fname = f"reviewed_{cfg['taxon']}.list" if key == "reviewed_list" \
            else _filename_for(url)
        path = os.path.join(raw, fname)
        if os.path.exists(path):
            all_files[url] = path
    all_files.update(fetched)

    versions = {
        "string_version": cfg.get("string_version", "?"),
        "taxon": cfg.get("taxon", "?"),
        "uniprot_release": fetch_uniprot_release(cfg),
        "build_threshold (config default)": cfg.get("threshold", "?"),
    }
    write_manifest(raw, all_files, versions)
    return fetched
