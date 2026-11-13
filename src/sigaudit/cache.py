"""Download-once file cache for public datasets."""
from __future__ import annotations

import hashlib
import os
import shutil
import sys
from pathlib import Path

import requests


def data_dir() -> Path:
    """Where downloads live. Override with the SIGAUDIT_DATA environment variable."""
    d = os.environ.get("SIGAUDIT_DATA")
    path = Path(d) if d else Path.home() / ".cache" / "sigaudit"
    path.mkdir(parents=True, exist_ok=True)
    return path


def fetch(url: str, name: str, subdir: str = "raw", quiet: bool = False) -> Path:
    """Download `url` to <data_dir>/<subdir>/<name> unless it is already there."""
    dest = data_dir() / subdir / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    if not quiet:
        print(f"[sigaudit] downloading {name} ...", file=sys.stderr)
    with requests.get(url, stream=True, timeout=120,
                      headers={"User-Agent": "sigaudit (research software)"}) as r:
        r.raise_for_status()
        r.raw.decode_content = True  # undo HTTP-level gzip (file-level .gz is kept as is)
        with open(tmp, "wb") as f:
            shutil.copyfileobj(r.raw, f, length=1 << 20)
    tmp.replace(dest)
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def processed(name: str) -> Path:
    p = data_dir() / "processed" / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
