"""Non-atlas reference data (gene constraint, bulk tissue expression, annotation).

Each file is listed in `config/config.yaml` under `external` with its URL and the
checksum published by the host; downloads are verified against it.
"""

from __future__ import annotations

import base64
import hashlib
import subprocess
from pathlib import Path


def fetch(url: str, dest: Path, md5_base64: str | None = None) -> str:
    """Download url to dest (resuming), verify the host's MD5 if given, return hex MD5."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(dest), url], check=True)
    h = hashlib.md5()
    with open(dest, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if md5_base64 is not None and base64.b64encode(h.digest()).decode() != md5_base64:
        raise RuntimeError(f"{dest.name}: MD5 mismatch with host")
    return h.hexdigest()
