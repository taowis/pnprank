"""Download, checksum, load atlases and harmonise cell-type labels.

The registry is `config/atlases.yaml`. Every entry must carry a verified accession,
the URL where it was verified, a licence statement and, once downloaded, checksums.
"""

from __future__ import annotations

from pathlib import Path

import yaml

REQUIRED_FIELDS = [
    "id",
    "title",
    "species",
    "tissue",
    "accession",
    "url",
    "source_url_verified",
    "licence",
    "publication",
    "access",
    "files",
]
SPECIES = {"human", "mouse", "rat", "macaque"}
ACCESS = {"open", "controlled", "on request"}


def load_registry(path: str | Path) -> list[dict]:
    """Atlas entries from the registry YAML."""
    return yaml.safe_load(Path(path).read_text()).get("atlases") or []


def validate_registry(path: str | Path) -> list[str]:
    """Return a list of problems in the atlas registry; empty means valid."""
    problems: list[str] = []
    seen: set[str] = set()
    for i, a in enumerate(load_registry(path)):
        label = a.get("id", f"entry {i}")
        missing = [f for f in REQUIRED_FIELDS if not a.get(f) and a.get(f) != []]
        if missing:
            problems.append(f"{label}: missing {missing}")
        if label in seen:
            problems.append(f"{label}: duplicate id")
        seen.add(label)
        if a.get("species") not in SPECIES:
            problems.append(f"{label}: unknown species {a.get('species')!r}")
        if a.get("access") not in ACCESS:
            problems.append(f"{label}: unknown access {a.get('access')!r}")
        if not str(a.get("source_url_verified", "")).startswith("http"):
            problems.append(f"{label}: source_url_verified is not a URL")
        for f in a.get("files") or []:
            if f.get("downloaded") and not f.get("sha256"):
                problems.append(f"{label}: downloaded file {f.get('name')} has no sha256")
    return problems
