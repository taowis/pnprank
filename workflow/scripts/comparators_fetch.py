"""Download the dated GO annotation, GO ontology and BioGRID releases listed in
config/comparators.yaml into <data_root>/external/comparators/ and write a manifest
with size, MD5 and SHA-256 of every file.

Usage: python workflow/scripts/comparators_fetch.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import yaml

cfg = yaml.safe_load(Path("config/comparators.yaml").read_text())
root = Path(yaml.safe_load(Path("config/paths.yaml").read_text())["data_root"])
out = root / "external" / "comparators"


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".part")
        subprocess.run(["curl", "-sSfL", "--retry", "5", "-A", "curl/8", "-o", str(tmp), url],
                       check=True)
        tmp.rename(dest)
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(dest, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return {"url": url, "path": str(dest.relative_to(root)), "bytes": dest.stat().st_size,
            "md5": md5.hexdigest(), "sha256": sha.hexdigest()}


manifest = []
for cutoff, r in cfg["go_annotation"]["releases"].items():
    rec = fetch(r["url"], out / "goa" / Path(r["url"]).name)
    manifest.append({"kind": "go_annotation", "cutoff": cutoff, "release": r["release"],
                     "date": r["date"], **rec})
for cutoff, r in cfg["go_ontology"]["releases"].items():
    rec = fetch(r["url"], out / "go" / r["release"] / Path(r["url"]).name)
    manifest.append({"kind": "go_ontology", "cutoff": cutoff, "release": r["release"],
                     "date": r["release"], **rec})
host = cfg["interactions"]["host"]
for cutoff, r in cfg["interactions"]["releases"].items():
    url = f"{host}BIOGRID-{r['release']}/{r['file']}"
    rec = fetch(url, out / "biogrid" / r["file"])
    manifest.append({"kind": "interactions", "cutoff": cutoff, "release": r["release"],
                     "date": r["date"], **rec})
    print(cutoff, r["release"], rec["bytes"], flush=True)

(out / "manifest.json").write_text(json.dumps(manifest, indent=1))
print(f"{len(manifest)} files; manifest at {out / 'manifest.json'}")
