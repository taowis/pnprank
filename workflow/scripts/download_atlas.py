"""Download one registered atlas file into <data_root>/raw/<atlas_id>/ and check it
against the SHA-256 recorded in config/atlases.yaml.

Only files with `downloaded: true` and a `url` in the registry can be fetched this way;
those are the files the pipeline was run on. The size and disk-budget checks of
`pnprank download-atlas` apply.

Usage: python workflow/scripts/download_atlas.py <atlas_id> <file name as in the registry>
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

from pnprank.atlases import download_atlas_files

atlas_id, name = sys.argv[1], sys.argv[2]
registry = Path("config/atlases.yaml")
paths = yaml.safe_load(Path("config/paths.yaml").read_text())
entry = next(a for a in yaml.safe_load(registry.read_text())["atlases"] if a["id"] == atlas_id)
rec = next(f for f in entry["files"] if f["name"] == name)
if not (rec.get("downloaded") and rec.get("url") and rec.get("sha256")):
    raise SystemExit(f"{atlas_id}/{name}: no recorded url and sha256 in {registry}")
if not rec["url"].endswith("/" + name):
    raise SystemExit(f"{atlas_id}/{name}: url does not end with the file name")

before = registry.read_text()
got = download_atlas_files(registry, atlas_id, [name], rec["url"][: -len(name)],
                           paths["data_root"], paths["disk_limit_gb"])[0]
if got["sha256"] != rec["sha256"]:
    registry.write_text(before)  # download_atlas_files records the new hash; undo that
    raise SystemExit(f"{atlas_id}/{name}: sha256 {got['sha256']} differs from the "
                     f"recorded {rec['sha256']}; the host file has changed")
print(f"ok  {atlas_id}/{name}")
