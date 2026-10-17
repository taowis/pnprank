"""Download the non-atlas reference files named in config/config.yaml, check each
against the host checksum, and build the two derived tables used downstream.

Files (all under data/external/):
- gnomAD v4.1 constraint metrics           (GCS object MD5)
- GTEx v10 gene median TPM                 (GCS object MD5)
- Ensembl 116 mouse GTF -> mouse_genes_e116.tsv (gene ID, name, biotype)  (BSD sum)
- Ensembl 116 human GTF (gene coordinates; table built by gene_coordinates.py)  (BSD sum)
- Ensembl Compara 116 mus_musculus homologies -> human_mouse_one2one_e116.tsv  (MD5SUM)
- UCSC hg38 cytoBand table                 (SHA-256 recorded at first download)

Existing files are kept if their checksum matches; a mismatch is an error, not a
silent re-download. Derived tables are rebuilt only when missing.

The truth-set sources (HGNC, PanelApp, ClinGen) are snapshotted by
`pnprank build-candidates`, the comparator releases by comparators_fetch.py and the
atlases by download_atlas.py.

Usage: python workflow/scripts/fetch_external.py
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import re
import subprocess
from pathlib import Path

import yaml

from pnprank.atlases import one_to_one_orthologues

cfg = yaml.safe_load(Path("config/config.yaml").read_text())
ext, orth = cfg["external"], cfg["reference"]["orthology"]


def digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 24):
            h.update(chunk)
    return h.hexdigest()


def bsd_sum(path: Path) -> str:
    # `sum` defaults to the BSD algorithm with 1 KB blocks on both macOS and GNU coreutils
    out = subprocess.run(["sum", str(path)], check=True, capture_output=True, text=True).stdout
    return " ".join(out.split()[:2])


def fetch(url: str, dest: Path, check) -> None:
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        subprocess.run(["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(tmp), url],
                       check=True)
        tmp.rename(dest)
    got, want = check(dest)
    if got != want:
        raise RuntimeError(f"{dest}: checksum {got} does not match {want}")
    print(f"ok  {dest}")


for key in ("gnomad_constraint", "gtex_median_tpm"):
    e = ext[key]
    fetch(e["url"], Path(e["dest"]), lambda p, e=e: (
        base64.b64encode(bytes.fromhex(digest(p, "md5"))).decode(), e["md5_base64"]))

g = ext["mouse_gtf"]
fetch(g["url"], Path(g["dest"]), lambda p: (bsd_sum(p), g["bsd_sum"]))

h = ext["human_gtf"]
fetch(h["url"], Path(h["dest"]), lambda p: (bsd_sum(p), h["bsd_sum"]))

c = ext["ucsc_cytoband"]
fetch(c["url"], Path(c["dest"]), lambda p: (digest(p, "sha256"), c["sha256"]))

fetch(orth["url"], Path(orth["dest"]), lambda p: (digest(p, "md5"), orth["md5"]))

# mouse gene ID -> name, biotype (gene_name is absent for some genes; kept as "")
genes = Path(g["genes_table"])
if not genes.exists():
    n = 0
    with gzip.open(g["dest"], "rt") as fh, open(genes, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            out.write(f"{a['gene_id']}\t{a.get('gene_name', '')}\t{a.get('gene_biotype', '')}\n")
            n += 1
    print(f"wrote {genes} ({n} genes)")

table = Path(orth["table"])
if not table.exists():
    o = one_to_one_orthologues(orth["dest"])
    with open(table, "w") as fh:
        fh.write("human_ensembl_gene_id\tmouse_ensembl_gene_id\n")
        for a, b in sorted(o.items()):
            fh.write(f"{a}\t{b}\n")
    print(f"wrote {table} ({len(o)} one-to-one pairs)")
