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


# --- Mouse-to-human orthology -------------------------------------------------

COMPARA_URL = (
    "https://ftp.ensembl.org/pub/release-{release}/tsv/ensembl-compara/homologies/"
    "{species}/Compara.{release}.protein_default.homologies.tsv.gz"
)


def one_to_one_orthologues(
    compara_tsv_gz: str | Path, species: str = "mus_musculus", high_confidence: bool = False
) -> dict[str, str]:
    """Human Ensembl gene ID -> other-species gene ID, one-to-one orthologues only.

    Reads an Ensembl Compara homologies dump. Row order is arbitrary (Compara README),
    so a human-mouse pair may appear with either species first; both are handled. The
    human-directory dump omits mouse in release 116, so pass the mus_musculus dump.
    A pair is kept only if Compara calls it `ortholog_one2one`, and (as a guard)
    neither gene appears in any other one2one pair with that species.
    """
    import csv
    import gzip
    from collections import Counter

    pairs = set()
    with gzip.open(compara_tsv_gz, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["homology_type"] != "ortholog_one2one":
                continue
            if high_confidence and r["is_high_confidence"] != "1":
                continue
            if r["species"] == "homo_sapiens" and r["homology_species"] == species:
                pairs.add((r["gene_stable_id"], r["homology_gene_stable_id"]))
            elif r["species"] == species and r["homology_species"] == "homo_sapiens":
                pairs.add((r["homology_gene_stable_id"], r["gene_stable_id"]))
    h = Counter(p[0] for p in pairs)
    m = Counter(p[1] for p in pairs)
    return {a: b for a, b in sorted(pairs) if h[a] == 1 and m[b] == 1}
