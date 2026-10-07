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


# --- Downloads ----------------------------------------------------------------

SINGLE_DOWNLOAD_LIMIT_BYTES = 5 * 10**9  # CLAUDE.md: ask before any single download over 5 GB


def dir_size_bytes(path: Path) -> int:
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())


def remote_size_bytes(url: str) -> int | None:
    """Content-Length from a HEAD request, or None if the server does not say."""
    import subprocess

    out = subprocess.run(
        ["curl", "-sIL", "-m", "60", url], capture_output=True, text=True, check=True
    ).stdout
    sizes = [int(line.split(":", 1)[1]) for line in out.splitlines()
             if line.lower().startswith("content-length:")]
    return sizes[-1] if sizes else None


def parallel_range_download(url: str, dest: Path, size: int, parts: int = 8) -> None:
    """Download `size` bytes of url as `parts` concurrent byte ranges, then concatenate.

    Some hosts (NCBI GEO) cap per-connection speed; parallel ranges work around that.
    Each part resumes if a partial file exists.
    """
    import subprocess

    bounds = [(i * size // parts, (i + 1) * size // parts - 1) for i in range(parts)]
    procs = []
    for i, (a, b) in enumerate(bounds):
        part = dest.with_name(f"{dest.name}.part{i}")
        have = part.stat().st_size if part.exists() else 0
        if have >= b - a + 1:
            continue
        procs.append(subprocess.Popen(
            ["curl", "-sSfL", "--retry", "5", "-r", f"{a + have}-{b}", "-o", "-", url],
            stdout=open(part, "ab"),
        ))
    if any(p.wait() != 0 for p in procs):
        raise RuntimeError(f"{url}: a range download failed; rerun to resume")
    with open(dest, "wb") as out:
        for i in range(parts):
            part = dest.with_name(f"{dest.name}.part{i}")
            with open(part, "rb") as fh:
                while chunk := fh.read(1 << 24):
                    out.write(chunk)
    for i in range(parts):
        dest.with_name(f"{dest.name}.part{i}").unlink()


def download_atlas_files(
    registry_path: str | Path,
    atlas_id: str,
    file_names: list[str],
    base_url: str,
    data_root: str | Path,
    disk_limit_gb: float,
) -> list[dict]:
    """Download named files of one atlas into <data_root>/raw/<atlas_id>/.

    Refuses any file over the single-download limit and any download that would push
    usage under data_root past disk_limit_gb. Records url, sha256 and size in the
    registry entry. Returns the updated file records.
    """
    import hashlib
    import subprocess

    registry_path = Path(registry_path)
    outdir = Path(data_root) / "raw" / atlas_id
    outdir.mkdir(parents=True, exist_ok=True)
    done = []
    for name in file_names:
        url = f"{base_url.rstrip('/')}/{name}"
        size = remote_size_bytes(url)
        if size is None:
            raise RuntimeError(f"{url}: server did not report a size; check manually")
        if size > SINGLE_DOWNLOAD_LIMIT_BYTES:
            raise RuntimeError(f"{name} is {size / 1e9:.1f} GB, over the single-download limit")
        used = dir_size_bytes(Path(data_root))
        if used + size > disk_limit_gb * 10**9:
            raise RuntimeError(f"{name} would exceed the {disk_limit_gb} GB disk budget")
        dest = outdir / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not (dest.exists() and dest.stat().st_size == size):
            if size > 200 * 10**6:
                parallel_range_download(url, dest, size)
            else:
                subprocess.run(
                    ["curl", "-sSfL", "-C", "-", "--retry", "5", "-o", str(dest), url], check=True
                )
        if dest.stat().st_size != size:
            raise RuntimeError(f"{name}: got {dest.stat().st_size} bytes, expected {size}")
        h = hashlib.sha256()
        with open(dest, "rb") as fh:
            while chunk := fh.read(1 << 24):
                h.update(chunk)
        done.append({"name": name, "url": url, "size_bytes": size, "sha256": h.hexdigest(),
                     "downloaded": True})
    _record_files(registry_path, atlas_id, done)
    return done


def _record_files(registry_path: Path, atlas_id: str, records: list[dict]) -> None:
    """Merge file records into one registry entry under an exclusive file lock, so
    concurrent downloads of different atlases cannot overwrite each other."""
    import fcntl

    with open(registry_path.with_suffix(".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        text = registry_path.read_text()
        header = text.split("atlases:")[0]
        reg = yaml.safe_load(text)
        entry = next(a for a in reg["atlases"] if a["id"] == atlas_id)
        for r in records:
            rec = next((f for f in entry["files"] if f["name"] == r["name"]), None)
            if rec is None:
                entry["files"].append(dict(r))
            else:
                rec.update(r)
        registry_path.write_text(
            header + yaml.safe_dump(reg, sort_keys=False, allow_unicode=True, width=100)
        )


# --- Gene identifiers ---------------------------------------------------------


def symbol_to_ensembl(hgnc_tsv: str | Path) -> tuple[dict[str, str], dict[str, str]]:
    """Map gene symbols to Ensembl gene IDs using an HGNC complete-set snapshot.

    Approved symbols win. Previous symbols and aliases are used only when they point to
    exactly one approved gene and do not collide with an approved symbol; ambiguous
    names are left unmapped. Returns (mapping, how) where how[symbol] is
    "approved" | "previous" | "alias".
    """
    import csv
    from collections import defaultdict

    approved: dict[str, str] = {}
    prev: dict[str, set[str]] = defaultdict(set)
    alias: dict[str, set[str]] = defaultdict(set)
    with open(hgnc_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ens = r["ensembl_gene_id"]
            if r["status"] != "Approved" or not ens:
                continue
            approved[r["symbol"]] = ens
            for s in filter(None, r["prev_symbol"].split("|")):
                prev[s].add(ens)
            for s in filter(None, r["alias_symbol"].split("|")):
                alias[s].add(ens)
    mapping, how = dict(approved), {s: "approved" for s in approved}
    for table, label in ((prev, "previous"), (alias, "alias")):
        for s, ids in table.items():
            if s not in mapping and len(ids) == 1:
                mapping[s] = next(iter(ids))
                how[s] = label
    return mapping, how


def load_mtx(counts_dir: str | Path):
    """Read matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz (genes x cells)."""
    import gzip

    import scipy.io
    import scipy.sparse as sp

    d = Path(counts_dir)
    m = sp.csr_matrix(scipy.io.mmread(d / "matrix.mtx.gz"))
    genes = [line.rstrip("\n") for line in gzip.open(d / "features.tsv.gz", "rt")]
    cells = [line.rstrip("\n") for line in gzip.open(d / "barcodes.tsv.gz", "rt")]
    if m.shape != (len(genes), len(cells)):
        raise ValueError(f"matrix {m.shape} vs {len(genes)} genes, {len(cells)} cells")
    return m, genes, cells


def collapse_to_ensembl(matrix, symbols: list[str], mapping: dict[str, str]):
    """Map rows from symbols to Ensembl IDs, summing rows that map to the same ID.

    Unmapped rows are dropped. Returns (matrix, ensembl_ids, unmapped_symbols).
    """
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    ids = pd.Series([mapping.get(s) for s in symbols])
    keep = np.flatnonzero(ids.notna().to_numpy())
    codes, uniq = pd.factorize(ids.iloc[keep])
    agg = sp.csr_matrix(
        (np.ones(len(keep)), (codes, keep)), shape=(len(uniq), len(symbols))
    )
    unmapped = [symbols[i] for i in np.flatnonzero(ids.isna().to_numpy())]
    return agg @ sp.csr_matrix(matrix), list(uniq), unmapped


# --- Standardised per-atlas cells and pseudobulk ----------------------------------

STANDARD_OBS = ["cell_type_author", "donor", "keep"]


def save_cells(path: str | Path, X, gene_ids: list[str], obs, species: str) -> None:
    """Write a standardised cells.h5ad: X = raw integer counts (cells x genes, CSR),
    var index = species-native Ensembl gene IDs, obs = cell_type_author, donor, keep."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    X = sp.csr_matrix(X)
    if X.nnz and not np.allclose(X.data, np.round(X.data)):
        raise ValueError("counts must be integers")
    missing = [c for c in STANDARD_OBS if c not in obs.columns]
    if missing:
        raise ValueError(f"obs missing {missing}")
    obs = obs[STANDARD_OBS].copy()
    obs["cell_type_author"] = obs["cell_type_author"].astype(str)
    obs["donor"] = obs["donor"].astype(str)
    obs["keep"] = obs["keep"].astype(bool)
    obs.index = obs.index.astype(str)
    a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=pd.DataFrame(index=pd.Index(gene_ids, dtype=str)))
    a.uns["species"] = species
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path, compression="gzip")


def atlas_pseudobulk(cells_h5ad: str | Path, human_of: dict[str, str] | None = None):
    """Sum kept cells' counts per (cell_type_author, donor), genes as human Ensembl IDs.

    For mouse atlases pass `human_of` (mouse ID -> human ID, one-to-one); unmapped genes
    are dropped. Returns (sums genes x groups, groups table).
    """
    import anndata as ad

    from pnprank.pseudobulk import pseudobulk

    a = ad.read_h5ad(cells_h5ad)
    a = a[a.obs["keep"].to_numpy()]
    genes = list(a.var_names)
    X = a.X.T.tocsr()
    if a.uns["species"] != "human":
        if human_of is None:
            raise ValueError("non-human atlas needs an orthology mapping")
        X, genes, _ = collapse_to_ensembl(X, genes, human_of)
    sums, n = pseudobulk(
        X, genes, a.obs["cell_type_author"].reset_index(drop=True),
        a.obs["donor"].reset_index(drop=True), min_cells=1,
    )
    groups = n.reset_index().rename(columns={"cell_type": "cell_type_author"})
    groups.insert(0, "group", [f"g{i}" for i in range(len(groups))])
    sums.columns = groups["group"]
    return sums, groups


def write_pseudobulk(outdir: str | Path, sums, groups) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sums.to_csv(outdir / "pseudobulk_sums.tsv.gz", sep="\t", compression="gzip")
    groups.to_csv(outdir / "pseudobulk_groups.tsv", sep="\t", index=False)


def read_pseudobulk(outdir: str | Path):
    import pandas as pd

    outdir = Path(outdir)
    sums = pd.read_csv(outdir / "pseudobulk_sums.tsv.gz", sep="\t", index_col=0)
    groups = pd.read_csv(outdir / "pseudobulk_groups.tsv", sep="\t", dtype={"donor": str})
    return sums, groups
