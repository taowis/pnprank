"""Standardise one downloaded atlas to data/interim/<id>/cells.h5ad, then pseudobulk it.

Each atlas has its own reader below (formats differ); all write the same standard
(raw counts, species-native Ensembl IDs, author cell type, donor, keep flag) and then
go through the shared `atlas_pseudobulk`. The `keep` flag encodes which cells may be
used (healthy/control/naive only).

Usage: python workflow/scripts/standardize.py <atlas_id>
"""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from pnprank.atlases import (
    atlas_pseudobulk,
    collapse_to_ensembl,
    load_mtx,
    save_cells,
    symbol_to_ensembl,
    write_pseudobulk,
)

CFG = yaml.safe_load(Path("config/config.yaml").read_text())
PATHS = yaml.safe_load(Path("config/paths.yaml").read_text())
ROOT = Path(PATHS["data_root"])
HGNC = Path(CFG["truth"]["sources_dir"]) / "hgnc_complete_set.txt"


def raw(atlas_id: str) -> Path:
    return ROOT / "raw" / atlas_id


def interim(atlas_id: str) -> Path:
    return ROOT / "interim" / atlas_id


def human_symbols_to_ids(m_genes_x_cells, symbols):
    mapping, _ = symbol_to_ensembl(HGNC)
    return collapse_to_ensembl(m_genes_x_cells, symbols, mapping)


def mouse_symbol_map() -> dict[str, str]:
    """Mouse gene symbol -> Ensembl ID from the Ensembl 116 GTF; ambiguous symbols dropped."""
    t = pd.read_csv(ROOT / "external" / "annotation" / "mouse_genes_e116.tsv", sep="\t",
                    names=["id", "symbol", "biotype"])
    counts = t.symbol.value_counts()
    t = t[t.symbol.map(counts) == 1]
    return dict(zip(t.symbol, t.id))


def rds_to_mtx(rds: Path, outdir: Path) -> None:
    if not (outdir / "matrix.mtx.gz").exists():
        subprocess.run(["Rscript", "workflow/scripts/rds_to_mtx.R", str(rds), str(outdir)], check=True)


# --- readers -------------------------------------------------------------------


def yang2022(atlas_id: str, species: str):
    tag = "human" if species == "human" else "mouse"
    counts_dir = interim(atlas_id) / "counts"
    rds_to_mtx(next(raw(atlas_id).glob(f"*snRNA-seq_{tag}_raw_counts.RDS.gz")), counts_dir)
    m, symbols, cells = load_mtx(counts_dir)
    meta = pd.read_csv(next(raw(atlas_id).glob(f"*snRNA-seq_{tag}_barcode_meta.csv.gz")))
    meta = meta.set_index(meta.columns[0]).loc[cells]
    if species == "human":
        m, ids, _ = human_symbols_to_ids(m, symbols)
    else:
        m, ids, _ = collapse_to_ensembl(m, symbols, mouse_symbol_map())
    keep = np.ones(len(meta), bool)
    if "model" in meta.columns:  # mouse TG: headache models vs naive
        keep = meta["model"].astype(str).str.lower().isin(["naive", "control", "ctrl", "sham"]).to_numpy()
    obs = pd.DataFrame({"cell_type_author": meta["subtype"], "donor": meta["donor"].astype(str),
                        "keep": keep}, index=cells)
    return m.T, ids, obs


def cellxgene(atlas_id: str, file_glob: str, label_col: str, donor_col: str, keep_fn=None):
    """CELLxGENE h5ad: raw counts in .raw.X (or .X), var index = Ensembl IDs."""
    import anndata as ad

    a = ad.read_h5ad(next(raw(atlas_id).glob(file_glob)))
    X = a.raw.X if a.raw is not None else a.X
    var_names = list(a.raw.var_names if a.raw is not None else a.var_names)
    keep = np.ones(a.n_obs, bool) if keep_fn is None else np.asarray(keep_fn(a.obs), bool)
    obs = pd.DataFrame({"cell_type_author": a.obs[label_col].astype(str),
                        "donor": a.obs[donor_col].astype(str), "keep": keep}, index=a.obs_names)
    return sp.csr_matrix(X), var_names, obs


def sparc465(atlas_id: str):
    """Multiome mtx with Gene Expression rows first, then ATAC Peaks; keep genes only.

    The matrix is filtered while streaming (awk) because the peak rows make it too large
    to load. Libraries are used as the replicate unit (library-to-donor mapping is not
    in meta.data.csv).
    """
    base = raw(atlas_id) / "primary" / "filtered_feature_bc_matrix"
    out = interim(atlas_id) / "gex"
    out.mkdir(parents=True, exist_ok=True)
    feats = pd.read_csv(base / "features.tsv.gz", sep="\t", header=None)
    gex = feats[feats[2] == "Gene Expression"]
    n_gex = len(gex)
    if not (gex.index == np.arange(n_gex)).all():
        raise ValueError("expected Gene Expression rows first")
    if not (out / "matrix.mtx.gz").exists():
        cmd = (f"gzcat {sh(base / 'matrix.mtx.gz')} | awk 'BEGIN{{h=0}} /^%/{{next}} "
               f"h==0{{h=1; ncol=$2; next}} $1<={n_gex}{{print; n++}} "
               f"END{{print n, ncol > \"/dev/stderr\"}}' > {sh(out / 'body.txt')}")
        res = subprocess.run(["bash", "-c", cmd], check=True, capture_output=True, text=True)
        nnz, ncol = res.stderr.split()
        with gzip.open(out / "matrix.mtx.gz", "wt") as fh:
            fh.write("%%MatrixMarket matrix coordinate integer general\n")
            fh.write(f"{n_gex} {ncol} {nnz}\n")
            with open(out / "body.txt") as body:
                for line in body:
                    fh.write(line)
        (out / "body.txt").unlink()
        gex[1].to_csv(out / "features.tsv.gz", index=False, header=False)
        subprocess.run(["cp", str(base / "barcodes.tsv.gz"), str(out / "barcodes.tsv.gz")], check=True)
    import scipy.io
    m = sp.csr_matrix(scipy.io.mmread(out / "matrix.mtx.gz"))
    cells = [x.strip() for x in gzip.open(out / "barcodes.tsv.gz", "rt")]
    meta = pd.read_csv(raw(atlas_id) / "primary" / "meta.data.csv").set_index("cell_name")
    keep_cells = [c for c in cells if c in meta.index]
    idx = [cells.index(c) for c in keep_cells] if len(keep_cells) < len(cells) else list(range(len(cells)))
    m = m[:, idx]
    meta = meta.loc[keep_cells]
    obs = pd.DataFrame({"cell_type_author": meta["Celltypes"], "donor": meta["sample id"],
                        "keep": True}, index=keep_cells)
    return m.T, list(gex[0]), obs


def sh(p) -> str:
    import shlex
    return shlex.quote(str(p))


def heming2025(atlas_id: str):
    """Per-sample CellBender h5 files inside GSE285983_RAW.tar; controls (level2 == CTRL) only.

    Counts are CellBender ambient-corrected integer counts, as deposited.
    """
    import tarfile

    import scanpy as sc

    meta = pd.read_csv(raw(atlas_id) / "GSE285983_metadata_all.csv.gz", low_memory=False)
    ctrl = meta[meta["level2"] == "CTRL"]
    samples = sorted(ctrl["sample"].unique())
    out = interim(atlas_id) / "h5"
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw(atlas_id) / "GSE285983_RAW.tar") as tf:
        members = {m.name.split("_")[1]: m for m in tf.getmembers()}
        for smp in samples:
            if not (out / members[smp].name).exists():
                tf.extract(members[smp], out)
    mats, ids, obs_parts = [], None, []
    for smp in samples:
        a = sc.read_10x_h5(out / members[smp].name, gex_only=True)
        a.obs_names = [f"{smp}_{b}" for b in a.obs_names]
        m = ctrl[ctrl["sample"] == smp].set_index("barcode")
        a = a[a.obs_names.isin(m.index)]
        gene_ids = list(a.var["gene_ids"])
        if ids is None:
            ids = gene_ids
        elif gene_ids != ids:
            raise ValueError("gene order differs between samples")
        mats.append(sp.csr_matrix(a.X))
        obs_parts.append(pd.DataFrame({"cell_type_author": m.loc[a.obs_names, "cluster"],
                                       "donor": smp, "keep": True}, index=a.obs_names))
    return sp.vstack(mats).tocsr(), ids, pd.concat(obs_parts)


READERS = {
    "yang2022_human_tg": lambda a: yang2022(a, "human"),
    "yang2022_mouse_tg": lambda a: yang2022(a, "mouse"),
    # pooled animals per replicate; donor_id identifies the sample (region x sex x rep)
    "alkaslasi2021_mouse_mn": lambda a: cellxgene(a, "*.h5ad", "cluster_assignment_fine", "donor_id"),
    # "Innervated only" dataset: all normal
    "ham2025_mouse_nmj": lambda a: cellxgene(a, "*.h5ad", "author_cell_type", "donor_id"),
    # nuclei only (cells and nuclei are mixed; nuclei capture myonuclei); all donors normal
    "kedlian2024_human_muscle": lambda a: cellxgene(
        a, "*.h5ad", "annotation_level1", "donor_id",
        keep_fn=lambda o: (o["suspension_type"] == "nucleus").to_numpy()),
    "sparc465_human_sympathetic": sparc465,
    "heming2025_human_sural": heming2025,
}
SPECIES = {"yang2022_human_tg": "human", "yang2022_mouse_tg": "mouse",
           "alkaslasi2021_mouse_mn": "mouse", "ham2025_mouse_nmj": "mouse",
           "kedlian2024_human_muscle": "human", "sparc465_human_sympathetic": "human",
           "heming2025_human_sural": "human"}


def main(atlas_id: str) -> None:
    X, ids, obs = READERS[atlas_id](atlas_id)
    species = SPECIES[atlas_id]
    path = interim(atlas_id) / "cells.h5ad"
    save_cells(path, X, ids, obs, species)
    human_of = None
    if species == "mouse":
        o = pd.read_csv(CFG["reference"]["orthology"]["table"], sep="\t")
        human_of = dict(zip(o.mouse_ensembl_gene_id, o.human_ensembl_gene_id))
    sums, groups = atlas_pseudobulk(path, human_of)
    write_pseudobulk(interim(atlas_id), sums, groups)
    print(atlas_id, "cells kept", int(obs.keep.sum()), "of", len(obs), "| genes", sums.shape[0],
          "| groups", len(groups))
    print(groups.groupby("cell_type_author")["n_cells"].sum().to_string())


if __name__ == "__main__":
    main(sys.argv[1])
