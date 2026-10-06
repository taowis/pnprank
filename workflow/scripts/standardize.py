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
}
SPECIES = {"yang2022_human_tg": "human", "yang2022_mouse_tg": "mouse",
           "alkaslasi2021_mouse_mn": "mouse", "ham2025_mouse_nmj": "mouse",
           "kedlian2024_human_muscle": "human"}


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
