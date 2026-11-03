"""Assemble the gene-by-feature matrix for a given cutoff year.

Feature families:
- atlas features (label-free): per atlas x harmonised cell type, log1p mean CPM and
  fraction specificity; computed once, independent of any cutoff.
- static features (label-free): gnomAD v4.1 constraint, GTEx v10 bulk expression.
- label-derived features: co-expression with genes known before the cutoff
  (`coexpression.py`); recomputed for every cutoff.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pnprank.atlases import read_pseudobulk
from pnprank.coexpression import known_gene_similarity
from pnprank.pseudobulk import cpm
from pnprank.specificity import fraction_specificity, tau


def load_celltype_map(path: str | Path) -> pd.DataFrame:
    """atlas_id, cell_type_author -> cell_type (harmonised) or 'exclude'."""
    m = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    need = {"atlas_id", "cell_type_author", "cell_type"}
    if not need <= set(m.columns):
        raise ValueError(f"cell-type map needs columns {need}")
    return m


def harmonised_profile(
    sums: pd.DataFrame, groups: pd.DataFrame, celltype_map: pd.DataFrame, atlas_id: str,
    min_cells: int, min_donors: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Re-aggregate author groups to harmonised cell types, then CPM and donor mean.

    Author labels missing from the map raise an error (every label must be decided).
    Returns (mean CPM genes x harmonised types, kept groups table).
    """
    m = celltype_map[celltype_map.atlas_id == atlas_id]
    lookup = dict(zip(m.cell_type_author, m.cell_type))
    unknown = sorted(set(groups.cell_type_author) - set(lookup))
    if unknown:
        raise ValueError(f"{atlas_id}: unmapped author labels {unknown}")
    g = groups.assign(cell_type=groups.cell_type_author.map(lookup))
    g = g[g.cell_type != "exclude"]
    key = g.cell_type + "||" + g.donor
    agg = sums[g.group].T.groupby(key.to_numpy()).sum().T
    n = g.groupby(key.to_numpy())["n_cells"].sum()
    agg = agg.loc[:, n[n >= min_cells].index]
    expr = cpm(agg)
    ct = pd.Series([c.split("||")[0] for c in expr.columns], index=expr.columns)
    n_donors = ct.value_counts()
    ok = n_donors[n_donors >= min_donors].index
    mean = expr.T.groupby(ct.to_numpy()).mean().T
    kept = pd.DataFrame({"cell_type": ct, "n_cells": n.reindex(expr.columns)})
    return mean[[c for c in mean.columns if c in ok]], kept[kept.cell_type.isin(ok)]


def atlas_profiles(
    atlas_ids: list[str], interim: Path, celltype_map: pd.DataFrame, min_cells: int, min_donors: int
) -> dict[str, pd.DataFrame]:
    out = {}
    for a in atlas_ids:
        sums, groups = read_pseudobulk(interim / a)
        prof, _ = harmonised_profile(sums, groups, celltype_map, a, min_cells, min_donors)
        out[a] = prof
    return out


def profiles_from_static(X_static: pd.DataFrame, atlas_ids: list[str]) -> dict[str, pd.DataFrame]:
    """Rebuild donor-mean CPM profiles (genes x cell types) from the static feature matrix.

    The `<atlas>|<cell type>|logcpm` columns are log1p of the donor-mean CPM of each
    harmonised cell type (`atlas_features`), so expm1 recovers the profiles that
    `harmonised_profile` returns, restricted to the gene universe. Genes absent from an
    atlas (all NaN) are dropped, as they are absent from the original profile. Used when
    the per-atlas pseudobulk files under interim/ have been removed to save disk space.
    """
    out = {}
    for a in atlas_ids:
        cols = [c for c in X_static.columns if c.startswith(a + "|") and c.endswith("|logcpm")]
        if not cols:
            raise ValueError(f"no logcpm columns for atlas {a}")
        prof = np.expm1(X_static[cols])
        prof.columns = [c.split("|")[1] for c in cols]
        out[a] = prof[prof.notna().all(axis=1)]
    return out


def load_profiles(atlas_ids: list[str], interim: Path, celltype_map: pd.DataFrame,
                  min_cells: int, min_donors: int, static_path: Path) -> dict[str, pd.DataFrame]:
    """Atlas profiles from interim/ pseudobulk if present, else from the static features."""
    if all((Path(interim) / a).exists() for a in atlas_ids):
        return atlas_profiles(atlas_ids, Path(interim), celltype_map, min_cells, min_donors)
    print(f"interim pseudobulk not found; profiles rebuilt from {static_path}", flush=True)
    return profiles_from_static(pd.read_csv(static_path, sep="\t", index_col=0), atlas_ids)


def atlas_features(profiles: dict[str, pd.DataFrame], universe: list[str]) -> pd.DataFrame:
    """log1p CPM and fraction specificity per atlas x cell type; absent genes are NaN.

    A gene present in the atlas but never detected gets 0 expression and NaN
    specificity.
    """
    cols = {}
    for a, prof in profiles.items():
        spec = fraction_specificity(prof)
        for c in prof.columns:
            cols[f"{a}|{c}|logcpm"] = np.log1p(prof[c]).reindex(universe)
            cols[f"{a}|{c}|spec"] = spec[c].reindex(universe)
        cols[f"{a}|tau"] = tau(prof).reindex(universe)
    return pd.DataFrame(cols, index=universe)


def gnomad_features(path: str | Path, universe: list[str]) -> pd.DataFrame:
    """LOEUF, pLI and missense Z from the MANE Select transcript (canonical if none),
    matched by Ensembl gene ID."""
    g = pd.read_csv(path, sep="\t", low_memory=False)
    g = g[g["gene_id"].astype(str).str.startswith("ENSG")]
    g = g.sort_values(["mane_select", "canonical"], ascending=False).drop_duplicates("gene_id")
    g = g.set_index("gene_id")
    return pd.DataFrame(
        {
            "gnomad|loeuf": g["lof.oe_ci.upper"],
            "gnomad|pli": g["lof.pLI"],
            "gnomad|mis_z": g["mis.z_score"],
        }
    ).reindex(universe)


def gtex_features(path: str | Path, universe: list[str], tissues: list[str]) -> pd.DataFrame:
    """log1p median TPM for selected tissues, mean over all tissues, and tau across tissues."""
    t = pd.read_csv(path, sep="\t", skiprows=2)
    t["gene"] = t["Name"].str.split(".").str[0]
    t = t.drop(columns=["Name", "Description"]).groupby("gene").sum()
    out = {f"gtex|{s}|logtpm": np.log1p(t[s]) for s in tissues}
    out["gtex|mean_all|logtpm"] = np.log1p(t.mean(axis=1))
    out["gtex|tau"] = tau(t)
    return pd.DataFrame(out).reindex(universe)


def known_before(truth: pd.DataFrame, cutoff_year: int, tiers: set[str]) -> set[str]:
    """Ensembl IDs of truth genes in `tiers` discovered strictly before cutoff_year."""
    t = truth[truth.tier.isin(tiers)]
    return set(t.ensembl_gene_id[t.discovery_year.astype(int) < cutoff_year])


def label_features(
    profiles: dict[str, pd.DataFrame], universe: list[str], known: set[str]
) -> pd.DataFrame:
    """Co-expression with known genes, one column per atlas."""
    cols = {}
    for a, prof in profiles.items():
        if prof.shape[1] >= 3:
            cols[f"{a}|coexp_known"] = known_gene_similarity(np.log1p(prof), known).reindex(universe)
    return pd.DataFrame(cols, index=universe)
