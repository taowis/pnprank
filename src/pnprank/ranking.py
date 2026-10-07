"""Released gene ranking (Stage 6) and lookups used by the `pnprank` CLI.

The ranking applies the procedure validated in the Stage 5 back-test to today's labels:
tier-1 genes are positives, tier-2 genes are left out of training, and every other gene
is scored out of fold by a model that never saw it as a training negative. Scores are
reported as percentiles within the fold model that scored the gene, because fold models
tune their regularisation separately (docs/decisions.md, Stage 5).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pnprank.backtest import fit_fold_models, training_labels, within_fold_percentile

STATUS_ORDER = ["known_tier1", "known_tier2", "listed_not_in_truth", "unlabelled"]
COLUMNS = [
    "rank", "hgnc_symbol", "ensembl_gene_id", "status", "percentile", "percentile_seed_min",
    "percentile_seed_max", "percentile_expr_matched", "percentile_non_atlas",
    "gtex_mean_decile", "top_cell_type", "n_atlases_detected", "aliases",
]


def _expr_bins(expr: pd.Series, genes, n_bins: int) -> pd.Series:
    """Quantile bin (0..n_bins-1) of `expr` over `genes`, as in within_fold_percentile."""
    e = expr.reindex(genes).dropna()
    return pd.qcut(e.rank(method="first"), n_bins, labels=False)


def _pct_against(ref: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Fraction of `ref` at or below each value (percentile of genes not in `ref`)."""
    ref = np.sort(ref[~np.isnan(ref)])
    return np.searchsorted(ref, values, side="right") / len(ref)


def _score_one_seed(X_static, y, held_out, expr, label_fn, seed, folds, inner_folds, n_bins,
                    n_jobs):
    """Within-fold percentiles for 0-genes, and fold-averaged percentiles for `held_out`
    genes (in no training set). Returns (pct, pct_expr_matched, best params)."""
    X, models, params, fold = fit_fold_models(X_static, y, "logistic", seed, folds,
                                              inner_folds, label_fn, n_jobs)
    zeros = fold.index
    raw = pd.concat([pd.Series(m.predict_proba(X.loc[zeros[fold == k]])[:, 1],
                               index=zeros[fold == k]) for k, m in enumerate(models)])
    raw = raw.reindex(zeros)
    pct = within_fold_percentile(raw, fold)
    pct_e = within_fold_percentile(raw, fold, expr, n_bins) if expr is not None else None

    if len(held_out):
        bins = _expr_bins(expr, zeros, n_bins) if expr is not None else None
        edges = None
        if bins is not None:
            e0 = expr.reindex(bins.index)
            edges = [e0[bins == b].max() for b in range(n_bins - 1)]
        hp, hpe = [], []
        for k, m in enumerate(models):
            s = m.predict_proba(X.loc[held_out])[:, 1]
            ref = raw[fold == k]
            hp.append(_pct_against(ref.to_numpy(), s))
            if bins is not None:
                he = expr.reindex(held_out).to_numpy()
                hb = np.searchsorted(edges, he, side="left")
                col = np.full(len(held_out), np.nan)
                for i, b in enumerate(hb):
                    if np.isnan(he[i]):
                        continue
                    r = ref[bins.reindex(ref.index) == b].to_numpy()
                    col[i] = _pct_against(r, s[i:i + 1])[0]
                hpe.append(col)
        pct = pd.concat([pct, pd.Series(np.mean(hp, axis=0), index=held_out)])
        if pct_e is not None:
            pct_e = pd.concat([pct_e, pd.Series(np.nanmean(hpe, axis=0), index=held_out)])
    return pct, pct_e, params


def top_cell_type(X_static: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Most specific atlas cell type per gene ("atlas|type (spec)") and number of atlases
    in which the gene was measured. Descriptive only; not a model output."""
    spec = X_static[[c for c in X_static.columns if c.endswith("|spec")]]
    best = spec.fillna(-np.inf).idxmax(axis=1)
    val = spec.max(axis=1)
    label = best.str.removesuffix("|spec") + " (" + val.round(2).astype(str) + ")"
    label[spec.isna().all(axis=1)] = np.nan
    atlas_tau = [c for c in X_static.columns if c.endswith("|tau") and not c.startswith("gtex|")]
    return label, X_static[atlas_tau].notna().sum(axis=1)


def build_ranking(X_static: pd.DataFrame, truth: pd.DataFrame, listed: set[str],
                  symbols: dict[str, str], aliases: dict[str, str], label_fn,
                  non_atlas_cols: list[str], expr: pd.Series | None, seeds: list[int],
                  folds: int, inner_folds: int, n_bins: int, n_jobs: int = 1):
    """Rank every gene in X_static with today's tier-1 labels.

    `seeds[0]` gives the reported percentiles (the back-tested procedure); the other seeds
    only measure how much the percentile depends on the random fold split and inner-CV
    choice of C. Returns (table, metadata dict).
    """
    universe = list(X_static.index)
    release = int(truth.discovery_year.astype(int).max()) + 1
    y = training_labels(universe, truth, release, {"1"})
    held_out = [g for g in universe if g not in y.index]   # tier-2 genes

    per_seed, params = [], {}
    for s in seeds:
        pct, pct_e, p = _score_one_seed(X_static, y, held_out, expr, label_fn, s, folds,
                                        inner_folds, n_bins, n_jobs)
        per_seed.append(pct)
        params[s] = [d["logisticregression__C"] for d in p]
        if s == seeds[0]:
            primary, primary_e = pct, pct_e
    non, _, p_non = _score_one_seed(X_static[non_atlas_cols], y, held_out, None, None,
                                    seeds[0], folds, inner_folds, n_bins, n_jobs)
    all_pct = pd.concat(per_seed, axis=1)

    t1 = set(truth.ensembl_gene_id[truth.tier == "1"])
    t2 = set(truth.ensembl_gene_id[truth.tier == "2"])
    status = pd.Series("unlabelled", index=universe)
    status[status.index.isin(listed)] = "listed_not_in_truth"
    status[status.index.isin(t2)] = "known_tier2"
    status[status.index.isin(t1)] = "known_tier1"
    top, n_atlas = top_cell_type(X_static)
    decile = (_expr_bins(expr, universe, 10) + 1) if expr is not None else None

    tab = pd.DataFrame({
        "hgnc_symbol": [symbols.get(g, "") for g in universe],
        "ensembl_gene_id": universe,
        "status": status.to_numpy(),
        "percentile": primary.reindex(universe).to_numpy(),
        "percentile_seed_min": all_pct.min(axis=1).reindex(universe).to_numpy(),
        "percentile_seed_max": all_pct.max(axis=1).reindex(universe).to_numpy(),
        "percentile_expr_matched": (primary_e.reindex(universe).to_numpy()
                                    if primary_e is not None else np.nan),
        "percentile_non_atlas": non.reindex(universe).to_numpy(),
        "gtex_mean_decile": decile.reindex(universe).astype("Int64").to_numpy()
        if decile is not None else pd.NA,
        "top_cell_type": top.reindex(universe).to_numpy(),
        "n_atlases_detected": n_atlas.reindex(universe).to_numpy(),
        "aliases": [aliases.get(g, "") for g in universe],
    })
    # rank genes that were scored (not tier-1 training positives); ties broken by the
    # mean percentile over seeds, then by Ensembl ID so the order is deterministic
    tab["_mean"] = all_pct.mean(axis=1).reindex(universe).to_numpy()
    scored = tab.percentile.notna()
    order = tab[scored].sort_values(["percentile", "_mean", "ensembl_gene_id"],
                                    ascending=[False, False, True]).index
    tab["rank"] = pd.Series(pd.NA, index=tab.index, dtype="Int64")
    tab.loc[order, "rank"] = np.arange(1, len(order) + 1)
    tab["_s"] = tab.status.map({s: i for i, s in enumerate(STATUS_ORDER)})
    tab = tab.sort_values(["rank", "_s", "hgnc_symbol"], na_position="last")
    for c in ["percentile", "percentile_seed_min", "percentile_seed_max",
              "percentile_expr_matched", "percentile_non_atlas"]:
        tab[c] = tab[c].round(4)
    tab = tab[COLUMNS].reset_index(drop=True)

    unl = [g for g in y.index[y == 0]]
    top100 = {s: set(per_seed[i].reindex(unl).nlargest(100).index) for i, s in enumerate(seeds)}
    meta = {
        "release_label_year": release,
        "n_positives_tier1": int(y.sum()),
        "n_held_out_tier2": len(held_out),
        "n_scored_out_of_fold": int((y == 0).sum()),
        "seeds": seeds,
        "scoring_folds": folds,
        "chosen_C_by_seed": {str(k): v for k, v in params.items()},
        "chosen_C_non_atlas": [d["logisticregression__C"] for d in p_non],
        "seed_stability": {
            "spearman_vs_primary": {
                str(s): round(float(all_pct.iloc[:, 0].corr(all_pct.iloc[:, i],
                                                             method="spearman")), 4)
                for i, s in enumerate(seeds[1:], 1)},
            "top100_overlap_vs_primary": {str(s): len(top100[seeds[0]] & top100[s])
                                          for s in seeds[1:]},
            "median_seed_range": round(float(
                (all_pct.max(axis=1) - all_pct.min(axis=1)).median()), 4),
        },
    }
    return tab, meta


# --- lookups ------------------------------------------------------------------------


def load_ranking(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"rank": "Int64", "gtex_mean_decile": "Int64",
                                             "aliases": str, "top_cell_type": str})


def resolve(tab: pd.DataFrame, queries: list[str]) -> tuple[pd.DataFrame, dict[str, str]]:
    """Match queries to rows by Ensembl ID (version ignored), then approved symbol, then
    previous/alias symbols (all case-insensitive). Returns (matched rows with a `query`
    column, {query: problem} for unmatched or ambiguous queries)."""
    by_id = {g.upper(): i for i, g in enumerate(tab.ensembl_gene_id)}
    by_sym = {s.upper(): i for i, s in enumerate(tab.hgnc_symbol.fillna("")) if s}
    by_alias: dict[str, set[int]] = {}
    for i, a in enumerate(tab.aliases.fillna("")):
        for s in filter(None, a.split("|")):
            by_alias.setdefault(s.upper(), set()).add(i)
    rows, problems = [], {}
    for q in queries:
        k = q.strip().upper()
        if not k:
            continue
        i = by_id.get(k.split(".")[0]) if k.startswith("ENSG") else by_sym.get(k)
        if i is None and k in by_alias:
            hits = by_alias[k]
            if len(hits) > 1:
                problems[q] = "ambiguous alias: " + ", ".join(
                    sorted(tab.hgnc_symbol.iloc[j] for j in hits))
                continue
            i = next(iter(hits))
        if i is None:
            problems[q] = "not found (not a protein-coding gene in the ranking?)"
            continue
        rows.append(tab.iloc[i].to_dict() | {"query": q.strip()})
    out = pd.DataFrame(rows, columns=["query"] + COLUMNS)
    return out, problems


def top_genes(tab: pd.DataFrame, n: int, include_known: bool = False,
              by: str = "percentile") -> pd.DataFrame:
    """The `n` highest-ranked genes by `by`, excluding known tier-1/2 genes by default."""
    d = tab[tab[by].notna()]
    if not include_known:
        d = d[~d.status.isin(["known_tier1", "known_tier2"])]
    if by == "percentile":
        return d.sort_values("rank").head(n)
    return d.sort_values([by, "percentile"], ascending=False).head(n)
