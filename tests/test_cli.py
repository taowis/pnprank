from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pnprank.cli import main
from pnprank.ranking import COLUMNS, build_ranking, resolve, top_genes

TRUTH = Path(__file__).parents[1] / "data" / "truth" / "neuropathy_genes.tsv"


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "pnprank" in capsys.readouterr().out


def test_validate_truth():
    assert main(["validate-truth", str(TRUTH)]) == 0


# --- ranking lookups -------------------------------------------------------------


def _table(**cols):
    return pd.DataFrame(cols).reindex(columns=COLUMNS)


def _ranking(tmp_path):
    tab = _table(
        rank=[1, 2, 3, 4, None],
        hgnc_symbol=["GENEA", "GENEB", "GENEC", "GENED", "GENEE"],
        ensembl_gene_id=["ENSG01", "ENSG02", "ENSG03", "ENSG04", "ENSG05"],
        status=["unlabelled", "listed_not_in_truth", "known_tier2", "unlabelled", "known_tier1"],
        percentile=[0.99, 0.95, 0.90, 0.50, None],
        percentile_expr_matched=[0.6, 0.97, 0.9, 0.8, None],
        aliases=["OLDA", "SHARED", "SHARED|OLDC", "", ""],
    )
    path = tmp_path / "ranking.tsv.gz"
    tab.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def test_resolve_by_id_symbol_and_alias():
    tab = _table(hgnc_symbol=["GENEA", "GENEB", "GENEC"],
                 ensembl_gene_id=["ENSG01", "ENSG02", "ENSG03"],
                 aliases=["OLDA", "SHARED", "SHARED|OLDC"])
    out, problems = resolve(tab, ["genea", "ENSG02.7", "OLDC", "SHARED", "NOPE", ""])
    assert list(out.hgnc_symbol) == ["GENEA", "GENEB", "GENEC"]
    assert list(out["query"]) == ["genea", "ENSG02.7", "OLDC"]
    assert problems["SHARED"].startswith("ambiguous")
    assert "NOPE" in problems and "" not in problems


def test_approved_symbol_wins_over_alias():
    tab = _table(hgnc_symbol=["GENEA", "GENEB"], ensembl_gene_id=["E1", "E2"],
                 aliases=["GENEB", ""])
    out, _ = resolve(tab, ["GENEB"])
    assert list(out.ensembl_gene_id) == ["E2"]


def test_top_excludes_known_by_default(tmp_path):
    from pnprank.ranking import load_ranking

    tab = load_ranking(_ranking(tmp_path))
    assert list(top_genes(tab, 10).hgnc_symbol) == ["GENEA", "GENEB", "GENED"]
    assert list(top_genes(tab, 2, include_known=True).hgnc_symbol) == ["GENEA", "GENEB"]
    assert list(top_genes(tab, 1, by="percentile_expr_matched").hgnc_symbol) == ["GENEB"]


def test_cli_query_and_top(tmp_path, capsys):
    path = str(_ranking(tmp_path))
    assert main(["query", "GENEA", "olda", "--ranking", path, "--tsv"]) == 0
    out = pd.read_csv(pd.io.common.StringIO(capsys.readouterr().out), sep="\t")
    assert list(out.ensembl_gene_id) == ["ENSG01", "ENSG01"]
    assert out["rank"].dtype.kind == "i"   # printed as 1, not 1.0
    assert main(["query", "NOPE", "--ranking", path]) == 1
    assert "not found" in capsys.readouterr().err
    genes = tmp_path / "genes.txt"
    genes.write_text("# header\nGENED\tvariant info\n")
    assert main(["query", "--file", str(genes), "--ranking", path]) == 0
    assert "GENED" in capsys.readouterr().out
    assert main(["top", "-n", "1", "--ranking", path]) == 0
    assert "GENEA" in capsys.readouterr().out
    assert main(["top", "--ranking", str(tmp_path / "missing.tsv.gz")]) == 2


# --- building the ranking ---------------------------------------------------------


def _synthetic(n=150, n_pos=15, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"ENSG{i:04d}" for i in range(n)]
    y = np.zeros(n)
    y[:n_pos] = 1
    X = pd.DataFrame({
        "gtex|mean_all|logtpm": rng.normal(size=n) + y,
        "gnomad|loeuf": rng.normal(size=n),
        "a1|t1|spec": rng.uniform(size=n) + 0.5 * y,
        "a1|t2|spec": rng.uniform(size=n),
        "a1|tau": rng.uniform(size=n),
    }, index=genes)
    truth = pd.DataFrame({
        "ensembl_gene_id": genes[:n_pos] + genes[n_pos:n_pos + 5],
        "discovery_year": ["2000"] * n_pos + ["2010"] * 5,
        "tier": ["1"] * n_pos + ["2"] * 5,
    })
    return X, truth


def _build(X, truth, seeds=(1, 2)):
    return build_ranking(X, truth, {X.index[30]}, {g: g.replace("ENSG", "S") for g in X.index},
                         {}, None, ["gtex|mean_all|logtpm", "gnomad|loeuf"],
                         X["gtex|mean_all|logtpm"], list(seeds), 3, 2, 4)


def test_build_ranking_statuses_and_scores():
    X, truth = _synthetic()
    tab, meta = _build(X, truth)
    assert list(tab.columns) == COLUMNS and len(tab) == len(X)
    st = tab.set_index("ensembl_gene_id").status
    assert (st == "known_tier1").sum() == 15 and (st == "known_tier2").sum() == 5
    assert st[X.index[30]] == "listed_not_in_truth"
    known1 = tab[tab.status == "known_tier1"]
    assert known1.percentile.isna().all() and known1["rank"].isna().all()
    rest = tab[tab.status != "known_tier1"]
    # every other gene, including held-out tier 2, gets a percentile and a unique rank
    assert rest.percentile.notna().all() and rest.percentile_expr_matched.notna().all()
    assert sorted(rest["rank"]) == list(range(1, len(rest) + 1))
    assert (rest.percentile_seed_min <= rest.percentile + 1e-9).all()
    assert (rest.percentile <= rest.percentile_seed_max + 1e-9).all()
    assert tab["rank"].dropna().is_monotonic_increasing
    assert meta["n_positives_tier1"] == 15 and meta["n_held_out_tier2"] == 5
    assert set(meta["chosen_C_by_seed"]) == {"1", "2"}


def test_build_ranking_is_reproducible_and_learns_from_tier1():
    X, truth = _synthetic()
    a, _ = _build(X, truth)
    b, _ = _build(X, truth)
    pd.testing.assert_frame_equal(a, b)
    # an unlabelled gene given the positives' feature profile should rank near the top
    X2 = X.copy()
    X2.loc["ENSG0100", ["gtex|mean_all|logtpm", "a1|t1|spec"]] = [3.0, 1.5]
    c, _ = _build(X2, truth)
    assert c.set_index("ensembl_gene_id").percentile["ENSG0100"] > 0.9
