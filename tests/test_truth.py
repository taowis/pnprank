from pathlib import Path

from pnprank.truth import REQUIRED_COLUMNS, validate_truth

TRUTH = Path(__file__).parents[1] / "data" / "truth" / "neuropathy_genes.tsv"


def test_truth_file_valid():
    assert validate_truth(TRUTH) == []


def test_truth_header():
    header = TRUTH.read_text().splitlines()[0].split("\t")
    assert header == REQUIRED_COLUMNS


def test_validator_catches_problems(tmp_path):
    bad = tmp_path / "bad.tsv"
    rows = [
        "\t".join(REQUIRED_COLUMNS),
        "GENE1\tENSG00000000001\tCMT\tAD\t1999\t1\thttps://example.org\t",
        "GENE1\tnot_an_id\tXYZ\tAD\t99\t\t\t",
    ]
    bad.write_text("\n".join(rows) + "\n")
    problems = validate_truth(bad)
    joined = " ".join(problems)
    for expected in ["duplicate", "ensembl_gene_id", "category", "discovery_year", "pmid", "source_url"]:
        assert expected in joined
