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
        "GENE1\tENSG00000000001\tCMT\t1\tno\tAD\t1999\t1\thttps://example.org\t",
        "GENE1\tnot_an_id\tXYZ\t3\tmaybe\tAD\t99\t\t\t",
    ]
    bad.write_text("\n".join(rows) + "\n")
    problems = validate_truth(bad)
    joined = " ".join(problems)
    for expected in ["duplicate", "ensembl_gene_id", "category", "tier", "also_hsp", "discovery_year", "pmid", "source_url"]:
        assert expected in joined


def test_discovery_year_takes_earliest():
    from pnprank.truth import discovery_year

    assert discovery_year({"pubdate": "2010 Jan", "epubdate": "2009 Dec 30"}) == 2009
    assert discovery_year({"pubdate": "2010 Jan", "epubdate": ""}) == 2010


def test_mentions_gene_whole_word():
    from pnprank.truth import mentions_gene

    rec = {"title": "Mutations in GARS cause CMT2D", "abstract": ""}
    assert mentions_gene(rec, ["GARS1", "GARS"])
    assert not mentions_gene(rec, ["ARS"])
    rec = {"title": "Mutation in the Mitofusin 2 gene", "abstract": ""}
    assert mentions_gene(rec, ["MFN2"], ("mitofusin 2",))


def test_clingen_parser(tmp_path):
    from pnprank.truth import clingen_genes

    csv_text = "\n".join(
        [
            '"CLINGEN GENE DISEASE VALIDITY CURATIONS"',
            (
                '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
                '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"'
            ),
            '"+++","+++"',
            (
                '"A","HGNC:1","cmt","M:1","AD","SOP10","Definitive","u","d",'
                '"Charcot-Marie-Tooth Disease Gene Curation Expert Panel"'
            ),
            (
                '"B","HGNC:2","cmt","M:2","AD","SOP10","Limited","u","d",'
                '"Charcot-Marie-Tooth Disease Gene Curation Expert Panel"'
            ),
            '"C","HGNC:3","x","M:3","AD","SOP10","Definitive","u","d","Other Panel"',
        ]
    )
    p = tmp_path / "c.csv"
    p.write_text(csv_text + "\n")
    assert set(clingen_genes(p, "Charcot-Marie-Tooth")) == {"HGNC:1"}
