"""No gene discovered at or after the cutoff may influence any feature.

Must pass before any back-test result is reported. Real checks are added in
Stage 4 with features.py; until then only the truth-set split is tested.
"""

from pnprank.truth import known_before


def test_known_before_is_strict():
    rows = [
        {"ensembl_gene_id": "ENSG_A", "discovery_year": "2009"},
        {"ensembl_gene_id": "ENSG_B", "discovery_year": "2010"},
        {"ensembl_gene_id": "ENSG_C", "discovery_year": "2011"},
    ]
    assert known_before(rows, 2010) == {"ENSG_A"}
