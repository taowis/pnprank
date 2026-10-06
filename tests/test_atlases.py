from pathlib import Path

from pnprank.atlases import validate_registry

REGISTRY = Path(__file__).parents[1] / "config" / "atlases.yaml"


def test_registry_valid():
    assert validate_registry(REGISTRY) == []


def test_validator_catches_problems(tmp_path):
    p = tmp_path / "a.yaml"
    p.write_text(
        "atlases:\n"
        "  - id: x\n"
        "    species: fish\n"
        "    access: open\n"
        "    source_url_verified: none\n"
        "    files:\n"
        "      - {name: a.h5ad, downloaded: true}\n"
    )
    joined = " ".join(validate_registry(p))
    for expected in ["missing", "species", "not a URL", "sha256"]:
        assert expected in joined


def test_one_to_one_orthologues_both_orders(tmp_path):
    import gzip

    from pnprank.atlases import one_to_one_orthologues

    cols = [
        "gene_stable_id", "protein_stable_id", "species", "identity", "homology_type",
        "homology_gene_stable_id", "homology_protein_stable_id", "homology_species",
        "homology_identity", "dn", "ds", "goc_score", "wga_coverage", "is_high_confidence",
        "homology_id",
    ]
    rows = [
        ["ENSMUSG1", "p", "mus_musculus", "90", "ortholog_one2one", "ENSG1", "p",
         "homo_sapiens", "90", "", "", "", "", "1", "1"],
        ["ENSG2", "p", "homo_sapiens", "90", "ortholog_one2one", "ENSMUSG2", "p",
         "mus_musculus", "90", "", "", "", "", "0", "2"],
        ["ENSG3", "p", "homo_sapiens", "90", "ortholog_one2many", "ENSMUSG3", "p",
         "mus_musculus", "90", "", "", "", "", "1", "3"],
    ]
    p = tmp_path / "h.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    assert one_to_one_orthologues(p) == {"ENSG1": "ENSMUSG1", "ENSG2": "ENSMUSG2"}
    assert one_to_one_orthologues(p, high_confidence=True) == {"ENSG1": "ENSMUSG1"}
