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
