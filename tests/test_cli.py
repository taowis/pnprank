from pathlib import Path

import pytest

from pnprank.cli import main

TRUTH = Path(__file__).parents[1] / "data" / "truth" / "neuropathy_genes.tsv"


def test_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "pnprank" in capsys.readouterr().out


def test_validate_truth():
    assert main(["validate-truth", str(TRUTH)]) == 0
