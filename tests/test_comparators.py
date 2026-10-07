import gzip
import zipfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from pnprank.comparators import (
    adjacency,
    ancestors_closure,
    go_gene_matrix,
    go_similarity_scores,
    parse_obo,
    read_biogrid,
    read_gaf,
    rwr_scores,
)

OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root

[Term]
id: GO:0000002
name: child
is_a: GO:0000001 ! root
alt_id: GO:0000099

[Term]
id: GO:0000003
name: part
relationship: part_of GO:0000002 ! child

[Term]
id: GO:0000004
name: old
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def test_parse_obo_and_closure(tmp_path):
    p = tmp_path / "go.obo"
    p.write_text(OBO)
    parents, alt = parse_obo(p)
    assert set(parents) == {"GO:0000001", "GO:0000002", "GO:0000003"}  # obsolete dropped
    assert alt == {"GO:0000099": "GO:0000002"}
    c = ancestors_closure(parents)
    assert c["GO:0000003"] == {"GO:0000001", "GO:0000002", "GO:0000003"}


def test_read_gaf_drops_not_nd_and_late(tmp_path):
    rows = [
        ["UniProt", "P1", "A", "", "GO:1", "ref", "IDA", "", "P", "", "", "protein", "taxon:9606", "20040101", "x"],
        ["UniProt", "P2", "B", "NOT", "GO:1", "ref", "IDA", "", "P", "", "", "protein", "taxon:9606", "20040101", "x"],
        ["UniProt", "P3", "C", "", "GO:1", "ref", "ND", "", "P", "", "", "protein", "taxon:9606", "20040101", "x"],
        ["UniProt", "P4", "D", "", "GO:1", "ref", "IEA", "", "P", "", "", "protein", "taxon:9606", "20050102", "x"],
    ]
    p = tmp_path / "a.gaf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("!gaf-version: 1.0\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    assert read_gaf(p).db_id.tolist() == ["P1", "P4"]
    assert read_gaf(p, before="20050101").db_id.tolist() == ["P1"]


def _go_fixture():
    parents = {"T0": set(), "T1": {"T0"}, "T2": {"T0"}, "T3": {"T1"}}
    ann = pd.DataFrame({"ensembl_gene_id": ["G0", "G1", "G2", "G3", "G4", "G5"],
                        "go_id": ["T3", "T3", "T1", "T2", "T2", "T2"]})
    universe = [f"G{i}" for i in range(7)]
    m, terms = go_gene_matrix(ann, parents, {}, universe, min_genes=1, max_genes=5)
    return m, terms, universe


def test_go_matrix_propagates_and_filters():
    m, terms, universe = _go_fixture()
    assert "T0" not in terms  # annotated to 6 genes > max_genes
    assert set(terms) == {"T1", "T2", "T3"}
    assert m[universe.index("G0"), terms.index("T1")] == 1  # propagated from T3


def test_go_similarity_ranks_functional_neighbours():
    m, _, universe = _go_fixture()
    s = go_similarity_scores(m, universe, {"G0"})
    assert s["G1"] > s["G2"] > s["G3"]   # same term > parent term > unrelated
    assert s["G6"] == 0                  # unannotated
    assert (go_similarity_scores(m, universe, set()) == 0).all()


def test_rwr_prefers_network_neighbours():
    universe = [f"G{i}" for i in range(6)]
    edges = pd.DataFrame({"a": ["G0", "G1", "G2"], "b": ["G1", "G2", "G3"]})
    adj = adjacency(edges, universe)
    assert (adj != adj.T).nnz == 0
    s = rwr_scores(adj, universe, {"G0"}, restart=0.5)
    assert s["G1"] > s["G2"] > s["G3"] > 0
    assert s["G4"] == 0 and s["G5"] == 0
    assert np.isclose(s.sum(), 1.0)


def _zip(tmp_path, text, name="BIOGRID-ORGANISM-Homo_sapiens-x.tab2.txt"):
    p = tmp_path / "b.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr(name, text)
    return p


def test_read_biogrid_tab2(tmp_path):
    cols = ["#BioGRID Interaction ID", "Entrez Gene Interactor A", "Entrez Gene Interactor B",
            "Official Symbol Interactor A", "Official Symbol Interactor B",
            "Experimental System Type", "Organism Interactor A", "Organism Interactor B"]
    rows = [["1", "10", "20", "A", "B", "physical", "9606", "9606"],
            ["2", "20", "10", "B", "A", "physical", "9606", "9606"],     # duplicate
            ["3", "10", "30", "A", "C", "genetic", "9606", "9606"],      # genetic
            ["4", "10", "40", "A", "D", "physical", "9606", "10090"],    # mouse
            ["5", "10", "10", "A", "A", "physical", "9606", "9606"],     # self
            ["6", "99", "20", "X", "B", "physical", "9606", "9606"]]     # symbol fallback
    text = "\n".join("\t".join(r) for r in [cols] + rows)
    e = read_biogrid(_zip(tmp_path, text), {"10": "E1", "20": "E2", "30": "E3", "40": "E4"},
                     {"X": "E9"})
    assert sorted(map(tuple, e.to_numpy())) == [("E1", "E2"), ("E2", "E9")]


def test_read_biogrid_tab_v2(tmp_path):
    cols = ["INTERACTOR_A", "INTERACTOR_B", "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B",
            "ALIASES_FOR_A", "ALIASES_FOR_B", "EXPERIMENTAL_SYSTEM",
            "SOURCE", "PUBMED_ID", "ORGANISM_A_ID", "ORGANISM_B_ID"]
    rows = [["HGNC:1", "HGNC:2", "A", "B", "-", "-", "Invivo", "x", "1", "9606", "9606"],
            ["HGNC:1", "HGNC:3", "A", "C", "-", "-", "Dosage Rescue", "x", "1", "9606", "9606"],
            ["HGNC:9", "HGNC:3", "Z", "C", "-", "-", "Two-hybrid", "x", "1", "9606", "9606"]]
    text = "header line\n" + "\n".join("\t".join(r) for r in [cols] + rows)
    e = read_biogrid(_zip(tmp_path, text, "BIOGRID-ORGANISM-Homo_sapiens-2.0.36.tab.txt"), {},
                     {"Z": "E9"}, hgnc_ids={"HGNC:1": "E1", "HGNC:2": "E2", "HGNC:3": "E3"})
    assert sorted(map(tuple, e.to_numpy())) == [("E1", "E2"), ("E3", "E9")]


def test_uniprot_to_ensembl_resolves_trembl_and_secondary(tmp_path):
    from pnprank.comparators import uniprot_to_ensembl

    hgnc = tmp_path / "hgnc.txt"
    hgnc.write_text("hgnc_id\tstatus\tensembl_gene_id\tuniprot_ids\tentrez_id\n"
                    "HGNC:1\tApproved\tE1\tP11111\t1\n"
                    "HGNC:2\tApproved\tE2\t\t2\n"
                    "HGNC:3\tApproved\tE3\tQ33333\t3\n")
    idm = tmp_path / "idm.dat.gz"
    with gzip.open(idm, "wt") as fh:
        fh.write("A0TREM\tHGNC\tHGNC:2\nA0TREM\tGeneID\t2\n"
                 "B0AMBI\tHGNC\tHGNC:1\nB0AMBI\tHGNC\tHGNC:3\n")
    sec = tmp_path / "sec_ac.txt"
    sec.write_text("header\nSecondary AC  Primary AC\n____________  __________\n"
                   "O00OLD        P11111\nO00TRM        A0TREM\nO00BAD        P11111\n"
                   "O00BAD        Q33333\n")
    m = uniprot_to_ensembl(hgnc, idm, sec)
    assert m["P11111"] == "E1" and m["A0TREM"] == "E2"
    assert m["O00OLD"] == "E1" and m["O00TRM"] == "E2"
    assert "B0AMBI" not in m and "O00BAD" not in m


def test_comparator_releases_predate_cutoffs():
    """Every dated release used at cutoff c is dated before 1 January of c."""
    cfg = yaml.safe_load(Path("config/comparators.yaml").read_text())
    n = 0
    for kind in ("go_annotation", "interactions"):
        for c, r in cfg[kind]["releases"].items():
            assert date.fromisoformat(r["date"]) < date(int(c), 1, 1), (kind, c)
            n += 1
    for c, r in cfg["go_ontology"]["releases"].items():
        assert date.fromisoformat(r["release"]) < date(int(c), 1, 1)
    assert date.fromisoformat(cfg["interactions"]["earliest_release"]["date"]) >= date(2005, 1, 1)
    assert 2005 not in cfg["interactions"]["releases"]
    assert n == 11
