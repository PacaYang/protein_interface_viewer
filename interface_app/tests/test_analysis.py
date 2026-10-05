import json
from pathlib import Path

from interface_app.analysis import analyze_loaded
from interface_app.structure import parse_structure


EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def test_existing_3bpn_pair_regression():
    loaded = parse_structure(EXAMPLES / "3BPN.pdb")
    result = analyze_loaded(loaded)
    counts = {tuple((pair["chain_a"], pair["chain_b"])): pair["n_contact_residue_pairs"] for pair in result["pairs"]}
    assert counts[("A", "B")] == 41
    assert counts[("A", "C")] == 55
    assert counts[("B", "C")] == 36
    ac = next(pair for pair in result["pairs"] if pair["chain_a"] == "A" and pair["chain_b"] == "C")
    assert 800 < ac["buried_surface_area_A2"] < 870
    assert ac["n_salt_bridges"] == 7
    assert {
        (pair["chain_a"], pair["chain_b"]): (
            len(pair["contact_map"]["rows"]),
            len(pair["contact_map"]["cols"]),
            len(pair["contact_map"]["cells"]),
        ) for pair in result["pairs"]
    } == {
        ("A", "B"): (25, 23, 67),
        ("A", "C"): (33, 27, 88),
        ("B", "C"): (23, 24, 64),
    }
    assert all("like_charges" in cell for cell in ac["contact_map"]["cells"])
    assert ac["n_like_charges"] == 4


def test_chain_and_residue_metadata_preserve_identity():
    loaded = parse_structure(EXAMPLES / "3BPO.pdb")
    result = analyze_loaded(loaded)
    assert {item["id"] for item in result["chains"]} >= {"A", "B", "C"}
    assert all("insertion_code" in row for pair in result["pairs"] for row in pair["residues"])
