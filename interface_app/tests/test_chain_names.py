import json
from pathlib import Path

import pytest
from fastapi import HTTPException

from interface_app.analysis import analyze_loaded
from interface_app.app import create_app
from interface_app.structure import (
    apply_source_chain_names, mmcif_chain_names, parse_structure, pdb_chain_names,
)


EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def endpoint(app, path, method):
    return next(route.endpoint for route in app.routes
                if route.path == path and method in getattr(route, "methods", ()))


def test_pdb_compnd_names_preserve_case_and_continuations():
    text = (
        "COMPND    MOL_ID: 1;\n"
        "COMPND   2 MOLECULE: Interleukin-4 receptor alpha\n"
        "COMPND   3 chain;\n"
        "COMPND   4 CHAIN: A, C;\n"
        "COMPND   5 MOL_ID: 2;\n"
        "COMPND   6 MOLECULE: IL-13;\n"
        "COMPND   7 CHAIN: B;\n"
    )
    assert pdb_chain_names(text) == {
        "A": "Interleukin-4 receptor alpha chain",
        "C": "Interleukin-4 receptor alpha chain",
        "B": "IL-13",
    }


def test_mmcif_entity_names_map_to_author_chains():
    names = mmcif_chain_names({
        "_entity.id": ["1", "2", "3"],
        "_entity.pdbx_description": ["Interleukin-4", "IL-4 receptor", "water"],
        "_entity_poly.entity_id": ["1", "2"],
        "_entity_poly.pdbx_strand_id": ["A", "B,C"],
        "_atom_site.auth_asym_id": ["A", "A-2"],
        "_atom_site.label_entity_id": ["1", "1"],
    })
    assert names == {"A": "Interleukin-4", "A-2": "Interleukin-4", "B": "IL-4 receptor", "C": "IL-4 receptor"}


def test_parsed_example_uses_molecule_names():
    loaded = parse_structure(EXAMPLES / "3BPN.pdb")
    names = {chain["id"]: chain["name"] for chain in loaded.metadata["chains"]}
    assert names["A"] == "INTERLEUKIN-4"
    assert names["B"] == "INTERLEUKIN-4 RECEPTOR ALPHA CHAIN"
    result = analyze_loaded(loaded)
    pair = next(item for item in result["pairs"] if (item["chain_a"], item["chain_b"]) == ("A", "B"))
    assert pair["name_a"] == "INTERLEUKIN-4"
    assert pair["contact_map"]["col_name"] == "INTERLEUKIN-4 RECEPTOR ALPHA CHAIN"


def test_missing_names_fall_back_to_chain_ids():
    assert pdb_chain_names("ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n") == {}


def test_saved_results_are_backfilled_once():
    source = (EXAMPLES / "1IAR.pdb").read_text()
    result = {
        "metadata": {"chains": [{"id": "A", "name": "Chain A"}, {"id": "B", "name": "Chain B"}]},
        "chains": [{"id": "A", "name": "Chain A"}],
        "pairs": [{"chain_a": "A", "chain_b": "B", "name_a": "Chain A", "name_b": "Chain B",
                   "contact_map": {"row_chain": "A", "col_chain": "B", "row_name": "Chain A", "col_name": "Chain B"}}],
    }
    assert apply_source_chain_names(result, source, "pdb")
    assert result["metadata"]["chains"][0]["name"] == "PROTEIN (INTERLEUKIN-4)"
    assert result["pairs"][0]["name_b"] == "PROTEIN (INTERLEUKIN-4 RECEPTOR ALPHA CHAIN)"
    assert result["pairs"][0]["contact_map"]["row_name"] == "PROTEIN (INTERLEUKIN-4)"
    assert not apply_source_chain_names(result, source, "pdb")


def test_delete_and_clear_history_keep_active_analyses(tmp_path):
    app = create_app(tmp_path)
    store = app.state.store
    try:
        made = [store.create_analysis(source_kind="upload", source_name=f"{i}.pdb", source_format="pdb",
                                      source_bytes=b"END\n") for i in range(3)]
        done, running, other = (item["id"] for item in made)
        store.update_job(store.create_job(done, "core")["id"], status="complete")
        store.update_job(store.create_job(running, "core")["id"], status="running")

        delete_one = endpoint(app, "/api/analyses/{analysis_id}", "DELETE")
        with pytest.raises(HTTPException) as error:
            delete_one(running)
        assert error.value.status_code == 409
        delete_one(done)
        assert store.get_analysis(done) is None
        assert not (tmp_path / "analyses" / done).exists()
        with pytest.raises(HTTPException) as error:
            delete_one(done)
        assert error.value.status_code == 404

        cleared = endpoint(app, "/api/analyses", "DELETE")()
        assert cleared == {"deleted": [other], "skipped": [running]}
        assert [item["id"] for item in store.list_analyses()] == [running]
    finally:
        app.state.jobs.shutdown()


def test_result_endpoint_backfills_saved_names(tmp_path):
    app = create_app(tmp_path)
    store = app.state.store
    try:
        analysis = store.create_analysis(source_kind="upload", source_name="1IAR.pdb", source_format="pdb",
                                         source_bytes=(EXAMPLES / "1IAR.pdb").read_bytes())
        path = store.write_json(tmp_path / "old.json", {
            "metadata": {"chains": [{"id": "A", "name": "Chain A"}]}, "chains": [], "pairs": [],
        })
        store.update_analysis(analysis["id"], result_path=path)
        response = endpoint(app, "/api/analyses/{analysis_id}/result", "GET")(analysis["id"])
        assert json.loads(response.body)["metadata"]["chains"][0]["name"] == "PROTEIN (INTERLEUKIN-4)"
        assert json.loads(Path(path).read_text())["metadata"]["chain_names_from_file"] is True
    finally:
        app.state.jobs.shutdown()


def test_history_list_returns_slim_chain_summaries(tmp_path):
    app = create_app(tmp_path)
    store = app.state.store
    try:
        analysis = store.create_analysis(source_kind="upload", source_name="1IAR.pdb", source_format="pdb",
                                         source_bytes=b"END\n")
        metadata = {"chains": [{"id": "A", "name": "IL-4", "residue_count": 2,
                                "residues": [{"number": 1}, {"number": 2}]}]}
        store.update_analysis(analysis["id"], metadata_json=json.dumps(metadata))
        item = endpoint(app, "/api/analyses", "GET")()["analyses"][0]
        assert item["chains"] == [{"id": "A", "name": "IL-4", "residue_count": 2}]
        assert "metadata" not in item
    finally:
        app.state.jobs.shutdown()
