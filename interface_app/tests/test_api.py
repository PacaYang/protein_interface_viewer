import json
from pathlib import Path

import pytest

from interface_app.app import create_app


ROOT = Path(__file__).resolve().parents[2]


def test_app_exposes_analysis_routes(tmp_path):
    app = create_app(tmp_path)
    paths = {route.path for route in app.routes}
    assert "/api/analyses/upload" in paths
    assert "/api/analyses/{analysis_id}/pockets" in paths
    assert "/api/jobs/{job_id}" in paths


def test_health(tmp_path):
    app = create_app(tmp_path)
    endpoint = next(route.endpoint for route in app.routes if route.path == "/api/health")
    assert endpoint()["status"] == "ok"


def test_curvature_endpoint_returns_latest_completed_result(tmp_path):
    app = create_app(tmp_path)
    store = app.state.store
    try:
        analysis = store.create_analysis(
            source_kind="upload", source_name="test.pdb", source_format="pdb", source_bytes=b"END\n",
        )
        for scale in ("old", "latest"):
            job = store.create_job(analysis["id"], "curvature")
            path = store.write_json(tmp_path / f"{scale}.json", {"result": scale})
            store.update_job(job["id"], status="complete", result_path=path)
        endpoint = next(route.endpoint for route in app.routes if route.path == "/api/analyses/{analysis_id}/curvature")
        response = endpoint(analysis["id"])
        assert json.loads(response.body) == {"result": "latest"}
    finally:
        app.state.jobs.shutdown()


def test_pocket_request_requires_exactly_one_definition():
    import pytest
    from pydantic import ValidationError

    from interface_app.app import PocketRequest

    base = {"target_chain": "A", "partner_chain": "B"}
    anchor = PocketRequest(**base, anchor_residue={"number": 67})
    assert anchor.model_dump(exclude_none=True)["anchor_residue"] == {"number": 67, "insertion_code": ""}
    residues = PocketRequest(**base, pocket_residues=[{"number": 67}, {"chain_id": "B", "number": 12, "insertion_code": "A"}])
    payload = residues.model_dump(exclude_none=True)
    assert "anchor_residue" not in payload
    assert payload["pocket_residues"] == [
        {"number": 67, "insertion_code": ""},
        {"chain_id": "B", "number": 12, "insertion_code": "A"},
    ]
    for invalid in ({}, {"anchor_residue": {"number": 1}, "pocket_residues": [{"number": 1}]},
                    {"pocket_residues": []}, {"pocket_residues": [{"number": 1}] * 61}):
        with pytest.raises(ValidationError):
            PocketRequest(**base, **invalid)


@pytest.mark.parametrize("definition", [
    {"anchor_residue": {"number": 67}},
    {"pocket_residues": [{"number": 67}, {"chain_id": "B", "number": 12}]},
])
def test_pocket_endpoint_serializes_only_the_selected_definition(tmp_path, monkeypatch, definition):
    from interface_app.app import PocketRequest

    app = create_app(tmp_path)
    store = app.state.store
    try:
        analysis = store.create_analysis(
            source_kind="upload", source_name="test.pdb", source_format="pdb", source_bytes=b"END\n",
        )
        store.update_analysis(analysis["id"], result_path="completed.json")
        submitted = []
        monkeypatch.setattr(app.state.jobs, "submit_pocket",
                            lambda analysis_id, request: submitted.append(request) or {"id": "pocket-job"})
        request = PocketRequest(target_chain="A", partner_chain="B", **definition)
        endpoint = next(route.endpoint for route in app.routes
                        if route.path == "/api/analyses/{analysis_id}/pockets" and "POST" in route.methods)
        assert endpoint(analysis["id"], request)["job_id"] == "pocket-job"
        assert submitted == [request.model_dump(exclude_none=True)]
        assert all(value is not None for value in submitted[0].values())
    finally:
        app.state.jobs.shutdown()


@pytest.mark.parametrize("definition", [
    {"anchor_residue": {"number": 67, "insertion_code": ""}},
    {"pocket_residues": [{"number": 67}, {"chain_id": "B", "number": 12}]},
])
def test_pocket_retry_preserves_anchor_and_residue_definitions(tmp_path, monkeypatch, definition):
    app = create_app(tmp_path)
    store = app.state.store
    try:
        analysis = store.create_analysis(
            source_kind="upload", source_name="test.pdb", source_format="pdb", source_bytes=b"END\n",
        )
        request = {"target_chain": "A", "partner_chain": "B", "radius_A": 8, "grid_A": 0.6, **definition}
        job = store.create_job(analysis["id"], "pocket", request=request)
        store.update_job(job["id"], status="interrupted")
        submitted = []
        monkeypatch.setattr(app.state.jobs, "submit_pocket",
                            lambda analysis_id, payload: submitted.append(payload) or {"id": "replacement"})
        endpoint = next(route.endpoint for route in app.routes if route.path == "/api/jobs/{job_id}/retry")
        assert endpoint(job["id"])["job_id"] == "replacement"
        assert submitted == [request]
    finally:
        app.state.jobs.shutdown()


@pytest.mark.parametrize("source_format", ["pdb", "mmcif"])
def test_face_view_reads_saved_results_without_starting_a_job(tmp_path, source_format):
    from Bio.PDB import MMCIFIO

    from interface_app.analysis import encode_pair
    from interface_app.structure import parse_structure

    loaded = parse_structure(ROOT / "examples" / "1IAR.pdb")
    source = ROOT / "examples" / "1IAR.pdb"
    if source_format == "mmcif":
        source = tmp_path / "1IAR.cif"
        writer = MMCIFIO()
        writer.set_structure(loaded.structure)
        writer.save(str(source))
    app = create_app(tmp_path / "data")
    try:
        store = app.state.store
        item = store.create_analysis(source_kind="upload", source_name=source.name,
                                     source_format=source_format, source_bytes=source.read_bytes())
        pair = {"id": encode_pair("A", "B"), "chain_a": "A", "chain_b": "B",
                "parameters": {"contact_cutoff_A": 5.0},
                "residues": [residue for chain in loaded.metadata["chains"]
                             for residue in chain["residues"][:4]]}
        path = store.write_json(store.analysis_dir(item["id"]) / "saved.json", {"pairs": [pair]})
        store.update_analysis(item["id"], status="complete", result_path=path)
        endpoint = next(route.endpoint for route in app.routes
                        if route.path == "/api/analyses/{analysis_id}/pairs/{pair_id}/face-view")
        response = endpoint(item["id"], pair["id"])
        assert response["available"]
        assert response["binding_cutoff_A"] == 5
        assert set(response["chains"]) == {"A", "B"}
        assert response["nearest"]["A"]["A:81"]["partner_key"] == "B:125"
        assert response["nearest"]["A"]["A:81"]["distance_A"] == pytest.approx(2.53, abs=0.01)
        assert store.list_jobs(item["id"]) == []
        assert json.loads(Path(path).read_text()) == {"pairs": [pair]}
        pair["residues"] = []
        store.write_json(path, {"pairs": [pair]})
        unavailable = endpoint(item["id"], pair["id"])
        assert not unavailable["available"]
        assert "No binding residues" in unavailable["reason"]
    finally:
        app.state.jobs.shutdown()


def test_face_view_missing_analysis_pair_and_pending_results(tmp_path):
    from fastapi import HTTPException

    app = create_app(tmp_path)
    try:
        store = app.state.store
        endpoint = next(route.endpoint for route in app.routes
                        if route.path == "/api/analyses/{analysis_id}/pairs/{pair_id}/face-view")
        with pytest.raises(HTTPException) as exc:
            endpoint("missing", "missing")
        assert exc.value.status_code == 404
        item = store.create_analysis(source_kind="upload", source_name="test.pdb",
                                     source_format="pdb", source_bytes=b"END\n")
        with pytest.raises(HTTPException) as exc:
            endpoint(item["id"], "missing")
        assert exc.value.status_code == 409
        path = store.write_json(tmp_path / "completed.json", {"pairs": []})
        store.update_analysis(item["id"], result_path=path)
        with pytest.raises(HTTPException) as exc:
            endpoint(item["id"], "missing")
        assert exc.value.status_code == 404
    finally:
        app.state.jobs.shutdown()
