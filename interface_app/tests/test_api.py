import json
from pathlib import Path

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
