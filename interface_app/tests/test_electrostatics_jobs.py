import json
from concurrent.futures import Future
from pathlib import Path

import pytest

from interface_app import jobs
from interface_app.storage import Store


ROOT = Path(__file__).resolve().parents[2]


class FakePool:
    def __init__(self, **_kwargs):
        self.tasks = []

    def submit(self, function, *args):
        future = Future()
        self.tasks.append((function, args, future))
        return future


def saved_analysis(tmp_path):
    store = Store(tmp_path)
    analysis = store.create_analysis(source_kind="upload", source_name="1IAR.pdb", source_format="pdb",
                                     source_bytes=(ROOT / "examples/1IAR.pdb").read_bytes())
    path = store.write_json(store.analysis_dir(analysis["id"]) / "analysis.json", {"pairs": []})
    store.update_analysis(analysis["id"], status="complete", result_path=path)
    return store, analysis["id"]


def test_electrostatics_is_visible_while_queued_and_runs_after_curvature(monkeypatch, tmp_path):
    monkeypatch.setattr(jobs, "ProcessPoolExecutor", FakePool)
    store, analysis_id = saved_analysis(tmp_path)
    manager = jobs.JobManager(store)
    core = manager.submit_analysis(analysis_id)
    initial = store.list_jobs(analysis_id)
    assert [job["kind"] for job in initial] == ["core", "electrostatics"]
    assert len(manager.pool.tasks) == 1
    store.update_job(core["id"], status="complete")
    manager.pool.tasks[0][2].set_result(True)
    assert manager.pool.tasks[1][0] is jobs._run_curvature
    assert len(manager.pool.tasks) == 2
    curvature = next(job for job in store.list_jobs(analysis_id) if job["kind"] == "curvature")
    path = store.write_json(tmp_path / "curvature.json", {"status": "complete", "meshes": {}})
    store.update_job(curvature["id"], status="complete", result_path=path)
    manager.pool.tasks[1][2].set_result(None)
    assert manager.pool.tasks[2][0] is jobs._run_electrostatics
    assert len([job for job in store.list_jobs(analysis_id) if job["kind"] == "electrostatics"]) == 1


@pytest.mark.parametrize("failure", ["core", "curvature", "surface unavailable", "unreadable surface"])
def test_failed_prerequisite_does_not_leave_electrostatics_queued(monkeypatch, tmp_path, failure):
    monkeypatch.setattr(jobs, "ProcessPoolExecutor", FakePool)
    store, analysis_id = saved_analysis(tmp_path)
    manager = jobs.JobManager(store)
    core = manager.submit_analysis(analysis_id)
    store.update_job(core["id"], status="failed" if failure == "core" else "complete")
    manager.pool.tasks[0][2].set_result(True)
    if failure != "core":
        curvature = next(job for job in store.list_jobs(analysis_id) if job["kind"] == "curvature")
        path = store.write_json(tmp_path / "curvature.json", {"status": "unavailable"})
        if failure == "unreadable surface":
            Path(path).write_text("incomplete JSON")
        store.update_job(curvature["id"], status="failed" if failure == "curvature" else "complete", result_path=path)
        manager.pool.tasks[1][2].set_result(None)
    pending = next(job for job in store.list_jobs(analysis_id) if job["kind"] == "electrostatics")
    assert pending["status"] == "cancelled"
    assert all(task[0] is not jobs._run_electrostatics for task in manager.pool.tasks)


@pytest.mark.parametrize("result_status", ["complete", "unavailable"])
def test_worker_persists_optional_artifacts_without_changing_analysis(monkeypatch, tmp_path, result_status):
    store, analysis_id = saved_analysis(tmp_path)
    curvature = store.create_job(analysis_id, "curvature")
    path = store.write_json(tmp_path / "curvature.json", {"status": "complete", "meshes": {}})
    store.update_job(curvature["id"], status="complete", result_path=path)

    def calculate(*_args, progress, **_kwargs):
        progress("APBS potential solve", 0.55)
        return {"version": "0.2.0", "status": result_status, "reason": "offline" if result_status == "unavailable" else None}

    monkeypatch.setattr(jobs, "compute_electrostatics", calculate)
    first = store.create_job(analysis_id, "electrostatics")
    jobs._run_electrostatics(str(store.root), analysis_id, first["id"])
    completed = store.get_job(first["id"])
    assert completed["status"] == ("complete" if result_status == "complete" else "failed")
    assert completed["progress"] == 1.0
    assert store.get_analysis(analysis_id)["status"] == "complete"
    artifact = json.loads(Path(completed["result_path"]).read_text())
    assert artifact["status"] == result_status
    assert artifact["curvature_job_id"] == curvature["id"]
    second = store.create_job(analysis_id, "electrostatics")
    jobs._run_electrostatics(str(store.root), analysis_id, second["id"])
    assert store.get_job(second["id"])["result_path"] != completed["result_path"]


def test_worker_honors_cancellation_and_does_not_publish_a_result(monkeypatch, tmp_path):
    store, analysis_id = saved_analysis(tmp_path)
    curvature = store.create_job(analysis_id, "curvature")
    path = store.write_json(tmp_path / "curvature.json", {"status": "complete"})
    store.update_job(curvature["id"], status="complete", result_path=path)
    job = store.create_job(analysis_id, "electrostatics")

    def calculate(*_args, cancelled, **_kwargs):
        store.update_job(job["id"], status="cancelled")
        assert cancelled()
        raise jobs.ElectrostaticsCancelled()

    monkeypatch.setattr(jobs, "compute_electrostatics", calculate)
    jobs._run_electrostatics(str(store.root), analysis_id, job["id"])
    assert store.get_job(job["id"])["status"] == "cancelled"
    assert store.get_job(job["id"])["result_path"] is None
    assert store.get_analysis(analysis_id)["status"] == "complete"


def test_worker_does_not_use_an_older_surface_during_a_retry(monkeypatch, tmp_path):
    store, analysis_id = saved_analysis(tmp_path)
    old_surface = store.create_job(analysis_id, "curvature")
    path = store.write_json(tmp_path / "curvature.json", {"status": "complete"})
    store.update_job(old_surface["id"], status="complete", result_path=path)
    store.create_job(analysis_id, "curvature")
    job = store.create_job(analysis_id, "electrostatics")
    monkeypatch.setattr(jobs, "compute_electrostatics", lambda *_a, **_kw: pytest.fail("Used an older surface"))
    jobs._run_electrostatics(str(store.root), analysis_id, job["id"])
    failed = store.get_job(job["id"])
    assert failed["status"] == "failed"
    assert "not ready" in failed["error"]
    assert store.get_analysis(analysis_id)["status"] == "complete"
