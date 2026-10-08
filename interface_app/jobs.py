"""Process-isolated analysis jobs and staged result persistence."""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .analysis import analyze_loaded
from .curvature import compute_curvature
from .electrostatics import compute_electrostatics, ElectrostaticsCancelled
from .pocket import analyze_pocket
from .storage import Store
from .structure import parse_structure


def _run_core(data_root: str, analysis_id: str, job_id: str):
    store = Store(data_root)
    analysis = store.get_analysis(analysis_id)
    if not analysis:
        return
    job = store.get_job(job_id)
    if not job:
        return
    if job["status"] == "cancelled":
        store.update_analysis(analysis_id, status="cancelled", error="Interface analysis was cancelled")
        return
    store.update_analysis(analysis_id, status="running", error=None)
    store.update_job(job["id"], status="running", stage="parsing", progress=0.05)
    try:
        loaded = parse_structure(analysis["source_path"], source_name=analysis["source_name"])
        store.update_job(job["id"], stage="interfaces", progress=0.25)
        result = analyze_loaded(loaded)
        result["analysis_id"] = analysis_id
        result["version"] = "0.2.0"
        result_path = store.write_json(store.analysis_dir(analysis_id) / "analysis.json", result)
        store.update_analysis(analysis_id, status="running", metadata_json=json.dumps(result["metadata"]), result_path=result_path)
        current = store.get_job(job["id"])
        if current and current["status"] == "cancelled":
            store.update_analysis(analysis_id, status="cancelled", error="Interface analysis was cancelled")
            return
        store.update_job(job["id"], status="complete", stage="complete", progress=1.0, result_path=result_path)
        return True
    except Exception as exc:
        current = store.get_job(job["id"])
        if current and current["status"] == "cancelled":
            store.update_analysis(analysis_id, status="cancelled", error="Interface analysis was cancelled")
        else:
            store.update_job(job["id"], status="failed", stage="failed", error=str(exc))
            store.update_analysis(analysis_id, status="failed", error=str(exc))


def _run_curvature(data_root: str, analysis_id: str, job_id: str):
    store = Store(data_root)
    analysis = store.get_analysis(analysis_id)
    if not analysis:
        return
    try:
        if store.get_job(job_id)["status"] == "cancelled":
            store.update_analysis(analysis_id, status="complete", error=None)
            return
        store.update_job(job_id, status="running", stage="parsing", progress=0.02)
        loaded = parse_structure(analysis["source_path"], source_name=analysis["source_name"])
        result = json.loads(Path(analysis["result_path"]).read_text())
        store.update_job(job_id, stage="surface curvature", progress=0.15)
        curvature = compute_curvature(loaded.model, result)
        result_path = store.write_json(store.analysis_dir(analysis_id) / "curvature.json", curvature)
        current = store.get_job(job_id)
        if current and current["status"] == "cancelled":
            store.update_analysis(analysis_id, status="complete", error=None)
            return
        store.update_job(job_id, status="complete" if curvature.get("status") != "failed" else "failed",
                         stage="complete", progress=1.0, result_path=result_path)
        # Core results remain usable even if curvature is unavailable.
        store.update_analysis(analysis_id, status="complete", error=None)
    except Exception as exc:
        current = store.get_job(job_id)
        if current and current["status"] == "cancelled":
            store.update_analysis(analysis_id, status="complete", error=None)
        else:
            store.update_job(job_id, status="failed", stage="failed", error=str(exc))
            store.update_analysis(analysis_id, status="complete", error=f"Curvature unavailable: {exc}")


def _run_electrostatics(data_root: str, analysis_id: str, job_id: str):
    """Run the optional APBS/PDB2PQR stage without affecting core results."""

    store = Store(data_root)
    analysis = store.get_analysis(analysis_id)
    if not analysis:
        return

    def cancelled():
        job = store.get_job(job_id)
        return not job or job["status"] not in {"queued", "running"}

    def progress(stage, amount):
        if not cancelled():
            store.update_job(job_id, stage=stage, progress=amount)

    try:
        current = store.get_job(job_id)
        if not current or current["status"] not in {"queued", "running"}:
            return
        store.update_job(job_id, status="running", stage="preparing electrostatics", progress=0.02)
        loaded = parse_structure(analysis["source_path"], source_name=analysis["source_name"])
        core = json.loads(Path(analysis["result_path"]).read_text()) if analysis.get("result_path") else {}
        curvature_job = next(
            (item for item in reversed(store.list_jobs(analysis_id))
             if item["kind"] == "curvature"),
            None,
        )
        if not curvature_job or curvature_job["status"] != "complete" or not curvature_job.get("result_path"):
            raise RuntimeError("Curvature results are not ready")
        curvature = json.loads(Path(curvature_job["result_path"]).read_text())
        if curvature.get("status") != "complete":
            raise RuntimeError("A complete surface mesh is required for electrostatics")
        store.update_job(job_id, stage="APBS/PDB2PQR", progress=0.12)
        electrostatics = compute_electrostatics(
            loaded.model,
            core,
            curvature,
            data_root=data_root,
            progress=progress,
            cancelled=cancelled,
        )
        if cancelled():
            return
        electrostatics["curvature_job_id"] = curvature_job["id"]
        result_path = store.write_json(
            store.analysis_dir(analysis_id) / f"electrostatics_{job_id}.json", electrostatics
        )
        current = store.get_job(job_id)
        if current and current["status"] == "cancelled":
            return
        if electrostatics.get("status") == "complete":
            store.update_job(job_id, status="complete", stage="complete", progress=1.0,
                             result_path=result_path, error=None)
        else:
            # Keep the diagnostic artifact available to the browser while
            # allowing the existing retry route to run after tools are fixed.
            store.update_job(job_id, status="failed", stage="unavailable", progress=1.0,
                             result_path=result_path, error=electrostatics.get("reason"))
    except ElectrostaticsCancelled:
        return
    except Exception as exc:
        current = store.get_job(job_id)
        if current and current["status"] == "cancelled":
            return
        artifact = {"version": "0.2.0", "status": "unavailable", "reason": str(exc),
                    "parameters": {}, "report": None, "meshes": {}}
        path = store.write_json(store.analysis_dir(analysis_id) / f"electrostatics_{job_id}.json", artifact)
        store.update_job(job_id, status="failed", stage="failed", progress=1.0, error=str(exc), result_path=path)


def _run_pocket(data_root: str, analysis_id: str, job_id: str, request: dict):
    store = Store(data_root)
    analysis = store.get_analysis(analysis_id)
    if not analysis:
        return
    try:
        if store.get_job(job_id)["status"] == "cancelled":
            return
        store.update_job(job_id, status="running", stage="parsing", progress=0.05)
        loaded = parse_structure(analysis["source_path"], source_name=analysis["source_name"])
        store.update_job(job_id, stage="grid flood-fill", progress=0.2)
        result = analyze_pocket(loaded.model, **request)
        result["analysis_id"] = analysis_id
        result["job_id"] = job_id
        result_path = store.write_json(store.analysis_dir(analysis_id) / f"pocket_{job_id}.json", result)
        if store.get_job(job_id)["status"] != "cancelled":
            store.update_job(job_id, status="complete", stage="complete", progress=1.0, result_path=result_path)
    except Exception as exc:
        store.update_job(job_id, status="failed", stage="failed", error=str(exc))


class JobManager:
    def __init__(self, store: Store, max_workers: int = 2):
        self.store = store
        self.pool = ProcessPoolExecutor(max_workers=max_workers)

    def submit_analysis(self, analysis_id: str):
        job = self.store.create_job(analysis_id, "core")
        electrostatics = self.store.create_job(analysis_id, "electrostatics")
        future = self.pool.submit(_run_core, str(self.store.root), analysis_id, job["id"])
        future.add_done_callback(lambda completed: self._queue_curvature(
            analysis_id, job["id"], electrostatics["id"], completed,
        ))
        return job

    def submit_pocket(self, analysis_id: str, request: dict):
        job = self.store.create_job(analysis_id, "pocket", request=request)
        self.pool.submit(_run_pocket, str(self.store.root), analysis_id, job["id"], request)
        return job

    def submit_curvature(self, analysis_id: str, electrostatics_job_id: str | None = None):
        job = self.store.create_job(analysis_id, "curvature")
        if electrostatics_job_id is None:
            electrostatics = self.store.create_job(analysis_id, "electrostatics")
            electrostatics_job_id = electrostatics["id"]
        future = self.pool.submit(_run_curvature, str(self.store.root), analysis_id, job["id"])
        future.add_done_callback(lambda completed: self._queue_electrostatics(
            analysis_id, job["id"], electrostatics_job_id, completed,
        ))
        return job

    def submit_electrostatics(self, analysis_id: str, job_id: str | None = None):
        job = self.store.get_job(job_id) if job_id else None
        if not job:
            job = self.store.create_job(analysis_id, "electrostatics")
        if job["status"] != "queued":
            return job
        self.pool.submit(_run_electrostatics, str(self.store.root), analysis_id, job["id"])
        return job

    def _queue_curvature(self, analysis_id: str, core_job_id: str, electrostatics_job_id: str, completed):
        try:
            completed.result()
        except Exception:
            self.store.update_job(electrostatics_job_id, status="cancelled", stage="cancelled",
                                  progress=1.0, error="Interface worker did not complete")
            return
        analysis = self.store.get_analysis(analysis_id)
        core = self.store.get_job(core_job_id)
        if not analysis or analysis["status"] == "failed" or not core or core["status"] != "complete":
            pending = self.store.get_job(electrostatics_job_id)
            if pending and pending["status"] in {"queued", "running"}:
                self.store.update_job(electrostatics_job_id, status="cancelled", stage="cancelled",
                                       progress=1.0, error="Interface analysis did not complete")
            return
        self.submit_curvature(analysis_id, electrostatics_job_id)

    def _queue_electrostatics(self, analysis_id: str, curvature_job_id: str,
                              electrostatics_job_id: str, completed):
        try:
            completed.result()
        except Exception:
            self.store.update_job(electrostatics_job_id, status="cancelled", stage="cancelled",
                                  progress=1.0, error="Curvature worker did not complete")
            return
        analysis = self.store.get_analysis(analysis_id)
        curvature = self.store.get_job(curvature_job_id)
        if not analysis or not curvature or curvature["status"] != "complete" or not curvature.get("result_path"):
            pending = self.store.get_job(electrostatics_job_id)
            if pending and pending["status"] in {"queued", "running"}:
                self.store.update_job(electrostatics_job_id, status="cancelled", stage="cancelled",
                                       progress=1.0, error="Curvature results were unavailable")
            return
        try:
            surface = json.loads(Path(curvature["result_path"]).read_text())
        except (OSError, ValueError) as exc:
            self.store.update_job(electrostatics_job_id, status="cancelled", stage="unavailable",
                                  progress=1.0, error=f"Could not read the surface mesh: {exc}")
            return
        if surface.get("status") != "complete":
            self.store.update_job(electrostatics_job_id, status="cancelled", stage="unavailable",
                                  progress=1.0, error="A complete surface mesh is required for electrostatics")
            return
        self.submit_electrostatics(analysis_id, electrostatics_job_id)

    def shutdown(self):
        self.pool.shutdown(wait=False, cancel_futures=True)
