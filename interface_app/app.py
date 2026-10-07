"""FastAPI entry point for the local interface-analysis browser app."""

from __future__ import annotations

import json
import os
import io
import zipfile
import gzip
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, model_validator

from .face_view import binding_face_view
from .jobs import JobManager
from .storage import Store
from .structure import apply_source_chain_names, fetch_rcsb, parse_text, rcsb_entry_metadata

MAX_UPLOAD_BYTES = 50 * 1024 * 1024
ACTIVE_JOB_STATES = {"queued", "running"}


class PDBRequest(BaseModel):
    pdb_id: str = Field(min_length=4, max_length=4)
    assembly_id: str | None = None


class AnchorRequest(BaseModel):
    number: int
    insertion_code: str = ""


class ResidueRef(AnchorRequest):
    chain_id: str | None = None


class PocketRequest(BaseModel):
    target_chain: str
    partner_chain: str
    anchor_residue: AnchorRequest | None = None
    pocket_residues: list[ResidueRef] | None = Field(default=None, min_length=1, max_length=60)
    radius_A: float = Field(default=8.0, ge=4.0, le=12.0)
    grid_A: float = Field(default=0.6, ge=0.4, le=1.2)

    @model_validator(mode="after")
    def require_one_definition(self):
        if (self.anchor_residue is None) == (self.pocket_residues is None):
            raise ValueError("Provide either an anchor residue or pocket residues, but not both")
        if self.target_chain == self.partner_chain:
            raise ValueError("Target and partner chains must be different")
        if self.pocket_residues is not None:
            for residue in self.pocket_residues:
                if residue.chain_id is not None and residue.chain_id not in {self.target_chain, self.partner_chain}:
                    raise ValueError("Pocket residues must belong to the target or partner chain")
        return self


def _format_from_name(name: str) -> str:
    lowered = name.lower()
    return "mmcif" if lowered.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")) else "pdb"


def create_app(data_dir: str | Path | None = None) -> FastAPI:
    store = Store(data_dir)
    manager = JobManager(store)
    static_dir = Path(__file__).parent / "static"
    app = FastAPI(title="Protein Interface Analyzer", version="0.2.0")
    app.state.store = store
    app.state.jobs = manager

    # A local server restart cannot resume a numerical process safely. Mark those
    # jobs explicitly so history never presents stale work as active.
    for item in store.list_analyses():
        for job in store.list_jobs(item["id"]):
            if job["status"] in {"queued", "running"}:
                store.update_job(job["id"], status="interrupted", stage="server restart",
                                 error="The server restarted before this job finished")
        if item["status"] in {"queued", "running"}:
            store.update_analysis(item["id"], status="interrupted",
                                  error="The server restarted before this analysis finished")

    @app.on_event("shutdown")
    def shutdown():
        manager.shutdown()

    @app.get("/", response_class=FileResponse)
    def index():
        return FileResponse(static_dir / "index.html")

    @app.get("/api/health")
    def health():
        return {"status": "ok", "version": "0.2.0"}

    def history_item(item: dict) -> dict:
        # The history list only needs chain composition; full metadata repeats
        # every residue of every chain and is served by the result endpoint.
        metadata = item.pop("metadata", {}) or {}
        item["chains"] = [
            {"id": chain.get("id"), "name": chain.get("name"), "residue_count": chain.get("residue_count", 0)}
            for chain in metadata.get("chains", [])
        ]
        if "rcsb" in metadata:
            item["rcsb"] = metadata["rcsb"]
        return {**item, "jobs": store.list_jobs(item["id"])}

    @app.get("/api/analyses")
    def analyses():
        return {"analyses": [history_item(item) for item in store.list_analyses()]}

    @app.post("/api/analyses/upload", status_code=202)
    async def upload_analysis(file: UploadFile = File(...)):
        source_name = file.filename or "uploaded.pdb"
        data = await file.read()
        if not data:
            raise HTTPException(400, "The uploaded file is empty")
        if len(data) > MAX_UPLOAD_BYTES:
            raise HTTPException(413, "Uploaded structures are limited to 50 MB")
        source_format = _format_from_name(source_name)
        if source_name.lower().endswith(".gz"):
            try:
                data = gzip.decompress(data)
                source_name = source_name[:-3]
            except OSError as exc:
                raise HTTPException(400, f"Could not decompress the uploaded file: {exc}") from exc
            if len(data) > MAX_UPLOAD_BYTES:
                raise HTTPException(413, "Uploaded structures are limited to 50 MB after decompression")
        try:
            parse_text(data.decode("utf-8", errors="replace"), source_name=source_name, source_format=source_format)
        except Exception as exc:
            raise HTTPException(400, f"Could not parse the coordinate file: {exc}") from exc
        analysis = store.create_analysis(source_kind="upload", source_name=source_name,
                                         source_format=source_format, source_bytes=data)
        job = manager.submit_analysis(analysis["id"])
        return {"analysis_id": analysis["id"], "job_id": job["id"], "status_url": f"/api/jobs/{job['id']}"}

    @app.post("/api/analyses/pdb", status_code=202)
    def pdb_analysis(request: PDBRequest):
        try:
            data, source_name, fetch_metadata = fetch_rcsb(request.pdb_id, assembly_id=request.assembly_id)
            source_format = "mmcif"
            parse_text(data.decode("utf-8", errors="replace"), source_name=source_name, source_format=source_format)
        except Exception as exc:
            raise HTTPException(400, str(exc)) from exc
        analysis = store.create_analysis(source_kind="pdb_id", source_name=source_name,
                                         source_format=source_format, source_bytes=data)
        store.update_analysis(analysis["id"], metadata_json=json.dumps({"rcsb": fetch_metadata}))
        job = manager.submit_analysis(analysis["id"])
        return {"analysis_id": analysis["id"], "job_id": job["id"], "status_url": f"/api/jobs/{job['id']}"}

    @app.get("/api/rcsb/{pdb_id}/metadata")
    def pdb_metadata(pdb_id: str):
        try:
            return rcsb_entry_metadata(pdb_id)
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc

    @app.get("/api/analyses/{analysis_id}")
    def analysis_detail(analysis_id: str):
        item = store.get_analysis(analysis_id)
        if not item:
            raise HTTPException(404, "Analysis not found")
        return {**item, "jobs": store.list_jobs(analysis_id)}

    @app.get("/api/analyses/{analysis_id}/result")
    def analysis_result(analysis_id: str):
        item = store.get_analysis(analysis_id)
        if not item or not item.get("result_path"):
            raise HTTPException(404, "Interface results are not ready")
        result = json.loads(Path(item["result_path"]).read_text())
        # Results saved before chain names were read from the file show
        # "Chain A"; backfill them once from the stored coordinates.
        source_text = Path(item["source_path"]).read_text(errors="replace")
        if apply_source_chain_names(result, source_text, item["source_format"]):
            store.write_json(item["result_path"], result)
            store.update_analysis(analysis_id, metadata_json=json.dumps({**item["metadata"], **result["metadata"]}))
        return JSONResponse(result)

    @app.get("/api/analyses/{analysis_id}/curvature")
    def curvature_result(analysis_id: str):
        jobs = store.list_jobs(analysis_id)
        job = next((x for x in reversed(jobs) if x["kind"] == "curvature" and x.get("result_path")), None)
        if not job:
            raise HTTPException(404, "Curvature results are not ready")
        return JSONResponse(json.loads(Path(job["result_path"]).read_text()))

    @app.get("/api/analyses/{analysis_id}/pairs/{pair_id}/face-view")
    def face_view(analysis_id: str, pair_id: str):
        item = store.get_analysis(analysis_id)
        if not item:
            raise HTTPException(404, "Analysis not found")
        if not item.get("result_path"):
            raise HTTPException(409, "Interface results are not ready")
        result = json.loads(Path(item["result_path"]).read_text())
        pair = next((pair for pair in result.get("pairs", []) if pair["id"] == pair_id), None)
        if pair is None:
            raise HTTPException(404, "Interface pair not found")
        try:
            loaded = parse_text(Path(item["source_path"]).read_text(errors="replace"),
                                source_name=item["source_name"], source_format=item["source_format"])
            return binding_face_view(loaded.model, pair)
        except ValueError as exc:
            raise HTTPException(422, str(exc)) from exc

    @app.get("/api/analyses/{analysis_id}/pockets")
    def pocket_history(analysis_id: str):
        item = store.get_analysis(analysis_id)
        if not item:
            raise HTTPException(404, "Analysis not found")
        jobs = [job for job in store.list_jobs(analysis_id) if job["kind"] == "pocket"]
        return {"jobs": jobs}

    @app.get("/api/analyses/{analysis_id}/structure")
    def structure(analysis_id: str):
        item = store.get_analysis(analysis_id)
        if not item:
            raise HTTPException(404, "Analysis not found")
        media = "chemical/x-pdb" if item["source_format"] == "pdb" else "chemical/x-mmcif"
        return PlainTextResponse(Path(item["source_path"]).read_text(errors="replace"), media_type=media)

    @app.post("/api/analyses/{analysis_id}/pockets", status_code=202)
    def pocket_job(analysis_id: str, request: PocketRequest):
        item = store.get_analysis(analysis_id)
        if not item:
            raise HTTPException(404, "Analysis not found")
        if not item.get("result_path"):
            raise HTTPException(409, "Wait for interface analysis to finish before measuring a pocket")
        payload = request.model_dump(exclude_none=True)
        job = manager.submit_pocket(analysis_id, payload)
        return {"analysis_id": analysis_id, "job_id": job["id"], "status_url": f"/api/jobs/{job['id']}"}

    @app.get("/api/jobs/{job_id}")
    def job_detail(job_id: str):
        job = store.get_job(job_id)
        if not job:
            raise HTTPException(404, "Job not found")
        return job

    @app.get("/api/jobs/{job_id}/result")
    def job_result(job_id: str):
        job = store.get_job(job_id)
        if not job or not job.get("result_path"):
            raise HTTPException(404, "Job result is not ready")
        return JSONResponse(json.loads(Path(job["result_path"]).read_text()))

    @app.post("/api/jobs/{job_id}/cancel")
    def cancel_job(job_id: str):
        job = store.get_job(job_id)
        if not job:
            raise HTTPException(404, "Job not found")
        if job["status"] in {"complete", "failed", "cancelled", "interrupted"}:
            return job
        store.update_job(job_id, status="cancelled", stage="cancelled", error="Cancelled by user")
        # Core cancellation changes the analysis state as well; otherwise a
        # queued job that never starts would leave history showing “running”.
        analysis = store.get_analysis(job["analysis_id"])
        if analysis and job["kind"] == "core" and not analysis.get("result_path"):
            store.update_analysis(job["analysis_id"], status="cancelled", error="Interface analysis was cancelled")
        elif analysis and job["kind"] == "curvature" and analysis.get("result_path"):
            # Core results are already usable when only the optional curvature
            # stage is cancelled.
            store.update_analysis(job["analysis_id"], status="complete", error=None)
        return store.get_job(job_id)

    @app.post("/api/jobs/{job_id}/retry", status_code=202)
    def retry_job(job_id: str):
        job = store.get_job(job_id)
        if not job:
            raise HTTPException(404, "Job not found")
        if job["status"] not in {"failed", "interrupted", "cancelled"}:
            raise HTTPException(409, "Only failed, interrupted, or cancelled jobs can be retried")
        if job["kind"] == "core":
            replacement = manager.submit_analysis(job["analysis_id"])
        elif job["kind"] == "pocket":
            request = json.loads(job.get("request_json") or "{}")
            if not request:
                raise HTTPException(409, "This pocket job did not retain its request")
            replacement = manager.submit_pocket(job["analysis_id"], request)
        elif job["kind"] == "curvature":
            replacement = manager.submit_curvature(job["analysis_id"])
        else:
            raise HTTPException(400, f"Unsupported job type {job['kind']!r}")
        return {"job_id": replacement["id"], "status_url": f"/api/jobs/{replacement['id']}"}

    @app.get("/api/analyses/{analysis_id}/export")
    def export_analysis(analysis_id: str):
        item = store.get_analysis(analysis_id)
        if not item:
            raise HTTPException(404, "Analysis not found")
        archive = io.BytesIO()
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as zipped:
            source = Path(item["source_path"])
            zipped.write(source, f"source/{source.name}")
            if item.get("result_path"):
                zipped.write(item["result_path"], "results/analysis.json")
            for job in store.list_jobs(analysis_id):
                if job.get("result_path"):
                    zipped.write(job["result_path"], f"results/{Path(job['result_path']).name}")
            zipped.writestr("README.txt", "Protein interface analyzer export.\n")
        archive.seek(0)
        return StreamingResponse(archive, media_type="application/zip",
                                 headers={"Content-Disposition": f'attachment; filename="{analysis_id}.zip"'})

    def has_active_jobs(analysis_id: str) -> bool:
        return any(job["status"] in ACTIVE_JOB_STATES for job in store.list_jobs(analysis_id))

    @app.delete("/api/analyses", status_code=200)
    def clear_history():
        # Workers write into each analysis directory, so active analyses are
        # kept rather than deleted out from under a running process.
        deleted, skipped = [], []
        for item in store.list_analyses(limit=-1):
            if has_active_jobs(item["id"]):
                skipped.append(item["id"])
            elif store.delete_analysis(item["id"]):
                deleted.append(item["id"])
        return {"deleted": deleted, "skipped": skipped}

    @app.delete("/api/analyses/{analysis_id}", status_code=204)
    def delete_analysis(analysis_id: str):
        if not store.get_analysis(analysis_id):
            raise HTTPException(404, "Analysis not found")
        if has_active_jobs(analysis_id):
            raise HTTPException(409, "Cancel or wait for running jobs before deleting this analysis")
        store.delete_analysis(analysis_id)
        return None

    app.mount("/static", StaticFiles(directory=static_dir), name="static")
    return app


app = create_app()
