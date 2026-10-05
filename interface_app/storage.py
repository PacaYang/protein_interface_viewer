"""Small SQLite-backed history and artifact store for local analyses."""

from __future__ import annotations

import json
import os
import sqlite3
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class Store:
    def __init__(self, root: str | Path | None = None):
        self.root = Path(root or os.environ.get("INTERFACE_APP_DATA", ".interface_app_data")).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        (self.root / "analyses").mkdir(exist_ok=True)
        self.db_path = self.root / "history.sqlite3"
        self._init_db()

    def _connect(self):
        connection = sqlite3.connect(self.db_path, timeout=30)
        connection.row_factory = sqlite3.Row
        return connection

    def _init_db(self):
        with self._connect() as db:
            db.executescript(
                """
                CREATE TABLE IF NOT EXISTS analyses (
                    id TEXT PRIMARY KEY,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    source_kind TEXT NOT NULL,
                    source_name TEXT NOT NULL,
                    source_path TEXT NOT NULL,
                    source_format TEXT NOT NULL,
                    status TEXT NOT NULL,
                    metadata_json TEXT NOT NULL DEFAULT '{}',
                    result_path TEXT,
                    error TEXT
                );
                CREATE TABLE IF NOT EXISTS jobs (
                    id TEXT PRIMARY KEY,
                    analysis_id TEXT NOT NULL,
                    kind TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    status TEXT NOT NULL,
                    progress REAL NOT NULL DEFAULT 0,
                    stage TEXT NOT NULL DEFAULT 'queued',
                    result_path TEXT,
                    error TEXT,
                    request_json TEXT NOT NULL DEFAULT '{}',
                    FOREIGN KEY(analysis_id) REFERENCES analyses(id) ON DELETE CASCADE
                );
                CREATE INDEX IF NOT EXISTS jobs_analysis_idx ON jobs(analysis_id);
                """
            )
            columns = {row["name"] for row in db.execute("PRAGMA table_info(jobs)").fetchall()}
            if "request_json" not in columns:
                db.execute("ALTER TABLE jobs ADD COLUMN request_json TEXT NOT NULL DEFAULT '{}'")

    def analysis_dir(self, analysis_id: str) -> Path:
        path = self.root / "analyses" / analysis_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    def create_analysis(self, *, source_kind: str, source_name: str, source_format: str, source_bytes: bytes) -> dict:
        analysis_id = uuid4().hex
        suffix = ".cif" if source_format == "mmcif" else ".pdb"
        directory = self.analysis_dir(analysis_id)
        source_path = directory / f"source{suffix}"
        source_path.write_bytes(source_bytes)
        now = utc_now()
        row = {
            "id": analysis_id, "created_at": now, "updated_at": now,
            "source_kind": source_kind, "source_name": source_name,
            "source_path": str(source_path), "source_format": source_format,
            "status": "queued", "metadata": {}, "result_path": None, "error": None,
        }
        with self._connect() as db:
            db.execute(
                "INSERT INTO analyses VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (analysis_id, now, now, source_kind, source_name, str(source_path), source_format,
                 "queued", "{}", None, None),
            )
        return row

    def get_analysis(self, analysis_id: str) -> dict | None:
        with self._connect() as db:
            row = db.execute("SELECT * FROM analyses WHERE id = ?", (analysis_id,)).fetchone()
        return self._analysis_row(row) if row else None

    def list_analyses(self, limit: int = 30) -> list[dict]:
        with self._connect() as db:
            rows = db.execute("SELECT * FROM analyses ORDER BY created_at DESC LIMIT ?", (int(limit),)).fetchall()
        return [self._analysis_row(row) for row in rows]

    @staticmethod
    def _analysis_row(row) -> dict:
        item = dict(row)
        item["metadata"] = json.loads(item.pop("metadata_json") or "{}")
        return item

    def update_analysis(self, analysis_id: str, **changes):
        allowed = {"status", "metadata_json", "result_path", "error", "source_name"}
        invalid = set(changes) - allowed
        if invalid:
            raise ValueError(f"Unsupported analysis fields: {sorted(invalid)}")
        changes["updated_at"] = utc_now()
        assignments = ", ".join(f"{key} = ?" for key in changes)
        with self._connect() as db:
            db.execute(
                f"UPDATE analyses SET {assignments} WHERE id = ?",
                (*changes.values(), analysis_id),
            )

    def create_job(self, analysis_id: str, kind: str, request: dict | None = None) -> dict:
        job_id = uuid4().hex
        now = utc_now()
        with self._connect() as db:
            db.execute(
                "INSERT INTO jobs (id, analysis_id, kind, created_at, updated_at, status, progress, stage, result_path, error, request_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (job_id, analysis_id, kind, now, now, "queued", 0.0, "queued", None, None,
                 json.dumps(request or {})),
            )
        return self.get_job(job_id)

    def get_job(self, job_id: str) -> dict | None:
        with self._connect() as db:
            row = db.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
        return dict(row) if row else None

    def list_jobs(self, analysis_id: str) -> list[dict]:
        with self._connect() as db:
            rows = db.execute("SELECT * FROM jobs WHERE analysis_id = ? ORDER BY created_at", (analysis_id,)).fetchall()
        return [dict(row) for row in rows]

    def update_job(self, job_id: str, **changes):
        allowed = {"status", "progress", "stage", "result_path", "error"}
        invalid = set(changes) - allowed
        if invalid:
            raise ValueError(f"Unsupported job fields: {sorted(invalid)}")
        changes["updated_at"] = utc_now()
        assignments = ", ".join(f"{key} = ?" for key in changes)
        with self._connect() as db:
            db.execute(f"UPDATE jobs SET {assignments} WHERE id = ?", (*changes.values(), job_id))

    def write_json(self, path: str | Path, value) -> str:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return str(target)

    def delete_analysis(self, analysis_id: str) -> bool:
        row = self.get_analysis(analysis_id)
        if not row:
            return False
        with self._connect() as db:
            db.execute("DELETE FROM jobs WHERE analysis_id = ?", (analysis_id,))
            db.execute("DELETE FROM analyses WHERE id = ?", (analysis_id,))
        directory = self.root / "analyses" / analysis_id
        if directory.exists():
            import shutil
            shutil.rmtree(directory)
        return True
