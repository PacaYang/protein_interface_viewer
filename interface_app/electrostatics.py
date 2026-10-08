"""APBS/PDB2PQR electrostatic potentials sampled on surface meshes.

The curvature worker produces the surface geometry used by the browser.  This
module keeps electrostatics as a separate, optional calculation: it prepares a
PQR with PDB2PQR, solves the linearized Poisson--Boltzmann equation with APBS,
and interpolates the resulting OpenDX map at the already serialized mesh
vertices.

The external programs are deliberately invoked through argument lists.  They
can be supplied with ``APBS_BIN`` and ``PDB2PQR_BIN`` or provisioned from the
checked-in tool manifest into the analysis data directory.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request
import zipfile
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from Bio.PDB import PDBIO, Model, Chain, Residue

from .constants import is_amino_acid, is_heavy_atom


PH = 7.4
FORCE_FIELD = "AMBER"
TEMPERATURE_K = 298.15
PROTEIN_DIELECTRIC = 2.0
SOLVENT_DIELECTRIC = 78.5
IONIC_STRENGTH_M = 0.15
SOLVENT_RADIUS_A = 1.4
GRID_MARGIN_A = 12.0
FINE_GRID_MARGIN_A = 8.0
GRID_DIME = (129, 129, 129)
COLOR_QUANTILE = 0.98
DEFAULT_TIMEOUT_S = 15 * 60
MANIFEST_PATH = Path(__file__).with_name("tool_manifest.json")


class ElectrostaticsUnavailable(RuntimeError):
    """Expected setup or numerical failure that should not hide convexity."""


class ElectrostaticsCancelled(RuntimeError):
    """The caller cancelled this optional job."""


def _check_cancel(cancelled: Callable | None):
    if cancelled and cancelled():
        raise ElectrostaticsCancelled("Electrostatics was cancelled")


def _platform_key() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "linux" and machine in {"x86_64", "amd64"}:
        return "linux-x86_64"
    return f"{system}-{machine}"


def _load_manifest(path: Path = MANIFEST_PATH) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ElectrostaticsUnavailable(f"Could not read the electrostatics tool manifest: {exc}") from exc


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download_artifact(artifact: dict, target_dir: Path, cancelled: Callable | None = None) -> Path:
    """Download and extract one manifest artifact atomically."""
    import fcntl  # Managed downloads currently support Linux only.

    target_dir.mkdir(parents=True, exist_ok=True)
    lock_path = target_dir / f"{artifact['name']}-{artifact['version']}.lock"
    with lock_path.open("w") as lock:
        while True:
            _check_cancel(cancelled)
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(0.2)
        return _install_artifact(artifact, target_dir, cancelled)


def _install_artifact(artifact: dict, target_dir: Path, cancelled: Callable | None) -> Path:
    name = str(artifact.get("name") or "tool")
    version = str(artifact.get("version") or "unknown")
    url = artifact.get("url")
    expected = str(artifact.get("sha256") or "").lower()
    if not str(url).startswith("https://") or not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ElectrostaticsUnavailable(
            f"No pinned download is configured for {name}; set the corresponding *_BIN override."
        )
    target_dir.mkdir(parents=True, exist_ok=True)
    install_dir = target_dir / f"{name}-{version}"
    executable = artifact.get("executable")
    executable_path = install_dir / executable if executable else None
    marker = install_dir / "installed.json"
    if executable_path and executable_path.is_file() and marker.is_file():
        try:
            if json.loads(marker.read_text()).get("sha256") == expected:
                return executable_path
        except (OSError, ValueError):
            pass

    with tempfile.TemporaryDirectory(prefix=f".{name}-", dir=target_dir) as raw_tmp:
        tmp = Path(raw_tmp)
        archive = tmp / "download"
        try:
            request = urllib.request.Request(str(url), headers={"User-Agent": "protein-interface-viewer/0.3"})
            with urllib.request.urlopen(request, timeout=120) as response, archive.open("wb") as handle:
                if not response.geturl().startswith("https://"):
                    raise ElectrostaticsUnavailable("Tool downloads must use HTTPS.")
                while True:
                    _check_cancel(cancelled)
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
        except ElectrostaticsCancelled:
            raise
        except Exception as exc:
            raise ElectrostaticsUnavailable(f"Could not download {name} {version}: {exc}") from exc
        actual = _sha256(archive)
        if actual != expected:
            raise ElectrostaticsUnavailable(
                f"Checksum mismatch while downloading {name} {version}: expected {expected}, got {actual}."
            )
        extracted = tmp / "extracted"
        extracted.mkdir()
        try:
            if str(url).lower().endswith((".zip", ".whl")):
                with zipfile.ZipFile(archive) as zipped:
                    root = extracted.resolve()
                    for member in zipped.infolist():
                        target = (extracted / member.filename).resolve()
                        if target != root and root not in target.parents:
                            raise ElectrostaticsUnavailable("The downloaded tool archive contains an unsafe path.")
                    zipped.extractall(extracted)
            else:
                with tarfile.open(archive, "r:*") as tar:
                    members = tar.getmembers()
                    root = extracted.resolve()
                    for member in members:
                        target = (extracted / member.name).resolve()
                        if (target != root and root not in target.parents) or not (member.isfile() or member.isdir()):
                            raise ElectrostaticsUnavailable("The downloaded tool archive contains an unsafe path.")
                    try:
                        tar.extractall(extracted, filter="data")
                    except TypeError:  # Python 3.10/3.11
                        tar.extractall(extracted)
        except Exception as exc:
            raise ElectrostaticsUnavailable(f"Could not unpack {name} {version}: {exc}") from exc
        if not executable or not (extracted / executable).is_file():
            raise ElectrostaticsUnavailable(f"The downloaded {name} archive does not contain {executable}.")
        (extracted / executable).chmod((extracted / executable).stat().st_mode | 0o111)
        (extracted / "installed.json").write_text(json.dumps({"sha256": expected, "version": version}))
        _check_cancel(cancelled)
        if install_dir.exists():
            shutil.rmtree(install_dir)
        os.replace(extracted, install_dir)
    if not executable_path or not executable_path.is_file():
        raise ElectrostaticsUnavailable(
            f"The downloaded {name} {version} archive did not contain {executable or 'the configured executable'}."
        )
    executable_path.chmod(executable_path.stat().st_mode | 0o111)
    return executable_path


def resolve_tools(data_root: str | Path | None = None, cancelled: Callable | None = None) -> tuple[list[str], list[str]]:
    """Return ``(pdb2pqr, apbs)`` executable paths.

    Explicit environment overrides take precedence.  Otherwise standard names
    are used when already available, followed by the pinned managed manifest.
    """

    def override(name: str, fallback_names: Iterable[str]) -> list[str] | None:
        value = os.environ.get(name)
        if value:
            candidate = Path(shutil.which(value) or value).expanduser().resolve()
            if not candidate.is_file() or not os.access(candidate, os.X_OK):
                raise ElectrostaticsUnavailable(f"{name} points to a non-executable file: {candidate}")
            return [str(candidate)]
        for executable in fallback_names:
            found = shutil.which(executable)
            if found:
                return [found]
        return None

    pdb2pqr = override("PDB2PQR_BIN", ("pdb2pqr30", "pdb2pqr"))
    apbs = override("APBS_BIN", ("apbs",))
    if pdb2pqr and apbs:
        return pdb2pqr, apbs

    manifest = _load_manifest()
    platform_entry = (manifest.get("platforms") or {}).get(_platform_key()) or {}
    if not platform_entry:
        raise ElectrostaticsUnavailable(
            f"Automatic APBS/PDB2PQR downloads support Linux x86_64; this host is {_platform_key()}. "
            "Supply both APBS_BIN and PDB2PQR_BIN to use installed tools."
        )
    root = (Path(data_root or os.environ.get("INTERFACE_APP_DATA", ".interface_app_data")) / "tools").resolve()
    missing = []
    if not pdb2pqr:
        artifact = platform_entry.get("pdb2pqr")
        if artifact:
            executable = _download_artifact({"name": "pdb2pqr", **artifact}, root, cancelled)
            paths = [str(executable.parent.parent)]
            for dependency in artifact.get("dependencies", []):
                _download_artifact(dependency, root, cancelled)
                paths.append(str(root / f"{dependency['name']}-{dependency['version']}"))
            pdb2pqr = [sys.executable, str(Path(__file__).with_name("_pdb2pqr_runner.py")), json.dumps(paths)]
        else:
            missing.append("PDB2PQR")
    if not apbs:
        artifact = platform_entry.get("apbs")
        if artifact:
            apbs = [str(_download_artifact({"name": "apbs", **artifact}, root, cancelled))]
        else:
            missing.append("APBS")
    if missing:
        raise ElectrostaticsUnavailable(
            f"{', '.join(missing)} is unavailable. Install the pinned tools or set "
            "PDB2PQR_BIN and APBS_BIN."
        )
    return pdb2pqr, apbs


def _write_pdb(model, path: Path) -> None:
    # PDB2PQR needs PDB identifiers. Renumber a filtered copy so mmCIF author
    # chain IDs, insertion codes and large residue numbers do not hit PDB's
    # field limits. Potential sampling uses coordinates, not these identifiers.
    chains = [chain for chain in model if any(is_amino_acid(r) for r in chain)]
    identifiers = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if len(chains) > len(identifiers):
        raise ElectrostaticsUnavailable("PDB2PQR preparation supports at most 62 protein chains.")
    prepared = Model.Model(0)
    for chain, ident in zip(chains, identifiers):
        new_chain = Chain.Chain(ident)
        for index, residue in enumerate((r for r in chain if is_amino_acid(r)), 1):
            if index > 9999:
                raise ElectrostaticsUnavailable("PDB2PQR preparation supports at most 9999 residues per chain.")
            # Copy atoms into a fresh residue. Deep-copying a residue also
            # follows Biopython's parent pointers through the entire assembly.
            copied = Residue.Residue((" ", index, " "), residue.resname, residue.segid)
            for atom in residue:
                if is_heavy_atom(atom):
                    selected = atom.selected_child if atom.is_disordered() == 2 else atom
                    copied.add(selected.copy())
            new_chain.add(copied)
        prepared.add(new_chain)
    writer = PDBIO()
    writer.set_structure(prepared)
    writer.save(str(path))


def _mesh_vertices(surface_result: dict) -> dict[str, np.ndarray]:
    meshes = surface_result.get("meshes") or {}
    result = {}
    for chain, mesh in meshes.items():
        position = np.asarray(mesh.get("position") or [], dtype=np.float64)
        if position.size == 0 or position.size % 3 or not np.isfinite(position).all():
            raise ElectrostaticsUnavailable(f"Surface mesh for chain {chain} has invalid vertex coordinates.")
        result[str(chain)] = position.reshape((-1, 3))
    if not result:
        raise ElectrostaticsUnavailable("The curvature result contains no surface meshes.")
    return result


def _bounds(vertices: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(list(vertices.values()), axis=0)
    low = points.min(axis=0) - GRID_MARGIN_A
    high = points.max(axis=0) + GRID_MARGIN_A
    return low, high


def _apbs_input(path: Path, pqr_name: str, low: np.ndarray, high: np.ndarray, output_name: str) -> None:
    extent = high - low
    center = (low + high) / 2
    fine_extent = np.maximum(extent - 2 * (GRID_MARGIN_A - FINE_GRID_MARGIN_A), 8.0)
    values = " ".join(str(value) for value in GRID_DIME)
    coarse = " ".join(f"{value:.3f}" for value in (extent + 4.0))
    fine = " ".join(f"{value:.3f}" for value in fine_extent)
    center_text = " ".join(f"{value:.3f}" for value in center)
    path.write_text(
        "\n".join([
            "read",
            f"  mol pqr {pqr_name}",
            "end",
            "elec",
            "  mg-auto",
            f"  dime {values}",
            f"  cglen {coarse}",
            f"  fglen {fine}",
            f"  cgcent {center_text}",
            f"  fgcent {center_text}",
            "  lpbe",
            "  bcfl mdh",
            f"  pdie {PROTEIN_DIELECTRIC:g}",
            f"  sdie {SOLVENT_DIELECTRIC:g}",
            f"  srad {SOLVENT_RADIUS_A:g}",
            "  swin 0.3",
            f"  temp {TEMPERATURE_K:g}",
            "  chgm spl2",
            "  srfm smol",
            "  sdens 10.0",
            f"  ion charge 1 conc {IONIC_STRENGTH_M:g} radius 2.0",
            f"  ion charge -1 conc {IONIC_STRENGTH_M:g} radius 1.8",
            "  mol 1",
            f"  write pot dx {output_name}",
            "end",
            "quit",
            "",
        ])
    )


_COUNTS_RE = re.compile(r"class\s+gridpositions\s+counts\s+(\d+)\s+(\d+)\s+(\d+)")
_ORIGIN_RE = re.compile(r"^origin\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$", re.MULTILINE)
_DELTA_RE = re.compile(r"^delta\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$", re.MULTILINE)
_ITEMS_RE = re.compile(r"items\s+(\d+)\s+data\s+follows", re.IGNORECASE)


def parse_open_dx(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read an APBS OpenDX scalar grid as ``(origin, spacing, values)``."""

    text = Path(path).read_text(errors="replace")
    counts_match = _COUNTS_RE.search(text)
    origin_match = _ORIGIN_RE.search(text)
    deltas = _DELTA_RE.findall(text)
    items_match = _ITEMS_RE.search(text)
    if not counts_match or not origin_match or len(deltas) < 3 or not items_match:
        raise ElectrostaticsUnavailable("APBS returned an unsupported or incomplete OpenDX potential map.")
    counts = np.asarray([int(value) for value in counts_match.groups()], dtype=int)
    origin = np.asarray([float(value) for value in origin_match.groups()], dtype=float)
    spacing = []
    for raw in deltas[:3]:
        vector = np.asarray([float(value) for value in raw], dtype=float)
        axis = int(np.argmax(np.abs(vector)))
        if np.count_nonzero(np.abs(vector) > 1e-10) != 1 or axis != len(spacing) or vector[axis] <= 0:
            raise ElectrostaticsUnavailable("APBS returned a rotated or non-positive OpenDX grid.")
        spacing.append(vector[axis])
    count = int(items_match.group(1))
    tokens = []
    for line in text[items_match.end():].splitlines():
        tokens.extend(line.split("#", 1)[0].split())
        if len(tokens) >= count:
            break
    try:
        values = np.asarray(tokens[:count], dtype=float)
    except ValueError as exc:
        raise ElectrostaticsUnavailable("APBS returned non-numeric OpenDX data.") from exc
    if len(values) != count or count != int(np.prod(counts)) or np.any(counts < 2):
        raise ElectrostaticsUnavailable("APBS OpenDX data length does not match its grid dimensions.")
    if not np.isfinite(values).all() or not np.isfinite(origin).all() or not np.isfinite(spacing).all():
        raise ElectrostaticsUnavailable("APBS returned non-finite potential or grid coordinates.")
    return origin, np.asarray(spacing, dtype=float), values[:count].reshape(tuple(counts))


def interpolate_grid(points: np.ndarray, origin: np.ndarray, spacing: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Trilinearly sample an axis-aligned OpenDX grid at ``points``."""

    coordinates = (np.asarray(points, dtype=float) - origin) / spacing
    if np.any(coordinates < -1e-6) or np.any(coordinates > (np.asarray(values.shape) - 1) + 1e-6):
        raise ElectrostaticsUnavailable("At least one surface vertex lies outside the APBS potential grid.")
    coordinates = np.clip(coordinates, 0, np.asarray(values.shape) - 1)
    lower = np.floor(coordinates).astype(int)
    upper = np.minimum(lower + 1, np.asarray(values.shape) - 1)
    fraction = coordinates - lower
    result = np.zeros(len(points), dtype=float)
    for bits in range(8):
        corner = np.asarray([(bits >> axis) & 1 for axis in range(3)], dtype=int)
        index = tuple(np.where(corner[axis], upper[:, axis], lower[:, axis]) for axis in range(3))
        weight = np.ones(len(points), dtype=float)
        for axis in range(3):
            weight *= fraction[:, axis] if corner[axis] else (1 - fraction[:, axis])
        result += weight * values[index]
    return result


def _color_limit(values: list[float]) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return 1.0
    limit = float(np.quantile(np.abs(finite), COLOR_QUANTILE))
    if limit <= 1e-12:
        limit = float(np.max(np.abs(finite)))
    return round(max(limit, 1e-6), 6)


def _stop_process(process: subprocess.Popen):
    def terminate(sig):
        try:
            if os.name == "posix":
                os.killpg(process.pid, sig)
            elif sig == signal.SIGTERM:
                process.terminate()
            else:
                process.kill()
        except ProcessLookupError:
            pass

    terminate(signal.SIGTERM)
    try:
        process.communicate(timeout=2)
    except subprocess.TimeoutExpired:
        terminate(signal.SIGKILL)
        process.communicate()


def _run(command: list[str], *, cwd: Path, timeout: float, runner: Callable | None = None,
         cancelled: Callable | None = None) -> subprocess.CompletedProcess:
    _check_cancel(cancelled)
    if runner is not None:
        return runner(command, cwd=cwd, timeout=timeout)
    try:
        env = {**os.environ, "OMP_NUM_THREADS": "1"}
        with subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, errors="replace", env=env, start_new_session=os.name == "posix") as process:
            deadline = time.monotonic() + timeout
            try:
                while True:
                    _check_cancel(cancelled)
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        raise ElectrostaticsUnavailable(
                            f"Electrostatics command timed out after {timeout:g} seconds: {command[0]}"
                        )
                    try:
                        stdout, stderr = process.communicate(timeout=min(0.25, remaining))
                        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
                    except subprocess.TimeoutExpired:
                        continue
            except BaseException:
                _stop_process(process)
                raise
    except OSError as exc:
        raise ElectrostaticsUnavailable(f"Could not run {command[0]}: {exc}") from exc


def compute_electrostatics(
    model,
    analysis: dict,
    surface_result: dict,
    *,
    data_root: str | Path | None = None,
    runner: Callable | None = None,
    timeout: float = DEFAULT_TIMEOUT_S,
    progress: Callable | None = None,
    cancelled: Callable | None = None,
) -> dict:
    """Calculate APBS potential values at every curvature mesh vertex."""

    base = {
        "version": "0.2.0",
        "status": "unavailable",
        "reason": None,
        "parameters": {
            "model": "APBS/PDB2PQR",
            "force_field": FORCE_FIELD,
            "pH": PH,
            "temperature_K": TEMPERATURE_K,
            "protein_dielectric": PROTEIN_DIELECTRIC,
            "solvent_dielectric": SOLVENT_DIELECTRIC,
            "ionic_strength_M": IONIC_STRENGTH_M,
            "units": "kT/e",
        },
        "report": None,
        "meshes": {},
    }
    try:
        def stage(name, amount):
            _check_cancel(cancelled)
            if progress:
                progress(name, amount)

        vertices = _mesh_vertices(surface_result)
        stage("provisioning APBS/PDB2PQR tools", 0.05)
        pdb2pqr, apbs = resolve_tools(data_root, cancelled)
        pdb2pqr = [pdb2pqr] if isinstance(pdb2pqr, str) else list(pdb2pqr)
        apbs = [apbs] if isinstance(apbs, str) else list(apbs)
        low, high = _bounds(vertices)
        temporary_parent = Path(data_root or os.environ.get("INTERFACE_APP_DATA", ".interface_app_data"))
        temporary_parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="electrostatics-", dir=temporary_parent) as raw_tmp:
            work = Path(raw_tmp)
            input_pdb = work / "structure.pdb"
            output_pqr = work / "structure.pqr"
            apbs_input = work / "apbs.in"
            output_base = work / "potential"
            output_dx = work / "potential.dx"
            stage("PDB2PQR preparation at pH 7.4", 0.20)
            _write_pdb(model, input_pdb)
            pqr_result = _run(
                [*pdb2pqr, f"--ff={FORCE_FIELD}", "--titration-state-method=propka", f"--with-ph={PH:g}", str(input_pdb), str(output_pqr)],
                cwd=work, timeout=timeout, runner=runner, cancelled=cancelled,
            )
            if pqr_result.returncode != 0 or not output_pqr.exists():
                details = (pqr_result.stderr or pqr_result.stdout or "no diagnostic output").strip()[-3000:]
                raise ElectrostaticsUnavailable(f"PDB2PQR failed: {details}")
            _apbs_input(apbs_input, output_pqr.name, low, high, output_base.name)
            stage("APBS potential solve", 0.55)
            apbs_result = _run([*apbs, str(apbs_input)], cwd=work, timeout=timeout, runner=runner, cancelled=cancelled)
            if apbs_result.returncode != 0 or not output_dx.exists():
                details = (apbs_result.stderr or apbs_result.stdout or "no diagnostic output").strip()[-3000:]
                raise ElectrostaticsUnavailable(f"APBS failed: {details}")
            stage("reading APBS OpenDX map", 0.85)
            origin, spacing, grid = parse_open_dx(output_dx)
            all_values = []
            sampled = {}
            stage("interpolating surface potentials", 0.92)
            for chain, points in vertices.items():
                values = interpolate_grid(points, origin, spacing, grid)
                if not np.isfinite(values).all():
                    raise ElectrostaticsUnavailable("APBS interpolation produced non-finite potentials.")
                sampled[chain] = np.round(values, 6).tolist()
                all_values.extend(values.tolist())
        finite = np.asarray(all_values, dtype=float)
        base["status"] = "complete"
        base["report"] = {
            "units": "kT/e",
            "color_limit_kT_e": _color_limit(all_values),
            "min_kT_e": round(float(np.min(finite)), 6),
            "max_kT_e": round(float(np.max(finite)), 6),
            "vertex_count": int(len(finite)),
        }
        base["meshes"] = {
            chain: {"potential_kT_e": values}
            for chain, values in sampled.items()
        }
        _check_cancel(cancelled)
        return base
    except ElectrostaticsCancelled:
        raise
    except ElectrostaticsUnavailable as exc:
        base["reason"] = str(exc)
        return base
    except Exception as exc:  # numerical adapters must never hide core results
        base["reason"] = f"Electrostatics calculation failed: {exc}"
        return base
