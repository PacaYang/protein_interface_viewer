from pathlib import Path
from subprocess import CompletedProcess
import hashlib
import io
import sys
import urllib.error
import zipfile

import numpy as np
import pytest

from interface_app import electrostatics
from interface_app.structure import parse_structure

ROOT = Path(__file__).resolve().parents[2]


def dx_text():
    return """object 1 class gridpositions counts 2 2 2
origin 0 0 0
delta 1 0 0
delta 0 1 0
delta 0 0 1
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
0 1 2 3 4 5 6 7
attribute \"dep\" string \"positions\"
"""


def test_parse_open_dx_and_trilinear_interpolation(tmp_path):
    path = tmp_path / "potential.dx"
    path.write_text(dx_text())
    origin, spacing, values = electrostatics.parse_open_dx(path)

    assert origin.tolist() == [0.0, 0.0, 0.0]
    assert spacing.tolist() == [1.0, 1.0, 1.0]
    assert values.shape == (2, 2, 2)
    sampled = electrostatics.interpolate_grid(
        np.asarray([[0.5, 0.5, 0.5], [1.0, 1.0, 1.0]]),
        origin, spacing, values,
    )
    assert sampled.tolist() == [3.5, 7.0]


def test_apbs_input_contains_documented_defaults(tmp_path):
    path = tmp_path / "apbs.in"
    electrostatics._apbs_input(
        path, "structure.pqr", np.asarray([-10.0, -11.0, -12.0]), np.asarray([10.0, 11.0, 12.0]), "potential",
    )
    text = path.read_text()
    assert "mol pqr structure.pqr" in text
    assert "lpbe" in text
    assert "pdie 2" in text
    assert "sdie 78.5" in text
    assert "temp 298.15" in text
    assert "ion charge 1 conc 0.15" in text
    assert "write pot dx potential" in text


def test_compute_electrostatics_samples_every_mesh_vertex(monkeypatch, tmp_path):
    monkeypatch.setattr(electrostatics, "resolve_tools", lambda _root, cancelled: (["pdb2pqr"], ["apbs"]))
    monkeypatch.setattr(electrostatics, "_write_pdb", lambda _model, path: path.write_text("ATOM\n"))

    commands = []
    stages = []

    def runner(command, cwd, timeout):
        commands.append(command)
        if command[0] == "pdb2pqr":
            (cwd / "structure.pqr").write_text("PQR\n")
        else:
            (cwd / "potential.dx").write_text(dx_text())
        return CompletedProcess(command, 0, "", "")

    surface = {
        "meshes": {
            "A": {"position": [0, 0, 0, 1, 0, 0, 0, 1, 1]},
            "B": {"position": [1, 1, 1]},
        }
    }
    result = electrostatics.compute_electrostatics(
        object(), {}, surface, data_root=tmp_path, runner=runner,
        progress=lambda name, amount: stages.append(name),
    )

    assert result["status"] == "complete"
    assert result["parameters"]["pH"] == 7.4
    assert result["report"]["units"] == "kT/e"
    assert len(result["meshes"]["A"]["potential_kT_e"]) == 3
    assert len(result["meshes"]["B"]["potential_kT_e"]) == 1
    assert result["report"]["vertex_count"] == 4
    assert "--ff=AMBER" in commands[0]
    assert "--titration-state-method=propka" in commands[0]
    assert "--with-ph=7.4" in commands[0]
    assert len(stages) == 5
    assert not list(tmp_path.glob("electrostatics-*"))


def test_compute_electrostatics_reports_missing_tools(monkeypatch, tmp_path):
    monkeypatch.setattr(
        electrostatics,
        "resolve_tools",
        lambda _root, cancelled: (_ for _ in ()).throw(electrostatics.ElectrostaticsUnavailable("APBS missing")),
    )
    result = electrostatics.compute_electrostatics(
        object(), {}, {"meshes": {"A": {"position": [0, 0, 0]}}}, data_root=tmp_path,
    )
    assert result["status"] == "unavailable"
    assert result["reason"] == "APBS missing"


@pytest.mark.parametrize("replace, message", [
    (("delta 1 0 0", "delta 1 1 0"), "rotated"),
    (("delta 1 0 0", "delta -1 0 0"), "non-positive"),
    (("delta 0 0 1", "delta 0 1 0"), "rotated"),
    (("items 8", "items 9"), "non-numeric"),
    (("0 1 2 3 4 5 6 7", "0 1 2 3 4 5 6 nan"), "non-finite"),
    (("0 1 2 3 4 5 6 7", "0 1 invalid"), "non-numeric"),
])
def test_open_dx_rejects_malformed_maps(tmp_path, replace, message):
    path = tmp_path / "potential.dx"
    path.write_text(dx_text().replace(*replace))
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match=message):
        electrostatics.parse_open_dx(path)


def test_interpolation_matches_linear_field_and_rejects_out_of_bounds():
    origin = np.asarray([-4, -2, 3])
    spacing = np.asarray([0.5, 2, 1])
    x, y, z = np.indices((5, 4, 3))
    values = 2 * x + 3 * y - z
    fractions = np.asarray([[1.3, 2.7, 0.2], [4, 3, 2], [0, 0, 0]])
    result = electrostatics.interpolate_grid(origin + fractions * spacing, origin, spacing, values)
    np.testing.assert_allclose(result, 2 * fractions[:, 0] + 3 * fractions[:, 1] - fractions[:, 2])
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match="outside"):
        electrostatics.interpolate_grid(np.asarray([[99, 0, 0]]), origin, spacing, values)


def test_global_color_limit_handles_neutral_and_uniform_maps():
    assert electrostatics._color_limit([0, 0, 0]) > 0
    assert electrostatics._color_limit([-3, -3, -3]) == 3
    assert electrostatics._color_limit([-1, 1, 2, -2]) == 2


@pytest.mark.parametrize("fmt", ["pdb", "mmcif"])
def test_preparation_preserves_coordinates_and_handles_mmcif_identifiers(tmp_path, fmt):
    from Bio.PDB import MMCIFIO

    loaded = parse_structure(ROOT / "examples/1IAR.pdb")
    # Large author identifiers cannot be written directly to a PDB.
    if fmt == "mmcif":
        loaded.model["A"].id = "alpha"
        writer = MMCIFIO()
        writer.set_structure(loaded.structure)
        source = tmp_path / "source.cif"
        writer.save(str(source))
        loaded = parse_structure(source)
    before = [atom.coord.copy() for chain in loaded.model for residue in chain
              if electrostatics.is_amino_acid(residue) for atom in residue
              if electrostatics.is_heavy_atom(atom)]
    target = tmp_path / "prepared.pdb"
    electrostatics._write_pdb(loaded.model, target)
    prepared = parse_structure(target)
    after = [atom.coord for atom in prepared.model.get_atoms()]
    np.testing.assert_allclose(after, before, atol=0.001)
    assert len(prepared.metadata["chains"]) == 2
    assert len(before) == len(after)
    assert ("alpha" in loaded.model.child_dict) == (fmt == "mmcif")


@pytest.mark.parametrize("program", ["pdb2pqr", "apbs"])
def test_solver_failure_is_an_unavailable_artifact(monkeypatch, tmp_path, program):
    monkeypatch.setattr(electrostatics, "resolve_tools", lambda *_args: (["pdb2pqr"], ["apbs"]))
    monkeypatch.setattr(electrostatics, "_write_pdb", lambda _model, path: path.write_text("ATOM\n"))

    def runner(command, cwd, timeout):
        if command[0] == program:
            return CompletedProcess(command, 1, "", "fixture solver failure")
        (cwd / "structure.pqr").write_text("PQR\n")
        return CompletedProcess(command, 0, "", "")

    result = electrostatics.compute_electrostatics(
        object(), {}, {"meshes": {"A": {"position": [0, 0, 0]}}}, data_root=tmp_path, runner=runner,
    )
    assert result["status"] == "unavailable"
    assert "fixture solver failure" in result["reason"]
    assert not list(tmp_path.glob("electrostatics-*"))


def test_process_timeout_and_cancellation_terminate_solver(tmp_path):
    command = [sys.executable, "-c", "import time; time.sleep(10)"]
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match="timed out"):
        electrostatics._run(command, cwd=tmp_path, timeout=0.05)
    checks = []

    def cancelled():
        checks.append(True)
        return len(checks) >= 3

    with pytest.raises(electrostatics.ElectrostaticsCancelled):
        electrostatics._run(command, cwd=tmp_path, timeout=5, cancelled=cancelled)


class Download(io.BytesIO):
    def geturl(self):
        return "https://example.org/tool.zip"


def archive_fixture(monkeypatch, filename="bin/apbs"):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as zipped:
        zipped.writestr(filename, b"binary fixture")
    raw = buffer.getvalue()
    monkeypatch.setattr(electrostatics.urllib.request, "urlopen", lambda *_args, **_kwargs: Download(raw))
    return {"name": "apbs", "version": "test", "url": "https://example.org/tool.zip",
            "executable": "bin/apbs", "sha256": hashlib.sha256(raw).hexdigest()}


def test_pinned_download_is_atomic_and_reuses_verified_cache(monkeypatch, tmp_path):
    artifact = archive_fixture(monkeypatch)
    executable = electrostatics._download_artifact(artifact, tmp_path)
    assert executable.read_bytes() == b"binary fixture"
    assert executable.stat().st_mode & 0o111
    monkeypatch.setattr(electrostatics.urllib.request, "urlopen", lambda *_a, **_kw: pytest.fail("Redownloaded cached tool"))
    assert electrostatics._download_artifact(artifact, tmp_path) == executable


def test_checksum_mismatch_and_unsafe_archives_do_not_install(monkeypatch, tmp_path):
    artifact = archive_fixture(monkeypatch)
    artifact["sha256"] = "a" * 64
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match="Checksum mismatch"):
        electrostatics._download_artifact(artifact, tmp_path)
    assert not (tmp_path / "apbs-test").exists()
    artifact = archive_fixture(monkeypatch, "../escaped")
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match="unsafe path"):
        electrostatics._download_artifact(artifact, tmp_path)
    assert not (tmp_path / "escaped").exists()


def test_failed_download_returns_setup_diagnostic(monkeypatch, tmp_path):
    artifact = archive_fixture(monkeypatch)

    def fail(*_args, **_kwargs):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr(electrostatics.urllib.request, "urlopen", fail)
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match="Could not download.*offline"):
        electrostatics._download_artifact(artifact, tmp_path)


def test_tool_discovery_supports_overrides_and_reports_unsupported_platform(monkeypatch, tmp_path):
    monkeypatch.delenv("PDB2PQR_BIN", raising=False)
    monkeypatch.delenv("APBS_BIN", raising=False)
    monkeypatch.setattr(electrostatics.shutil, "which", lambda _name: None)
    monkeypatch.setattr(electrostatics, "_platform_key", lambda: "unsupported-arm64")
    with pytest.raises(electrostatics.ElectrostaticsUnavailable, match="Linux x86_64"):
        electrostatics.resolve_tools(tmp_path)
    for name in ["APBS_BIN", "PDB2PQR_BIN"]:
        tool = tmp_path / name
        tool.write_text("#!/bin/sh\n")
        tool.chmod(0o755)
        monkeypatch.setenv(name, str(tool))
    pdb2pqr, apbs = electrostatics.resolve_tools(tmp_path)
    assert pdb2pqr == [str(tmp_path / "PDB2PQR_BIN")]
    assert apbs == [str(tmp_path / "APBS_BIN")]
