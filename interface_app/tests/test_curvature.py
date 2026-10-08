import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from interface_app import curvature
from interface_app.curvature import compute_curvature
from interface_app.structure import parse_structure


@pytest.mark.parametrize("mode", ["cpu", "cuda", "unavailable", "failure", "mixed"])
def test_curvature_payload_serializes_scales_meshes_and_interface_metadata(monkeypatch, mode):
    vertices = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    faces = np.asarray([[0, 1, 2], [0, 1, 3]], dtype=np.int32)

    def atoms(_chain):
        return (
            vertices.copy(),
            np.ones(len(vertices), dtype=float),
            ["ALA:1", "GLY:2", "SER:3", "THR:4"],
        )

    def stats(weights, principal, residues, mesh_vertices):
        good = np.isfinite(principal).all(axis=1) & (weights > 0)
        return {
            "area_A2": float(weights.sum()),
            "valid_area_A2": float(weights[good].sum()),
            "mean_H_Ainv": 0.1 if good.any() else None,
            "convex_fraction": 1.0 if good.any() else None,
            "concave_fraction": 0.0 if good.any() else None,
            "mixed_flat_fraction": 0.0 if good.any() else None,
            "residues": [
                {"label": label, "area_A2": 1.0, "centroid": mesh_vertices[index].tolist()}
                for index, label in enumerate(sorted(set(residues)))
                if weights[index] > 0
            ],
        }

    cpu_calls = []

    def cpu_curvatures(mesh_vertices, selected, field, origin, radius):
        cpu_calls.append(radius)
        return np.full((len(mesh_vertices), 2), 0.1)

    fake = SimpleNamespace(
        RADII_A=(4.0, 6.0, 8.0),
        PROBE_A=1.4,
        GRID_A=0.6,
        PATCH_RADIUS_A=6.0,
        MAX_CONTACT_PADDING_A=8.0,
        ZERNIKE_ORDER=6,
        SAMPLES=np.zeros((3, 3)),
        CURVATURE_THRESHOLD=0.02,
        _atoms=atoms,
        _surface=lambda coords, radii: (
            np.zeros(3),
            np.zeros((2, 2, 2)),
            vertices.copy(),
            faces.copy(),
        ),
        _buried_by_partner=lambda centers, coords, radii: np.asarray([True, False]),
        _mesh_edge_graph=lambda mesh_vertices, mesh_faces: None,
        _contact_distances=lambda graph, mesh_faces, mask, count: np.zeros(count, dtype=np.uint8),
        _vertex_curvatures=cpu_curvatures,
        _stats=stats,
    )
    monkeypatch.setitem(sys.modules, "zernike_convexity", fake)

    class FakeCuda:
        device_name = "Test CUDA device"

        def __init__(self, _zernike):
            if mode == "cpu":
                pytest.fail("CPU mode initialized CUDA")
            if mode == "unavailable":
                raise RuntimeError("CUDA device unavailable")
            self.calls = 0

        def calculate(self, mesh_vertices, field, origin, radii):
            self.calls += 1
            if mode == "failure" or (mode == "mixed" and self.calls == 2):
                raise RuntimeError("CUDA allocation failed")
            return {str(int(radius)): np.full((len(mesh_vertices), 2), 0.1) for radius in radii}

    monkeypatch.setattr(curvature, "CudaCurvature", FakeCuda)

    loaded = parse_structure(Path(__file__).resolve().parents[2] / "examples" / "3BPN.pdb")
    analysis = {
        "source": {"name": "synthetic.pdb"},
        "chains": [{"id": "A", "name": "Chain A"}, {"id": "B", "name": "Chain B"}],
        "pairs": [{
            "id": "QQBC",
            "chain_a": "A",
            "chain_b": "B",
            "residues": [
                {"chain_id": "A", "dSASA_A2": 12.0},
                {"chain_id": "B", "dSASA_A2": 10.0},
            ],
        }],
    }
    result = compute_curvature(loaded.model, analysis, backend="cpu" if mode == "cpu" else "cuda")

    assert result["status"] == "complete"
    assert result["version"] == "0.2.0"
    assert set(result["report"]["scales"]) == {"4", "6", "8"}
    assert set(result["meshes"]) == {"A", "B"}
    assert len(result["meshes"]["A"]["position"]) == 12
    assert len(result["meshes"]["A"]["h"]["6"]) == 4
    assert result["report"]["scales"]["6"]["interfaces"]["QQBC"]["A"]["reference_dSASA_A2"] == 12.0
    assert len(result["report"]["pair_directions"]["QQBC"]) == 3
    expected = {"cpu": "cpu", "cuda": "cuda", "unavailable": "cpu", "failure": "cpu", "mixed": "mixed"}
    acceleration = result["parameters"]["acceleration"]
    assert acceleration["backend"] == expected[mode]
    assert acceleration["precision"] == "float64"
    assert acceleration["chains"] == {
        "A": "cuda" if mode in {"cuda", "mixed"} else "cpu",
        "B": "cuda" if mode == "cuda" else "cpu",
    }
    assert len(cpu_calls) == {"cpu": 6, "cuda": 0, "unavailable": 6, "failure": 6, "mixed": 3}[mode]
    assert bool(acceleration["fallback_reason"]) == (mode in {"unavailable", "failure", "mixed"})
