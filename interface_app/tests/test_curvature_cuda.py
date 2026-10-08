"""Backend selection/fallback and optional numerical checks on a real GPU.

Run the CUDA checks with INTERFACE_APP_TEST_CUDA=1 on a CUDA-capable host.
They exercise the production backend and compare against the SciPy reference.
"""

import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from interface_app import curvature
from interface_app.curvature_cuda import CudaCurvature, CudaUnavailable


def test_cpu_selection_does_not_load_the_optional_backend(monkeypatch):
    monkeypatch.setenv("INTERFACE_APP_CURVATURE_BACKEND", "cpu")
    monkeypatch.setattr(curvature, "CudaCurvature", lambda _z: pytest.fail("Initialized CUDA in CPU mode"))
    engine, metadata = curvature._acceleration(object(), None)
    assert engine is None
    assert metadata["requested"] == "cpu"
    assert metadata["fallback_reason"] is None


def test_missing_pytorch_falls_back_to_cpu(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)
    engine, metadata = curvature._acceleration(object(), "auto")
    assert engine is None
    assert "PyTorch" in metadata["fallback_reason"]


def test_unavailable_cuda_falls_back_to_cpu(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: False)))
    engine, metadata = curvature._acceleration(object(), "cuda")
    assert engine is None
    assert "CUDA device" in metadata["fallback_reason"]


def test_backend_override_and_validation(monkeypatch):
    monkeypatch.setenv("INTERFACE_APP_CURVATURE_BACKEND", "invalid")
    engine, metadata = curvature._acceleration(object(), "cpu")
    assert engine is None
    assert metadata["requested"] == "cpu"
    with pytest.raises(ValueError, match="auto, cpu, or cuda"):
        curvature._acceleration(object(), None)


@pytest.fixture(scope="module")
def cuda_backend():
    if os.environ.get("INTERFACE_APP_TEST_CUDA") != "1":
        pytest.skip("Set INTERFACE_APP_TEST_CUDA=1 to run numerical checks on a real CUDA GPU")
    import zernike_convexity as zernike

    try:
        return CudaCurvature(zernike, batch_size=17)
    except (CudaUnavailable, ImportError) as exc:
        pytest.fail(f"CUDA numerical checks were requested but CUDA is unavailable: {exc}")


@pytest.mark.parametrize("radius", [4., 6., 8.])
def test_cuda_analytic_curvatures_and_invalid_gradients_match_reference(cuda_backend, radius):
    import zernike_convexity as z

    p = radius * z.SAMPLES
    sphere = np.linalg.norm(p + [0, 0, 10.], axis=1) - 10.
    samples = np.column_stack([
        p[:, 2], sphere, -sphere,
        p[:, 2] - 0.06 * (p[:, 0] ** 2 - p[:, 1] ** 2),
        np.zeros(len(p)), 0.1 * p[:, 2],
    ])
    expected = z.curvature_from_samples(samples, radius)
    actual = cuda_backend.curvature_from_samples(samples, radius)
    np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=2e-11)
    assert np.all(actual[1] > 0), "Sphere must remain convex"
    assert np.all(actual[2] < 0), "Cavity must remain concave"
    assert actual[3, 0] < 0 < actual[3, 1], "Saddle must retain opposite principal signs"
    assert np.isnan(actual[4:]).all(), "Zero and sub-threshold gradients remain unavailable"


def test_cuda_voxel_surface_all_scales_and_batch_tail(cuda_backend):
    import zernike_convexity as z

    origin, field, vertices, _ = z._surface(np.asarray([[2., -7., 4.]]), np.asarray([10.]))
    # A translated surface covers every grid axis and a final incomplete batch.
    vertices = vertices[::137][:43]
    assert len(vertices) % cuda_backend.batch_size != 0
    actual = cuda_backend.calculate(vertices, field, origin, z.RADII_A)
    assert set(actual) == {"4", "6", "8"}
    for radius in z.RADII_A:
        expected = z._vertex_curvatures(vertices, np.arange(len(vertices)), field, origin, radius)
        np.testing.assert_array_equal(np.isfinite(actual[str(int(radius))]), np.isfinite(expected))
        np.testing.assert_allclose(actual[str(int(radius))], expected, rtol=1e-9, atol=2e-11)


def test_cuda_grid_axes_and_border_padding_match_scipy(cuda_backend):
    import zernike_convexity as z

    origin = np.asarray([-12.3, 10.7, -6.4])
    shape = (31, 35, 41)
    xyz = np.indices(shape).astype(float) * z.GRID_A + origin[:, None, None, None]
    x, y, height = xyz
    field = height + 0.01 * x ** 2 - 0.02 * y ** 2 + 0.03 * x * y
    indices = np.asarray([[15, 17, 20], [1, 2, 3], [29, 33, 39], [-5, 10, 20], [36, 10, 20]])
    vertices = origin + indices * z.GRID_A
    actual = cuda_backend.calculate(vertices, field, origin, z.RADII_A)
    for radius in z.RADII_A:
        expected = z._vertex_curvatures(vertices, np.arange(len(vertices)), field, origin, radius)
        np.testing.assert_array_equal(np.isfinite(actual[str(int(radius))]), np.isfinite(expected))
        np.testing.assert_allclose(actual[str(int(radius))], expected, rtol=1e-9, atol=2e-11)
