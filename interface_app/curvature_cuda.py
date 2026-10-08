"""Optional double-precision CUDA implementation of the Zernike vertex fit.

PyTorch is imported only when a curvature worker requests acceleration. The
surface geometry and sampling stencil still come from zernike_convexity; this
backend moves trilinear interpolation, stencil evaluation, and the principal
curvature calculation to CUDA without changing their definitions.
"""

from __future__ import annotations

import numpy as np


class CudaUnavailable(RuntimeError):
    """CUDA could not be initialized; the caller can use the CPU reference."""


class CudaCurvature:
    def __init__(self, zernike, *, batch_size: int = 2048):
        try:
            import torch
        except ImportError as exc:
            raise CudaUnavailable("Install CUDA-enabled PyTorch to accelerate curvature") from exc
        if not torch.cuda.is_available():
            raise CudaUnavailable("PyTorch cannot access a CUDA device")
        if batch_size < 1:
            raise ValueError("CUDA curvature batch size must be positive")

        self.torch = torch
        self.device = torch.device("cuda", torch.cuda.current_device())
        self.device_name = torch.cuda.get_device_name(self.device)
        self.batch_size = batch_size
        self.grid_A = zernike.GRID_A
        self.derivative_step = zernike.DERIVATIVE_STEP
        self.min_gradient = zernike.MIN_GRADIENT
        # Keep the CPU reference's precision. Float32 loses precision when
        # finite differences subtract nearly equal stencil values.
        self.samples = torch.as_tensor(zernike.SAMPLES, dtype=torch.float64, device=self.device)
        self.weights = torch.as_tensor(zernike.ZERNIKE_WEIGHTS, dtype=torch.float64, device=self.device)
        self.axes = torch.eye(3, dtype=torch.float64, device=self.device)

    def _principal_curvatures(self, value, radius):
        """Equivalent to curvature_from_samples after evaluating the stencil."""
        torch = self.torch
        e = self.derivative_step
        grad = torch.stack([
            (value[1 + 2 * i] - value[2 + 2 * i]) / (2 * e * radius)
            for i in range(3)
        ], dim=1)
        hessian = torch.zeros((value.shape[1], 3, 3), dtype=torch.float64, device=self.device)
        for i in range(3):
            hessian[:, i, i] = (
                value[1 + 2 * i] - 2 * value[0] + value[2 + 2 * i]
            ) / (e * radius) ** 2
        offset = 7
        for i in range(3):
            for j in range(i + 1, 3):
                hessian[:, i, j] = hessian[:, j, i] = (
                    value[offset] - value[offset + 1] - value[offset + 2] + value[offset + 3]
                ) / (4 * (e * radius) ** 2)
                offset += 4

        norm = torch.linalg.vector_norm(grad, dim=1)
        valid = norm >= self.min_gradient
        # Avoid NaNs in the eigensolver for invalid, zero-gradient samples.
        # The final mask retains exactly the reference's unavailable values.
        safe_norm = norm.clamp_min(1e-12)
        normal = grad / safe_norm[:, None]
        axis = self.axes[torch.argmin(normal.abs(), dim=1)]
        tangent_1 = torch.linalg.cross(normal, axis)
        tangent_1 /= torch.linalg.vector_norm(tangent_1, dim=1).clamp_min(1e-12)[:, None]
        tangent_2 = torch.linalg.cross(normal, tangent_1)
        tangents = torch.stack((tangent_1, tangent_2), dim=1)
        shape = torch.einsum("npi,nij,nqj->npq", tangents, hessian, tangents)
        shape /= safe_norm[:, None, None]
        principal = torch.linalg.eigvalsh(shape)
        return torch.where(valid[:, None], principal, torch.nan)

    def curvature_from_samples(self, field_samples, radius):
        """Evaluate an existing sampling stencil, useful for reference checks."""
        torch = self.torch
        with torch.inference_mode():
            values = torch.as_tensor(field_samples, dtype=torch.float64, device=self.device)
            return self._principal_curvatures(self.weights @ values, radius).cpu().numpy()

    def calculate(self, vertices, signed_distance, origin, radii) -> dict[str, np.ndarray]:
        """Calculate every vertex at all radii, uploading each grid only once.

        Coordinate buffers are bounded by batch_size rather than the full mesh
        size. Only the final two curvatures per vertex are copied back to CPU.
        """
        torch = self.torch
        if signed_distance.ndim != 3 or min(signed_distance.shape) < 2:
            raise ValueError("CUDA curvature requires a three-dimensional grid with at least two nodes per axis")
        with torch.inference_mode():
            field = torch.as_tensor(signed_distance, dtype=torch.float64, device=self.device)[None, None]
            points = torch.as_tensor(vertices, dtype=torch.float64, device=self.device)
            grid_origin = torch.as_tensor(origin, dtype=torch.float64, device=self.device)
            dimensions = torch.as_tensor(signed_distance.shape, dtype=torch.float64, device=self.device) - 1
            result = {}
            for radius in radii:
                principal = torch.empty((len(vertices), 2), dtype=torch.float64, device=self.device)
                for start in range(0, len(vertices), self.batch_size):
                    part = points[start:start + self.batch_size]
                    coordinates = (part[:, None, :] + radius * self.samples[None, :, :] - grid_origin) / self.grid_A
                    # grid_sample expects coordinates in width/height/depth
                    # order, the reverse of SciPy's grid array axes. Border
                    # padding reproduces map_coordinates(mode="nearest").
                    normalized = (2 * coordinates / dimensions - 1).flip(-1)
                    sampled = torch.nn.functional.grid_sample(
                        field, normalized[None, None], mode="bilinear",
                        padding_mode="border", align_corners=True,
                    ).reshape(len(part), len(self.samples)).T
                    principal[start:start + len(part)] = self._principal_curvatures(
                        self.weights @ sampled, radius,
                    )
                values = principal.cpu().numpy()
                if np.isinf(values).any():
                    raise RuntimeError("CUDA curvature returned non-finite principal curvatures")
                result[str(int(radius))] = values
            return result
