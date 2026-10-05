"""Geometry checks for the signed 3D Zernike curvature calculation."""

import unittest

import numpy as np

from zernike_convexity import (
    RADII_A, SAMPLES, _contact_distances, _mesh_edge_graph, _surface,
    _vertex_curvatures, curvature_from_samples,
)


class ZernikeCurvatureTests(unittest.TestCase):
    def principal(self, signed_distance, radius):
        points = radius * SAMPLES
        return curvature_from_samples(signed_distance(points)[:, None], radius)[0]

    def test_plane(self):
        for radius in RADII_A:
            with self.subTest(radius=radius):
                k = self.principal(lambda p: p[:, 2], radius)
                np.testing.assert_allclose(k, [0, 0], atol=0.002)

    def test_sphere_and_cavity(self):
        radius = 10.0
        sphere = lambda p: np.linalg.norm(p + [0, 0, radius], axis=1) - radius
        for scale in RADII_A:
            with self.subTest(radius=scale):
                np.testing.assert_allclose(self.principal(sphere, scale), [0.1, 0.1], atol=0.02)
                np.testing.assert_allclose(self.principal(lambda p: -sphere(p), scale), [-0.1, -0.1], atol=0.02)

    def test_saddle(self):
        for radius in RADII_A:
            with self.subTest(radius=radius):
                k = self.principal(lambda p: p[:, 2] - 0.06 * (p[:, 0] ** 2 - p[:, 1] ** 2), radius)
                np.testing.assert_allclose(k, [-0.12, 0.12], atol=0.003)

    def test_voxel_surface_sphere(self):
        origin, field, vertices, _ = _surface(np.array([[0., 0., 0.]]), np.array([10.]))
        vertex_id = int(np.argmin(np.linalg.norm(vertices - [0, 0, 10], axis=1)))
        principal = _vertex_curvatures(vertices, np.array([vertex_id]), field, origin, 6.0)[vertex_id]
        self.assertTrue(np.all(principal > 0))
        self.assertAlmostEqual(float(np.mean(principal)), 0.1, delta=0.025)

    def test_contact_padding_follows_connected_mesh(self):
        vertices = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                             [1., 1., 0.], [0., 0., 0.1], [1., 0., 0.1],
                             [0., 1., 0.1]])
        faces = np.array([[0, 1, 2], [1, 3, 2], [4, 5, 6]])
        graph = _mesh_edge_graph(vertices, faces)
        distances = _contact_distances(graph, faces, np.array([True, False, False]), len(vertices))
        np.testing.assert_array_equal(distances[:3], [0, 0, 0])
        self.assertEqual(int(distances[3]), 10)
        np.testing.assert_array_equal(distances[4:], [255, 255, 255])


if __name__ == "__main__":
    unittest.main()
