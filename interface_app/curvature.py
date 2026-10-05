"""Background-friendly generic Zernike surface data.

The repository's detailed Zernike implementation remains the numerical source
of truth. This adapter keeps its surface and curvature calculations, but
serializes the complete mesh/report shape needed by the browser app and works
with arbitrary protein-chain identifiers.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .constants import AA3TO1


def _residue_labels(vertices, coords, radii, atom_labels):
    """Assign each surface vertex to the nearest atom's residue."""

    tree = cKDTree(coords)
    labels = []
    for start in range(0, len(vertices), 2048):
        part = vertices[start:start + 2048]
        candidates = tree.query_ball_point(part, r=float(radii.max() + 1))
        for point, atoms in zip(part, candidates):
            if not atoms:
                _, nearest = tree.query(point)
                atoms = [int(nearest)]
            atoms = np.asarray(atoms, dtype=int)
            delta = coords[atoms] - point
            score = np.linalg.norm(delta, axis=1) - radii[atoms]
            name = atom_labels[int(atoms[int(np.argmin(score))])]
            residue_name, number = name.split(":", 1)
            labels.append(f"{AA3TO1.get(residue_name.upper(), 'X')}{number}")
    return np.asarray(labels)


def _chain_names(analysis: dict) -> dict[str, str]:
    items = analysis.get("chains", [])
    if not items:
        items = analysis.get("metadata", {}).get("chains", [])
    if isinstance(items, dict):
        items = [
            {"id": chain_id, **(value if isinstance(value, dict) else {})}
            for chain_id, value in items.items()
        ]
    names = {}
    for item in items:
        if isinstance(item, str):
            names[item] = f"Chain {item}"
            continue
        chain_id = str(item["id"])
        names[chain_id] = item.get("name", f"Chain {chain_id}")
    return names


def _stats(zernike, weights, principal, residues, vertices):
    return zernike._stats(weights, principal, residues, vertices)


def _pair_direction(report, pair_id, chain_a, chain_b, atom_cache):
    sides = report["scales"]["6"]["interfaces"].get(pair_id, {})

    def centroid(chain_id):
        points = [
            residue["centroid"]
            for residue in sides.get(chain_id, {}).get("residues", [])
            if residue.get("centroid")
        ]
        if points:
            return np.mean(np.asarray(points, dtype=float), axis=0)
        return np.mean(atom_cache[chain_id][0], axis=0)

    direction = centroid(chain_b) - centroid(chain_a)
    length = float(np.linalg.norm(direction))
    if length <= 1e-12:
        direction = np.asarray([0.0, 0.0, 1.0])
    else:
        direction = direction / length
    return np.round(direction, 5).tolist()


def _unavailable(exc: Exception) -> dict:
    return {
        "version": "0.2.0",
        "status": "unavailable",
        "reason": "Install scikit-image to enable Zernike surface curvature.",
        "error": str(exc),
        "parameters": {},
        "chains": {},
        "interfaces": {},
        "report": None,
        "meshes": {},
    }


def compute_curvature(model, analysis: dict) -> dict:
    try:
        import zernike_convexity as zernike
    except Exception as exc:
        return _unavailable(exc)

    try:
        chain_names = _chain_names(analysis)
        chain_ids = list(chain_names)
        if not chain_ids:
            raise ValueError("No protein chains were available for surface calculation")

        radii = tuple(float(radius) for radius in zernike.RADII_A)
        report = {
            "source": analysis.get("source", {}).get("name"),
            "parameters": {
                "probe_A": zernike.PROBE_A,
                "surface_grid_A": zernike.GRID_A,
                "zernike_radius_A": zernike.PATCH_RADIUS_A,
                "zernike_radii_A": list(radii),
                "max_contact_padding_A": zernike.MAX_CONTACT_PADDING_A,
                "zernike_order": zernike.ZERNIKE_ORDER,
                "sample_nodes": len(zernike.SAMPLES),
                "curvature_threshold_Ainv": zernike.CURVATURE_THRESHOLD,
                "curvature_sign": "positive = convex/outward",
            },
            "proteins": {},
            "interfaces": {},
            "scales": {
                str(int(radius)): {"proteins": {}, "interfaces": {}}
                for radius in radii
            },
            "pair_directions": {},
        }
        atom_cache = {}
        surface_cache = {}
        meshes = {}
        all_h = []

        for chain_id in chain_ids:
            atom_cache[chain_id] = zernike._atoms(model[chain_id])

        for chain_id in chain_ids:
            coords, atom_radii, atom_labels = atom_cache[chain_id]
            origin, field, vertices, faces = zernike._surface(coords, atom_radii)
            triangle_points = vertices[faces]
            triangle_area = np.linalg.norm(
                np.cross(
                    triangle_points[:, 1] - triangle_points[:, 0],
                    triangle_points[:, 2] - triangle_points[:, 0],
                ),
                axis=1,
            ) / 2
            centers = triangle_points.mean(axis=1)
            del triangle_points

            masks = {}
            for partner_id in chain_ids:
                if partner_id == chain_id:
                    continue
                partner_coords, partner_radii, _ = atom_cache[partner_id]
                masks[partner_id] = zernike._buried_by_partner(
                    centers, partner_coords, partner_radii
                )
            union = (
                np.logical_or.reduce(list(masks.values()))
                if masks else np.zeros(len(centers), dtype=bool)
            )

            edge_graph = zernike._mesh_edge_graph(vertices, faces)
            contact_distance = {
                partner_id: zernike._contact_distances(
                    edge_graph, faces, mask, len(vertices)
                ).tolist()
                for partner_id, mask in masks.items()
            }
            del edge_graph

            principals = {}
            for radius in radii:
                key = str(int(radius))
                principals[key] = zernike._vertex_curvatures(
                    vertices,
                    np.arange(len(vertices)),
                    field,
                    origin,
                    radius,
                )

            residues = _residue_labels(vertices, coords, atom_radii, atom_labels)
            whole_weights = np.zeros(len(vertices))
            np.add.at(whole_weights, faces.ravel(), np.repeat(triangle_area / 3, 3))
            union_weights = np.zeros(len(vertices))
            if np.any(union):
                np.add.at(union_weights, faces[union].ravel(), np.repeat(triangle_area[union] / 3, 3))

            surface_cache[chain_id] = {
                "coords": coords,
                "vertices": vertices,
                "faces": faces,
                "triangle_area": triangle_area,
                "weights": whole_weights,
                "principals": principals,
                "masks": masks,
                "residues": residues,
            }

            vertex_h = {}
            for scale, principal in principals.items():
                whole = _stats(zernike, whole_weights, principal, residues, vertices)
                whole["name"] = chain_names[chain_id]
                report["scales"][scale]["proteins"][chain_id] = whole

                patch = _stats(zernike, union_weights, principal, residues, vertices)
                patch["name"] = chain_names[chain_id]
                if scale == "6":
                    report["proteins"][chain_id] = patch

                h = np.mean(principal, axis=1)
                vertex_h[scale] = np.round(
                    np.where(np.isfinite(h), h, 0), 4
                ).tolist()
                all_h.extend(h[np.isfinite(h)].tolist())

            meshes[chain_id] = {
                "position": np.round(vertices, 2).ravel().tolist(),
                "index": faces.ravel().tolist(),
                "h": vertex_h,
                "residue": residues.tolist(),
                "contact_distance_dA": contact_distance,
            }
            del field

        if not all_h:
            raise ValueError("No valid Zernike curvature samples on any surface")

        for pair in analysis.get("pairs", []):
            pair_id = pair["id"]
            chain_a, chain_b = pair["chain_a"], pair["chain_b"]
            if chain_a not in surface_cache or chain_b not in surface_cache:
                continue
            for chain_id, partner_id in ((chain_a, chain_b), (chain_b, chain_a)):
                cached = surface_cache[chain_id]
                mask = cached["masks"].get(partner_id)
                interface_weights = np.zeros(len(cached["weights"]))
                if mask is not None and np.any(mask):
                    np.add.at(
                        interface_weights,
                        cached["faces"][mask].ravel(),
                        np.repeat(cached["triangle_area"][mask] / 3, 3),
                    )
                for scale, principal in cached["principals"].items():
                    side = _stats(
                        zernike,
                        interface_weights,
                        principal,
                        cached["residues"],
                        cached["vertices"],
                    )
                    side["name"] = chain_names[chain_id]
                    # Keep the independent interface-area reference used by
                    # contact_map.html alongside the mesh-derived summary.
                    # This is useful in the browser when checking that the
                    # visual contact patch and the SASA calculation describe
                    # the same interface, and preserves the historical
                    # surface-convexity payload shape.
                    side["reference_dSASA_A2"] = round(
                        float(sum(
                            row.get("dSASA_A2", 0.0) or 0.0
                            for row in pair.get("residues", [])
                            if row.get("chain_id") == chain_id
                        )),
                        1,
                    )
                    report["scales"][scale]["interfaces"].setdefault(pair_id, {})[
                        chain_id
                    ] = side
                    if scale == "6":
                        report["interfaces"].setdefault(pair_id, {})[chain_id] = side

        report["color_limit_Ainv"] = round(
            float(np.quantile(np.abs(np.asarray(all_h)), 0.98)), 4
        )
        for pair in analysis.get("pairs", []):
            pair_id = pair["id"]
            report["pair_directions"][pair_id] = _pair_direction(
                report, pair_id, pair["chain_a"], pair["chain_b"], atom_cache
            )

        scale_six = report["scales"]["6"]
        return {
            "version": "0.2.0",
            "status": "complete",
            "parameters": report["parameters"],
            "chains": scale_six["proteins"],
            "interfaces": scale_six["interfaces"],
            "report": report,
            "meshes": meshes,
        }
    except Exception as exc:
        return _unavailable(exc)
