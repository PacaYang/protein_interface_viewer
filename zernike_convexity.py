#!/usr/bin/env python3
"""Local 3D Zernike curvature of ternary cytokine binding surfaces.

The molecular surface is the boundary of the union of protein-atom spheres
expanded by a 1.4 Å water probe (the solvent-accessible surface). For each
surface point buried by another chain, a positive-outside signed distance
field is fitted in a 6 Å ball with real 3D Zernike polynomials. Derivatives of
that fitted field give the two principal curvatures and signed mean curvature.
This retains the sign of curvature; rotationally invariant descriptor norms
alone would discard it.
"""

from __future__ import annotations

import warnings

import numpy as np
from Bio.PDB.SASA import ATOMIC_RADII
from scipy.ndimage import distance_transform_edt, map_coordinates
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.special import eval_jacobi, sph_harm
from skimage.measure import marching_cubes


PROBE_A = 1.4
GRID_A = 0.6
PATCH_RADIUS_A = 6.0
RADII_A = (4.0, 6.0, 8.0)
MAX_CONTACT_PADDING_A = 8.0
ZERNIKE_ORDER = 6
SAMPLE_AXIS = 13
DERIVATIVE_STEP = 0.02  # in unit-ball coordinates
CURVATURE_THRESHOLD = 0.02  # Å^-1; suppresses near-flat numerical noise
MIN_GRADIENT = 0.2
CHAIN_NAMES = {"A": "IL-4", "B": "IL-4Rα", "C": "IL-13Rα1"}
PAIRS = ("AB", "AC", "BC")
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def _atoms(chain):
    coords, radii, labels = [], [], []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        for atom in residue:
            if atom.element == "H" or atom.get_altloc() not in (" ", "A"):
                continue
            if atom.element not in ATOMIC_RADII:
                raise ValueError(f"No radius for element {atom.element!r}")
            coords.append(atom.coord)
            radii.append(ATOMIC_RADII[atom.element] + PROBE_A)
            labels.append(f"{residue.resname}:{residue.id[1]}{residue.id[2].strip()}")
    if not coords:
        raise ValueError(f"No protein heavy atoms in chain {chain.id}")
    return np.asarray(coords, dtype=np.float64), np.asarray(radii), labels


def _residue_label(atom_label):
    name, number = atom_label.split(":", 1)
    return AA3TO1[name] + number


def _real_zernike(points, modes):
    """Orthonormal real 3D Zernike basis on the unit ball."""
    r = np.linalg.norm(points, axis=1)
    theta = np.arccos(np.divide(points[:, 2], r, out=np.zeros_like(r), where=r > 0).clip(-1, 1))
    azimuth = np.arctan2(points[:, 1], points[:, 0])
    columns = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning, module="scipy")
        for n, ell, m in modes:
            radial = (np.sqrt(2 * n + 3) * r ** ell
                      * eval_jacobi((n - ell) // 2, 0, ell + 0.5, 2 * r * r - 1))
            if m == 0:
                harmonic = sph_harm(0, ell, azimuth, theta).real
            elif m > 0:
                harmonic = np.sqrt(2) * (-1) ** m * sph_harm(m, ell, azimuth, theta).real
            else:
                harmonic = np.sqrt(2) * (-1) ** (-m) * sph_harm(-m, ell, azimuth, theta).imag
            columns.append(radial * harmonic)
    return np.column_stack(columns)


def _zernike_stencil():
    axis = np.linspace(-1, 1, SAMPLE_AXIS)
    samples = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    samples = samples[np.sum(samples * samples, axis=1) <= 1 + 1e-10]
    modes = [(n, ell, m) for n in range(ZERNIKE_ORDER + 1)
             for ell in range(n % 2, n + 1, 2) for m in range(-ell, ell + 1)]

    e = DERIVATIVE_STEP
    evaluation = [np.zeros(3)]
    for i in range(3):
        v = np.zeros(3)
        v[i] = e
        evaluation.extend((v.copy(), -v.copy()))
    for i in range(3):
        for j in range(i + 1, 3):
            for a, b in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                v = np.zeros(3)
                v[i], v[j] = a * e, b * e
                evaluation.append(v)

    basis = _real_zernike(samples, modes)
    if np.linalg.matrix_rank(basis) != len(modes):
        raise ValueError("3D Zernike sampling stencil is rank deficient")
    # W evaluates the least-squares Zernike fit at the finite-difference stencil.
    weights = _real_zernike(np.asarray(evaluation), modes) @ np.linalg.pinv(basis)
    return samples, weights


SAMPLES, ZERNIKE_WEIGHTS = _zernike_stencil()


def curvature_from_samples(field_samples, radius=PATCH_RADIUS_A):
    """Return principal curvatures for fields sampled at SAMPLES around a point.

    `field_samples` has shape (925, n_points). Coordinates are physical Å. For
    a sphere, positive-outside signed distance yields k1=k2=+1/radius.
    """
    value = ZERNIKE_WEIGHTS @ field_samples
    e = DERIVATIVE_STEP
    grad = np.stack([(value[1 + 2 * i] - value[2 + 2 * i]) / (2 * e * radius)
                     for i in range(3)], axis=1)
    hessian = np.zeros((value.shape[1], 3, 3))
    for i in range(3):
        hessian[:, i, i] = (value[1 + 2 * i] - 2 * value[0] + value[2 + 2 * i]) / (e * radius) ** 2
    offset = 7
    for i in range(3):
        for j in range(i + 1, 3):
            hessian[:, i, j] = hessian[:, j, i] = (
                value[offset] - value[offset + 1] - value[offset + 2] + value[offset + 3]
            ) / (4 * (e * radius) ** 2)
            offset += 4

    norm = np.linalg.norm(grad, axis=1)
    valid = norm >= MIN_GRADIENT
    result = np.full((len(norm), 2), np.nan)
    if np.any(valid):
        normal = grad[valid] / norm[valid, None]
        axis = np.eye(3)[np.argmin(np.abs(normal), axis=1)]
        tangent_1 = np.cross(normal, axis)
        tangent_1 /= np.linalg.norm(tangent_1, axis=1)[:, None]
        tangent_2 = np.cross(normal, tangent_1)
        tangents = np.stack((tangent_1, tangent_2), axis=1)
        shape = np.einsum("npi,nij,nqj->npq", tangents, hessian[valid], tangents)
        shape /= norm[valid, None, None]
        result[valid] = np.linalg.eigvalsh(shape)
    return result


def _surface(coords, radii):
    # Retain the original 6 Å grid alignment for comparable interface results.
    # The 8 Å stencil still fits between the molecular surface and grid edge.
    pad = float(radii.max() + PATCH_RADIUS_A + 2)
    origin = coords.min(axis=0) - pad
    high = coords.max(axis=0) + pad
    shape = np.ceil((high - origin) / GRID_A).astype(int) + 1
    inside = np.zeros(tuple(shape), dtype=bool)
    for center, radius in zip(coords, radii):
        lo = np.maximum(0, np.floor((center - radius - origin) / GRID_A).astype(int))
        hi = np.minimum(shape, np.ceil((center + radius - origin) / GRID_A).astype(int) + 1)
        x = origin[0] + np.arange(lo[0], hi[0]) * GRID_A - center[0]
        y = origin[1] + np.arange(lo[1], hi[1]) * GRID_A - center[1]
        z = origin[2] + np.arange(lo[2], hi[2]) * GRID_A - center[2]
        inside[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] |= (
            x[:, None, None] ** 2 + y[None, :, None] ** 2 + z[None, None, :] ** 2 <= radius ** 2
        )
    signed_distance = (distance_transform_edt(~inside) - distance_transform_edt(inside)) * GRID_A
    del inside
    vertices, faces, _, _ = marching_cubes(signed_distance, level=0, spacing=(GRID_A,) * 3)
    vertices = vertices.astype(np.float64) + origin
    faces = faces.astype(np.int32)
    return origin, signed_distance, vertices, faces


def _buried_by_partner(points, partner_coords, partner_radii):
    tree = cKDTree(partner_coords)
    mask = np.zeros(len(points), dtype=bool)
    search_radius = float(partner_radii.max())
    for start in range(0, len(points), 2048):
        part = points[start:start + 2048]
        candidates = tree.query_ball_point(part, r=search_radius)
        for offset, atoms in enumerate(candidates):
            if atoms:
                delta = partner_coords[atoms] - part[offset]
                mask[start + offset] = np.any(np.sum(delta * delta, axis=1) < partner_radii[atoms] ** 2)
    return mask


def _residue_for_vertices(vertices, coords, radii, labels):
    tree = cKDTree(coords)
    res = []
    for start in range(0, len(vertices), 2048):
        part = vertices[start:start + 2048]
        candidates = tree.query_ball_point(part, r=float(radii.max() + 1))
        for point, atoms in zip(part, candidates):
            if not atoms:
                _, nearest = tree.query(point)
                atoms = [int(nearest)]
            delta = coords[atoms] - point
            at = atoms[int(np.argmin(np.linalg.norm(delta, axis=1) - radii[atoms]))]
            res.append(_residue_label(labels[at]))
    return np.asarray(res)


def _vertex_curvatures(vertices, selected, signed_distance, origin, radius=PATCH_RADIUS_A):
    principal = np.full((len(vertices), 2), np.nan)
    for start in range(0, len(selected), 128):
        ids = selected[start:start + 128]
        sample_xyz = vertices[ids, None, :] + radius * SAMPLES[None, :, :]
        ijk = ((sample_xyz - origin) / GRID_A).reshape(-1, 3).T
        values = map_coordinates(signed_distance, ijk, order=1, mode="nearest")
        values = values.reshape(len(ids), len(SAMPLES)).T
        principal[ids] = curvature_from_samples(values, radius)
    return principal


def _mesh_edge_graph(vertices, faces):
    """Weighted surface edges; paths cannot jump between nearby mesh folds."""
    edges = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    lengths = np.linalg.norm(vertices[edges[:, 0]] - vertices[edges[:, 1]], axis=1)
    n = len(vertices)
    graph = coo_matrix((lengths, (edges[:, 0], edges[:, 1])), shape=(n, n)).tocsr()
    return graph.maximum(graph.T)


def _contact_distances(graph, faces, mask, n_vertices):
    """Deciångström distance from contact-face vertices; 255 means beyond 8 Å."""
    seeds = np.unique(faces[mask])
    if not len(seeds):
        return np.full(n_vertices, 255, dtype=np.uint8)
    distance = dijkstra(graph, directed=False, indices=seeds, min_only=True,
                        limit=MAX_CONTACT_PADDING_A)
    return np.where(np.isfinite(distance), np.rint(distance * 10), 255).astype(np.uint8)


def _stats(weights, principal, residues, vertices):
    good = (weights > 0) & np.all(np.isfinite(principal), axis=1)
    area = float(np.sum(weights))
    valid_area = float(np.sum(weights[good]))
    mean_h = float(np.average(np.mean(principal[good], axis=1), weights=weights[good])) if valid_area else None
    convex = np.all(principal[good] > CURVATURE_THRESHOLD, axis=1)
    concave = np.all(principal[good] < -CURVATURE_THRESHOLD, axis=1)
    summary = {
        "area_A2": round(area, 2),
        "valid_area_A2": round(valid_area, 2),
        "mean_H_Ainv": round(mean_h, 5) if mean_h is not None else None,
        "convex_fraction": round(float(np.sum(weights[good][convex]) / valid_area), 4) if valid_area else None,
        "concave_fraction": round(float(np.sum(weights[good][concave]) / valid_area), 4) if valid_area else None,
    }
    summary["mixed_flat_fraction"] = (round(1 - summary["convex_fraction"] - summary["concave_fraction"], 4)
                                      if valid_area else None)
    per_residue = []
    for label in sorted(set(residues[weights > 0]), key=lambda x: (int(''.join(c for c in x if c in '-0123456789')), x)):
        own = (residues == label) & (weights > 0)
        item = _stats_summary_only(weights[own], principal[own])
        w = weights[own]
        item.update({"label": label, "centroid": np.average(vertices[own], weights=w, axis=0).round(2).tolist()})
        per_residue.append(item)
    summary["residues"] = per_residue
    return summary


def _stats_summary_only(weights, principal):
    good = (weights > 0) & np.all(np.isfinite(principal), axis=1)
    area = float(np.sum(weights))
    valid_area = float(np.sum(weights[good]))
    if not valid_area:
        return {"area_A2": round(area, 2), "valid_area_A2": 0.0,
                "mean_H_Ainv": None, "convex_fraction": None,
                "concave_fraction": None, "mixed_flat_fraction": None}
    w = weights[good]
    k = principal[good]
    convex = float(np.sum(w[np.all(k > CURVATURE_THRESHOLD, axis=1)]) / valid_area)
    concave = float(np.sum(w[np.all(k < -CURVATURE_THRESHOLD, axis=1)]) / valid_area)
    return {"area_A2": round(area, 2), "valid_area_A2": round(valid_area, 2),
            "mean_H_Ainv": round(float(np.average(np.mean(k, axis=1), weights=w)), 5),
            "convex_fraction": round(convex, 4), "concave_fraction": round(concave, 4),
            "mixed_flat_fraction": round(1 - convex - concave, 4)}


def _patch_mesh(vertices, faces, mask, principal, residues):
    chosen_faces = faces[mask]
    ids, inverse = np.unique(chosen_faces, return_inverse=True)
    local_faces = inverse.reshape(-1, 3)
    local_h = np.mean(principal[ids], axis=1)
    valid = np.isfinite(local_h)
    local_h = np.where(valid, local_h, 0)
    return {
        "x": np.round(vertices[ids, 0], 2).tolist(),
        "y": np.round(vertices[ids, 1], 2).tolist(),
        "z": np.round(vertices[ids, 2], 2).tolist(),
        "i": local_faces[:, 0].tolist(),
        "j": local_faces[:, 1].tolist(),
        "k": local_faces[:, 2].tolist(),
        "h": np.round(local_h, 4).tolist(),
        "valid": valid.tolist(),
        "residue": residues[ids].tolist(),
    }


def compute(model, analysis, *, chain_names=CHAIN_NAMES, source="3BPN.pdb",
            source_resolution_A=3.02):
    """Return scale-keyed whole surfaces while retaining the 6 Å interface report."""
    atoms = {chain: _atoms(model[chain]) for chain in chain_names}
    report = {
        "source": source, "source_resolution_A": source_resolution_A,
        "parameters": {"probe_A": PROBE_A, "surface_grid_A": GRID_A,
                       "zernike_radius_A": PATCH_RADIUS_A, "zernike_radii_A": RADII_A,
                       "max_contact_padding_A": MAX_CONTACT_PADDING_A,
                       "zernike_order": ZERNIKE_ORDER,
                       "sample_nodes": len(SAMPLES), "curvature_threshold_Ainv": CURVATURE_THRESHOLD,
                       "curvature_sign": "positive = convex/outward"},
        "proteins": {}, "interfaces": {},
        "scales": {str(int(radius)): {"proteins": {}, "interfaces": {}} for radius in RADII_A},
    }
    meshes = {}
    all_h = []

    for chain in chain_names:
        coords, radii, atom_labels = atoms[chain]
        print(f"  surface {chain} ({chain_names[chain]}): {len(coords)} heavy atoms", flush=True)
        origin, field, vertices, faces = _surface(coords, radii)
        triangle_points = vertices[faces]
        triangle_area = np.linalg.norm(np.cross(triangle_points[:, 1] - triangle_points[:, 0],
                                                triangle_points[:, 2] - triangle_points[:, 0]), axis=1) / 2
        centers = triangle_points.mean(axis=1)
        del triangle_points
        masks = {}
        for partner in chain_names:
            if partner == chain:
                continue
            pc, pr, _ = atoms[partner]
            masks[partner] = _buried_by_partner(centers, pc, pr)
        union = np.logical_or.reduce(list(masks.values()))
        edge_graph = _mesh_edge_graph(vertices, faces)
        contact_distance = {
            partner: _contact_distances(edge_graph, faces, mask, len(vertices)).tolist()
            for partner, mask in masks.items()
        }
        del edge_graph
        print(f"    {len(vertices)} vertices, {len(faces)} faces; fitting all vertices at {RADII_A} Å", flush=True)
        principals = {}
        selected = np.arange(len(vertices))
        for radius in RADII_A:
            principals[str(int(radius))] = _vertex_curvatures(vertices, selected, field, origin, radius)
            print(f"    completed {radius:g} Å", flush=True)
        del field
        residues = _residue_for_vertices(vertices, coords, radii, atom_labels)
        whole_weights = np.zeros(len(vertices))
        np.add.at(whole_weights, faces.ravel(), np.repeat(triangle_area / 3, 3))
        union_weights = np.zeros(len(vertices))
        np.add.at(union_weights, faces[union].ravel(), np.repeat(triangle_area[union] / 3, 3))
        vertex_h = {}
        for scale, principal in principals.items():
            whole = _stats(whole_weights, principal, residues, vertices)
            whole["name"] = chain_names[chain]
            report["scales"][scale]["proteins"][chain] = whole
            patch_union = _stats(union_weights, principal, residues, vertices)
            patch_union["name"] = chain_names[chain]
            if scale == "6":
                report["proteins"][chain] = patch_union
            h = np.mean(principal, axis=1)
            vertex_h[scale] = np.round(np.where(np.isfinite(h), h, 0), 4).tolist()
            all_h.extend(h[np.isfinite(h)].tolist())

        meshes[chain] = {
            "position": np.round(vertices, 2).ravel().tolist(),
            "index": faces.ravel().tolist(),
            "h": vertex_h,
            "residue": residues.tolist(),
            "contact_distance_dA": contact_distance,
        }

        for partner, mask in masks.items():
            weights = np.zeros(len(vertices))
            np.add.at(weights, faces[mask].ravel(), np.repeat(triangle_area[mask] / 3, 3))
            pair = "".join(sorted((chain, partner)))
            for scale, principal in principals.items():
                side = _stats(weights, principal, residues, vertices)
                side["name"] = chain_names[chain]
                side["reference_dSASA_A2"] = round(sum(r["dSASA"] for r in analysis[pair]["residues"]
                                                       if r["chain"] == chain), 1)
                if abs(side["area_A2"] / side["reference_dSASA_A2"] - 1) > 0.25:
                    raise ValueError(f"{pair} chain {chain} buried mesh area disagrees with ΔSASA")
                report["scales"][scale]["interfaces"].setdefault(pair, {})[chain] = side
                if scale == "6":
                    report["interfaces"].setdefault(pair, {})[chain] = side
                    print(f"    {pair} {chain}: {side['area_A2']:.0f} Å², H={side['mean_H_Ainv']}", flush=True)

    if not all_h:
        raise ValueError("No valid Zernike curvature samples on any interface")
    limit = float(np.quantile(np.abs(all_h), 0.98))
    report["color_limit_Ainv"] = round(limit, 4)
    for pair in PAIRS:
        a, b = pair
        ac = np.asarray([r["centroid"] for r in report["interfaces"][pair][a]["residues"]])
        bc = np.asarray([r["centroid"] for r in report["interfaces"][pair][b]["residues"]])
        direction = bc.mean(axis=0) - ac.mean(axis=0)
        direction /= np.linalg.norm(direction)
        report.setdefault("pair_directions", {})[pair] = np.round(direction, 5).tolist()
    return {"report": report, "meshes": meshes}
