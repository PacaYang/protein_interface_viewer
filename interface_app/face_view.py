"""Binding-face frames and unrestricted residue matches in complex coordinates."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .analysis import CONTACT_CUTOFF_A, _atom_records


def _unit(vector):
    vector = np.asarray(vector, dtype=np.float64)
    length = float(np.linalg.norm(vector))
    return vector / length if length > 1e-7 else None


def _canonical_axis(vector):
    """Fix an SVD vector's arbitrary sign using its largest component."""
    return -vector if vector[int(np.argmax(np.abs(vector)))] < 0 else vector


def _in_plane(reference, normal):
    axis = _unit(reference - np.dot(reference, normal) * normal)
    if axis is None:
        reference = np.eye(3)[int(np.argmin(np.abs(normal)))]
        axis = _unit(reference - np.dot(reference, normal) * normal)
    return axis


def _plane(points, toward):
    center = points.mean(axis=0)
    _, singular, axes = np.linalg.svd(points - center, full_matrices=len(points) < 3)
    # A line cannot define a plane; near-isotropic clouds have no unique
    # least-variance direction. In either case use the separation direction.
    stable = (len(points) >= 3 and singular[1] > max(1e-7, singular[0] * 1e-6)
              and singular[2] < singular[1] * 0.95)
    normal = axes[-1] if stable else toward
    method = "residue_centroid_plane" if stable else "centroid_direction"
    if normal is not None and toward is not None:
        alignment = float(np.dot(normal, toward))
        if stable and abs(alignment) < 1e-7:
            normal, method = toward, "centroid_direction"
        elif alignment < 0:
            normal = -normal
    else:
        # A plane's sign has no physical meaning without a direction to the
        # partner. Do not invent a binding face for coincident centroids.
        normal = None
    principal = axes[0] if stable and singular[0] > singular[1] * 1.01 else None
    return center, normal, principal, method


def nearest_residues(source: list[dict], partner: list[dict]) -> dict:
    """Map each residue to its nearest partner by minimum heavy-atom distance.

    KD-tree queries avoid a dense atom-distance matrix. Expand only the nearest
    distance shell to resolve ties deterministically in partner sequence order.
    Distances are compared and returned at full precision, without a cutoff.
    """
    coordinates = np.asarray([atom["coord"] for atom in partner])
    tree = cKDTree(coordinates)
    source_coords = np.asarray([atom["coord"] for atom in source])
    _, indices = tree.query(source_coords)
    order = {key: index for index, key in enumerate(dict.fromkeys(
        atom["residue_key"] for atom in partner
    ))}
    groups = {}
    for index, atom in enumerate(source):
        groups.setdefault(atom["residue_key"], []).append(index)
    result = {}
    for key, atom_indices in groups.items():
        best = min(float(np.linalg.norm(source_coords[i] - coordinates[indices[i]]))
                   for i in atom_indices)
        radius = np.nextafter(best, np.inf) + np.spacing(best) * 8
        best_key = None
        for i in atom_indices:
            for j in tree.query_ball_point(source_coords[i], radius):
                distance = float(np.linalg.norm(source_coords[i] - coordinates[j]))
                partner_key = partner[j]["residue_key"]
                if distance < best or (distance == best and (
                        best_key is None or order[partner_key] < order[best_key])):
                    best, best_key = distance, partner_key
        # Include the original nearest candidate even if a tree shell is
        # empty due to a platform's floating point boundary behavior.
        if best_key is None:
            i = min(atom_indices, key=lambda i: (
                float(np.linalg.norm(source_coords[i] - coordinates[indices[i]])),
                order[partner[indices[i]]["residue_key"]],
            ))
            best_key = partner[indices[i]]["residue_key"]
        result[key] = {"partner_key": best_key, "distance_A": best}
    return result


def binding_face_view(model, pair: dict) -> dict:
    """Build two outward frames from the stored <=reported-cutoff residues."""
    chains = [pair["chain_a"], pair["chain_b"]]
    cutoff = pair.get("parameters", {}).get("contact_cutoff_A", CONTACT_CUTOFF_A)
    result = {"pair_id": pair["id"], "chain_a": chains[0], "chain_b": chains[1],
              "binding_cutoff_A": cutoff, "available": False, "reason": None,
              "chains": {}, "nearest": {}}
    atoms = {chain: _atom_records(model, chain) for chain in chains}
    groups = {}
    binding = {}
    for chain in chains:
        groups[chain] = {}
        for atom in atoms[chain]:
            groups[chain].setdefault(atom["residue_key"], []).append(atom["coord"])
        keys = {row["key"] for row in pair.get("residues", []) if row["chain_id"] == chain}
        # Keep sequence order rather than the interface table's dSASA order.
        binding[chain] = {key: np.mean(coords, axis=0) for key, coords in groups[chain].items()
                          if key in keys}
    if any(not binding[chain] for chain in chains):
        result["reason"] = f"No binding residues were reported within {cutoff:g} Å for this interface."
        return result
    points = {chain: np.asarray(list(binding[chain].values())) for chain in chains}
    separation = _unit(points[chains[1]].mean(axis=0) - points[chains[0]].mean(axis=0))
    direction_source = "binding_centroids"
    if separation is None:
        separation = _unit(np.mean([atom["coord"] for atom in atoms[chains[1]]], axis=0)
                           - np.mean([atom["coord"] for atom in atoms[chains[0]]], axis=0))
        direction_source = "chain_centroids"
    planes = [_plane(points[chain], separation if index == 0 else (
        -separation if separation is not None else None)) for index, chain in enumerate(chains)]
    if any(plane[1] is None for plane in planes):
        result["reason"] = "The binding residues do not define a stable direction toward the partner."
        return result
    reference = planes[0][2]
    if reference is None:
        reference = np.eye(3)[int(np.argmin(np.abs(planes[0][1])))]
    reference = _canonical_axis(_in_plane(reference, planes[0][1]))
    for chain, (center, normal, _, method) in zip(chains, planes):
        result["chains"][chain] = {
            "center": center.tolist(), "normal": normal.tolist(),
            "up": _in_plane(reference, normal).tolist(),
            "binding_residue_keys": list(binding[chain]), "method": method,
            "estimated": method != "residue_centroid_plane", "direction_source": direction_source,
        }
    result["nearest"] = {
        chain: nearest_residues(atoms[chain], atoms[partner])
        for chain, partner in (chains, chains[::-1])
    }
    result["available"] = True
    return result
