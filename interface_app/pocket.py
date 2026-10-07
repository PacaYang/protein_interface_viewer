"""Grid-based binding-pocket geometry and free/bound comparison."""

from __future__ import annotations

from collections import deque
from pathlib import Path

import numpy as np
from Bio.PDB.SASA import ATOMIC_RADII
from scipy.ndimage import binary_propagation, distance_transform_edt, find_objects, generate_binary_structure, label
from scipy.spatial import cKDTree

from .analysis import residue_sasa, residues_for_chain
from .constants import atom_element, is_heavy_atom, residue_identity

GRID_A = 0.6
SOLVENT_PROBE_A = 1.4
BULK_PROBE_A = 4.0
MAX_VOXELS = 14_000_000
MAX_COMPONENTS = 5
MAX_POCKET_RESIDUES = 60
NEIGHBORHOOD = generate_binary_structure(3, 1)


def _atom_records(model, chains: list[str]) -> list[dict]:
    result = []
    for chain_id in chains:
        for residue in residues_for_chain(model, chain_id):
            ident = residue_identity(chain_id, residue)
            for atom in residue:
                if not is_heavy_atom(atom):
                    continue
                element = atom_element(atom)
                radius = float(ATOMIC_RADII.get(element, 1.70))
                result.append({
                    "coord": np.asarray(atom.coord, dtype=np.float64),
                    "radius_A": radius,
                    "chain_id": chain_id,
                    "residue": ident,
                })
    if not result:
        raise ValueError("The selected pocket state has no supported heavy atoms")
    return result


def _grid_for(records: list[dict], grid_A: float, padding_A: float):
    coords = np.asarray([record["coord"] for record in records])
    origin = coords.min(axis=0) - padding_A
    high = coords.max(axis=0) + padding_A
    shape = np.ceil((high - origin) / grid_A).astype(int) + 1
    voxels = int(np.prod(shape))
    if voxels > MAX_VOXELS:
        raise ValueError(
            f"Pocket grid would contain {voxels:,} voxels; reduce the structure size or use a coarser grid"
        )
    return origin.astype(float), tuple(int(x) for x in shape)


def _voxelize(records: list[dict], origin: np.ndarray, shape: tuple[int, ...], grid_A: float) -> np.ndarray:
    solid = np.zeros(shape, dtype=bool)
    for record in records:
        center = record["coord"]
        radius = record["radius_A"]
        lo = np.maximum(0, np.floor((center - radius - origin) / grid_A).astype(int))
        hi = np.minimum(shape, np.ceil((center + radius - origin) / grid_A).astype(int) + 1)
        if np.any(hi <= lo):
            continue
        x = origin[0] + np.arange(lo[0], hi[0]) * grid_A - center[0]
        y = origin[1] + np.arange(lo[1], hi[1]) * grid_A - center[1]
        z = origin[2] + np.arange(lo[2], hi[2]) * grid_A - center[2]
        solid[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] |= (
            x[:, None, None] ** 2 + y[None, :, None] ** 2 + z[None, None, :] ** 2 <= radius ** 2
        )
    return solid


def _border_mask(mask: np.ndarray) -> np.ndarray:
    result = np.zeros_like(mask, dtype=bool)
    result[0, :, :] = result[-1, :, :] = True
    result[:, 0, :] = result[:, -1, :] = True
    result[:, :, 0] = result[:, :, -1] = True
    return result & mask


def _exterior(free: np.ndarray) -> np.ndarray:
    seeds = _border_mask(free)
    return binary_propagation(seeds, structure=NEIGHBORHOOD, mask=free)


def _neighbor_masks(mask: np.ndarray):
    """Yield (axis, direction, neighbor) without periodic wrap-around."""

    for axis in range(3):
        for direction in (-1, 1):
            neighbor = np.zeros_like(mask, dtype=bool)
            source = [slice(None)] * 3
            target = [slice(None)] * 3
            if direction == -1:
                source[axis] = slice(1, None)
                target[axis] = slice(None, -1)
            else:
                source[axis] = slice(None, -1)
                target[axis] = slice(1, None)
            neighbor[tuple(target)] = mask[tuple(source)]
            yield axis, direction, neighbor


def _shifted(mask: np.ndarray, axis: int, direction: int) -> np.ndarray:
    """Return a mask shifted into the neighbor position without wrapping."""

    neighbor = np.zeros_like(mask, dtype=bool)
    source = [slice(None)] * 3
    target = [slice(None)] * 3
    if direction == -1:
        source[axis] = slice(1, None)
        target[axis] = slice(None, -1)
    else:
        source[axis] = slice(None, -1)
        target[axis] = slice(1, None)
    neighbor[tuple(target)] = mask[tuple(source)]
    return neighbor


def _component_geometry(
    component: np.ndarray,
    bulk_external: np.ndarray,
    solid: np.ndarray,
    grid_A: float,
    clearance: np.ndarray | None = None,
    lining_clearance_A: float | None = None,
):
    boundary_faces = 0
    opening_faces = 0
    mouth = np.zeros_like(component, dtype=bool)
    lining_voxels = np.zeros_like(component, dtype=bool)
    for axis, direction, neighbor in _neighbor_masks(component):
        face = component & ~neighbor
        neighbor_bulk = _shifted(bulk_external, axis, direction)
        neighbor_solid = _shifted(solid, axis, direction)
        boundary_faces += int(face.sum())
        opening_faces += int((face & neighbor_bulk).sum())
        mouth |= face & neighbor_bulk
        lining_voxels |= face & neighbor_solid
        # A solvent-center voxel is separated from the atom sphere by the
        # probe radius, so exact solid adjacency is uncommon at 0.6 Å. Use
        # the clearance field to capture the first solvent shell as well.
        if clearance is not None and lining_clearance_A is not None:
            lining_voxels |= face & ~neighbor_bulk & (clearance <= lining_clearance_A)
    boundary_area = boundary_faces * grid_A * grid_A
    opening_area = opening_faces * grid_A * grid_A
    return boundary_area, opening_area, mouth, lining_voxels


def _geodesic_depth(component: np.ndarray, mouth: np.ndarray, grid_A: float) -> float | None:
    if not np.any(mouth):
        return None
    distance = np.full(component.shape, -1, dtype=np.int32)
    queue = deque(int(x) for x in np.flatnonzero(mouth))
    distance[mouth] = 0
    shape = component.shape
    while queue:
        flat = queue.popleft()
        point = np.unravel_index(flat, shape)
        for axis in range(3):
            for direction in (-1, 1):
                nxt = list(point)
                nxt[axis] += direction
                if not (0 <= nxt[axis] < shape[axis]):
                    continue
                nxt = tuple(nxt)
                if component[nxt] and distance[nxt] < 0:
                    distance[nxt] = distance[point] + 1
                    queue.append(np.ravel_multi_index(nxt, shape))
    valid = distance[component]
    return round(float(valid.max()) * grid_A, 2) if len(valid) else None


def _mesh(component: np.ndarray, origin: np.ndarray, grid_A: float) -> dict | None:
    try:
        from skimage.measure import marching_cubes
    except ImportError:
        return None
    if int(component.sum()) < 4:
        return None
    try:
        vertices, faces, _, _ = marching_cubes(component.astype(np.float32), level=0.5, spacing=(grid_A,) * 3)
    except (RuntimeError, ValueError):
        # A component that is entirely clipped by the selected search window
        # may have no level-crossing surface. Geometry metrics remain valid;
        # omit only the optional visualization artifact.
        return None
    vertices = vertices + origin
    # Meshes are visualization artifacts. Keep them bounded for a responsive API.
    if len(vertices) > 20_000:
        stride = int(np.ceil(len(vertices) / 20_000))
        # Sample faces and remap their referenced vertices. Returning vertices
        # without faces would make the browser silently disable the requested
        # pocket highlight for larger components.
        chosen_faces = faces[::stride]
        if len(chosen_faces):
            vertex_ids, inverse = np.unique(chosen_faces, return_inverse=True)
            vertices = vertices[vertex_ids]
            faces = inverse.reshape(-1, 3)
        return {
            "vertices": np.round(vertices, 3).tolist(),
            "faces": faces.astype(int).tolist(),
            "decimated": True,
        }
    return {
        "vertices": np.round(vertices, 3).tolist(),
        "faces": faces.astype(int).tolist(),
        "decimated": False,
    }


def _lining(records: list[dict], lining_voxels: np.ndarray, origin: np.ndarray, grid_A: float, sasa: dict):
    if not np.any(lining_voxels):
        return [], 0.0
    coordinates = np.asarray([record["coord"] for record in records])
    tree = cKDTree(coordinates)
    voxel_indices = np.argwhere(lining_voxels)
    points = origin + voxel_indices * grid_A
    _, nearest = tree.query(points)
    keys = {}
    for index in np.atleast_1d(nearest):
        record = records[int(index)]
        ident = record["residue"]
        keys[ident["key"]] = ident
    output = []
    total = 0.0
    for key, ident in sorted(keys.items(), key=lambda item: (item[1]["chain_id"], item[1]["number"], item[1]["insertion_code"])):
        value = float(sasa.get((ident["chain_id"], key), 0.0))
        total += value
        output.append({**ident, "sasa_A2": round(value, 2)})
    return output, total


def _search_mask(
    origin: np.ndarray,
    shape: tuple[int, ...],
    grid_A: float,
    centers: np.ndarray,
    radius_A: float,
) -> np.ndarray:
    """Return the union of spheres without allocating a whole-grid coordinate mesh."""

    mask = np.zeros(shape, dtype=bool)
    centers = np.asarray(centers, dtype=float).reshape(-1, 3)
    if not len(centers):
        return mask
    lo = np.maximum(0, np.floor((centers.min(axis=0) - radius_A - origin) / grid_A).astype(int))
    hi = np.minimum(shape, np.ceil((centers.max(axis=0) + radius_A - origin) / grid_A).astype(int) + 1)
    if np.any(hi <= lo):
        return mask
    region = mask[tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))]
    tree = cKDTree(centers)
    # Flat chunks also bound allocations when listed residues are far apart.
    for start in range(0, region.size, 65_536):
        stop = min(start + 65_536, region.size)
        indices = np.column_stack(np.unravel_index(np.arange(start, stop), region.shape))
        points = origin + (indices + lo) * grid_A
        distances, _ = tree.query(points)
        region.flat[start:stop] = distances <= radius_A
    return mask


def _state_result(
    model,
    state_name: str,
    state_chains: list[str],
    centers: np.ndarray,
    search_mask: np.ndarray,
    origin: np.ndarray,
    shape: tuple[int, ...],
    grid_A: float,
    solvent_probe_A: float,
    bulk_probe_A: float,
    defined_residues: list[dict] | None = None,
):
    records = _atom_records(model, state_chains)
    solid = _voxelize(records, origin, shape, grid_A)
    clearance = distance_transform_edt(~solid, sampling=grid_A)
    small_free = clearance >= solvent_probe_A
    large_free = clearance >= bulk_probe_A
    bulk_external = _exterior(large_free)
    candidates = small_free & ~bulk_external
    full_labels, _ = label(candidates, structure=NEIGHBORHOOD)
    full_sizes = np.bincount(full_labels.ravel())
    inside_sizes = np.bincount(full_labels[search_mask], minlength=len(full_sizes))
    if defined_residues is None:
        # Retain the historical component identities and largest-first ranking.
        labels = np.where(search_mask, full_labels, 0)
    else:
        # A union of windows can split one global component into multiple pockets.
        labels, _ = label(candidates & search_mask, structure=NEIGHBORHOOD)
    sizes = np.bincount(labels.ravel())
    candidate_ids = [int(x) for x in np.flatnonzero(sizes[1:]) + 1]
    bounds = find_objects(labels)
    sasa = residue_sasa(model, state_chains)

    def region(component_id):
        slices = tuple(
            slice(max(0, part.start - 1), min(shape[axis], part.stop + 1))
            for axis, part in enumerate(bounds[component_id - 1])
        )
        local_origin = origin + np.array([part.start for part in slices]) * grid_A
        return slices, labels[slices] == component_id, local_origin

    present_keys = [
        item["key"] for item in (defined_residues or [])
        if item["chain_id"] in state_chains
    ]
    ranked = []
    if defined_residues is None:
        candidate_ids.sort(key=lambda item: (-int(full_sizes[item]), item))
        ranked = [{"component_id": item} for item in candidate_ids[:MAX_COMPONENTS]]
    else:
        # Coverage must be evaluated before limiting the returned components.
        # Only retained candidates need depth searches and visualization meshes.
        for component_id in candidate_ids:
            slices, component, local_origin = region(component_id)
            _, _, _, lining_voxels = _component_geometry(
                component, bulk_external[slices], solid[slices], grid_A,
                clearance[slices], solvent_probe_A + 2.0 * grid_A,
            )
            lining, lining_sasa = _lining(records, lining_voxels, local_origin, grid_A, sasa)
            lining_keys = {item["key"] for item in lining}
            defined_lining = [key for key in present_keys if key in lining_keys]
            coverage = len(defined_lining) / len(present_keys) if present_keys else None
            ranked.append({
                "component_id": component_id,
                "lining_residues": lining,
                "lining_sasa_A2": lining_sasa,
                "defined_residues_lining": defined_lining,
                "defined_residues_present": len(present_keys),
                "defined_residue_coverage": coverage,
            })
        ranked.sort(key=lambda item: (
            -(item["defined_residue_coverage"] or 0.0),
            -int(sizes[item["component_id"]]),
            item["component_id"],
        ))
        ranked = ranked[:MAX_COMPONENTS]

    center_tree = cKDTree(centers)
    pockets = []
    for entry in ranked:
        component_id = entry["component_id"]
        slices, component, local_origin = region(component_id)
        volume = float(component.sum()) * grid_A ** 3
        indices = np.argwhere(component)
        nearest_anchor = float(center_tree.query(local_origin + indices * grid_A)[0].min())
        boundary_area, opening_area, mouth, lining_voxels = _component_geometry(
            component, bulk_external[slices], solid[slices], grid_A,
            clearance[slices], solvent_probe_A + 2.0 * grid_A,
        )
        if defined_residues is None:
            lining, lining_sasa = _lining(records, lining_voxels, local_origin, grid_A, sasa)
        else:
            lining, lining_sasa = entry["lining_residues"], entry["lining_sasa_A2"]
        depth = _geodesic_depth(component, mouth, grid_A)
        touches_grid_boundary = any(
            (part.start == 0 and np.any(np.take(component, 0, axis=axis)))
            or (part.stop == shape[axis] and np.any(np.take(component, -1, axis=axis)))
            for axis, part in enumerate(slices)
        )
        first_index = tuple(indices[0] + np.array([part.start for part in slices]))
        parent_id = int(full_labels[first_index])
        touches_search_boundary = bool(full_sizes[parent_id] > inside_sizes[parent_id])
        touches_boundary = touches_grid_boundary or touches_search_boundary
        enclosure = 1.0 if boundary_area <= 0 else max(0.0, min(1.0, 1.0 - opening_area / boundary_area))
        class_name = "truncated" if touches_boundary else ("solvent_accessible" if np.any(mouth) else "enclosed")
        deepest = None
        if depth is not None:
            distances_from_mouth = np.full(component.shape, -1, dtype=np.int32)
            # Reuse the small BFS's result cheaply for the deepest point marker.
            queue = deque(int(x) for x in np.flatnonzero(mouth))
            distances_from_mouth[mouth] = 0
            while queue:
                flat = queue.popleft()
                point = np.unravel_index(flat, component.shape)
                for axis in range(3):
                    for direction in (-1, 1):
                        nxt = list(point); nxt[axis] += direction
                        if 0 <= nxt[axis] < component.shape[axis]:
                            nxt = tuple(nxt)
                            if component[nxt] and distances_from_mouth[nxt] < 0:
                                distances_from_mouth[nxt] = distances_from_mouth[point] + 1
                                queue.append(np.ravel_multi_index(nxt, component.shape))
            deepest_index = np.argwhere(distances_from_mouth == distances_from_mouth[component].max())[0]
            deepest = (local_origin + deepest_index * grid_A).round(3).tolist()
        pocket = {
            "component_id": component_id,
            "state": state_name,
            "class": class_name,
            "volume_A3": round(volume, 2),
            "depth_A": depth,
            "opening_area_A2": round(opening_area, 2),
            "boundary_area_A2": round(boundary_area, 2),
            "enclosure_fraction": round(enclosure, 4),
            "lining_sasa_A2": round(lining_sasa, 2),
            "solvent_exposure_A2": round(opening_area + lining_sasa, 2),
            "anchor_distance_A": round(nearest_anchor, 2),
            "lining_residues": lining,
            "deepest_point": deepest,
            "mouth_voxel_count": int(mouth.sum()),
            "truncated": touches_boundary,
            "search_region_truncated": touches_search_boundary,
            "mesh": _mesh(component, local_origin, grid_A),
        }
        if defined_residues is not None:
            for field in ("defined_residues_lining", "defined_residues_present", "defined_residue_coverage"):
                pocket[field] = entry[field]
        pockets.append(pocket)
    return {
        "state": state_name,
        "parameters": {
            "grid_A": grid_A,
            "solvent_probe_A": solvent_probe_A,
            "bulk_probe_A": bulk_probe_A,
            "definition": "standard-probe space not connected to bulk-probe exterior",
        },
        "pockets": pockets,
        "primary": pockets[0] if pockets else None,
    }


def _anchor_residue(model, target_chain: str, anchor: dict):
    residues = residues_for_chain(model, target_chain)
    number = int(anchor["number"])
    insertion = str(anchor.get("insertion_code", "") or "")
    for residue in residues:
        ident = residue_identity(target_chain, residue)
        if ident["number"] == number and ident["insertion_code"] == insertion:
            return residue, ident
    raise ValueError(f"Residue {target_chain}:{number}{insertion} was not found")


def _residue_set(model, target_chain: str, partner_chain: str, references: list[dict]):
    if not references:
        raise ValueError("Define a pocket with at least one residue")
    if len(references) > MAX_POCKET_RESIDUES:
        raise ValueError(f"Pocket definitions are limited to {MAX_POCKET_RESIDUES} residues")
    identities = {}
    coordinates = []
    for reference in references:
        chain_id = reference.get("chain_id")
        if chain_id is None:
            chain_id = target_chain
        if chain_id not in {target_chain, partner_chain}:
            raise ValueError("Pocket residues must belong to the target or partner chain")
        residue, ident = _anchor_residue(model, chain_id, reference)
        if ident["key"] in identities:
            continue
        atoms = [atom.coord for atom in residue if is_heavy_atom(atom)]
        if not atoms:
            raise ValueError(f"Residue {ident['key']} has no supported heavy atoms")
        identities[ident["key"]] = ident
        coordinates.extend(atoms)
    return list(identities.values()), np.asarray(coordinates, dtype=float)


def analyze_pocket(
    model,
    target_chain: str,
    partner_chain: str,
    anchor_residue: dict | None = None,
    radius_A: float = 8.0,
    *,
    pocket_residues: list[dict] | None = None,
    grid_A: float = GRID_A,
    solvent_probe_A: float = SOLVENT_PROBE_A,
    bulk_probe_A: float = BULK_PROBE_A,
) -> dict:
    if not (4.0 <= float(radius_A) <= 12.0):
        raise ValueError("Pocket search radius must be between 4 and 12 Å")
    if target_chain == partner_chain:
        raise ValueError("Target and partner chains must be different")
    if (anchor_residue is None) == (pocket_residues is None):
        raise ValueError("Provide either an anchor residue or pocket residues, but not both")
    defined_residues = None
    anchor_ident = None
    if pocket_residues is not None:
        defined_residues, centers = _residue_set(model, target_chain, partner_chain, pocket_residues)
        anchor = centers.mean(axis=0)
    else:
        residue, anchor_ident = _anchor_residue(model, target_chain, anchor_residue)
        anchor_atoms = [a.coord for a in residue if is_heavy_atom(a)]
        if not anchor_atoms:
            raise ValueError("The anchor residue has no supported heavy atoms")
        anchor = np.mean(anchor_atoms, axis=0).astype(float)
        centers = anchor[None, :]
    all_records = _atom_records(model, [target_chain, partner_chain])
    origin, shape = _grid_for(all_records, grid_A, max(12.0, float(radius_A) + bulk_probe_A + 1.0))
    search_mask = _search_mask(origin, shape, grid_A, centers, float(radius_A))
    free = _state_result(
        model, "free", [target_chain], centers, search_mask, origin, shape,
        grid_A, solvent_probe_A, bulk_probe_A, defined_residues,
    )
    bound = _state_result(
        model, "bound", [target_chain, partner_chain], centers, search_mask, origin, shape,
        grid_A, solvent_probe_A, bulk_probe_A, defined_residues,
    )
    free_primary = free["primary"]
    bound_primary = bound["primary"]
    comparison = {}
    if free_primary and bound_primary:
        for field in ("volume_A3", "depth_A", "opening_area_A2", "enclosure_fraction", "solvent_exposure_A2"):
            a, b = free_primary.get(field), bound_primary.get(field)
            comparison[f"delta_{field}"] = None if a is None or b is None else round(float(b - a), 4)
        if bound_primary["class"] == "enclosed" and free_primary["class"] != "enclosed":
            comparison["state_change"] = "occluded"
        elif bound_primary["volume_A3"] < 0.75 * free_primary["volume_A3"]:
            comparison["state_change"] = "collapsed"
        elif bound_primary["volume_A3"] > 1.25 * free_primary["volume_A3"]:
            comparison["state_change"] = "opened"
        else:
            comparison["state_change"] = "unchanged"
    else:
        comparison["state_change"] = "not_comparable"
    return {
        "target_chain": target_chain,
        "partner_chain": partner_chain,
        "mode": "residues" if defined_residues is not None else "anchor",
        "pocket_residues": defined_residues or [],
        "anchor_residue": anchor_ident,
        "anchor_coord": anchor.round(3).tolist(),
        "radius_A": float(radius_A),
        "free": free,
        "bound": bound,
        "comparison": comparison,
        "warnings": [
            "Volume is estimated on a voxel grid; rerun with a finer grid for convergence checks.",
            "The free state removes the partner without relaxing the target coordinates.",
        ],
    }
