import numpy as np
import pytest
from pathlib import Path

from interface_app.analysis import analyze_loaded
from interface_app import pocket
from interface_app.pocket import _component_geometry, _geodesic_depth, _residue_set, _search_mask, _state_result
from interface_app.pocket import analyze_pocket
from interface_app.structure import parse_structure

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def test_open_component_has_mouth_and_depth():
    component = np.zeros((9, 9, 9), dtype=bool)
    component[2:7, 2:7, 1:7] = True
    bulk = np.zeros_like(component)
    bulk[2:7, 2:7, 7] = True
    solid = np.zeros_like(component)
    boundary, opening, mouth, _ = _component_geometry(component, bulk, solid, 0.6)
    assert boundary > opening > 0
    assert mouth.any()
    assert _geodesic_depth(component, mouth, 0.6) is not None


def test_closed_component_depth_is_unavailable():
    component = np.zeros((7, 7, 7), dtype=bool)
    component[2:5, 2:5, 2:5] = True
    bulk = np.zeros_like(component)
    solid = np.zeros_like(component)
    _, opening, mouth, _ = _component_geometry(component, bulk, solid, 0.6)
    assert opening == 0
    assert not mouth.any()
    assert _geodesic_depth(component, mouth, 0.6) is None


def test_real_complex_returns_free_bound_measurements():
    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    result = analyze_pocket(loaded.model, "A", "B", {"number": 67}, radius_A=8.0)
    assert result["free"]["primary"] is not None
    assert result["bound"]["primary"] is not None
    assert result["free"]["primary"]["volume_A3"] > 0
    assert "delta_volume_A3" in result["comparison"]


def _legacy_sphere_mask(origin, shape, grid_A, anchor, radius_A):
    mask = np.zeros(shape, dtype=bool)
    anchor_index = np.rint((anchor - origin) / grid_A).astype(int)
    radius_voxels = int(np.ceil(radius_A / grid_A))
    lo = np.maximum(0, anchor_index - radius_voxels)
    hi = np.minimum(np.asarray(shape), anchor_index + radius_voxels + 1)
    xs, ys, zs = (origin[axis] + np.arange(lo[axis], hi[axis]) * grid_A for axis in range(3))
    distances = ((xs[:, None, None] - anchor[0]) ** 2 + (ys[None, :, None] - anchor[1]) ** 2
                 + (zs[None, None, :] - anchor[2]) ** 2)
    mask[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = distances <= radius_A ** 2
    return mask


def test_single_center_search_mask_matches_anchor_sphere():
    origin = np.array([-3.1, 0.2, 5.0])
    shape = (60, 55, 50)
    anchor = np.array([10.37, 12.0, 19.91])
    expected = _legacy_sphere_mask(origin, shape, 0.6, anchor, 8.0)
    assert np.array_equal(_search_mask(origin, shape, 0.6, anchor[None, :], 8.0), expected)


def test_search_mask_is_union_of_spheres():
    origin = np.zeros(3)
    shape = (40, 40, 40)
    centers = np.array([[5.0, 5.0, 5.0], [15.0, 5.0, 5.0]])
    union = _search_mask(origin, shape, 0.6, centers, 3.0)
    expected = _search_mask(origin, shape, 0.6, centers[:1], 3.0) | _search_mask(origin, shape, 0.6, centers[1:], 3.0)
    assert np.array_equal(union, expected)


def test_residue_defined_pocket_reports_coverage():
    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    residues = [{"number": number} for number in (63, 64, 66, 67, 70)]
    result = analyze_pocket(loaded.model, "A", "B", radius_A=6.0, pocket_residues=residues)
    assert result["mode"] == "residues"
    assert result["anchor_residue"] is None
    assert [item["key"] for item in result["pocket_residues"]] == ["A:63", "A:64", "A:66", "A:67", "A:70"]
    for state in ("free", "bound"):
        primary = result[state]["primary"]
        assert primary is not None
        assert primary["defined_residues_present"] == 5
        assert 0.0 <= primary["defined_residue_coverage"] <= 1.0
        lining = {item["key"] for item in primary["lining_residues"]}
        assert set(primary["defined_residues_lining"]) <= lining
        pockets = result[state]["pockets"]
        assert primary["defined_residue_coverage"] == max(item["defined_residue_coverage"] for item in pockets)


def test_partner_residues_are_excluded_from_free_state_coverage():
    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    contact = next(iter(analyze_loaded(loaded)["pairs"][0]["contacts"]))
    residues = [
        {"chain_id": "A", "number": contact["residue_a"]["number"], "insertion_code": contact["residue_a"]["insertion_code"]},
        {"chain_id": "B", "number": contact["residue_b"]["number"], "insertion_code": contact["residue_b"]["insertion_code"]},
    ]
    result = analyze_pocket(loaded.model, "A", "B", radius_A=8.0, pocket_residues=residues)
    for item in result["free"]["pockets"]:
        assert item["defined_residues_present"] == 1
    for item in result["bound"]["pockets"]:
        assert item["defined_residues_present"] == 2


@pytest.mark.parametrize("kwargs, message", [
    ({}, "either"),
    ({"anchor_residue": {"number": 67}, "pocket_residues": [{"number": 67}]}, "either"),
    ({"pocket_residues": []}, "at least one"),
    ({"pocket_residues": [{"number": 99999}]}, "not found"),
    ({"pocket_residues": [{"chain_id": "Z", "number": 67}]}, "target or partner"),
    ({"pocket_residues": [{"number": 67}] * 61}, "limited"),
])
def test_invalid_pocket_definitions_are_rejected(kwargs, message):
    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    with pytest.raises(ValueError, match=message):
        analyze_pocket(loaded.model, "A", "B", radius_A=8.0, **kwargs)


def test_search_mask_handles_strided_chunks_and_outside_centers():
    origin = np.zeros(3)
    shape = (90, 95, 100)
    center = np.array([45.0, 47.0, 50.0])
    x, y, z = np.ogrid[:shape[0], :shape[1], :shape[2]]
    expected = ((x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2) <= 38.0 ** 2
    assert np.array_equal(_search_mask(origin, shape, 1.0, center[None, :], 38.0), expected)
    assert not _search_mask(origin, shape, 1.0, np.array([[-100.0] * 3]), 4.0).any()


def test_residue_resolution_deduplicates_and_defaults_null_chain():
    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    identities, coordinates = _residue_set(
        loaded.model, "A", "B",
        [{"number": 67}, {"chain_id": None, "number": 67}, {"chain_id": "B", "number": 12}],
    )
    assert [item["key"] for item in identities] == ["A:67", "B:12"]
    _, expected = _residue_set(loaded.model, "A", "B", [{"number": 67}, {"chain_id": "B", "number": 12}])
    assert np.array_equal(coordinates, expected)


def test_residue_without_heavy_atoms_is_rejected():
    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    residue = loaded.model["A"][(" ", 67, " ")]
    for atom in list(residue):
        residue.detach_child(atom.id)
    with pytest.raises(ValueError, match="no supported heavy atoms"):
        _residue_set(loaded.model, "A", "B", [{"number": 67}])


def synthetic_state(monkeypatch, candidates, mask, defined_residues):
    """Supply a known clearance field to isolate region grouping and ranking."""
    clearance = np.where(candidates, 2.0, 5.0)
    ident = {"chain_id": "A", "number": 99, "insertion_code": "", "key": "A:99"}
    records = [{"coord": np.array([0.0, 0.0, 0.0]), "radius_A": 1.7, "residue": ident}]
    monkeypatch.setattr(pocket, "_atom_records", lambda *_: records)
    monkeypatch.setattr(pocket, "_voxelize", lambda *_: np.zeros(candidates.shape, dtype=bool))
    monkeypatch.setattr(pocket, "distance_transform_edt", lambda *_args, **_kwargs: clearance)
    monkeypatch.setattr(pocket, "residue_sasa", lambda *_: {})
    return _state_result(
        None, "free", ["A"], np.array([[3.0, 3.0, 3.0]]), mask,
        np.zeros(3), candidates.shape, 1.0, 1.4, 4.0, defined_residues,
    )


def test_coverage_is_ranked_before_five_component_limit(monkeypatch):
    candidates = np.zeros((25, 7, 7), dtype=bool)
    candidates[2:23:4, 3, 3] = True
    candidates[2, 4, 3] = True  # The first component is larger than the matching sixth one.
    ident = {"chain_id": "A", "number": 99, "insertion_code": "", "key": "A:99"}
    monkeypatch.setattr(pocket, "_lining", lambda records, voxels, origin, *_:
                        ([ident], 1.0) if origin[0] > 20 else ([], 0.0))
    mesh_calls = []
    monkeypatch.setattr(pocket, "_mesh", lambda component, *_: mesh_calls.append(int(component.sum())))
    result = synthetic_state(monkeypatch, candidates, np.ones_like(candidates), [ident])
    assert len(result["pockets"]) == 5
    assert result["primary"]["component_id"] == 6
    assert result["primary"]["defined_residue_coverage"] == 1.0
    assert result["primary"]["volume_A3"] == 1.0
    assert len(mesh_calls) == 5


def test_disconnected_search_regions_are_separate_pockets(monkeypatch):
    candidates = np.zeros((15, 7, 7), dtype=bool)
    candidates[2:13, 3, 3] = True  # One global component joins outside the search windows.
    mask = np.zeros_like(candidates)
    mask[2:5, 2:5, 2:5] = True
    mask[10:13, 2:5, 2:5] = True
    result = synthetic_state(monkeypatch, candidates, mask, [])
    assert len(result["pockets"]) == 2
    assert all(item["volume_A3"] == 3.0 for item in result["pockets"])
    assert all(item["search_region_truncated"] for item in result["pockets"])


def test_partner_only_definition_has_no_free_state_coverage(monkeypatch):
    candidates = np.zeros((7, 7, 7), dtype=bool)
    candidates[3, 3, 3] = True
    partner = {"chain_id": "B", "number": 12, "key": "B:12"}
    result = synthetic_state(monkeypatch, candidates, np.ones_like(candidates), [partner])
    assert result["primary"] is not None
    assert result["primary"]["defined_residues_present"] == 0
    assert result["primary"]["defined_residues_lining"] == []
    assert result["primary"]["defined_residue_coverage"] is None
