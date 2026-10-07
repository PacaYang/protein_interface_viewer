"""Plane orientation and nearest-residue matching use original coordinates."""

from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue
from scipy.spatial.transform import Rotation

from interface_app.analysis import _atom_records, encode_pair
from interface_app.constants import residue_identity
from interface_app.face_view import binding_face_view, nearest_residues
from interface_app.structure import parse_structure


def fixture(points_a, points_b, *, chain_ids=("A", "B"), transform=None):
    model = Model.Model(0)
    rows = []
    for chain_id, points in zip(chain_ids, (points_a, points_b)):
        chain = Chain.Chain(chain_id)
        model.add(chain)
        for index, point in enumerate(points):
            residue = Residue.Residue((" ", index + 1, "A" if index == 1 else " "), "ALA", "")
            chain.add(residue)
            coords = np.asarray(point, dtype=float)
            if transform:
                coords = transform(coords)
            residue.add(Atom.Atom("CA", coords, 1, 1, " ", " CA ", index + 1, element="C"))
            rows.append(residue_identity(chain_id, residue))
    return model, {"id": encode_pair(*chain_ids), "chain_a": chain_ids[0], "chain_b": chain_ids[1],
                   "parameters": {"contact_cutoff_A": 5.0}, "residues": rows}


FLAT = np.array([[-4, -1, 0], [4, -1, 0], [3, 1, 0], [-3, 1, 0]])


def assert_frames(view):
    assert view["available"], view
    for data in view["chains"].values():
        normal, up = np.asarray(data["normal"]), np.asarray(data["up"])
        assert np.linalg.norm(normal) == pytest.approx(1)
        assert np.linalg.norm(up) == pytest.approx(1)
        assert np.dot(normal, up) == pytest.approx(0, abs=1e-12)
        # Frontend uses this proper rotation with normal toward camera -Z.
        frame = np.array([np.cross(normal, up), up, -normal])
        assert frame @ frame.T == pytest.approx(np.eye(3), abs=1e-12)
        assert np.linalg.det(frame) == pytest.approx(1)
        assert frame @ normal == pytest.approx([0, 0, -1], abs=1e-12)


@pytest.mark.parametrize("rotation", [np.eye(3), Rotation.from_euler("xyz", [42, -18, 61], degrees=True).as_matrix()])
def test_flat_and_rotated_planes_face_each_other(rotation):
    model, pair = fixture(FLAT, FLAT + [0, 0, 3], transform=lambda p: rotation @ p + [11, -5, 32])
    view = binding_face_view(model, pair)
    assert_frames(view)
    assert view["chains"]["A"]["normal"] == pytest.approx(rotation @ [0, 0, 1])
    assert view["chains"]["B"]["normal"] == pytest.approx(rotation @ [0, 0, -1])
    assert view["chains"]["A"]["method"] == "residue_centroid_plane"
    assert np.dot(view["chains"]["A"]["up"], view["chains"]["B"]["up"]) == pytest.approx(1)
    assert binding_face_view(model, pair) == view


def test_residue_centroids_have_equal_weight_despite_atom_counts():
    model, pair = fixture(FLAT, FLAT + [0, 0, 3])
    residue = model["A"][1]
    for i, name in enumerate(["N", "C", "O", "CB"]):
        residue.add(Atom.Atom(name, FLAT[0], 1, 1, " ", name, 50+i, element=name[0]))
    residue.add(Atom.Atom("H", FLAT[0] + [0, 0, 3], 1, 1, " ", " H ", 90, element="H"))
    residue.add(Atom.Atom("OXT", FLAT[0] + [0, 0, 3], 1, 1, "B", " OXT", 91, element="O"))
    view = binding_face_view(model, pair)
    assert view["chains"]["A"]["center"] == pytest.approx(FLAT.mean(axis=0))
    assert view["nearest"]["A"]["A:1"]["distance_A"] == pytest.approx(3)


@pytest.mark.parametrize("points", [
    [[0, 0, 0]], [[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
    [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]],
])
def test_sparse_collinear_and_ambiguous_planes_use_estimated_direction(points):
    points = np.array(points)
    model, pair = fixture(points, points + [0, 0, 3])
    view = binding_face_view(model, pair)
    assert_frames(view)
    assert all(data["estimated"] for data in view["chains"].values())
    assert view["chains"]["A"]["normal"] == pytest.approx([0, 0, 1])


def test_roll_fallback_handles_circular_plane_and_parallel_partner_normal():
    square = np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]])
    # A's deterministic in-plane X is parallel to B's normal.
    second = np.array([[3, -1, 2], [3, 1, 2], [3, 1, 4], [3, -1, 4]])
    model, pair = fixture(square, second)
    view = binding_face_view(model, pair)
    assert_frames(view)
    assert view["chains"]["A"]["up"] == pytest.approx([1, 0, 0])
    assert binding_face_view(model, pair) == view


def test_missing_binding_residues_and_coincident_centroids_are_unavailable():
    model, pair = fixture(FLAT, FLAT)
    assert not binding_face_view(model, pair)["available"]
    pair["residues"] = []
    view = binding_face_view(model, pair)
    assert not view["available"]
    assert "No binding residues" in view["reason"]
    assert view["chains"] == view["nearest"] == {}


def test_whole_chain_centroids_supply_direction_when_binding_centroids_coincide():
    model, pair = fixture(np.vstack([FLAT, [0, 0, -50]]), np.vstack([FLAT, [0, 0, 50]]))
    pair["residues"] = [row for row in pair["residues"] if row["number"] != 5]
    view = binding_face_view(model, pair)
    assert_frames(view)
    assert view["chains"]["A"]["direction_source"] == "chain_centroids"


def test_unrestricted_nearest_matches_brute_force_in_both_directions():
    loaded = parse_structure(Path(__file__).resolve().parents[2] / "examples" / "1IAR.pdb")
    atoms = [_atom_records(loaded.model, chain) for chain in ("A", "B")]
    distant = 0
    for source, partner in (atoms, atoms[::-1]):
        partner_order = list(dict.fromkeys(atom["residue_key"] for atom in partner))
        mapping = nearest_residues(source, partner)
        assert len(mapping) == len(set(atom["residue_key"] for atom in source))
        for key, match in mapping.items():
            coords = np.array([atom["coord"] for atom in source if atom["residue_key"] == key])
            candidates = []
            for partner_key in partner_order:
                target = np.array([atom["coord"] for atom in partner if atom["residue_key"] == partner_key])
                distance = float(np.linalg.norm(coords[:, None] - target, axis=2).min())
                candidates.append((distance, partner_key))
            expected_distance, expected_key = min(candidates, key=lambda item: (item[0], partner_order.index(item[1])))
            assert match["partner_key"] == expected_key
            assert match["distance_A"] == pytest.approx(expected_distance, abs=1e-12)
            distant += match["distance_A"] > 6
    assert distant > 50


def test_nearest_ties_use_sequence_order_and_preserve_chain_and_insertion_ids():
    model, pair = fixture([[0, 0, 0], [100, 0, 0]], [[1, 0, 0], [-1, 0, 0]], chain_ids=("alpha", "beta"))
    view = binding_face_view(model, pair)
    assert_frames(view)
    assert view["nearest"]["alpha"]["alpha:1"] == {"partner_key": "beta:1", "distance_A": 1.0}
    assert view["nearest"]["beta"]["beta:2A"]["partner_key"] == "alpha:1"
    # Distances differing below display precision must not be considered tied.
    model["beta"][(" ", 2, "A")]["CA"].coord = np.array([-0.99999, 0, 0])
    match = binding_face_view(model, pair)["nearest"]["alpha"]["alpha:1"]
    assert match["partner_key"] == "beta:2A"
    assert match["distance_A"] == 0.99999
