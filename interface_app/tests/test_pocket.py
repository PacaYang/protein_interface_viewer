import numpy as np
from pathlib import Path

from interface_app.pocket import _component_geometry, _geodesic_depth
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
