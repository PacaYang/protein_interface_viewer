from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "static"


def test_workspace_contains_requested_sections_and_controls():
    html = (STATIC / "index.html").read_text()
    for element_id in (
        "interfaces-sidebar",
        "interfaces-toggle",
        "pair-workspace",
        "residue-tab",
        "contact-tab",
        "contact-grid",
        "surface-resolution",
        "surface-opacity",
        "surface-padding",
        "surface-separation",
        "chain-toggles",
        "pocket-result",
    ):
        assert f'id="{element_id}"' in html
    assert 'id="residue-panel"' in html
    assert 'id="contact-panel"' in html
    assert "Residue interface" in html
    assert "Closest contact pairs" in html


def test_browser_orchestration_contains_linked_selection_and_pocket_actions():
    javascript = (STATIC / "app.js").read_text()
    for function_name in (
        "renderContactMap",
        "selectResidues",
        "renderSurface",
        "showPocketHighlight",
        "clearPocketHighlight",
        "stopPocketPolling",
    ):
        assert f"function {function_name}" in javascript
    assert "contact_distance_dA" in javascript
    assert "reference_dSASA_A2" in javascript
    assert "Highlight bound pocket" in javascript
    assert "Cancel highlight" in javascript
