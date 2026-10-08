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
        "surface-mode-control",
        "surface-mode-convexity",
        "surface-mode-electrostatics",
        "electrostatics-status",
        "electrostatics-retry",
        "pocket-result",
    ):
        assert f'id="{element_id}"' in html
    assert 'id="residue-panel"' in html
    assert 'id="contact-panel"' in html
    assert "Residue interface" in html
    assert "Closest contact pairs" in html
    for element_id in ("source-drawer", "drawer-toggle", "drawer-open", "drawer-backdrop", "main-column",
                       "upload-form", "upload-submit", "pdb-form", "history", "clear-history"):
        assert f'id="{element_id}"' in html
    # The 3D surface is the first panel in the workspace.
    assert html.index('id="surface-panel"') < html.index('id="pair-workspace"')


def test_browser_orchestration_contains_linked_selection_and_pocket_actions():
    javascript = (STATIC / "app.js").read_text()
    for function_name in (
        "renderContactMap",
        "selectResidues",
        "renderSurface",
        "setSurfaceMode",
        "loadElectrostatics",
        "meshColors",
        "showPocketHighlight",
        "clearPocketHighlight",
        "stopPocketPolling",
    ):
        assert f"function {function_name}" in javascript
    assert "contact_distance_dA" in javascript
    assert "reference_dSASA_A2" in javascript
    assert "Highlight bound pocket" in javascript
    assert "Cancel highlight" in javascript
    assert "potential_kT_e" in javascript
    assert "electrostaticsColors" in javascript
    assert "convexitySurfaceColors" in javascript


def test_history_renders_chain_strips_and_overlay_drawer():
    javascript = (STATIC / "app.js").read_text()
    for function_name in ("chainStrip", "historyRow", "openDrawer", "closeDrawer", "bindDropZone"):
        assert f"function {function_name}" in javascript
    assert "max-width: 1100px" in javascript


def test_residue_highlight_panel_and_residue_defined_pockets():
    html = (STATIC / "index.html").read_text()
    for element_id in ("surface-stage", "residue-highlight-panel", "residue-highlight-toggle",
                       "residue-highlight-body", "highlight-chain", "highlight-residues",
                       "highlight-feedback", "pocket-residues", "pocket-use-selection"):
        assert f'id="{element_id}"' in html
    # The side panel sits beside the NGL viewport inside the same stage.
    stage = html.index('id="surface-stage"')
    assert stage < html.index('id="surface-view"') < html.index('id="residue-highlight-panel"')
    assert 'name="pocket-mode" value="residues"' in html
    javascript = (STATIC / "app.js").read_text()
    for function_name in ("parseResidueSpec", "applyResidueHighlight", "setHighlightPanelCollapsed",
                          "setPocketMode", "pocketResidueRequest"):
        assert f"function {function_name}" in javascript
    assert "pocket_residues" in javascript
