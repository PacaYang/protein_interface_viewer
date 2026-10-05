"""Exercise the workspace with real coordinates, API jobs, and NGL/WebGL.

Run separately from pytest; requires Playwright and a local Chrome browser.
All fixtures, server history, logs, and screenshots stay in the artifact folder.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
import urllib.request


ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"
sys.path.insert(0, str(ROOT))


def prepare_fixtures(artifacts: Path, reuse: bool) -> dict:
    fixture_path = artifacts / "fixtures.json"
    if reuse and fixture_path.exists():
        return json.loads(fixture_path.read_text())

    from Bio.PDB import MMCIFIO

    from interface_app.analysis import analyze_loaded
    from interface_app.curvature import compute_curvature
    from interface_app.storage import Store
    from interface_app.structure import parse_structure

    loaded = parse_structure(EXAMPLES / "1IAR.pdb")
    result = analyze_loaded(loaded)
    result["version"] = "0.2.0"
    surface = compute_curvature(loaded.model, result)
    assert surface["status"] == "complete", surface
    cif_path = artifacts / "1IAR.cif"
    writer = MMCIFIO()
    writer.set_structure(loaded.structure)
    writer.save(str(cif_path))
    store = Store(artifacts / "data")
    fixtures = {"data_dir": str(store.root), "analyses": {}, "surface": surface}
    for fmt, source in (("pdb", EXAMPLES / "1IAR.pdb"), ("mmcif", cif_path)):
        analysis = store.create_analysis(
            source_kind="upload", source_name=source.name,
            source_format=fmt, source_bytes=source.read_bytes(),
        )
        analysis_id = analysis["id"]
        payload = copy.deepcopy(result)
        payload["analysis_id"] = analysis_id
        payload["source"] = {"name": source.name, "format": fmt}
        payload["metadata"]["source_name"] = source.name
        payload["metadata"]["source_format"] = fmt
        path = store.write_json(store.analysis_dir(analysis_id) / "analysis.json", payload)
        store.update_analysis(analysis_id, status="complete", result_path=path)
        for kind, data in (("core", payload), ("curvature", surface)):
            job = store.create_job(analysis_id, kind)
            job_path = store.write_json(store.analysis_dir(analysis_id) / f"{kind}.json", data)
            store.update_job(job["id"], status="complete", stage="complete", progress=1, result_path=job_path)
        fixtures["analyses"][fmt] = analysis_id
    fixture_path.write_text(json.dumps(fixtures))
    return fixtures


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def exercise_browser(url: str, browser_path: str, artifacts: Path, fixtures: dict):
    from playwright.sync_api import expect, sync_playwright

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=browser_path, headless=True, timeout=30_000,
            args=["--no-sandbox", "--disable-dev-shm-usage", "--enable-unsafe-swiftshader"],
        )
        try:
            page = browser.new_page(viewport={"width": 1440, "height": 1000})
            errors = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(url)

            def open_analysis(fmt):
                analysis_id = fixtures["analyses"][fmt]
                page.locator(f'.history-row[data-id="{analysis_id}"]').click()
                page.wait_for_function(
                    "id => state.analysisId === id && state.surface?.status === 'complete' "
                    "&& Object.keys(state.meshComponents).length === 2", arg=analysis_id,
                    timeout=60_000,
                )
                return analysis_id

            def slider(element_id, value):
                page.locator(f"#{element_id}").evaluate(
                    "(el, value) => { el.value = value; el.dispatchEvent(new Event('input', {bubbles:true})); }",
                    str(value),
                )

            def check_structure_selection():
                assert page.evaluate("""() => state.selectedResidues.length > 0
                    && state.selectedResidues.every(key => {
                      const part = selectionParts(key);
                      return part && state.structureComponents[part.chain].structure
                        .getView(new NGL.Selection(selectionString(key))).atomCount > 0;
                    }) && state.selectionReps.length === state.selectedResidues.length""")

            def check_surface_selection():
                stats = page.evaluate("""() => {
                  const require = (ok, message) => { if (!ok) throw Error(message); };
                  let yellowFaces = 0;
                  let yellowOutsidePadding = 0;
                  for (const [chain, component] of Object.entries(state.meshComponents)) {
                    const mesh = state.surface.meshes[chain];
                    const display = state.meshGeometry[chain];
                    const geometry = component.object.bufferList[0].geometry;
                    const colors = geometry.attributes.color.array;
                    const positions = geometry.attributes.position.array;
                    const normals = geometry.attributes.normal.array;
                    const primitiveIds = geometry.attributes.primitiveId.array;
                    const index = geometry.index.array;
                    const source = state.meshSourceData.get(mesh);
                    const labels = new Set(state.selectedResidues
                      .filter(key => selectionParts(key)?.chain === chain).map(residueLabelFromKey));
                    require(index.length === mesh.index.length, `${chain}: no extra or missing faces`);
                    require(positions.length === display.vertexMap.length * 3, `${chain}: vertex count`);
                    const pair = state.result.pairs.find(pair => pair.id === state.pairId);
                    const partner = chain === pair.chain_a ? pair.chain_b : pair.chain_a;
                    const distances = mesh.contact_distance_dA[partner];
                    for (let vertex = 0; vertex < display.vertexMap.length; vertex++) {
                      const original = display.vertexMap[vertex];
                      require(primitiveIds[vertex] === vertex, `${chain}: displayed picking ID`);
                      for (let axis = 0; axis < 3; axis++) {
                        require(positions[vertex * 3 + axis] === source.position[original * 3 + axis],
                          `${chain}: original position preserved`);
                        require(normals[vertex * 3 + axis] === source.normal[original * 3 + axis],
                          `${chain}: original smooth normal preserved`);
                      }
                    }
                    for (let face = 0; face < mesh.index.length; face += 3) {
                      const selected = [0, 1, 2].some(corner => labels.has(mesh.residue[mesh.index[face + corner]]));
                      if (selected) yellowFaces++;
                      for (let corner = 0; corner < 3; corner++) {
                        const vertex = index[face + corner];
                        const original = mesh.index[face + corner];
                        require(display.vertexMap[vertex] === original, `${chain}: face topology preserved`);
                        const yellow = colors[vertex * 3] === 1 && colors[vertex * 3 + 1] === 1
                          && colors[vertex * 3 + 2] === 0;
                        require(yellow === selected, `${chain}: full yellow face without neighboring color bleed`);
                        require(selected ? vertex >= display.highlightedVertexStart : vertex === original,
                          `${chain}: only selected faces use duplicated vertices`);
                        if (selected && distances?.[original] > state.surfacePadding * 10) yellowOutsidePadding++;
                      }
                    }
                    const isPair = chain === pair.chain_a || chain === pair.chain_b;
                    require(state.meshReps[chain].getParameters().opacity ===
                      (isPair ? state.surfaceOpacity : state.surfaceOpacity * 0.28), `${chain}: highlight opacity`);
                    require(component.visible === (state.visibleChains[chain] !== false), `${chain}: visibility`);
                  }
                  require(!state.stage.compList.some(component => component.name.startsWith('selected-')),
                    'Centroid sphere markers must be absent');
                  return {yellowFaces, yellowOutsidePadding};
                }""")
                assert (stats["yellowFaces"] > 0) == bool(page.evaluate("state.selectedResidues.length")), stats
                return stats

            def check_boundary_mesh():
                # Non-coplanar triangles share vertices across a residue boundary.
                # Explicit expected faces also cover insertion codes and two-residue selection.
                assert page.evaluate("""() => {
                  const mesh = {
                    position: [0,0,0, 1,0,0, 0,1,0, 1,1,1, 2,1,0],
                    index: [0,1,2, 1,3,2, 3,4,2],
                    residue: ['R12', 'G13B', 'G13B', 'S14', 'T15'],
                    h: {'6': [-0.8,-0.3,0,0.3,0.8]},
                    contact_distance_dA: {B: [0,10,20,30,40]},
                  };
                  const saved = {surface: state.surface, result: state.result, geometry: state.meshGeometry,
                    scale: state.surfaceScale, padding: state.surfacePadding};
                  const require = (ok, message) => { if (!ok) throw Error(message); };
                  try {
                    state.surface = {report: {color_limit_Ainv: 1}, meshes: {A: mesh}};
                    state.result = {pairs: [{id: state.pairId, chain_a: 'A', chain_b: 'B'}]};
                    state.surfaceScale = '6';
                    state.surfacePadding = 0;
                    const original = surfaceMeshGeometry(mesh, new Set());
                    state.meshGeometry = {A: original};
                    const originalColors = meshColors('A');
                    for (const [labels, expectedFaces] of [
                      [['R12'], [true, false, false]],
                      [['R12', 'S14'], [true, true, true]],
                      [['G13B'], [true, true, true]],
                      [['T15'], [false, false, true]],
                      [[], [false, false, false]],
                    ]) {
                      const display = surfaceMeshGeometry(mesh, new Set(labels));
                      state.meshGeometry.A = display;
                      const shape = new NGL.Shape('boundary-test');
                      shape.addMesh(display.position, meshColors('A'), display.index, display.normal, 'boundary');
                      try {
                        const geometry = shape.bufferList[0].geometry;
                        const colors = geometry.attributes.color.array;
                        for (let face = 0; face < 3; face++) {
                          for (let corner = 0; corner < 3; corner++) {
                            const vertex = geometry.index.array[face * 3 + corner];
                            const source = mesh.index[face * 3 + corner];
                            require(display.vertexMap[vertex] === source, 'Boundary picking mapping');
                            for (let axis = 0; axis < 3; axis++) {
                              const expected = expectedFaces[face] ? [1,1,0][axis] : originalColors[source * 3 + axis];
                              require(colors[vertex * 3 + axis] === expected, 'Boundary face color');
                              require(geometry.attributes.normal.array[vertex * 3 + axis] === original.normal[source * 3 + axis],
                                'Boundary smooth normal');
                            }
                          }
                        }
                      } finally { shape.dispose(); }
                    }
                    return true;
                  } finally {
                    state.surface = saved.surface;
                    state.result = saved.result;
                    state.meshGeometry = saved.geometry;
                    state.surfaceScale = saved.scale;
                    state.surfacePadding = saved.padding;
                  }
                }""")

            pdb_id = open_analysis("pdb")
            check_boundary_mesh()
            check_surface_selection()
            expect(page.locator("#residue-panel")).to_be_visible()
            page.locator("#residue-table .focus-row").first.click()
            check_structure_selection()
            check_surface_selection()
            page.locator("#contact-tab").click()
            expect(page.locator("#contact-panel")).to_be_visible()
            expect(page.locator("#residue-panel")).to_be_hidden()
            page.locator("#contact-table .focus-row").first.click()
            assert page.evaluate("state.selectedResidues.length") == 2
            check_structure_selection()
            check_surface_selection()
            page.locator("#contact-grid .map-cell.hit").first.click()
            check_structure_selection()
            check_surface_selection()
            page.locator("#contact-grid .map-row-label").first.click()
            assert page.evaluate("state.selectedResidues.length") == 1
            check_structure_selection()
            check_surface_selection()
            page.locator("#curvature-table .focus-row").first.click()
            check_structure_selection()
            check_surface_selection()

            # Exercise hover/click IDs added by the duplicated yellow vertices.
            assert page.evaluate("""() => {
              const key = state.selectedResidues[0];
              const {chain} = selectionParts(key);
              const mesh = state.surface.meshes[chain];
              const display = state.meshGeometry[chain];
              const pid = display.vertexMap.findIndex((original, vertex) =>
                vertex >= display.highlightedVertexStart && mesh.residue[original] === residueLabelFromKey(key));
              const pick = {type: 'mesh', component: state.meshComponents[chain], pid};
              const hit = focusSurfaceHit(pick);
              if (!hit || hit.vertex !== display.vertexMap[pid] || surfaceResidueKey(chain, hit.label) !== key
                || hit.h !== mesh.h[state.surfaceScale][hit.vertex]) return false;
              state.stage.signals.hovered.dispatch(pick);
              if (!document.getElementById('surface-tip').textContent.includes(hit.label)) return false;
              state.stage.signals.clicked.dispatch(pick);
              hideSurfaceTip();
              return state.selectedResidues.length === 1 && state.selectedResidues[0] === key
                && focusSurfaceHit({...pick, pid: -1}) === null
                && focusSurfaceHit({...pick, pid: display.vertexMap.length}) === null;
            }""")
            check_surface_selection()

            page.locator("#contact-cutoff").select_option("6")
            for kind in ("all", "tight", "hb", "sb", "like"):
                page.locator("#contact-type").select_option(kind)
                expected = page.evaluate("""() => contactMapSlice(state.result.pairs
                    .find(pair => pair.id === state.pairId).contact_map).cells.length""")
                assert page.locator("#contact-grid .map-cell.hit").count() == expected
                assert page.locator("#contact-map-table .focus-row").count() == expected
            page.locator("#contact-type").select_option("all")
            assert page.locator("#contact-grid .map-lane.q-0").first.inner_text() in "●▲◆○"
            assert page.locator("#contact-grid .map-gap-line").count() > 0
            page.locator("#contact-sasa").uncheck()
            assert page.locator("#contact-grid .map-sasa").count() == 0
            page.locator("#contact-sasa").check()
            page.locator("#contact-values").uncheck()
            assert not page.locator("#contact-grid .map-cell.hit").first.inner_text().replace("+", "").replace("×", "").replace("H", "").strip()
            page.locator("#contact-values").check()
            page.locator("#surface-resolution").select_option("4")
            check_surface_selection()
            slider("surface-opacity", 0.3)
            check_surface_selection()
            slider("surface-padding", 0)
            assert check_surface_selection()["yellowOutsidePadding"] > 0
            slider("surface-padding", 5)
            check_surface_selection()
            slider("surface-separation", 6)
            check_surface_selection()
            assert page.evaluate("state.surfaceScale === '4' && state.surfaceOpacity === 0.3 && state.surfacePadding === 5 && state.surfaceSeparation === 6")
            assert page.evaluate("Object.values(state.meshComponents).every(c => c.object.bufferList[0].geometry.attributes.position.count > 0)")
            chain_toggle = page.locator("#chain-toggles input").first
            chain_toggle.uncheck()
            assert page.evaluate("!state.structureComponents.A.visible && !state.meshComponents.A.visible")
            check_surface_selection()
            chain_toggle.check()
            check_surface_selection()
            selection = page.evaluate("state.selectedResidues")
            page.evaluate("selectResidues([])")
            check_surface_selection()
            assert page.evaluate("state.selectionReps.length") == 0
            page.locator("#curvature-table .focus-row").nth(1).click()
            assert page.evaluate("state.selectedResidues") != selection
            check_surface_selection()
            page.evaluate("keys => selectResidues(keys)", selection)
            check_surface_selection()

            # The viewer compatibility fallback must reinstall the yellow mesh.
            page.evaluate("""() => {
              const buffer = state.meshComponents.A.object.bufferList[0];
              buffer.setAttributes = () => { throw Error('Exercise color-update fallback'); };
              state.meshColorKey = null;
              renderSurface();
            }""")
            check_surface_selection()
            page.locator("#interfaces-toggle").click()
            expect(page.locator("#interfaces-body")).to_be_hidden()
            page.locator("#interfaces-toggle").click()
            expect(page.locator("#interfaces-body")).to_be_visible()

            page.locator("#pocket-anchor").select_option("67:")
            page.locator("#pocket-form button").click()
            expect(page.locator("#cancel-pocket-highlight")).to_be_enabled(timeout=120_000)
            assert page.evaluate("state.pocketHighlighted && !!state.pocketComponent")
            assert page.evaluate("state.pocketComponent.reprList[0].getParameters().opacity") == 0.2
            check_surface_selection()
            metrics = page.locator("#pocket-result").inner_text()
            page.locator("#surface-panel").screenshot(path=str(artifacts / "surface-and-pocket.png"))
            page.locator("#cancel-pocket-highlight").click()
            assert page.evaluate("!state.pocketHighlighted && !state.pocketComponent && state.pocketHighlightCancelled")
            assert page.locator("#pocket-result").inner_text() == metrics
            slider("surface-padding", 2)
            page.locator("#surface-resolution").select_option("8")
            assert page.evaluate("!state.pocketHighlighted && !state.pocketComponent")
            check_surface_selection()
            page.locator("#show-pocket-highlight").click()
            assert page.evaluate("state.pocketHighlighted && !!state.pocketComponent")
            assert page.evaluate("state.pocketComponent.reprList[0].getParameters().opacity") == 0.2
            chain_toggle.uncheck()
            assert page.evaluate("!state.pocketComponent.visible")
            chain_toggle.check()

            open_analysis("mmcif")
            page.locator("#residue-table .focus-row").first.click()
            check_structure_selection()
            check_surface_selection()
            assert page.evaluate("state.structureFormat") == "mmcif"

            # An older surface response must not overwrite the newly opened analysis.
            held = []
            page.route(f"**/api/analyses/{pdb_id}/curvature", lambda route: held.append(route))
            page.locator(f'.history-row[data-id="{pdb_id}"]').click()
            page.wait_for_function("state.curvatureLoaded && state.surface === null")
            assert held
            open_analysis("mmcif")
            stale = copy.deepcopy(fixtures["surface"])
            stale["report"]["source"] = "stale response"
            held[0].fulfill(status=200, content_type="application/json", body=json.dumps(stale))
            page.wait_for_timeout(250)
            assert page.evaluate("state.surface.report.source") != "stale response"
            page.unroute(f"**/api/analyses/{pdb_id}/curvature")

            page.locator("#pair-workspace").screenshot(path=str(artifacts / "interface-and-contact-map.png"))
            page.set_viewport_size({"width": 390, "height": 900})
            page.locator("#workspace").screenshot(path=str(artifacts / "mobile-workspace.png"))
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
            assert not errors, errors
            print("Browser smoke passed: yellow residue faces and boundaries, smooth normals, linked selections, picking, surface controls, transparent pockets, tabs, map filters, PDB/mmCIF, stale responses, and mobile layout.", flush=True)
        finally:
            browser.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--browser", default="/usr/bin/google-chrome")
    parser.add_argument("--artifacts", type=Path, default=Path("/tmp/interface-app-browser"))
    parser.add_argument("--reuse-fixtures", action="store_true")
    args = parser.parse_args()
    artifacts = args.artifacts.resolve()
    artifacts.mkdir(parents=True, exist_ok=True)
    print("Preparing browser fixtures…", flush=True)
    fixtures = prepare_fixtures(artifacts, args.reuse_fixtures)
    port = free_port()
    url = f"http://127.0.0.1:{port}"
    env = {**os.environ, "INTERFACE_APP_DATA": fixtures["data_dir"]}
    with (artifacts / "server.log").open("w") as log:
        server = subprocess.Popen(
            [sys.executable, str(ROOT / "run_interface_app.py"), "--port", str(port), "--data-dir", fixtures["data_dir"]],
            cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT,
        )
        try:
            deadline = time.monotonic() + 30
            while time.monotonic() < deadline:
                if server.poll() is not None:
                    raise RuntimeError(f"Server exited; see {artifacts / 'server.log'}")
                try:
                    with urllib.request.urlopen(f"{url}/api/health", timeout=1):
                        break
                except OSError:
                    time.sleep(0.1)
            else:
                raise RuntimeError("Server did not become ready")
            exercise_browser(url, args.browser, artifacts, fixtures)
        finally:
            server.terminate()
            try:
                server.wait(timeout=15)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()


if __name__ == "__main__":
    main()
