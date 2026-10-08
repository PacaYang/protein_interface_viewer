"""Exercise the workspace with real coordinates, API jobs, and NGL/WebGL.

Run separately from pytest; requires Playwright and a local Chrome browser.
All fixtures, server history, logs, and screenshots stay in the artifact folder.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
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
        cached = json.loads(fixture_path.read_text())
        if cached.get("fixture_version") == 2:
            return cached

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
    electrostatics = {
        "version": "0.2.0",
        "status": "complete",
        "parameters": {"model": "synthetic browser fixture", "units": "kT/e"},
        "report": {"units": "kT/e", "color_limit_kT_e": 1.0, "min_kT_e": -1.0, "max_kT_e": 1.0},
        "meshes": {
            chain: {"potential_kT_e": [round(((index % 9) - 4) / 4, 3)
                                      for index in range(len(mesh["position"]) // 3)]}
            for chain, mesh in surface["meshes"].items()
        },
    }
    cif_path = artifacts / "1IAR.cif"
    writer = MMCIFIO()
    writer.set_structure(loaded.structure)
    writer.save(str(cif_path))
    store = Store(artifacts / "data")
    fixtures = {"fixture_version": 2, "data_dir": str(store.root), "analyses": {}, "surface": surface}
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
        for kind, data in (("core", payload), ("curvature", surface), ("electrostatics", electrostatics)):
            job = store.create_job(analysis_id, kind)
            if kind == "curvature":
                curvature_id = job["id"]
            if kind == "electrostatics":
                data = {**data, "curvature_job_id": curvature_id}
            job_path = store.write_json(store.analysis_dir(analysis_id) / f"{kind}.json", data)
            store.update_job(job["id"], status="complete", stage="complete", progress=1, result_path=job_path)
        fixtures["analyses"][fmt] = analysis_id
    fixture_path.write_text(json.dumps(fixtures))
    return fixtures


def free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def exercise_browser(url: str, browser_path: str, artifacts: Path, fixtures: dict, *, faces_only: bool = False, picking_only: bool = False):
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
                    const index = geometry.index.array;
                    const pickingGeometry = component.object.bufferList[0].getPickingMesh().geometry;
                    const pickingIds = pickingGeometry.attributes.primitiveId.array;
                    const pickingPositions = pickingGeometry.attributes.position.array;
                    const source = state.meshSourceData.get(mesh);
                    const labels = new Set(state.selectedResidues
                      .filter(key => selectionParts(key)?.chain === chain).map(residueLabelFromKey));
                    require(index.length === mesh.index.length, `${chain}: no extra or missing faces`);
                    require(positions.length === display.vertexMap.length * 3, `${chain}: vertex count`);
                    require(pickingIds.length === mesh.index.length, `${chain}: one picking ID per triangle corner`);
                    const pair = state.result.pairs.find(pair => pair.id === state.pairId);
                    const partner = chain === pair.chain_a ? pair.chain_b : pair.chain_a;
                    const distances = mesh.contact_distance_dA[partner];
                    for (let vertex = 0; vertex < display.vertexMap.length; vertex++) {
                      const original = display.vertexMap[vertex];
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
                        require(pickingIds[face + corner] === face / 3, `${chain}: constant triangle picking ID`);
                        for (let axis = 0; axis < 3; axis++) {
                          require(pickingPositions[(face + corner) * 3 + axis] === positions[vertex * 3 + axis],
                            `${chain}: picking and display triangles aligned`);
                        }
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

            def check_surface_modes():
                page.wait_for_function("state.electrostatics?.status === 'complete' && !document.getElementById('surface-mode-electrostatics').disabled")
                assert page.evaluate("state.surfaceMode") == "convexity"
                page.locator("#residue-table .focus-row").first.click()
                page.mouse.move(1, 1)
                check_surface_selection()
                before = page.evaluate("() => Object.fromEntries(Object.entries(state.meshComponents).map(([chain, component]) => [chain, Array.from(component.object.bufferList[0].geometry.attributes.color.array)]))")
                page.locator("#surface-mode-electrostatics").click()
                page.wait_for_function("state.surfaceMode === 'electrostatics'")
                page.mouse.move(1, 1)
                check_surface_selection()
                after = page.evaluate("() => Object.fromEntries(Object.entries(state.meshComponents).map(([chain, component]) => [chain, Array.from(component.object.bufferList[0].geometry.attributes.color.array)]))")
                assert any(before[chain] != after[chain] for chain in before), "Electrostatic mode did not recolor the interactive mesh"
                expect(page.locator("#surface-color-label")).to_contain_text("kT/e")
                check_interactive_picking()
                page.locator("#surface-faces-tab").click()
                wait_faces()
                page.mouse.move(1, 1)
                check_face_highlights()
                # With no temporary hover, the same selected residues give
                # identical vertex colors in the interactive and both panes.
                page.wait_for_function("state.faceView.hoverKey === null")
                face_colors = page.evaluate("() => Object.fromEntries(state.faceView.panels.map(panel => [panel.chain, Array.from(panel.meshComponent.object.bufferList[0].geometry.attributes.color.array)]))")
                assert face_colors == {chain: after[chain] for chain in face_colors}
                hover_face("a")
                expect(page.locator("#surface-tip")).to_contain_text("kT/e")
                page.mouse.move(1, 1)
                page.wait_for_function("state.faceView.hoverKey === null")
                page.locator("#surface-panel").screenshot(path=str(artifacts / "electrostatics-binding-faces.png"))
                page.locator("#surface-mode-convexity").click()
                page.wait_for_function("state.surfaceMode === 'convexity'")
                page.mouse.move(1, 1)
                expect(page.locator("#surface-color-label")).to_contain_text("H (Å")
                face_colors = page.evaluate("() => Object.fromEntries(state.faceView.panels.map(panel => [panel.chain, Array.from(panel.meshComponent.object.bufferList[0].geometry.attributes.color.array)]))")
                assert face_colors == {chain: before[chain] for chain in face_colors}
                page.locator("#surface-interactive-tab").click()
                check_surface_selection()
                page.evaluate("selectResidues([])")
                print("Surface mode checks passed: electrostatic and convexity colors match in all views; selections and picking are preserved.", flush=True)

            def check_electrostatics_loading(pdb_id):
                detail = page.evaluate("state.detail")
                potential = page.evaluate("state.electrostatics")
                detail_url = f"**/api/analyses/{pdb_id}"
                potential_url = f"**/api/analyses/{pdb_id}/electrostatics"
                phase = {"ready": False, "requests": 0}

                # A complete interface analysis must keep polling an active
                # optional calculation until the electrostatic map arrives.
                def pending_detail(route):
                    phase["requests"] += 1
                    response = copy.deepcopy(detail)
                    if not phase["ready"]:
                        job = next(job for job in response["jobs"] if job["kind"] == "electrostatics")
                        job.update(status="running", stage="Fixture APBS solve",
                                   progress=min(phase["requests"], 9) / 10, result_path=None)
                    route.fulfill(status=200, content_type="application/json", body=json.dumps(response))

                page.route(detail_url, pending_detail)
                open_analysis("pdb")
                page.wait_for_function("state.detail.status === 'complete' && state.timer !== null "
                                       "&& state.detail.jobs.some(job => job.kind === 'electrostatics' "
                                       "&& job.status === 'running' && job.progress >= .2)")
                expect(page.locator("#surface-mode-electrostatics")).to_be_disabled()
                expect(page.locator("#electrostatics-status")).to_contain_text("Fixture APBS solve")
                assert page.evaluate("state.surfaceMode") == "convexity"
                phase["ready"] = True
                page.wait_for_function("electrostaticsReady() && state.timer === null")
                page.unroute(detail_url)

                # A failed optional artifact returns both viewers to convexity
                # and offers a local reload while the interface stays usable.
                open_faces()
                page.mouse.move(1, 1)
                page.wait_for_function("state.faceView.hoverKey === null")
                page.locator("#surface-mode-electrostatics").click()
                unavailable = {"status": "unavailable", "reason": "APBS fixture unavailable",
                               "report": None, "meshes": {}}
                page.route(potential_url, lambda route: route.fulfill(
                    status=200, content_type="application/json", body=json.dumps(unavailable)))
                page.evaluate("() => { state.electrostaticsLoaded = false; void loadElectrostatics(); }")
                expect(page.locator("#electrostatics-status")).to_contain_text("APBS fixture unavailable")
                expect(page.locator("#surface-mode-electrostatics")).to_be_disabled()
                expect(page.locator("#surface-mode-convexity")).to_have_attribute("aria-pressed", "true")
                assert page.evaluate("""() => state.detail.status === 'complete'
                  && Object.entries(state.meshComponents).every(([chain, component]) => {
                    const actual = component.object.bufferList[0].geometry.attributes.color.array;
                    const expected = meshColors(chain);
                    return actual.every((value,index) => value === expected[index]);
                  }) && state.faceView.panels.every(panel => {
                    const actual = panel.meshComponent.object.bufferList[0].geometry.attributes.color.array;
                    const expected = meshColors(panel.chain, panel.meshGeometry);
                    return actual.every((value,index) => value === expected[index]);
                  })""")
                expect(page.locator("#electrostatics-retry")).to_have_text("Reload")
                page.unroute(potential_url)
                page.locator("#electrostatics-retry").click()
                page.wait_for_function("electrostaticsReady()")

                # A response for a previous analysis cannot overwrite the map
                # belonging to the analysis opened while it was in flight.
                held = []
                page.route(potential_url, lambda route: held.append(route))
                with page.expect_request(f"{url}/api/analyses/{pdb_id}/electrostatics"):
                    page.evaluate("() => { state.electrostaticsLoaded = false; void loadElectrostatics(); }")
                page.wait_for_timeout(100)
                assert len(held) == 1
                open_analysis("mmcif")
                page.wait_for_function("electrostaticsReady()")
                stale = {**potential, "version": "stale electrostatics response"}
                held[0].fulfill(status=200, content_type="application/json", body=json.dumps(stale))
                page.wait_for_timeout(100)
                assert page.evaluate("state.electrostatics.version") != stale["version"]
                page.unroute(potential_url)
                open_analysis("pdb")
                page.wait_for_function("electrostaticsReady()")
                print("Electrostatic loading checks passed: polling after core completion, unavailable fallback, reload, and stale responses.", flush=True)

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

            def check_surface_picking():
                # Only three widely separated vertex IDs form the visible
                # triangle. An interpolated ID points to a remote unused vertex.
                # Two more triangles cross the picking color's 255/256 boundary.
                stats = page.evaluate("""async () => {
                  const require = (ok, message) => { if (!ok) throw Error(message); };
                  const host = document.createElement('div');
                  Object.assign(host.style, {position:'fixed', left:'0', top:'0', width:'480px', height:'400px', zIndex:100});
                  document.body.appendChild(host);
                  const stage = new NGL.Stage(host, {cameraType:'orthographic', sampleLevel:0,
                    clipMode:'camera', clipNear:0.1, clipFar:2000, tooltip:false});
                  try {
                    const positions = new Float32Array(1024 * 3).fill(200);
                    positions.set([-6,-3,0], 0);
                    positions.set([6,-3,0], 511 * 3);
                    positions.set([0,6,0], 1023 * 3);
                    positions.set([-16,-2,0, -10,-2,0, -13,4,0], 255 * 3);
                    positions.set([10,-2,0, 16,-2,0, 13,4,0], 767 * 3);
                    const index = new Uint32Array(257 * 3).fill(1);
                    index.set([0,511,1023], 0);
                    index.set([255,256,257], 255 * 3);
                    index.set([767,768,769], 256 * 3);
                    const mesh = {position:positions, index, residue:Array(1024).fill('G999')};
                    mesh.residue[0] = 'R1'; mesh.residue[511] = 'G2'; mesh.residue[1023] = 'L3';
                    let samples = 0;
                    for (const selected of [new Set(), new Set(['R1', 'L3'])]) {
                      const geometry = surfaceMeshGeometry(mesh, selected);
                      const shape = new NGL.Shape('picking-regression');
                      shape.addMesh(geometry.position, new Float32Array(geometry.position.length).fill(.5),
                        geometry.index, geometry.normal);
                      installSurfacePicking(shape, geometry);
                      const component = stage.addComponentFromObject(shape);
                      component.addRepresentation('buffer', {opacity:.3, side:'double'});
                      await new Promise(resolve => stage.tasks.onZeroOnce(resolve));
                      let disposed = false;
                      shape.bufferList[0].getPickingMesh().geometry.addEventListener('dispose', () => { disposed = true; });
                      for (const mode of ['front', 'back', 'transformed', 'perspective', 'zoom']) {
                        stage.setParameters({cameraType:mode === 'perspective' ? 'perspective' : 'orthographic'});
                        component.setTransform(new NGL.Matrix4().makeRotationY(mode === 'back' ? Math.PI :
                          mode === 'transformed' ? .4 : 0).setPosition(new NGL.Vector3(2,-1,0)));
                        stage.viewerControls.rotate(new NGL.Quaternion().setFromAxisAngle(new NGL.Vector3(0,0,1), .2));
                        stage.viewerControls.distance(mode === 'zoom' ? -38 : -48);
                        stage.viewer.render();
                        for (const triangle of [0,255,256]) {
                          const source = Array.from(index.subarray(triangle * 3, triangle * 3 + 3));
                          for (let target = 0; target < 3; target++) {
                            const point = new NGL.Vector3();
                            for (let corner = 0; corner < 3; corner++) {
                              point.addScaledVector(new NGL.Vector3().fromArray(positions, source[corner] * 3),
                                corner === target ? .8 : .1);
                            }
                            point.applyMatrix4(component.matrix);
                            const cursor = stage.viewerControls.getPositionOnCanvas(point);
                            const x = Math.round(cursor.x), y = Math.round(cursor.y);
                            stage.mouseObserver.canvasPosition.set(x, y);
                            const pick = stage.pickingControls.pick(x, y);
                            require(pick?.type === 'mesh' && pick.pid === triangle,
                              `${mode}: triangle ${triangle} under cursor, got ${pick?.type}/${pick?.pid}`);
                            require(surfaceHitVertex(pick, mesh, geometry) === source[target],
                              `${mode}: residue belongs to hovered corner of triangle ${triangle}`);
                            samples++;
                          }
                        }
                        // At the center of triangle 0, the top corner is closest.
                        const center = stage.viewerControls.getPositionOnCanvas(new NGL.Vector3().applyMatrix4(component.matrix));
                        stage.mouseObserver.canvasPosition.copy(center);
                        const pick = stage.pickingControls.pick(center.x, center.y);
                        require(pick?.pid === 0 && surfaceHitVertex(pick, mesh, geometry) === 1023,
                          `${mode}: interior picks stay on the visible triangle`);
                        require(surfaceHitVertex({...pick, pid:-1}, mesh, geometry) === null, 'Negative triangle ID');
                        require(surfaceHitVertex({...pick, pid:257}, mesh, geometry) === null, 'Out of bounds triangle ID');
                        require(surfaceHitVertex({...pick, canvasPosition:null}, mesh, geometry) === null, 'Missing cursor');
                        samples++;
                      }
                      stage.removeComponent(component);
                      require(disposed, 'Picking geometry disposed with the displayed component');
                    }
                    return {samples};
                  } finally { disposeNglStage(stage); host.remove(); }
                }""")
                print(f"Surface picking regression: {stats}", flush=True)
                assert stats["samples"] == 100, stats

            def check_interactive_picking():
                selection = page.evaluate("state.selectedResidues")
                page.locator("#surface-view").scroll_into_view_if_needed()
                for _ in range(2):
                    page.wait_for_function("state.stage.tasks.count === 0")
                    hit = page.evaluate("""() => {
                      const stage = state.stage;
                      const offsets = [0, .08, -.08, .16, -.16, .24, -.24, .32, -.32];
                      const edgeDistance = (p, a, b) => {
                        const dx = b.x-a.x, dy = b.y-a.y;
                        const t = Math.max(0, Math.min(1, ((p.x-a.x)*dx + (p.y-a.y)*dy) / (dx*dx+dy*dy || 1)));
                        return Math.hypot(p.x-a.x-t*dx, p.y-a.y-t*dy);
                      };
                      for (const dy of offsets) for (const dx of offsets) {
                        const x = Math.round(stage.viewer.width * (.5 + dx));
                        const y = Math.round(stage.viewer.height * (.5 + dy));
                        stage.mouseObserver.canvasPosition.set(x, y);
                        const pick = stage.pickingControls.pick(x, y);
                        const hit = focusSurfaceHit(pick);
                        if (!hit) continue;
                        const mesh = state.surface.meshes[hit.chain];
                        const vertices = Array.from(mesh.index.slice(pick.pid * 3, pick.pid * 3 + 3));
                        if (!vertices.includes(hit.vertex)) throw Error('Residue outside the picked triangle');
                        const triangle = vertices.map(vertex => stage.viewerControls.getPositionOnCanvas(
                          new NGL.Vector3().fromArray(mesh.position, vertex * 3).applyMatrix4(pick.component.matrix)));
                        const sides = triangle.map((a, i) => {
                          const b = triangle[(i+1)%3];
                          return (b.x-a.x)*(y-a.y) - (b.y-a.y)*(x-a.x);
                        });
                        const inside = sides.every(side => side >= 0) || sides.every(side => side <= 0);
                        const distance = inside ? 0 : Math.min(...triangle.map((a,i) =>
                          edgeDistance({x,y}, a, triangle[(i+1)%3])));
                        // NGL samples a 5x5 pixel region around the cursor.
                        if (distance > 3) throw Error(`Picked triangle is ${distance} px from the cursor`);
                        return {x, y, key:surfaceResidueKey(hit.chain, hit.label), label:hit.label};
                      }
                      throw Error('No real interactive surface pick');
                    }""")
                    box = page.locator("#surface-view canvas").bounding_box()
                    x, y = box["x"] + hit["x"], box["y"] + box["height"] - hit["y"]
                    page.mouse.move(x, y)
                    expect(page.locator("#surface-tip")).to_be_visible()
                    expect(page.locator("#surface-tip")).to_contain_text(hit["label"])
                    expect(page.locator("#surface-tip")).to_contain_text(
                        "kT/e" if page.evaluate("state.surfaceMode") == "electrostatics" else "local H")
                    page.mouse.click(x, y)
                    page.wait_for_function("key => state.selectedResidues.length === 1 && state.selectedResidues[0] === key",
                                           arg=hit["key"])
                    check_surface_selection()
                    page.mouse.move(1, 1)
                page.evaluate("keys => selectResidues(keys)", selection)
                print("Interactive surface checks passed: real cursor hover/click and highlighted triangles.", flush=True)

            def wait_faces(surfaces=True):
                condition = "panel.meshComponent && " if surfaces else ""
                page.wait_for_function("state.faceView?.ready && state.faceView.panels.length === 2 "
                                       f"&& state.faceView.panels.every(panel => {condition}panel.stage.tasks.count === 0)",
                                       timeout=30_000)
                expect(page.locator("#face-view-grid")).to_be_visible()
                page.locator("#face-view-a").scroll_into_view_if_needed()

            def open_faces():
                page.locator("#surface-faces-tab").click()
                wait_faces()
                expect(page.locator("#surface-separation-control")).to_be_hidden()
                expect(page.locator("#chain-toggles")).to_be_hidden()
                page.locator("#face-view-a").scroll_into_view_if_needed()

            def check_face_frames():
                assert page.evaluate("""() => {
                  const require = (ok, message) => { if (!ok) throw Error(message); };
                  const face = state.faceView;
                  for (const panel of face.panels) {
                    const frame = face.geometry.chains[panel.chain];
                    const n = new NGL.Vector3(...frame.normal).transformDirection(panel.transform);
                    const up = new NGL.Vector3(...frame.up).transformDirection(panel.transform);
                    const camera = panel.stage.viewer.camera;
                    require(camera.isOrthographicCamera, 'Orthographic binding view');
                    require(n.distanceTo(new NGL.Vector3(0,0,-1)) < 1e-6, 'Normal faces NGL camera');
                    require(up.distanceTo(new NGL.Vector3(0,1,0)) < 1e-6, 'Shared up axis');
                    require(Math.abs(panel.transform.determinant() - 1) < 1e-6, 'Proper rotation without reflection');
                    require(panel.component.matrix.equals(panel.meshComponent.matrix), 'Atoms and mesh aligned');
                    require(panel.component.visible && panel.meshComponent.visible, 'Both partners visible');
                  }
                  const scales = face.panels.map(panel => panel.stage.viewerControls.getCanvasScaleFactor());
                  require(Math.abs(scales[0] - scales[1]) < 1e-6, 'Same initial angstroms per pixel');
                  return true;
                }""")
                check_face_rotation_alignment()

            def face_snapshot():
                return page.evaluate("""() => state.faceView.panels.map(panel => ({
                  orientation: panel.stage.viewerControls.getOrientation().elements,
                  rotation: panel.stage.viewerControls.rotation.toArray(),
                  position: panel.stage.viewerControls.position.toArray(),
                  transform: panel.component.matrix.elements,
                  meshTransform: panel.meshComponent?.matrix.elements || null,
                  zoom: panel.stage.viewer.camera.position.z,
                }))""")

            def check_face_rotation_alignment():
                # Corresponding binding coordinates are horizontally mirrored.
                # After rotation they must retain equal height and depth, with
                # opposite horizontal positions. This checks the resulting
                # geometry independently of the quaternion mapping helper.
                assert page.evaluate("""() => {
                  const [a,b] = state.faceView.panels.map(panel => panel.stage.viewerControls.rotation);
                  for (const [x,y,z] of [[1,0,0], [0,1,0], [0,0,1], [2,3,-4]]) {
                    const first = new NGL.Vector3(x,y,z).applyQuaternion(a);
                    const second = new NGL.Vector3(-x,y,z).applyQuaternion(b);
                    first.x = -first.x;
                    if (first.distanceTo(second) > 1e-6) throw Error('Corresponding sites lose mirrored alignment');
                  }
                  return true;
                }""")

            def check_mirrored_rotation_axes():
                saved = face_snapshot()[0]["rotation"]
                for side in ("a", "b"):
                    for axis in ([1,0,0], [0,1,0], [0,0,1], [1,2,3]):
                        page.evaluate("""({side,axis}) => {
                          const rotation = new NGL.Quaternion().setFromAxisAngle(new NGL.Vector3(...axis).normalize(), .37);
                          state.faceView.panels.find(panel => panel.side === side).stage.viewerControls.rotate(rotation);
                        }""", {"side": side, "axis": axis})
                        check_face_rotation_alignment()
                # Real mouse gestures must produce opposite visible Z rolls.
                def screen_angles():
                    return page.evaluate("""() => state.faceView.panels.map(panel => {
                      panel.stage.viewer.render();
                      const origin = panel.stage.viewerControls.getPositionOnCanvas(new NGL.Vector3());
                      const point = panel.stage.viewerControls.getPositionOnCanvas(new NGL.Vector3(1,0,0));
                      return Math.atan2(point.y-origin.y, point.x-origin.x);
                    })""")

                for side in ("a", "b"):
                    for dx, dy in ((30,20), (-30,-20)):
                        page.evaluate("state.faceView.panels[0].stage.viewerControls.rotate(new NGL.Quaternion())")
                        before = screen_angles()
                        rotate_faces(side, button="right", modifier="Control", dx=dx, dy=dy)
                        after = screen_angles()
                        angles = [math.atan2(math.sin(new-old), math.cos(new-old)) for old,new in zip(before,after)]
                        assert abs(angles[0]) > 1e-4, "The Z gesture must rotate the visible structure"
                        assert abs(angles[0]+angles[1]) < 1e-6, "Clockwise roll must produce equal counterclockwise partner roll"
                page.evaluate("rotation => state.faceView.panels[0].stage.viewerControls.rotate(new NGL.Quaternion(...rotation))", saved)
                print("Mirrored rotation passed: X/Y/Z and mixed axes from either pane, plus clockwise/counterclockwise mouse rolls.", flush=True)

            def check_face_highlights():
                page.wait_for_function("state.faceView.panels.every(panel => panel.stage.tasks.count === 0)")
                assert page.evaluate("""() => {
                  const require = (ok, message) => { if (!ok) throw Error(message); };
                  const face = state.faceView;
                  for (const panel of face.panels) {
                    const expected = [...new Set([...state.selectedResidues, ...face.hoverKeys])]
                      .filter(key => selectionParts(key)?.chain === panel.chain).sort();
                    const actual = new Set();
                    const view = panel.component.structure.getView(panel.highlightRep.repr.selection);
                    try {
                      view.eachAtom(atom => actual.add(`${atom.chainname}:${atom.resno}${atom.inscode || ''}`));
                    } finally { view.dispose(); }
                    require(JSON.stringify([...actual].sort()) === JSON.stringify(expected),
                      `${panel.chain}: selected and hovered atoms are highlighted`);
                    if (!panel.meshComponent) continue;
                    const labels = new Set(expected.map(residueLabelFromKey));
                    const geometry = panel.meshComponent.object.bufferList[0].geometry;
                    const colors = geometry.attributes.color.array;
                    const index = geometry.index.array;
                    for (let face = 0; face < index.length; face += 3) {
                      const yellow = [0,1,2].some(corner => labels.has(panel.mesh.residue[panel.mesh.index[face+corner]]));
                      for (let corner = 0; corner < 3; corner++) {
                        const offset = index[face+corner] * 3;
                        require((colors[offset] === 1 && colors[offset+1] === 1 && colors[offset+2] === 0) === yellow,
                          `${panel.chain}: each selected/hovered residue has its full yellow surface patch`);
                      }
                    }
                  }
                  return true;
                }""")

            def rotate_faces(side, *, button="left", modifier=None, dx=30, dy=20):
                page.locator(f"#face-view-{side}").scroll_into_view_if_needed()
                before = face_snapshot()
                box = page.locator(f"#face-view-{side} canvas").bounding_box()
                x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
                if modifier:
                    page.keyboard.down(modifier)
                page.mouse.move(x, y)
                page.mouse.down(button=button)
                page.mouse.move(x + dx, y + dy, steps=4)
                page.mouse.up(button=button)
                if modifier:
                    page.keyboard.up(modifier)
                source = 0 if side == "a" else 1
                page.wait_for_function("({previous, source}) => state.faceView.panels[source].stage.viewerControls.rotation.toArray() "
                                       ".some((value, index) => Math.abs(value-previous[index]) > 1e-6)",
                                       arg={"previous": before[source]["rotation"], "source": source})
                after = face_snapshot()
                if page.evaluate("state.faceRotationCoupled"):
                    check_face_rotation_alignment()
                else:
                    assert after[1-source] == before[1-source], "Uncoupled rotation must leave the other pane unchanged"
                for old, new in zip(before, after):
                    for key in ("zoom", "position", "transform", "meshTransform"):
                        assert old[key] == new[key], f"Rotation must preserve {key}"
                return after

            def check_face_rotation_toggle():
                toggle = page.locator("#face-rotation-coupled")
                expect(toggle).to_be_checked()
                original = face_snapshot()
                selection = page.evaluate("state.selectedResidues")
                toggle.uncheck()
                expect(page.locator("#face-view-status")).to_contain_text("Independent rotation")
                expect(page.locator("#surface-view-help")).to_contain_text("independently")
                assert face_snapshot() == original, "Disabling coupling preserves both orientations"
                rotate_faces("a")
                rotate_faces("b")
                hover_face("a")
                page.mouse.move(1, 1)
                page.wait_for_function("state.faceView.hoverKey === null")
                check_face_highlights()
                held = face_snapshot()
                page.locator("#surface-interactive-tab").click()
                open_faces()
                expect(toggle).not_to_be_checked()
                assert face_snapshot() == held, "Tab switches preserve independent orientations"
                page.evaluate("cleanupFaceView(); void ensureFaceView()")
                wait_faces()
                expect(toggle).not_to_be_checked()
                rotate_faces("b")
                # Re-coupling keeps the last rotated pane as reference, even
                # if the other pane was zoomed most recently.
                for reference in ("a", "b"):
                    rotate_faces(reference)
                    other = "b" if reference == "a" else "a"
                    page.locator(f'[data-face-zoom="{other}"][data-zoom="in"]').click()
                    before = face_snapshot()
                    toggle.check()
                    after = face_snapshot()
                    index = 0 if reference == "a" else 1
                    assert after[index]["rotation"] == before[index]["rotation"], "Re-coupling preserves the last rotated pane"
                    for old, new in zip(before, after):
                        for key in ("zoom", "position", "transform", "meshTransform"):
                            assert old[key] == new[key], f"The toggle must preserve {key}"
                    check_face_rotation_alignment()
                    expect(page.locator("#face-view-status")).to_contain_text("Mirrored rotation")
                    rotate_faces(other)
                    toggle.uncheck()
                toggle.check()
                page.evaluate("""saved => {
                  state.faceView.panels[0].stage.viewerControls.rotate(new NGL.Quaternion(...saved[0].rotation));
                  state.faceView.panels.forEach((panel, index) => panel.stage.viewerControls.distance(saved[index].zoom));
                }""", original)
                page.mouse.move(1, 1)
                page.wait_for_function("state.faceView.hoverKey === null")
                assert page.evaluate("state.selectedResidues") == selection
                check_face_highlights()
                check_face_rotation_alignment()
                print("Rotation toggle passed: independent drags in both panes, linked hover, tab/rebuild persistence, re-coupling from either pane, and independent zoom.", flush=True)

            def check_face_table_selection():
                selection = page.evaluate("state.selectedResidues")
                active_tab = page.evaluate("state.activeTab")
                camera = face_snapshot()
                for tab, selector in (("residue", "#residue-table"), ("contact", "#contact-table"),
                                      ("contact", "#contact-map-table"), (None, "#curvature-table")):
                    if tab:
                        page.locator(f"#{tab}-tab").click()
                    row = page.locator(f"{selector} .focus-row").first
                    expected = (row.get_attribute("data-keys") or row.get_attribute("data-key")).split("|")
                    row.click()
                    page.mouse.move(1, 1)
                    page.wait_for_function("state.faceView.hoverKey === null")
                    assert page.evaluate("state.selectedResidues") == expected
                    check_face_highlights()
                    assert face_snapshot() == camera, "Table selections must preserve rotation and zoom"
                selected = page.evaluate("state.selectedResidues")
                page.locator("#face-view-a").scroll_into_view_if_needed()
                hover_face("a")
                page.mouse.move(1, 1)
                page.wait_for_function("state.faceView.hoverKey === null")
                assert page.evaluate("state.selectedResidues") == selected
                check_face_highlights()
                # Selections made while the tab is hidden also appear on return.
                page.locator("#surface-interactive-tab").click()
                page.locator("#highlight-residues").fill("A:63-67, B:12")
                page.locator('#residue-highlight-form button[type="submit"]').click()
                assert page.evaluate("state.selectedResidues.length") == 6
                open_faces()
                check_face_highlights()
                assert face_snapshot() == camera
                page.locator("#surface-faces-panel").screenshot(path=str(artifacts / "binding-faces-table-selection.png"))
                page.locator("#surface-interactive-tab").click()
                page.locator("#highlight-clear").click()
                open_faces()
                check_face_highlights()
                page.evaluate("keys => selectResidues(keys)", selection)
                page.locator(f"#{active_tab}-tab").click()
                check_face_highlights()
                page.locator("#face-view-a").scroll_into_view_if_needed()
                print("Binding-face selections passed: residue/contact/map/curvature tables, multiple residues, hover persistence, hidden tab, and clear.", flush=True)

            def find_face_pick(side, surfaces=True):
                return page.evaluate("""({side, surfaces}) => {
                  const panel = state.faceView.panels.find(item => item.side === side);
                  const offsets = [0, .08, -.08, .16, -.16, .24, -.24, .32, -.32];
                  const candidates = [];
                  for (const dy of offsets) for (const dx of offsets) {
                    candidates.push({x:Math.round(panel.stage.viewer.width * (.5 + dx)),
                      y:Math.round(panel.stage.viewer.height * (.5 + dy))});
                  }
                  if (!surfaces) {
                    // A rotated cartoon can fall between a coarse screen grid.
                    // Sample actual backbone positions as well as grid points.
                    panel.component.structure.eachAtom(atom => {
                      if (atom.atomname !== 'CA') return;
                      const point = atom.positionToVector3(new NGL.Vector3()).applyMatrix4(panel.component.matrix);
                      const screen = panel.stage.viewerControls.getPositionOnCanvas(point);
                      candidates.push({x:Math.round(screen.x), y:Math.round(screen.y)});
                    });
                  }
                  let first = null;
                  const misses = {};
                  const examples = [];
                  for (const {x,y} of candidates) {
                    if (x < 3 || y < 3 || x > panel.stage.viewer.width-3 || y > panel.stage.viewer.height-3) continue;
                    panel.stage.mouseObserver.canvasPosition.set(x, y);
                    const pick = panel.stage.pickingControls.pick(x, y);
                    const key = facePickedResidue(panel, pick);
                    if (!key || (surfaces && pick.type !== 'mesh')) {
                      const type = pick?.type || 'empty';
                      misses[type] = (misses[type] || 0) + 1;
                      if (pick && examples.length < 4) examples.push({x,y,type,pid:pick.pid,
                        atom:pick.atom && {chain:pick.atom.chainname,number:pick.atom.resno,name:pick.atom.resname},
                        component:pick.component?.name});
                      continue;
                    }
                    if (pick.type === 'mesh') {
                      const vertices = panel.mesh.index.slice(pick.pid * 3, pick.pid * 3 + 3);
                      if (!vertices.some(vertex => surfaceResidueKey(panel.chain, panel.mesh.residue[vertex]) === key)) {
                        throw Error('Picked residue does not belong to the visible triangle');
                      }
                    }
                    if (!first) first = {x, y, key};
                    const match = state.faceView.geometry.nearest[panel.chain][key];
                    if (match.distance_A < 5) return {x, y, key};
                  }
                  if (first) return first;
                  throw Error('No real binding-face pick: ' + JSON.stringify({side,surfaces,misses,examples,
                    atoms:panel.component.structure.atomCount, visible:panel.component.visible,
                    bounds:panel.host.getBoundingClientRect().toJSON(),
                    rotation:panel.stage.viewerControls.rotation.toArray(),
                    camera:{near:panel.stage.viewer.camera.near,far:panel.stage.viewer.camera.far,z:panel.stage.viewer.camera.position.z},
                    representations:panel.component.reprList.map(rep => ({type:rep.repr.type,visible:rep.visible,parameters:rep.getParameters()}))}));
                }""", {"side": side, "surfaces": surfaces})

            def hover_face(side, surfaces=True):
                page.wait_for_function("state.faceView.panels.every(panel => panel.stage.tasks.count === 0)")
                try:
                    hit = find_face_pick(side, surfaces)
                except Exception:
                    page.locator("#surface-faces-panel").screenshot(path=str(artifacts / "binding-faces-picking-failure.png"))
                    raise
                box = page.locator(f"#face-view-{side}").bounding_box()
                page.mouse.move(box["x"] + hit["x"], box["y"] + box["height"] - hit["y"])
                try:
                    page.wait_for_function("key => state.faceView.hoverKey === key", arg=hit["key"], timeout=15_000)
                except Exception:
                    details = page.evaluate("""side => {
                      const face = state.faceView;
                      const panel = face.panels.find(panel => panel.side === side);
                      const mouse = panel.stage.mouseObserver;
                      const pick = panel.stage.pickingControls.pick(mouse.canvasPosition.x, mouse.canvasPosition.y);
                      return {hover:face.hoverKey, pending:face.pendingHover, hoverInside:panel.hoverInside,
                        canvas:mouse.canvasPosition.toArray(), mouse:mouse.position.toArray(),
                        overElement:mouse.overElement, moving:mouse.moving, hovering:mouse.hovering,
                        tasks:panel.stage.tasks.count, pid:pick?.pid, type:pick?.type,
                        actual:facePickedResidue(panel,pick), bounds:panel.host.getBoundingClientRect().toJSON()};
                    }""", side)
                    print(f"Binding-face hover mismatch: expected {hit}; actual {details}", flush=True)
                    page.screenshot(path=str(artifacts / "binding-faces-hover-failure.png"))
                    raise
                page.wait_for_function("state.faceView.panels.every(panel => panel.stage.tasks.count === 0)")
                assert page.evaluate("""() => {
                  const face = state.faceView;
                  if (face.hoverKeys.length !== 2) return false;
                  for (const panel of face.panels) {
                    const key = face.hoverKeys.find(key => selectionParts(key).chain === panel.chain);
                    if (!panel.meshComponent) continue;
                    const colors = panel.meshComponent.object.bufferList[0].geometry.attributes.color.array;
                    if (!colors.some((_, index) => index % 3 === 0 && colors[index] === 1
                      && colors[index+1] === 1 && colors[index+2] === 0)) return false;
                    const geometry = panel.meshGeometry;
                    const vertex = geometry.vertexMap.findIndex((original, vertex) => vertex >= geometry.highlightedVertexStart
                      && panel.mesh.residue[original] === residueLabelFromKey(key));
                    const pid = Math.floor(geometry.index.indexOf(vertex) / 3);
                    const point = new NGL.Vector3().fromArray(geometry.position, vertex * 3).applyMatrix4(panel.meshComponent.matrix);
                    const canvasPosition = panel.stage.viewerControls.getPositionOnCanvas(point);
                    if (facePickedResidue(panel, {type:'mesh', component: panel.meshComponent, pid, canvasPosition}) !== key) return false;
                  }
                  return true;
                }""")
                check_face_highlights()
                expect(page.locator("#face-hover-readout")).to_contain_text("minimum heavy-atom distance")
                return hit

            def check_face_interaction():
                original_selection = page.evaluate("state.selectedResidues")
                original_camera = page.evaluate("state.stage.viewerControls.getOrientation().elements")
                open_faces()
                check_face_frames()
                check_mirrored_rotation_axes()
                check_face_rotation_toggle()
                # Rebuilding the tab beneath a stationary pointer must not
                # require leaving the pane before hover starts working again.
                hover_face("b")
                page.evaluate("cleanupFaceView(); void ensureFaceView()")
                wait_faces()
                hover_face("b")
                # A completed hover pick during a mesh rebuild must recover
                # when the new surface arrives, even without a mouse move.
                hovered = page.evaluate("state.faceView.hoverKey")
                page.evaluate("""() => {
                  const face = state.faceView;
                  applyFaceHover(face, null);
                  const panel = face.panels.find(panel => panel.side === 'b');
                  panel.stage.mouseObserver.moving = false;
                  panel.stage.mouseObserver.hovering = true;
                }""")
                page.wait_for_function("key => state.faceView.hoverKey === key", arg=hovered, timeout=15_000)
                check_face_highlights()
                hover_face("a")
                assert page.evaluate("state.selectedResidues") == original_selection
                page.locator("#surface-faces-panel").screenshot(path=str(artifacts / "binding-faces-linked-hover.png"))
                hover_face("b")
                assert page.evaluate("state.selectedResidues") == original_selection
                # Blank-space picking clears hover while table selections persist.
                box = page.locator("#face-view-b").bounding_box()
                page.mouse.move(box["x"] + 6, box["y"] + 6)
                page.wait_for_function("state.faceView.hoverKey === null")
                hover_face("a")
                page.mouse.move(1, 1)
                page.wait_for_function("state.faceView.hoverKey === null")
                check_face_highlights()
                rotate_faces("a")
                hover_face("a")
                rotate_faces("b")
                hover_face("b")
                rotate_faces("b", button="right", modifier="Control")
                check_face_table_selection()
                locked = face_snapshot()
                box = page.locator("#face-view-a").bounding_box()
                x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
                for button in ("left", "middle", "right"):
                    for modifier in (None, "Shift", "Control", "Alt"):
                        if (button, modifier) in (("left", None), ("right", "Control")):
                            continue
                        if modifier:
                            page.keyboard.down(modifier)
                        page.mouse.move(x, y)
                        page.mouse.down(button=button)
                        page.mouse.move(x + 25, y + 18, steps=4)
                        page.mouse.up(button=button)
                        if modifier:
                            page.keyboard.up(modifier)
                page.mouse.click(x, y)
                page.mouse.dblclick(x, y)
                page.locator("#face-view-a canvas").focus()
                for key in ("i", "k", "r", "ArrowLeft", "ArrowRight"):
                    page.keyboard.press(key)
                page.wait_for_timeout(350)
                assert face_snapshot() == locked, "Movement, recenter, and keyboard actions must remain locked"
                page.locator('[data-face-zoom="a"][data-zoom="in"]').click()
                zoomed = face_snapshot()
                assert abs(zoomed[0]["zoom"]) < abs(locked[0]["zoom"])
                assert zoomed[1] == locked[1], "Partner zooms independently"
                page.mouse.move(x, y)
                page.mouse.wheel(0, 120)
                page.wait_for_function("z => state.faceView.panels[0].stage.viewer.camera.position.z !== z", arg=zoomed[0]["zoom"])
                zoomed = face_snapshot()
                slider("surface-opacity", 0.3)
                slider("surface-padding", 4)
                page.locator("#surface-resolution").select_option("4")
                hover_face("a")
                assert face_snapshot() == zoomed, "Surface colors and hover must preserve rotation and zoom"
                page.locator("#surface-interactive-tab").click()
                print("Coupled binding-face interaction checks passed.", flush=True)
                assert page.evaluate("state.faceView.hoverKey") is None
                assert page.evaluate("state.selectedResidues") == original_selection
                assert page.evaluate("state.stage.viewerControls.getOrientation().elements") == original_camera
                check_structure_selection()
                check_surface_selection()
                open_faces()
                assert face_snapshot() == zoomed, "Tab switches preserve rotation and zoom"
                page.locator(".interface-card.active").click()
                page.wait_for_function("state.faceView?.ready && state.faceView.panels.every(panel => panel.stage.tasks.count === 0)")
                assert page.evaluate("state.faceView.hoverKey") is None
                check_face_frames()
                page.locator("#surface-interactive-tab").click()

            def check_face_loading_and_mobile(pdb_id):
                open_analysis("mmcif")
                expect(page.locator("#highlight-residues")).to_have_value("")
                expect(page.locator("#pocket-residues")).to_have_value("")
                expect(page.locator("#pocket-mode-anchor")).to_be_checked()
                page.locator("#residue-table .focus-row").first.click()
                check_structure_selection()
                check_surface_selection()
                assert page.evaluate("state.structureFormat") == "mmcif"
                open_faces()
                check_face_frames()
                mmcif_selection = page.evaluate("state.selectedResidues")
                hover_face("a")
                hover_face("b")
                assert page.evaluate("state.selectedResidues") == mmcif_selection
                page.locator("#surface-faces-panel").screenshot(path=str(artifacts / "binding-faces-mmcif.png"))
                page.locator("#surface-interactive-tab").click()

                # An older surface response must not overwrite the newly opened analysis.
                held = []
                page.route(f"**/api/analyses/{pdb_id}/curvature", lambda route: held.append(route))
                page.locator(f'.history-row[data-id="{pdb_id}"]').click()
                page.wait_for_function("state.curvatureLoaded && state.surface === null")
                assert held
                page.locator("#highlight-residues").fill("A:67")
                page.locator('#residue-highlight-form button[type="submit"]').click()
                check_structure_selection()
                expect(page.locator("#highlight-feedback")).to_contain_text("atom highlighting")
                page.locator("#surface-faces-tab").click()
                wait_faces(surfaces=False)
                assert page.evaluate("state.faceView.panels.every(panel => panel.meshComponent === null)")
                check_face_highlights()
                rotate_faces("a")
                fallback_zoom = page.evaluate("state.faceView.panels.map(panel => panel.stage.viewer.camera.position.z)")
                fallback_rotation = page.evaluate("state.faceView.panels.map(panel => panel.stage.viewerControls.rotation.toArray())")
                hover_face("a", surfaces=False)
                # Adding the late surface keeps rotation, zoom, and selections.
                held[0].fulfill(status=200, content_type="application/json", body=json.dumps(fixtures["surface"]))
                wait_faces()
                assert page.evaluate("state.faceView.panels.map(panel => panel.stage.viewer.camera.position.z)") == fallback_zoom
                assert page.evaluate("state.faceView.panels.map(panel => panel.stage.viewerControls.rotation.toArray())") == fallback_rotation
                check_face_frames()
                check_face_highlights()
                hover_face("a")
                # An unavailable optional mesh still leaves cartoon picking usable.
                held.clear()
                page.evaluate("() => { state.curvatureLoaded = false; state.surface = null; void loadCurvature(); }")
                page.wait_for_function("state.curvatureLoaded && state.surface === null")
                held[0].fulfill(status=200, content_type="application/json", body=json.dumps({
                    "status": "unavailable", "reason": "Surface fixture unavailable", "report": None, "meshes": {},
                }))
                expect(page.locator("#face-view-status")).to_contain_text("surfaces unavailable")
                assert page.evaluate("state.faceView.panels.every(panel => panel.meshComponent === null)")
                check_face_highlights()
                hover_face("a", surfaces=False)
                page.locator("#surface-interactive-tab").click()
                # Hold a fresh request to test an older response across analyses.
                held.clear()
                page.evaluate("() => { state.curvatureLoaded = false; state.surface = null; void loadCurvature(); }")
                page.wait_for_function("state.curvatureLoaded && state.surface === null")
                assert held
                open_analysis("mmcif")
                stale = copy.deepcopy(fixtures["surface"])
                stale["report"]["source"] = "stale response"
                held[0].fulfill(status=200, content_type="application/json", body=json.dumps(stale))
                page.wait_for_timeout(250)
                assert page.evaluate("state.surface.report.source") != "stale response"
                page.unroute(f"**/api/analyses/{pdb_id}/curvature")

                # Errors have a local retry, and stale pair/analysis responses cannot
                # build a viewport or populate the current geometry cache.
                face_url = f"**/api/analyses/{pdb_id}/pairs/*/face-view"
                held_faces = []
                page.route(face_url, lambda route: held_faces.append(route))
                open_analysis("pdb")
                page.locator("#surface-faces-tab").click()
                page.wait_for_function("state.faceLoading")
                assert held_faces
                held_faces.pop(0).fulfill(status=503, content_type="application/json", body=json.dumps({"detail": "Temporary geometry error"}))
                expect(page.locator("#face-view-retry")).to_be_visible()
                expect(page.locator("#face-view-status")).to_contain_text("Temporary geometry error")
                assert page.evaluate("state.faceView") is None
                page.unroute(face_url)
                page.locator("#face-view-retry").click()
                wait_faces()
                geometry = page.evaluate("state.faceView.geometry")
                page.route(face_url, lambda route: held_faces.append(route))
                page.evaluate("state.faceCache.clear()")
                page.locator(".interface-card.active").click()
                page.wait_for_function("state.faceLoading")
                assert len(held_faces) == 1
                page.locator(".interface-card.active").click()
                page.wait_for_function("state.faceLoading")
                page.wait_for_timeout(100)
                assert len(held_faces) == 2
                stale_geometry = {**geometry, "pair_id": "stale pair response"}
                held_faces.pop(0).fulfill(status=200, content_type="application/json", body=json.dumps(stale_geometry))
                page.wait_for_timeout(100)
                assert page.evaluate("state.faceView") is None
                held_faces.pop(0).fulfill(status=200, content_type="application/json", body=json.dumps(geometry))
                wait_faces()
                assert page.evaluate("state.faceView.geometry.pair_id") == geometry["pair_id"]
                page.evaluate("state.faceCache.clear()")
                page.locator(".interface-card.active").click()
                page.wait_for_function("state.faceLoading")
                assert len(held_faces) == 1
                open_analysis("mmcif")
                open_faces()
                held_faces.pop(0).fulfill(status=200, content_type="application/json", body=json.dumps(stale_geometry))
                page.wait_for_timeout(100)
                assert page.evaluate("state.faceView.geometry.pair_id") == geometry["pair_id"]
                assert page.evaluate("state.faceCache.get(state.pairId).pair_id") == geometry["pair_id"]
                page.unroute(face_url)
                page.locator("#surface-interactive-tab").click()

                page.locator("#pair-workspace").screenshot(path=str(artifacts / "interface-and-contact-map.png"))
                page.set_viewport_size({"width": 390, "height": 900})
                panel_box = page.locator("#residue-highlight-panel").bounding_box()
                viewer_box = page.locator("#surface-view").bounding_box()
                assert panel_box["y"] >= viewer_box["y"] + viewer_box["height"]
                page.locator("#residue-highlight-toggle").click()
                expect(page.locator("#residue-highlight-body")).to_be_hidden()
                page.locator("#residue-highlight-toggle").click()
                expect(page.locator("#residue-highlight-body")).to_be_visible()
                page.locator("#workspace").screenshot(path=str(artifacts / "mobile-workspace.png"))
                assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
                open_faces()
                page.wait_for_function("state.faceView.panels.every(panel => Math.abs(panel.lastWidth "
                                       "- panel.host.getBoundingClientRect().width) < 1)")
                check_face_frames()
                boxes = [page.locator(f"#face-view-{side}").bounding_box() for side in ("a", "b")]
                assert abs(boxes[0]["y"] - boxes[1]["y"]) < 2
                assert boxes[0]["x"] + boxes[0]["width"] <= boxes[1]["x"]
                assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
                hover_face("a")
                # Chrome sends a real one-finger gesture to the mobile canvas.
                before_touch = face_snapshot()
                canvas = page.locator("#face-view-b canvas").bounding_box()
                x, y = canvas["x"] + canvas["width"] / 2, canvas["y"] + canvas["height"] / 2
                touch = page.context.new_cdp_session(page)
                try:
                    touch.send("Input.dispatchTouchEvent", {"type": "touchStart", "touchPoints": [{"x": x, "y": y}]})
                    touch.send("Input.dispatchTouchEvent", {"type": "touchMove", "touchPoints": [{"x": x + 16, "y": y + 12}]})
                    touch.send("Input.dispatchTouchEvent", {"type": "touchEnd", "touchPoints": []})
                finally:
                    touch.detach()
                page.wait_for_function("previous => state.faceView.panels[0].stage.viewerControls.rotation.toArray() "
                                       ".some((value, index) => Math.abs(value-previous[index]) > 1e-6)",
                                       arg=before_touch[0]["rotation"])
                after_touch = face_snapshot()
                check_face_rotation_alignment()
                assert [panel["zoom"] for panel in after_touch] == [panel["zoom"] for panel in before_touch]
                check_face_highlights()
                page.locator("#surface-faces-panel").screenshot(path=str(artifacts / "binding-faces-mobile.png"))
                transforms = page.evaluate("state.faceView.panels.map(panel => panel.transform.elements)")
                page.set_viewport_size({"width": 1200, "height": 950})
                page.wait_for_function("state.faceView.panels.every(panel => Math.abs(panel.lastWidth "
                                       "- panel.host.getBoundingClientRect().width) < 1)")
                check_face_frames()
                assert page.evaluate("state.faceView.panels.map(panel => panel.stage.viewerControls.rotation.toArray())") == [panel["rotation"] for panel in after_touch]
                assert page.evaluate("state.faceView.panels.map(panel => panel.transform.elements)") == transforms
                page.evaluate("window.disposedFacePanels = state.faceView.panels; closeWorkspace()")
                assert page.evaluate("state.faceView === null && state.faceCache.size === 0 "
                                     "&& disposedFacePanels.every(panel => panel.stage.compList.length === 0 "
                                     "&& panel.host.childElementCount === 0 "
                                     "&& !panel.stage.viewerControls.signals.changed.has(panel.rotationChanged))")
                assert not errors, errors
                print("Binding-face browser checks passed: mirrored mouse/touch rotation, persistent table selections, linked hover, movement locks, independent zoom, loading/retry, stale geometry, resize/disposal, PDB/mmCIF, and mobile layout.", flush=True)

            pdb_id = open_analysis("pdb")
            check_surface_picking()
            check_interactive_picking()
            check_surface_modes()
            check_electrostatics_loading(pdb_id)
            if picking_only:
                assert not errors, errors
                print("Surface picking browser checks passed.", flush=True)
                return
            if faces_only:
                page.locator("#residue-table .focus-row").first.click()
                check_face_interaction()
                check_face_loading_and_mobile(pdb_id)
                return
            check_boundary_mesh()
            check_surface_selection()
            contact_key = page.evaluate("state.result.pairs.find(pair => pair.id === state.pairId).contacts[0].residue_a.key")
            page.locator("#highlight-residues").fill(f"A:63-67, {contact_key}, 99999")
            page.locator('#residue-highlight-form button[type="submit"]').click()
            expected_keys = list(dict.fromkeys([f"A:{number}" for number in range(63, 68)] + [contact_key]))
            assert page.evaluate("state.selectedResidues") == expected_keys
            expect(page.locator("#highlight-feedback")).to_contain_text("99999")
            assert page.locator("#residue-table .focus-row.selected").count() > 0
            check_structure_selection()
            check_surface_selection()
            selected = page.evaluate("state.selectedResidues")
            page.locator("#highlight-residues").fill("Z:99999")
            page.locator('#residue-highlight-form button[type="submit"]').click()
            assert page.evaluate("state.selectedResidues") == selected
            expect(page.locator("#highlight-feedback")).to_contain_text("retained")
            page.locator("#highlight-clear").click()
            assert page.evaluate("state.selectedResidues") == []
            expect(page.locator("#highlight-residues")).to_have_value("Z:99999")
            check_surface_selection()
            page.locator("#highlight-residues").fill("A:63-67")
            page.locator('#residue-highlight-form button[type="submit"]').click()
            width = page.locator("#surface-view").bounding_box()["width"]
            page.locator("#residue-highlight-toggle").click()
            expect(page.locator("#residue-highlight-body")).to_be_hidden()
            expect(page.locator("#residue-highlight-toggle")).to_have_attribute("aria-expanded", "false")
            page.wait_for_function("width => document.getElementById('surface-view').getBoundingClientRect().width > width + 150", arg=width)
            assert page.evaluate("localStorage.getItem('highlightPanelCollapsed')") == "1"
            page.locator("#residue-highlight-toggle").click()
            expect(page.locator("#residue-highlight-body")).to_be_visible()
            page.wait_for_function("width => Math.abs(document.getElementById('surface-view').getBoundingClientRect().width - width) < 2", arg=width)
            # Changing interfaces retains typed input but resets the shared selection.
            page.locator(".interface-card.active").click()
            assert page.evaluate("state.selectedResidues") == []
            expect(page.locator("#highlight-residues")).to_have_value("A:63-67")
            expect(page.locator("#highlight-feedback")).to_contain_text("not been applied")
            page.locator('#residue-highlight-form button[type="submit"]').click()
            expect(page.locator("#residue-panel")).to_be_visible()
            page.locator("#residue-table .focus-row").first.click()
            assert page.evaluate("state.selectedResidues.length") == 1
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
            check_face_interaction()
            # The pair-reset check intentionally clears shared selection.
            page.locator("#curvature-table .focus-row").first.click()

            # Highlighted triangles must map their duplicated corners back to
            # the original residue after chain separation and camera transforms.
            assert page.evaluate("""() => {
              const key = state.selectedResidues[0];
              const {chain} = selectionParts(key);
              const mesh = state.surface.meshes[chain];
              const display = state.meshGeometry[chain];
              const vertex = display.vertexMap.findIndex((original, vertex) =>
                vertex >= display.highlightedVertexStart && mesh.residue[original] === residueLabelFromKey(key));
              const pid = Math.floor(display.index.indexOf(vertex) / 3);
              const component = state.meshComponents[chain];
              const point = new NGL.Vector3().fromArray(display.position, vertex * 3).applyMatrix4(component.matrix);
              const canvasPosition = state.stage.viewerControls.getPositionOnCanvas(point);
              const pick = {type:'mesh', component, pid, canvasPosition};
              const hit = focusSurfaceHit(pick);
              if (!hit || hit.vertex !== display.vertexMap[vertex] || surfaceResidueKey(chain, hit.label) !== key
                || hit.h !== mesh.h[state.surfaceScale][hit.vertex]) return false;
              state.stage.signals.hovered.dispatch(pick);
              if (!document.getElementById('surface-tip').textContent.includes(hit.label)) return false;
              state.stage.signals.clicked.dispatch(pick);
              hideSurfaceTip();
              return state.selectedResidues.length === 1 && state.selectedResidues[0] === key
                && focusSurfaceHit({...pick, pid: -1}) === null
                && focusSurfaceHit({...pick, pid: display.index.length / 3}) === null
                && focusSurfaceHit({...pick, canvasPosition:null}) === null;
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
            check_interactive_picking()
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
            page.locator("#pocket-submit").click()
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

            page.locator("#highlight-residues").fill("A:63, 64, 66, 67, 70, B:12")
            page.locator('#residue-highlight-form button[type="submit"]').click()
            check_structure_selection()
            check_surface_selection()
            page.locator("#pocket-radius").fill("9")
            page.locator("#pocket-mode-residues").check()
            expect(page.locator("#pocket-anchor-field")).to_be_hidden()
            expect(page.locator("#pocket-residue-field")).to_be_visible()
            expect(page.locator("#pocket-radius")).to_have_value("6")
            page.locator("#pocket-radius").fill("6.5")
            page.locator("#pocket-mode-anchor").check()
            expect(page.locator("#pocket-radius")).to_have_value("9")
            page.locator("#pocket-mode-residues").check()
            expect(page.locator("#pocket-radius")).to_have_value("6.5")
            page.locator("#pocket-radius").fill("6")
            page.locator("#pocket-use-selection").click()
            expect(page.locator("#pocket-residues")).to_have_value("A:63, A:64, A:66, A:67, A:70, B:12")
            page.locator("#pocket-submit").click()
            page.wait_for_function("state.pocketResult?.mode === 'residues'", timeout=120_000)
            assert page.evaluate("state.pocketResult.free.primary.defined_residues_present") == 5
            assert page.evaluate("state.pocketResult.bound.primary.defined_residues_present") == 6
            expect(page.locator("#pocket-result")).to_contain_text("Defined residues lining")
            expect(page.locator("#cancel-pocket-highlight")).to_be_enabled()
            assert page.evaluate("state.pocketComponent.reprList[0].getParameters().opacity") == 0.2
            metrics = page.locator("#pocket-result").inner_text()
            page.locator("#pocket-residues").fill("67, 99999")
            page.locator("#pocket-submit").click()
            expect(page.locator("#pocket-feedback")).to_contain_text("Unknown or invalid")
            assert page.locator("#pocket-result").inner_text() == metrics
            assert page.evaluate("state.pocketHighlighted")
            page.locator("#pocket-use-selection").click()
            page.locator("#surface-panel").screenshot(path=str(artifacts / "residue-highlight-panel.png"))
            page.locator("#pocket-panel").screenshot(path=str(artifacts / "residue-defined-pocket.png"))

            check_face_loading_and_mobile(pdb_id)
            print("Browser smoke passed: residue input, linked selections, anchor and residue-defined pockets, and coupled binding faces.", flush=True)
        finally:
            browser.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--browser", default="/usr/bin/google-chrome")
    parser.add_argument("--artifacts", type=Path, default=Path("/tmp/interface-app-browser"))
    parser.add_argument("--reuse-fixtures", action="store_true")
    parser.add_argument("--binding-faces-only", action="store_true", help="Run side-by-side rotation, selection, and hover checks without repeating pocket jobs")
    parser.add_argument("--picking-only", action="store_true", help="Run surface picking regression checks")
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
            exercise_browser(url, args.browser, artifacts, fixtures, faces_only=args.binding_faces_only, picking_only=args.picking_only)
        finally:
            server.terminate()
            try:
                server.wait(timeout=15)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()


if __name__ == "__main__":
    main()
