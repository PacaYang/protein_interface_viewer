"""Exercise cached surface highlights with the bundled NGL buffer API."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest


NODE = shutil.which("node")
APP_JS = Path(__file__).resolve().parents[1] / "static/app.js"
pytestmark = pytest.mark.skipif(NODE is None, reason="Node is needed to exercise surface geometry")

DRIVER = r"""
const fs = require('fs'), vm = require('vm'), assert = require('assert');
global.window = {performance, navigator:{userAgent:'Node geometry test'}, document:{}, location:{search:''}};
global.self = global;
global.navigator = window.navigator;
try {
  const NGL = require(process.argv[2]);
  const source = fs.readFileSync(process.argv[1], 'utf8');
  const context = {state:{meshSourceData:new WeakMap()}, NGL};
  vm.createContext(context);
  for (const name of ['surfaceMeshGeometry', 'updateSurfaceMeshSelection', 'installSurfacePicking']) {
    const start = source.indexOf('function ' + name + '(');
    vm.runInContext(source.slice(start, source.indexOf('\n}', start) + 2), context);
  }
  // Shared vertices, a repeated residue in a face, an insertion code, and
  // non-coplanar faces expose color bleeding and incorrect selection counts.
  const mesh = {
    position:[0,0,0, 1,0,0, 0,1,0, 1,1,1, 2,1,0],
    index:[0,1,2, 1,3,2, 3,4,2],
    residue:['R12','G13B','G13B','S14','T15'],
  };
  const display = context.surfaceMeshGeometry(mesh, new Set());
  const independent = context.surfaceMeshGeometry(mesh, new Set(['T15']));
  assert.strictEqual(display.source, independent.source);
  assert.strictEqual(display.position, independent.position);
  assert.strictEqual(display.normal, independent.normal);
  assert.notStrictEqual(display.index, independent.index);
  const shape = new NGL.Shape('selection-regression');
  const colors = new Float32Array(display.position.length).fill(.5);
  for (let v = 5; v < 10; v++) colors.set([1,1,0], v * 3);
  shape.addMesh(display.position, colors, display.index, display.normal);
  context.installSurfacePicking(shape, display);
  const buffer = shape.bufferList[0], geometry = buffer.geometry;
  const positions = geometry.attributes.position.array, normals = geometry.attributes.normal.array;
  const index = geometry.index.array, picking = buffer.getPickingMesh().geometry;
  const snapshots = [];
  for (const [labels, expected] of [
    [['R12'],[true,false,false]],
    [['R12','S14'],[true,true,true]],
    [['S14'],[false,true,true]],
    [['G13B','S14'],[true,true,true]],
    [['G13B'],[true,true,true]],
    [['T15'],[false,false,true]],
    [[],[false,false,false]],
    [['missing'],[false,false,false]],
  ]) {
    if (context.updateSurfaceMeshSelection(display, new Set(labels))) buffer.setAttributes({index:display.index});
    assert.strictEqual(buffer.geometry, geometry);
    assert.strictEqual(geometry.attributes.position.array, positions);
    assert.strictEqual(geometry.attributes.normal.array, normals);
    assert.strictEqual(geometry.index.array, index);
    assert.strictEqual(buffer.getPickingMesh().geometry, picking);
    assert.strictEqual(context.updateSurfaceMeshSelection(display, new Set(labels)), false);
    for (let f = 0; f < 3; f++) for (let c = 0; c < 3; c++) {
      const original = mesh.index[f * 3 + c], vertex = index[f * 3 + c];
      assert.strictEqual(vertex, original + (expected[f] ? 5 : 0));
      assert.strictEqual(display.vertexMap[vertex], original);
      assert.strictEqual(picking.attributes.primitiveId.array[f * 3 + c], f);
      for (let axis = 0; axis < 3; axis++) {
        assert.strictEqual(positions[vertex * 3 + axis], positions[original * 3 + axis]);
        assert.strictEqual(normals[vertex * 3 + axis], normals[original * 3 + axis]);
        assert.strictEqual(picking.attributes.position.array[(f * 3 + c) * 3 + axis], positions[vertex * 3 + axis]);
        assert.strictEqual(geometry.attributes.color.array[vertex * 3 + axis], expected[f] ? [1,1,0][axis] : .5);
      }
    }
    snapshots.push(Array.from(index));
  }
  assert.deepStrictEqual(Array.from(display.selectedFaceCounts), [0,0,0]);
  assert.deepStrictEqual(Array.from(independent.index), [0,1,2,1,3,2,8,9,7]);
  shape.dispose();
  process.stdout.write(JSON.stringify({snapshots, vertices:positions.length / 3}));
} catch (error) { console.error(error.stack); process.exitCode = 1; }
"""


def test_selection_updates_preserve_geometry_picking_normals_and_other_viewers():
    response = subprocess.run(
        [NODE, "-e", DRIVER, str(APP_JS), str(APP_JS.with_name("ngl.js"))],
        text=True, capture_output=True, check=True, timeout=10,
    )
    result = json.loads(response.stdout)
    assert result["vertices"] == 10
    assert len(result["snapshots"]) == 8
