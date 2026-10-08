"""Exercise sign conventions and masks in the browser's shared color function."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest


NODE = shutil.which("node")
APP_JS = Path(__file__).resolve().parents[1] / "static/app.js"
pytestmark = pytest.mark.skipif(NODE is None, reason="Node is needed to exercise surface colors")

DRIVER = r"""
const fs = require('fs');
const vm = require('vm');
const input = JSON.parse(fs.readFileSync(0, 'utf8'));
const source = fs.readFileSync(process.argv[1], 'utf8');
const context = {state: {
  surfaceMode: input.mode, surfaceScale:'6', surfacePadding:0, pairId:'AB', curvatureJobId:'surface',
  surface:{report:{color_limit_Ainv:1}, meshes:{A:{
    position:Array(12).fill(0), h:{'6':[-1,0,1,1]}, contact_distance_dA:{B:[0,0,0,255]}
  }}},
  electrostatics:{status:'complete', curvature_job_id:input.curvatureId || 'surface',
    report:{color_limit_kT_e:1}, meshes:{A:{potential_kT_e:[-1,0,1,1]}}},
  result:{pairs:[{id:'AB',chain_a:'A',chain_b:'B'}]},
  meshGeometry:{A:{vertexMap:[0,1,2,3,2],highlightedVertexStart:4}},
}};
vm.createContext(context);
vm.runInContext(source.slice(source.indexOf('const curvatureColors ='), source.indexOf('const cartoonColors =')), context);
for (const name of ['curvatureReport','electrostaticsReady','interpolateColor','meshColors']) {
  const start = source.indexOf('function ' + name + '(');
  vm.runInContext(source.slice(start, source.indexOf('\n}', start) + 2), context);
}
process.stdout.write(JSON.stringify({ready:context.electrostaticsReady(), colors:Array.from(context.meshColors('A'))}));
"""


def colors(mode, curvature_id=None):
    response = subprocess.run([NODE, "-e", DRIVER, str(APP_JS)],
                              input=json.dumps({"mode": mode, "curvatureId": curvature_id}),
                              text=True, capture_output=True, check=True, timeout=5)
    return json.loads(response.stdout)


@pytest.mark.parametrize("mode, negative, neutral, positive", [
    ("convexity", [217, 95, 138], [250, 250, 250], [27, 158, 119]),
    ("electrostatics", [33, 102, 172], [247, 247, 247], [178, 24, 43]),
])
def test_signs_neutral_padding_and_yellow_selection(mode, negative, neutral, positive):
    result = colors(mode)
    assert result["ready"]
    expected = [component / 255 for rgb in [negative, neutral, positive, [157, 164, 172]] for component in rgb]
    expected += [1, 1, 0]
    assert result["colors"] == pytest.approx(expected)


def test_potentials_from_another_surface_are_unavailable():
    result = colors("electrostatics", "old-surface")
    assert not result["ready"]
    assert result["colors"] == []
