"""Exercise the browser's residue parsing and pocket input validation with Node."""

import json
import shutil
import subprocess
from pathlib import Path

import pytest


APP_JS = Path(__file__).resolve().parents[1] / "static" / "app.js"
NODE = shutil.which("node")
pytestmark = pytest.mark.skipif(NODE is None, reason="Node is needed to exercise browser input helpers")

DRIVER = r"""
const fs = require("fs");
const vm = require("vm");
const input = JSON.parse(fs.readFileSync(0, "utf8"));
const source = fs.readFileSync(process.argv[1], "utf8");
const context = {
  state: {result: {metadata: {chains: input.chains}}},
  fields: {
    "pocket-residues": {value: input.text},
    "pocket-target": {value: input.defaultChain},
    "pocket-partner": {value: "B"},
  },
};
context.$ = id => context.fields[id];
vm.createContext(context);
// Load input helpers without bootstrapping the viewer or a DOM.
for (const name of ["residueKey", "parseResidueSpec", "selectionParts", "pocketResidueRequest"]) {
  const start = source.indexOf("function " + name + "(");
  if (start === -1) throw Error("Missing helper " + name);
  vm.runInContext(source.slice(start, source.indexOf("\n}", start) + 2), context);
}
try {
  const result = input.operation === "pocket"
    ? context.pocketResidueRequest()
    : context.parseResidueSpec(input.text, input.defaultChain);
  process.stdout.write(JSON.stringify(result));
} catch (error) {
  process.stdout.write(JSON.stringify({error: error.message}));
}
"""


def identity(chain, number, insertion=""):
    return {"chain_id": chain, "number": number, "insertion_code": insertion,
            "key": f"{chain}:{number}{insertion}"}


CHAINS = [
    {"id": "A", "residues": [
        identity("A", -3), identity("A", -2), identity("A", 1),
        identity("A", 3), identity("A", 3, "A"), identity("A", 3, "B"),
        identity("A", 5), identity("A", 63), identity("A", 63, "A"),
        identity("A", 64), identity("A", 66), identity("A", 67), identity("A", 70),
    ]},
    {"id": "B", "residues": [identity("B", 12, "A"), identity("B", 14)]},
    {"id": "alpha", "residues": [identity("alpha", 100)]},
]


def run_input(text, default_chain="A", operation="parse", chains=None):
    result = subprocess.run(
        [NODE, "-e", DRIVER, str(APP_JS)],
        input=json.dumps({"text": text, "defaultChain": default_chain,
                         "operation": operation, "chains": chains or CHAINS}),
        text=True, capture_output=True, check=True, timeout=5,
    )
    return json.loads(result.stdout)


@pytest.mark.parametrize("text, default_chain, keys", [
    ("63, 63A 64\n66", "A", ["A:63", "A:63A", "A:64", "A:66"]),
    ("A:-3-3B, 5 B:12A, 14 alpha:100", "A",
     ["A:-3", "A:-2", "A:1", "A:3", "A:3A", "A:3B", "A:5", "B:12A", "B:14", "alpha:100"]),
    ("-3--2", "A", ["A:-3", "A:-2"]),
    ("63-67", "A", ["A:63", "A:63A", "A:64", "A:66", "A:67"]),
    ("A: 63, 63, 63A B:12A 12A", "A", ["A:63", "A:63A", "B:12A"]),
    ("14 A:63", "B", ["B:14", "A:63"]),
])
def test_residue_syntax_resolves_existing_identities(text, default_chain, keys):
    assert run_input(text, default_chain) == {"keys": keys, "unmatched": []}


def test_invalid_residues_and_ranges_are_reported_without_guessing():
    result = run_input("A:70-63 63-69 3a G63 999 B:12 Z:63 64 A:67")
    assert result == {
        "keys": ["A:67"],
        "unmatched": ["A:70-63", "63-69", "3a", "G63", "999", "B:12", "Z:63", "64"],
    }


def test_pocket_residue_request_preserves_chains_and_insertion_codes():
    assert run_input("A:-3 3B B:12A 12A", operation="pocket") == [
        {"chain_id": "A", "number": -3, "insertion_code": ""},
        {"chain_id": "A", "number": 3, "insertion_code": "B"},
        {"chain_id": "B", "number": 12, "insertion_code": "A"},
    ]


@pytest.mark.parametrize("text, message", [
    ("", "at least one"),
    ("63, 999", "Unknown or invalid"),
    ("alpha:100", "target or partner"),
])
def test_pocket_input_rejects_invalid_definitions(text, message):
    assert message in run_input(text, operation="pocket")["error"]


def test_pocket_range_cannot_exceed_sixty_residues():
    chains = [{"id": "A", "residues": [identity("A", number) for number in range(1, 66)]}]
    assert "limited to 60" in run_input("1-65", operation="pocket", chains=chains)["error"]
