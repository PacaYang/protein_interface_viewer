# Protein interface analyzer

This repository contains the complete local browser application and its
scientific calculations. It accepts PDB/mmCIF uploads and RCSB PDB IDs, analyzes every
protein-chain pair, and lets a user measure selected binding pockets in free and
partner-bound states.

## Run

Use Python 3.10 or newer. Python 3.12 was used for validation. From this
repository's root, create an environment and install the runtime dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_interface_app.py --data-dir .interface_app_data
```

Open <http://127.0.0.1:8000>. The default data directory contains only the local
SQLite history, uploaded/downloaded coordinates, and JSON result artifacts.
The three structures in `examples/` can be uploaded directly for an initial
analysis. The viewer's NGL JavaScript and license are bundled locally.

All runtime source files are contained in this repository. The default app
history is created on first launch and is excluded by `.gitignore`. To reuse an
existing history, pass its location with `--data-dir`.

## Repository contents

| Path | Purpose |
| --- | --- |
| `interface_app/` | FastAPI endpoints, analysis workers, storage, browser interface, and app tests |
| `run_interface_app.py` | Server launcher |
| `zernike_convexity.py` | Shared surface and Zernike curvature calculation used by the app and standalone tools |
| `requirements.txt` | Runtime dependencies; SciPy retains the spherical-harmonic API used by the calculation |
| `requirements-dev.txt` | Runtime dependencies plus pytest and Playwright |
| `pytest.ini` and `tests/` | Test discovery configuration and scientific geometry tests |
| `examples/` | PDB fixtures and upload examples: 1IAR, 3BPN, and 3BPO |
| `standalone/` | Original contact-map generator, HTML template, analysis scripts, and reference data |
| `standalone/reference_results/` | Packaged HTML pages, numerical reports, and preview; ignored by Git as generated outputs |

See `standalone/README.md` for regeneration commands. The app and standalone
tools resolve their assets and test fixtures within this repository; the
original research workspace is not required.

## Start a new Git repository

Run these commands from this folder:

```bash
git init -b main
git add .
git status --short
git commit -m "Initial protein interface analyzer"
```

The ignore rules keep local history, uploaded structures, calculation outputs,
Python environments, caches, and generated reference pages out of commits.
The example structures and small reference JSON inputs are included in commits
so tests and standalone regeneration remain reproducible. The bundled NGL MIT
license is retained in `interface_app/static/License_for_NGL.txt`.

## Workspace behavior

The app starts the core interface calculation and curvature calculation in
separate worker processes. Contact tables become available before curvature
finishes. Pocket jobs are submitted explicitly from a selected interface after
choosing a target chain, residue anchor, partner chain, and 4–12 Å radius.

The workspace mirrors the standalone `contact_map.html` view while keeping the
analysis data local to the browser app:

- `INTERFACES` is a collapsible chain-pair list on the left. Selecting a pair
  opens its `Residue interface` and `Closest contact pairs` tabs in the main
  pane.
- The contact map retains every residue pair with a heavy-atom contact through
  6 Å. Its distance cutoff, interaction type (all, tight, H-bond, salt, or
  like-charge), numeric cell labels, and ΔSASA lane can be changed without
  rerunning the calculation. Clicking a row, map cell, axis label, or surface
  vertex selects the corresponding residue(s) in all views.
  Charge lanes use the same glyphs as the standalone view for polar,
  hydrophobic, aromatic, and Gly/Pro/Cys residues, and the contact table
  includes each residue's class, charge, and ΔSASA. Sequence gaps are spaced
  apart, axis labels carry curvature marks, and missing ΔSASA is shown as a
  hollow lane.
- `3D Zernike surface` displays the signed mean curvature `H` from the local
  Zernike fit. Positive `H` is convex/outward and negative `H` is
  concave/inward. The 4, 6, and 8 Å neighborhoods are retained in the result;
  the resolution selector changes the displayed scale. Surface opacity,
  distance padding from the selected contact face (0–8 Å), separation of the
  opposing chains, and per-chain visibility are view-only controls.
  Selected residues appear as yellow surface patches: every triangle touching
  a vertex assigned to a selected residue is colored yellow, including faces
  outside the contact padding. The patches follow the surface opacity control
  and preserve neighboring faces' curvature colors and smooth shading.
  Atom-level ball-and-stick selection remains available alongside the patches.
- Pocket analysis reports free and partner-bound geometry. When a bound primary
  pocket has a mesh, it is highlighted automatically in the 3D view at 20%
  opacity. The `Cancel highlight` action removes that mesh while leaving the calculated
  pocket metrics visible; `Highlight bound pocket` restores it.

Curvature results use a versioned payload (`0.2.0`) with a `report.scales`
mapping, whole-chain and pair-interface summaries, per-vertex mesh positions,
faces, residue labels, signed `H` values, and deci-Å contact-distance arrays.
The adapter accepts any number of protein chains and arbitrary PDB/mmCIF chain
identifiers. If the optional numerical dependencies are unavailable, core
contact and pocket results remain usable and the surface is reported as
unavailable.

Saved results from earlier versions do not contain the contact grid or surface
meshes. Upload those coordinates again or retrieve the PDB entry again to
generate the new views.

Leaving the assembly field blank retrieves the deposited asymmetric unit. Some
crystal entries contain more than one copy of the biological complex in that
file; for example, 5VI4 contains protein chains A–F, while assembly 1 selects
A–C and assembly 2 selects D–F.

## Analysis definitions

- Contacts use heavy atoms. The contact map retains contacts through 6 Å; the
  reported interface contact-pair count uses the historical 5 Å cutoff.
- Buried surface area is calculated from all retained amino-acid residues:
  `(SASA(A) + SASA(B) - SASA(A:B)) / 2`, using a 1.4 Å probe and 960 points.
- A pocket is standard-probe-accessible space that is not connected to exterior
  space accessible to a larger 4 Å probe. Its volume is a 0.6 Å voxel estimate.
- Pocket depth is a six-neighbor shortest path from a detected mouth to the
  deepest voxel. Fully enclosed cavities have no mouth and report depth as
  unavailable. Enclosure is `1 - opening area / boundary area`.
- Pocket solvent exposure reports the opening area and the SASA of residues
  lining the selected component; both are retained in the JSON result.

The app does not relax coordinates after removing a partner. “Free” means the
target chain is analyzed alone using its deposited coordinates.

## Validation

```bash
python -m pip install -r requirements-dev.txt
pytest -q
node --check interface_app/static/app.js
python -m py_compile interface_app/*.py
```

The optional browser smoke test requires Playwright and a Chrome browser. It
uses `examples/1IAR.pdb`, a converted mmCIF, a temporary API server, and real pocket jobs
to check yellow residue faces and boundaries, preserved normals, linked
selections and picking, surface controls, pocket transparency and cancellation,
delayed responses, and mobile layout. It saves screenshots and logs under the
artifact directory:

```bash
python interface_app/tests/browser_smoke.py --browser /usr/bin/google-chrome --artifacts /tmp/interface-app-browser
```

Add `--reuse-fixtures` when rerunning browser checks without numerical changes.
An existing Chrome/Chromium installation can be selected with `--browser`.
