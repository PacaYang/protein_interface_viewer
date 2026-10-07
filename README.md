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
choosing a target chain, partner chain, residue anchor or residue list, and
4–12 Å search distance.

The workspace mirrors the standalone `contact_map.html` view while keeping the
analysis data local to the browser app:

- A collapsible side menu holds `Upload coordinates`, `Retrieve from RCSB`,
  and `Local history`. History entries can be deleted individually or with
  `Clear all`. Analyses with queued or running jobs are kept until they finish.
- Chains are labeled with the molecule name declared in the file (PDB
  `COMPND MOLECULE` or mmCIF `_entity.pdbx_description`) plus the chain ID.
  Analyses saved before this change are relabeled when they are opened.
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
  Surface hover and clicks identify the triangle beneath the cursor, then use
  its nearest corner on screen to resolve the residue. A separate picking mesh
  assigns one ID to each triangle, preventing blended vertex IDs from selecting
  a distant residue. This applies to both surface viewer tabs and accounts for
  camera orientation, zoom, and chain separation.
  Atom-level ball-and-stick selection remains available alongside the patches.
- The surface viewer's **Binding faces** tab shows the selected interface's
  two partner chains side by side. Each chain initially faces the viewer along
  its own binding-plane normal, using an orthographic projection and a shared up
  reference. The initial views use the same Å-per-pixel scale. Drag either
  partner to rotate both views with mirrored motion; a one-finger touch drag
  also rotates both, and Ctrl + right drag rotates around the viewing axis.
  A clockwise Z rotation in one pane produces a counterclockwise rotation in
  the other. Turns about Y are also reversed, while tilts about X move together,
  so corresponding sites keep matching heights and depths across the two
  outward binding faces. The coupling works from either pane and across mixed
  rotations. The **Couple rotation** checkbox above the views is enabled by
  default. Turn it off to rotate each chain independently; turning it back on
  aligns the views using the last pane you rotated as the reference. Zooming
  a pane does not change that reference. The setting persists across tabs and
  interfaces within the open page. Scroll or use the **+ / −** buttons to zoom
  each partner independently. Resizing refits both partners to their current
  rotation. Movement, recentering, and keyboard
  camera actions are disabled. Curvature resolution, opacity, and contact
  padding apply to both views.
  Clicking a residue, closest-contact, contact-map, or curvature table highlights
  the selected residues in this tab as yellow surface patches and atoms. Typed
  residue selections appear here too, including selections made in Interactive
  3D before opening Binding faces. Multiple selected residues can appear on
  either chain and remain highlighted until the shared selection changes.
  Hover a surface or cartoon residue to highlight it and its nearest residue
  on the partner chain in yellow. The readout includes their minimum heavy-atom
  distance in the original complex, without a contact cutoff. These temporary
  highlights clear on leaving the structure, switching tabs, or changing the
  interface; shared residue/table/contact-map selections stay highlighted.
  Binding planes are fitted by SVD to one equally weighted heavy-atom centroid
  per reported binding residue (normally within 5 Å). Sparse, collinear, or
  ambiguous interfaces use a labeled estimated direction toward the partner's
  centroid. Interfaces with no binding residues or no stable partner direction
  show an explanatory message. Cartoons remain usable while surfaces load.
  This view also works with saved analyses, using
  `GET /api/analyses/{analysis_id}/pairs/{pair_id}/face-view` without rerunning
  analysis or changing stored results.
- A collapsible **Highlight residues** panel sits beside the surface viewer
  (below it on mobile). Choose a default chain and enter numbers, insertion codes,
  or inclusive ranges, such as `A:45-52, 60, 63A, B:101`. Commas and spaces
  separate entries; a chain prefix applies until the next prefix. Ranges follow
  actual residues in the coordinate file, including insertion codes and sequence
  gaps, and require both endpoints to exist. Negative numbers are supported,
  such as `A:-3--1`. **Highlight** applies matched residues to the shared yellow
  surface, atom, table, and contact-map selection and reports unmatched input.
  An entirely invalid entry retains the previous selection. **Clear** clears
  the selection and keeps the text. Later residue clicks replace the selection.
  The panel remembers whether it is collapsed.
- Pocket analysis offers **Anchor residue** and **Residue list** modes.
  **Residue list** accepts the same syntax for up to 60 residues from either
  chain of the selected interface; unqualified numbers use the target chain.
  **Use current selection** copies the highlighted residues into the definition.
  Invalid pocket input is reported before a job starts, preserving the displayed
  result. The default search distance is 8 Å for an anchor and 6 Å for a residue
  list, and each mode remembers its adjusted distance within the open analysis.
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
- An anchor searches within a sphere around that residue's heavy-atom centroid.
  A residue list searches the union of spheres around every heavy atom of the
  listed residues. The free and bound states use the same spatial region,
  including the deposited coordinates of any listed partner residues.
  Disconnected regions are measured separately in residue-list mode.
- Residue-defined pockets report which listed residues line the component and
  their coverage fraction. The denominator includes only listed residues present
  in that state, so partner residues are excluded in the free state. A
  partner-only definition has unavailable free-state coverage. Candidates are
  ranked by coverage, then volume, before retaining the best five; free-state
  candidates with unavailable coverage are ranked by volume. The selected
  primary pocket is reported for each state.
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
to check typed residue highlighting, validation, collapse/expand, yellow residue
faces and boundaries, preserved normals, linked selections and picking,
surface controls, anchor and residue-list pocket jobs, pocket transparency and
cancellation, delayed responses, and mobile layout. It saves screenshots and logs under the
artifact directory:

```bash
python interface_app/tests/browser_smoke.py --browser /usr/bin/google-chrome --artifacts /tmp/interface-app-browser
```

Add `--reuse-fixtures` when rerunning browser checks without numerical changes.
Add `--binding-faces-only` to check mirrored coupled and independent rotation,
the rotation toggle and re-coupling, persistent table selections,
linked hover, movement locks, independent zoom, loading and retry, stale geometry,
and mobile resizing without repeating the pocket jobs. The geometry tests also
compare unrestricted nearest residue matches against brute force and cover
rotated, sparse, and ambiguous
binding planes, insertion codes, and chain identifiers with multiple characters.
Add `--picking-only` to check GPU triangle picking, highlighted vertices,
both surface sides, transforms, zoom, picking ID byte boundaries, and real
interactive hover/clicks without repeating the pocket jobs.
An existing Chrome/Chromium installation can be selected with `--browser`.
