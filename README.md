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

Open <http://127.0.0.1:8000>. The default data directory contains the local
SQLite history, uploaded/downloaded coordinates, JSON result artifacts, and
cached electrostatics tools downloaded on first use.
The three structures in `examples/` can be uploaded directly for an initial
analysis. The viewer's NGL JavaScript and license are bundled locally.

All runtime source files are contained in this repository. The default app
history is created on first launch and is excluded by `.gitignore`. To reuse an
existing history, pass its location with `--data-dir`.

## Repository contents

| Path | Purpose |
| --- | --- |
| `interface_app/` | FastAPI endpoints, analysis workers, storage, browser interface, electrostatics adapter, and app tests |
| `run_interface_app.py` | Server launcher |
| `zernike_convexity.py` | Shared surface and Zernike curvature calculation used by the app and standalone tools |
| `requirements.txt` | Runtime dependencies; SciPy retains the spherical-harmonic API used by the calculation |
| `requirements-cuda.txt` | Optional CUDA-enabled PyTorch backend for curvature calculations |
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

The app runs the core interface calculation, surface curvature, and optional
electrostatics in successive background stages. Contact tables become available
before curvature finishes, and convexity becomes available while electrostatics
is still calculating. Pocket jobs are submitted explicitly from a selected interface after
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
- `3D protein surface` defaults to the signed mean curvature `H` from the local
  Zernike fit. Positive `H` is convex/outward and negative `H` is
  concave/inward. The 4, 6, and 8 Å neighborhoods are retained in the result;
  the resolution selector changes the displayed scale. The **Surface coloring**
  toggle switches between the pink/green convexity map and an optional
  blue/white/red APBS electrostatic-potential map. Concave regions are pink and
  convex regions are green; negative potential is blue and positive potential
  is red, with neutral values near white. Electrostatics is reported in `kT/e`.
  The toggle is shared by Interactive 3D and both
  Binding faces panes. Surface opacity,
  distance padding from the selected contact face (0–8 Å), separation of the
  opposing chains, and per-chain visibility are view-only controls.
  Selected residues appear as yellow surface patches: every triangle touching
  a vertex assigned to a selected residue is colored yellow, including faces
  outside the contact padding. The patches follow the surface opacity control
  and preserve neighboring faces' colors and smooth shading in either mode.
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
  camera actions are disabled. Surface coloring, curvature resolution, opacity,
  and contact padding apply to both views. Hover tooltips show the local
  curvature or electrostatic potential for the selected coloring mode.
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

Curvature automatically uses CUDA when CUDA-enabled PyTorch and a working
NVIDIA GPU are available to the Python worker. The accelerated stage samples
the existing distance field, evaluates the same Zernike stencil, and calculates
principal curvatures in double precision (`float64`) at all three neighborhood
scales. Each chain's field is uploaded once and processed in bounded batches.
Surface generation, contact padding, residue assignment, and summary statistics
keep their existing CPU calculations. The surface geometry and browser rendering
use the same result format.

Install the optional backend in the environment used to run the server:

```bash
python -m pip install -r requirements-cuda.txt
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

The availability check must print `True` to use CUDA. An existing CUDA-enabled
PyTorch installation is reused. If the installed wheel lacks CUDA support,
choose the appropriate CUDA wheel using the
[PyTorch installation instructions](https://pytorch.org/get-started/locally/).
Restart the server after installing it, and analyze the coordinates again to
calculate new results; opening saved results does not recalculate curvature.

`INTERFACE_APP_CURVATURE_BACKEND` selects `auto` (the default), `cpu`, or `cuda`:

```bash
INTERFACE_APP_CURVATURE_BACKEND=cuda python run_interface_app.py
```

Both `auto` and `cuda` fall back to the CPU reference if CUDA cannot initialize
or a GPU calculation fails. An explicit `cuda` request logs initialization
failures; GPU calculation failures are also logged. The saved curvature JSON's
`parameters.acceleration` records the requested backend, actual backend for
each chain, device name, precision, and any fallback reason. `backend: "mixed"`
means an earlier chain completed on CUDA before a later chain fell back to CPU.
The base runtime dependencies remain sufficient for CPU calculation.

Electrostatics uses PDB2PQR to prepare the full retained protein assembly and
APBS to solve the linearized Poisson–Boltzmann equation. The potential is
interpolated onto the same vertices used by both surface viewers. The fixed
defaults are AMBER charges with PROPKA protonation at pH 7.4, protein dielectric
2, solvent dielectric 78.5, 0.15 M monovalent ions, 298.15 K, and a 1.4 Å solvent
radius. A single symmetric color range covers all chains and clips display
colors at the 98th percentile of absolute potential. Hover values remain
unclipped. Curvature summaries, tables, and contact-map shape marks continue
to show curvature when the surface is colored by potential.

On Linux x86_64, missing tools are downloaded automatically on the first new
analysis from the pinned, SHA256-verified artifacts in
`interface_app/tool_manifest.json`: APBS 3.4.1, PDB2PQR 3.6.2, and PDB2PQR's
Python dependencies. The app caches them under `<data-dir>/tools` and reuses
them for later analyses. Existing `apbs`, `pdb2pqr30`, or `pdb2pqr` executables
on `PATH` take precedence over downloads; `APBS_BIN` and `PDB2PQR_BIN` override
discovery with explicit executable paths. Other platforms require installed
tools. For example:

```bash
APBS_BIN=/path/to/apbs PDB2PQR_BIN=/path/to/pdb2pqr python run_interface_app.py
```

The Electrostatics option stays disabled until a matching potential map is
ready. If tool installation or the solver fails, its status explains the
failure and offers **Retry** while convexity and the core results remain
usable. A failed result fetch offers **Reload**. Cancelling or retrying the
optional job leaves the interface analysis complete. Electrostatic values and
calculation parameters are saved in a separate versioned JSON artifact,
available through `GET /api/analyses/{analysis_id}/electrostatics` and included
in the analysis export. Saved analyses without an electrostatics result explain
that the coordinates need to be analyzed again.

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

CUDA numerical checks are opt-in and require access to a real GPU. They compare
the production CUDA backend against the CPU reference for planes, spheres,
cavities, saddles, invalid gradients, every curvature scale, grid boundaries,
and batches with incomplete tails:

```bash
INTERFACE_APP_TEST_CUDA=1 python -m pytest -q interface_app/tests/test_curvature_cuda.py
```

The optional browser smoke test requires Playwright and a Chrome browser. It
uses `examples/1IAR.pdb`, a converted mmCIF, a temporary API server, and real pocket jobs
to check typed residue highlighting, validation, collapse/expand, yellow residue
faces and boundaries, preserved normals, linked selections and picking,
surface controls, shared convexity/electrostatics colors and tooltips in all
three panes, optional-job polling, unavailable results and reload, anchor and
residue-list pocket jobs, pocket transparency and cancellation, delayed responses,
and mobile layout. Browser potential values are synthetic fixtures; numerical
tests cover OpenDX parsing, interpolation, tool provisioning, failure handling,
and job lifecycle independently. It saves screenshots and logs under the
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
