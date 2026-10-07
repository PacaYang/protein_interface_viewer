# Plan: residue highlight panel and residue-defined pockets

## Decisions

- Typed residues reuse the existing yellow selection (`selectResidues`), so they stay linked to the tables, contact map, and ball+stick reps. Any later click replaces them.
- A residue-defined pocket searches the union of R-Å spheres around every heavy atom of the listed residues. Each result reports how many of those residues line the detected pocket.

## Shared: residue list parser (`app.js`)

`parseResidueSpec(text, defaultChain)` returns `{keys, unmatched}`.

- Tokens are separated by commas or whitespace, e.g. `A:45-52, 60 63A, B:101`.
- A `X:` chain prefix applies to its own token and every later token until the next prefix. Tokens without a prefix use `defaultChain`.
- Single residues match `-?\d+[A-Za-z]?`. Ranges use `start-end` and are expanded by position in `chainResidues(chain)`, so insertion codes and sequence gaps are handled. Negative numbers work (`-3-5`).
- Each token is validated against `state.result.metadata.chains`. Unknown chains or residues go into `unmatched`, and nothing is guessed.
- Returned keys use `residueKey()` format (`A:63A`), which `selectResidues`, `residueLabelFromKey`, and the mesh code already accept.

## Feature 1: collapsible highlight panel beside `surface-view`

### `index.html`
Wrap `#surface-view` in a `div.surface-stage` that holds:
- `#surface-view` (unchanged id, so the NGL stage and tooltip code still work)
- `<aside id="residue-highlight-panel" class="surface-side">`, containing:
  - toggle `#residue-highlight-toggle` (`aria-expanded`, `aria-controls="residue-highlight-body"`), using the same `−`/`+` pattern as `#interfaces-toggle`
  - `#residue-highlight-body`: a `<form id="residue-highlight-form">` with a default-chain select `#highlight-chain`, a text input `#highlight-residues` (placeholder `A:45-52, 60, B:101`), and **Highlight** / **Clear** buttons
  - `#highlight-feedback` (`aria-live="polite"`): matched count, unmatched tokens, and residues with no surface vertices at the current scale

### `app.js`
- `initResidueHighlightPanel()`: bind submit, clear, and toggle. Store the collapsed state in `localStorage` (`highlightPanelCollapsed`), as the drawer does. After a toggle, call `state.stage.handleResize()` so the canvas resizes.
- `applyResidueHighlight()`: parse with `#highlight-chain`, call `selectResidues(keys)`, then `renderHighlightFeedback()`.
- `renderHighlightFeedback(keys, unmatched)`: use `wholeSurface(chain).residues` to list matched residues that have no surface patch (buried). Before curvature finishes, say that only ball+stick is shown.
- **Clear** calls `selectResidues([])` and keeps the text so it can be applied again.
- Fill `#highlight-chain` with every chain from `renderStructureMetadata()`/`renderSurfaceControls()`, defaulting to the pair's `chain_a`.
- `selectPair()` already clears `selectedResidues`. Leave the input text and reset the feedback to "Not applied to this interface".
- The panel is a view-only control and needs no backend change. The mesh rebuild in `installSurfaceMeshes()` already reacts to `selectedResidues`.

### `style.css`
- `.surface-stage { display: grid; grid-template-columns: minmax(0,1fr) 250px; gap: 12px; }`
- `.surface-stage.side-collapsed { grid-template-columns: minmax(0,1fr) 36px; }`, which hides the body and leaves the toggle rail
- At `max-width: 850px`, use a single column with the panel below the view
- Reuse the existing panel, field, and button tokens

## Feature 2: residue-defined pockets

### `pocket.py`
- `_search_mask(origin, shape, grid_A, centers, radius_A)`: builds a union-of-spheres mask by running a `cKDTree(centers).query(points, distance_upper_bound=radius_A)` over the grid points in the centers' bounding box ± R. The current anchor mode passes one center, so its results stay the same. The inline anchor-mask code in `_state_result` goes away.
- `_state_result(...)` takes `centers` and `radius_A` instead of `anchor` and `anchor_radius_A`. `anchor_distance_A` becomes the minimum distance to any center.
- `_residue_set(model, chains, residues)`: resolves each `{chain_id, number, insertion_code}` with the logic in `_anchor_residue`. It raises `ValueError` for missing residues, an empty list, or more than 60 residues. Residues may come from the target or partner chain, because the region is spatial and interface pockets can be lined by both.
- Coverage per pocket, when residues are defined:
  `defined_residues_lining` (keys found in `lining_residues`), `defined_residue_coverage` (`n_lining / n_defined_present_in_state`).
  In the free state, partner residues are left out of the denominator.
- Primary pocket in residue mode is ranked by highest coverage, then volume. Anchor mode keeps the largest-first ranking.
- `analyze_pocket(model, target_chain, partner_chain, anchor_residue=None, radius_A=8.0, *, pocket_residues=None, ...)`:
  exactly one of `anchor_residue` or `pocket_residues` must be given. The output adds `mode` (`"anchor"`/`"residues"`), `pocket_residues` (identities), and `anchor_coord` (centroid in residue mode). The 4–12 Å bound stays the same.

### `app.py`
- `ResidueRef(BaseModel)`: `chain_id: str | None = None`, `number: int`, `insertion_code: str = ""`. `chain_id` defaults to `target_chain`.
- `PocketRequest`: `anchor_residue: AnchorRequest | None = None` and `pocket_residues: list[ResidueRef] | None = Field(None, min_length=1, max_length=60)`. A `model_validator` enforces exactly one. Dump with `exclude_none=True`.
- `jobs._run_pocket` passes `**request`, so it works as is. Stored anchor-only requests still retry correctly.

### Pocket panel (`index.html` / `app.js`)
- Mode radios `name="pocket-mode"`: **Anchor residue** (current) and **Residue list**.
- Residue-list mode shows `#pocket-residues` (same parser, default chain = target) and a **Use current selection** button (`#pocket-use-selection`). That button fills the field from `state.selectedResidues`, so feature 1 feeds into feature 2.
- The radius label changes to "Distance from residues". The default is 6 Å in residue mode and 8 Å in anchor mode.
- `submitPocket()` builds `pocket_residues` from the parsed keys. Unmatched tokens block submission and show an inline error.
- `renderPocketCard()` adds a "Defined residues lining" row (`n / m`). The result header lists the defined residues. Mesh highlight, cancel, and polling behavior stay the same.
- Update the `#pocket-result` placeholder text in `selectPair()` to cover both modes.

## Tests
- `test_pocket.py`
  - `_search_mask` with one center equals the old single-sphere mask.
  - Residue mode on `examples/1IAR.pdb` (A residues around 67) returns a free/bound primary, `mode == "residues"`, and coverage fields in [0, 1].
  - `ValueError` when both or neither of anchor/residues are given, for a missing residue, and for more than 60 residues. A partner-chain residue is accepted.
- `test_api.py`: the `PocketRequest` validator accepts each mode alone and rejects both or none.
- `test_ui_contract.py`: new ids (`residue-highlight-panel`, `residue-highlight-toggle`, `highlight-residues`, `pocket-residues`, `pocket-use-selection`) and functions (`parseResidueSpec`, `applyResidueHighlight`).
- `browser_smoke.py`: type a range into the panel, confirm yellow faces appear and the table rows are selected, collapse and expand, and run one residue-list pocket job.

## Docs
In the README "Workspace behavior" section, add the highlight panel and its syntax. In the pocket definitions, describe residue-list mode, union-of-spheres search, coverage, and primary ranking.

## Verification
```bash
pytest -q
node --check interface_app/static/app.js
python -m py_compile interface_app/*.py
python interface_app/tests/browser_smoke.py --browser /usr/bin/google-chrome --artifacts /tmp/interface-app-browser
```
