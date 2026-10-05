# Standalone contact maps

This folder retains the contact-map implementation that informed the browser
app. Its scripts use the shared Zernike module in the repository root, the PDB
files in `examples/`, and the NGL library and license in `interface_app/static/`.
They do not depend on a local conda installation or the original workspace.

Install `requirements.txt`, then run from the repository root:

```bash
python standalone/contact_map.py
python standalone/contact_map.py --complex il13
```

These commands generate `contact_map.html` or `contact_map_IL13.html` and their
surface reports under `outputs/standalone/`. Fitting all surface vertices at
three neighborhood sizes can take several minutes. Each HTML page includes
its coordinates, data, and NGL viewer and can be opened directly in a browser.
Use `--output-dir PATH` to choose another output location.

The input interface-analysis snapshots are in `reference_data/`. To recompute
those inputs with the original analysis scripts:

```bash
python standalone/analyze_interfaces.py
python standalone/analyze_3bpo.py
python standalone/contact_map.py --analysis-dir outputs/standalone
python standalone/contact_map.py --complex il13 --analysis-dir outputs/standalone
```

The analysis scripts write their JSON outputs to `outputs/standalone/`, keeping
the committed reference inputs available for regression comparisons.

`reference_results/` contains the existing generated IL-4 and IL-13 contact-map
pages, Zernike JSON reports, contact-map preview, and text reports copied from
the research workspace. They are packaged for local reference and ignored by
Git by default because they are generated artifacts. The scripts, template,
and `reference_data/` inputs are included in a normal source commit.
