#!/usr/bin/env python3
"""
Interactive contact maps and 3D Zernike interface convexity for PDB 3BPN/3BPO.

Rows/columns = residues of the selected chain pair (AB, AC, or BC)
Cell    = minimum heavy-atom distance between the two residues

Three variables are encoded:
  * AA type  -> one-letter code + number in the axis label, plus a class shape
  * distance -> the cell fill, a 5-bin sequential blue ramp
  * charge   -> a diverging orange/aqua band on each axis (per-residue formal charge)
                plus per-cell glyphs for salt bridges and like-charge proximity

Writes a self-contained HTML page (data, PDB, and NGL inlined, opens under
file://), a numeric surface-convexity JSON report, and verification to stdout.

Regenerate from the repository root with the app's Python dependencies:
  python standalone/contact_map.py
  python standalone/contact_map.py --complex il13

Generated pages and surface reports go to outputs/standalone/ by default.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = PROJECT_ROOT / "examples"
REFERENCE_DATA = Path(__file__).resolve().parent / "reference_data"
sys.path.insert(0, str(PROJECT_ROOT))

PAIRS = ("AB", "AC", "BC")

CONTACT_CUTOFF = 6.0   # grid inclusion; the house scripts use 5.0
TIGHT_CUTOFF = 4.0     # packing contacts
HBOND_CUTOFF = 3.5     # N/O ... N/O
SALT_CUTOFF = 4.0      # charged group atom .. charged group atom
LIKE_CHARGE_CUTOFF = 6.0   # like-charge functional groups in proximity

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

ANION_ATOMS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
CATION_ATOMS = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "HIS": ["ND1", "NE2"]}

# Formal side-chain charge at pH 7.4. HIS is 0 here -- side-chain pKa ~6.0, so only
# ~4% protonated at 7.4 -- even though CATION_ATOMS above includes HIS for atom-level
# salt-bridge geometry. That asymmetry is deliberate: the atom-level convention is what
# produced the salt-bridge counts in interface_analysis.json that this script verifies
# against, while the band and the +/- glyph report the actual charge state. Do not
# "fix" one half without the other.
CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}

# Class for residues that carry no formal charge. Charged residues are shown by the
# charge band instead, so they are absent here by design.
RES_CLASS = {
    "ASN": "polar", "GLN": "polar", "SER": "polar", "THR": "polar", "HIS": "polar",
    "ALA": "hyd", "VAL": "hyd", "LEU": "hyd", "ILE": "hyd", "MET": "hyd",
    "PHE": "arom", "TRP": "arom", "TYR": "arom",
    "GLY": "spec", "PRO": "spec", "CYS": "spec",
}

# IL-4 secondary structure, from the 3BPN HELIX records (validate_target.py:33-44).
IL4_ELEMENTS = [
    ("helix A", 9, 21),
    ("loop AB", 22, 44),
    ("helix B", 45, 64),
    ("loop BC", 65, 73),
    ("helix C", 74, 96),
    ("loop CD", 97, 108),
    ("helix D", 109, 127),
]

# From the 3BPO HELIX records, using the same element grouping as analyze_3bpo.py.
IL13_ELEMENTS = [
    ("helix A", 6, 21), ("loop AB", 22, 27), ("mini-helix A'", 28, 32),
    ("loop", 33, 43), ("helix B", 44, 50), ("loop BC", 51, 56),
    ("helix C", 57, 68), ("loop CD", 69, 90), ("helix D", 91, 110),
    ("C-term", 111, 112),
]

COMPLEXES = {
    "il4": {
        "pdb": "3BPN.pdb", "analysis": "interface_analysis.json",
        "html": "contact_map.html", "surface_json": "surface_convexity.json",
        "title": "IL-4 ternary complex", "resolution": 3.02,
        "chain_names": {"A": "IL-4", "B": "IL-4Rα", "C": "IL-13Rα1"},
        "elements": IL4_ELEMENTS,
        "expected": {"AB": (41, 67, 25, 23), "AC": (55, 88, 33, 27),
                     "BC": (36, 64, 23, 24)},
    },
    "il13": {
        "pdb": "3BPO.pdb", "analysis": "interface_analysis_3BPO.json",
        "html": "contact_map_IL13.html", "surface_json": "surface_convexity_IL13.json",
        "title": "IL-13 ternary complex", "resolution": 3.00,
        "chain_names": {"A": "IL-13", "B": "IL-4Rα", "C": "IL-13Rα1"},
        "elements": IL13_ELEMENTS,
        "expected": {"AB": (39, 60, 23, 19), "AC": (49, 72, 25, 23),
                     "BC": (36, 61, 23, 23)},
    },
}

# IL-13Ra1 target sub-patches (design_peptide.py:25-26).
D3_PATCH = [223, 254, 256, 257, 318, 319, 320, 321, 322]
D1_PATCH = [65, 74, 75, 76, 77, 78, 79, 80, 81, 104]


def clean(p):
    m = PDBParser(QUIET=True).get_structure("x", p)[0]
    for ch in m:
        for r in list(ch):
            if r.get_id()[0] != " " or r.get_resname() not in AA3TO1:
                ch.detach_child(r.get_id())
    for ch in list(m):
        if not len(list(ch)):
            m.detach_child(ch.get_id())
    return m


def heavy(r):
    return [a for a in r if a.element != "H" and a.get_altloc() in (" ", "A")]


def res_list(model, chain):
    """Residues of `chain` in sequence order.

    Sorted on the full (resid, icode) tuple rather than resid alone -- 3BPN has no
    insertion codes in A or C, but this keeps the ordering right elsewhere.
    """
    return sorted(model[chain], key=lambda r: (r.get_id()[1], r.get_id()[2]))


def label(res):
    return f"{AA3TO1[res.get_resname()]}{res.get_id()[1]}"


def min_dist_matrix(rows, cols):
    """(len(rows), len(cols)) array of minimum heavy-atom distances, np.inf where none.

    Initialised to inf rather than the house 99.0 sentinel so a non-contact can never
    leak into the colour scale as a real value.
    """
    a1 = [(i, a) for i, r in enumerate(rows) for a in heavy(r)]
    a2 = [(j, a) for j, r in enumerate(cols) for a in heavy(r)]

    # A residue contributing zero heavy atoms (e.g. all atoms at altloc B) would leave
    # a row/column silently untouched below. Catch it rather than trust the structure.
    assert len({i for i, _ in a1}) == len(rows), "row residue with no heavy atoms"
    assert len({j for j, _ in a2}) == len(cols), "column residue with no heavy atoms"

    x1 = np.array([a.coord for _, a in a1])
    x2 = np.array([a.coord for _, a in a2])
    ci = np.array([j for j, _ in a2])

    D = np.full((len(rows), len(cols)), np.inf)
    for k, (i, _) in enumerate(a1):
        d = np.linalg.norm(x2 - x1[k], axis=-1)
        np.minimum.at(D[i], ci, d)
    return D


def select_grid(D, cutoff):
    """Row/column indices of residues making at least one contact, in sequence order."""
    return np.where((D <= cutoff).any(axis=1))[0], np.where((D <= cutoff).any(axis=0))[0]


def seq_breaks(residues, idx):
    """Positions in the filtered axis where the PDB numbering jumps.

    Chain C is modelled 30-338 but only 295 of 309 positions are present, and the
    contacting columns jump 106 -> 199. Marking the breaks stops the axis from
    implying sequence continuity it does not have.
    """
    nums = [residues[i].get_id()[1] for i in idx]
    return [k for k in range(1, len(nums)) if nums[k] - nums[k - 1] > 1]


def pair_features(rows, cols, row_idx, col_idx):
    """dict[(i, j)] -> atom-level features for every pair in the selected sub-grid.

    Adapted from pairwise_contacts (analyze_interfaces.py:85-122), with the cutoff
    raised to CONTACT_CUTOFF and a like-charge test added.

    Both charge tests are deliberately atom-level. Gating them on the residue
    min-distance overcounts badly: two lysines whose backbones pass within 6 A are not
    repelling if their NZ atoms point apart. Residue-level logic gives 8 attractive and
    6 repulsive pairs here, against the 3 real salt-bridge pairs.
    """
    sub_r = [(i, rows[i]) for i in row_idx]
    sub_c = [(j, cols[j]) for j in col_idx]
    a1 = [(i, res, a) for i, res in sub_r for a in heavy(res)]
    a2 = [(j, res, a) for j, res in sub_c for a in heavy(res)]

    x1 = np.array([a.coord for _, _, a in a1])
    x2 = np.array([a.coord for _, _, a in a2])
    d = np.linalg.norm(x1[:, None, :] - x2[None, :, :], axis=-1)

    feats = defaultdict(lambda: {
        "n_contacts": 0, "n_tight": 0, "hbonds": [], "saltbridges": [], "like": [],
    })

    for k, l in np.argwhere(d <= CONTACT_CUTOFF):
        i, res1, at1 = a1[k]
        j, res2, at2 = a2[l]
        dist = float(d[k, l])
        f = feats[(i, j)]
        f["n_contacts"] += 1
        if dist <= TIGHT_CUTOFF:
            f["n_tight"] += 1

        if at1.element in ("N", "O") and at2.element in ("N", "O") and dist <= HBOND_CUTOFF:
            f["hbonds"].append((at1.get_id(), at2.get_id(), round(dist, 2)))

        n1, n2 = res1.get_resname(), res2.get_resname()
        a1_anion = at1.get_id() in ANION_ATOMS.get(n1, [])
        a1_cation = at1.get_id() in CATION_ATOMS.get(n1, [])
        a2_anion = at2.get_id() in ANION_ATOMS.get(n2, [])
        a2_cation = at2.get_id() in CATION_ATOMS.get(n2, [])

        if dist <= SALT_CUTOFF and ((a1_anion and a2_cation) or (a1_cation and a2_anion)):
            f["saltbridges"].append((at1.get_id(), at2.get_id(), round(dist, 2)))

        if dist <= LIKE_CHARGE_CUTOFF and ((a1_anion and a2_anion) or (a1_cation and a2_cation)):
            f["like"].append((at1.get_id(), at2.get_id(), round(dist, 2)))

    return feats


def normalize_analysis(raw):
    """Give both historical interface reports one small, stable page schema."""
    report = {}
    for pair in PAIRS:
        source = raw[pair]
        residues = []
        for row in source.get("residues", source.get("rows", [])):
            residues.append({
                "chain": row["chain"], "label": row["label"], "dSASA": row["dSASA"],
                "partners": row["partners"],
                "n_contacts": row.get("n_contacts", row.get("n")),
                "n_tight": row.get("n_tight", row.get("tight")),
                "min_dist": row.get("min_dist", row.get("min")),
                "n_hbond": row.get("n_hbond", row.get("hb")),
                "n_saltbridge": row.get("n_saltbridge", row.get("sb")),
            })
        report[pair] = {
            "buried_surface_area_A2": source.get("buried_surface_area_A2", source.get("bsa")),
            "n_residue_pairs": source.get("n_residue_pairs"),
            "n_interface_residues": len(residues),
            "n_hbond_like": source.get("n_hbond_like", source.get("nhb")),
            "n_saltbridge": source.get("n_saltbridge", source.get("nsb")),
            "residues": residues,
        }
    return report


def load_dsasa(report, pair):
    """Pairwise dSASA and residues also buried against the third chain.

    Read rather than recomputed: ShrakeRupley on deep-copied sub-models is the slow
    path. The JSON was computed at a 5.0 A cutoff, so a few residues in the 6.0 A grid
    have no dSASA -- they are rendered as an open mark, never as a filled zero.
    """
    ds = {(r["chain"], r["label"]): r["dSASA"] for r in report[pair]["residues"]}
    third = ({"A", "B", "C"} - set(pair)).pop()
    other = {}
    for chain in pair:
        other_pair = "".join(sorted((chain, third)))
        other[chain] = {r["label"] for r in report[other_pair]["residues"]
                        if r["chain"] == chain and r["dSASA"] > 5}
    return ds, other, third


def element_of(resid, elements):
    for name, lo, hi in elements:
        if lo <= resid <= hi:
            return name
    return ""


def patch_of(resid):
    if resid in D3_PATCH:
        return "D3"
    if resid in D1_PATCH:
        return "D1"
    return ""


def axis_entry(res, chain, dsasa, other_epitope, third, config):
    name = res.get_resname()
    resid = res.get_id()[1]
    one = AA3TO1[name]
    q = CHARGE.get(name, 0)
    return {
        "label": f"{one}{resid}",
        "resid": resid,
        "aa1": one,
        "aa3": name,
        "charge": q,
        "cls": "pos" if q > 0 else "neg" if q < 0 else RES_CLASS[name],
        "dsasa": dsasa.get((chain, f"{one}{resid}")),
        "element": element_of(resid, config["elements"]) if chain == "A" else "",
        "patch": patch_of(resid) if chain == "C" else "",
        "also_other": config["chain_names"][third] if f"{one}{resid}" in other_epitope[chain] else "",
    }


def build_payload(rows, cols, D, pair, dsasa, other_epitope, third, config):
    """The dict inlined into the HTML.

    Ships the full 6.0 A grid plus every cell's atom-level detail, so the in-page
    cutoff filter re-slices client-side without needing the structure again.
    """
    row_idx, col_idx = select_grid(D, CONTACT_CUTOFF)
    feats = pair_features(rows, cols, row_idx, col_idx)

    cells = []
    for (i, j), f in feats.items():
        cells.append({
            "r": int(np.where(row_idx == i)[0][0]),
            "c": int(np.where(col_idx == j)[0][0]),
            "d": round(float(D[i, j]), 2),
            "n": f["n_contacts"],
            "tight": f["n_tight"],
            "hb": [list(h) for h in f["hbonds"]],
            "sb": [list(s) for s in f["saltbridges"]],
            "like": [list(x) for x in f["like"]],
        })
    cells.sort(key=lambda c: c["d"])

    row_chain, col_chain = pair
    return {
        "pdb": config["pdb"],
        "pair": pair,
        "cutoff": CONTACT_CUTOFF,
        "tight": TIGHT_CUTOFF,
        "rowChain": row_chain,
        "colChain": col_chain,
        "rowName": config["chain_names"][row_chain],
        "colName": config["chain_names"][col_chain],
        "rows": [axis_entry(rows[i], row_chain, dsasa, other_epitope, third, config) for i in row_idx],
        "cols": [axis_entry(cols[j], col_chain, dsasa, other_epitope, third, config) for j in col_idx],
        "rowBreaks": seq_breaks(rows, row_idx),
        "colBreaks": seq_breaks(cols, col_idx),
        "cells": cells,
    }


def verify(rows, cols, D, payload, analysis, config):
    """Assert the values established during planning, so a regression is visible."""
    pair = payload["pair"]
    lab_r = {label(r): i for i, r in enumerate(rows)}
    lab_c = {label(c): j for j, c in enumerate(cols)}

    print("=" * 78)
    print(f"VERIFICATION {pair}")
    print("=" * 78)

    n5 = int((D <= 5.0).sum())
    n6 = int((D <= CONTACT_CUTOFF).sum())
    ri, ci = select_grid(D, CONTACT_CUTOFF)
    reference_pairs = analysis["n_residue_pairs"]
    reference_note = (f"   ({config['analysis']}: {reference_pairs})"
                      if reference_pairs is not None else "")
    print(f"  pairs at 5.0 A                  : {n5:4d}{reference_note}")
    print(f"  pairs at {CONTACT_CUTOFF} A                  : {n6:4d}   grid "
          f"{len(ri)} x {len(ci)}")
    expected = config["expected"][pair]
    assert n5 == expected[0], f"{pair} 5.0 A pair count changed"
    if reference_pairs is not None:
        assert n5 == reference_pairs, f"{pair} report pair count disagrees"
    assert (n6, len(ri), len(ci)) == expected[1:], f"{pair} 6.0 A grid shape changed"
    at5 = [x for x in payload["cells"] if x["d"] <= 5.0]
    assert len(at5) == n5, f"{pair} payload is missing 5.0 A contacts"
    assert sum(len(x["hb"]) for x in at5) == analysis["n_hbond_like"], \
        f"{pair} H-bond-like count disagrees with the interface report"
    assert sum(len(x["sb"]) for x in at5) == analysis["n_saltbridge"], \
        f"{pair} salt-bridge atom count disagrees with the interface report"

    if config["pdb"] == "3BPO.pdb":
        expected_salt = {
            "AB": {("R65", "D67"), ("R65", "D72")},
            "AC": {("E91", "K76")}, "BC": set(),
        }[pair]
        got_salt = {(payload["rows"][x["r"]]["label"], payload["cols"][x["c"]]["label"])
                    for x in at5 if x["sb"]}
        assert got_salt == expected_salt, f"{pair} salt-bridge residue pairs changed"
        if pair == "AC":
            for a, b, want in (("R108", "C320", 2.64), ("E91", "K76", 2.74),
                               ("F107", "R256", 2.84), ("K104", "K318", 3.03)):
                got = round(float(D[lab_r[a], lab_c[b]]), 2)
                assert got == want, f"{a}-{b} is {got}, expected {want}"
        print(f"  BSA (from JSON)                 : {analysis['buried_surface_area_A2']:.0f} A^2")
        print(f"  H-bond-like / salt bridges      : {analysis['n_hbond_like']} / "
              f"{analysis['n_saltbridge']}   (atom pairs, at 5.0 A)")
        print("  all checks passed")
        return

    if pair != "AC":
        print("  all checks passed")
        return

    spot = {
        ("C127", "R256"): 2.47, ("Y124", "R256"): 2.53, ("S125", "R256"): 3.17,
        ("R121", "E322"): 2.77, ("R121", "K318"): 3.16, ("R121", "Y321"): 3.32,
        ("R121", "C320"): 3.98,
    }
    print("\n  spot distances")
    for (a, b), want in spot.items():
        got = round(float(D[lab_r[a], lab_c[b]]), 2)
        print(f"    {a:>5} .. {b:<5} {got:5.2f} A   expected {want:.2f}")
        assert got == want, f"{a}-{b} is {got}, expected {want}"

    gmin = float(D.min())
    gi, gj = np.unravel_index(D.argmin(), D.shape)
    print(f"\n  global minimum                  : {gmin:.2f} A at "
          f"{label(rows[gi])} .. {label(cols[gj])}")
    assert round(gmin, 2) == 2.47 and label(rows[gi]) == "C127" and label(cols[gj]) == "R256"

    # R88 is a known IL-4Ra hotspot. If it shows up here the chains are mixed up.
    r88_min = float(D[lab_r["R88"]].min())
    in_grid = "R88" in {r["label"] for r in payload["rows"]}
    print(f"  R88 nearest chain-C residue     : {r88_min:.2f} A   in grid: {in_grid}")
    assert not in_grid and round(r88_min, 2) == 16.61

    # R256 buries 118.9 A^2 and makes 6 H-bonds, but all seven of its partners are
    # uncharged -- a salt-bridge glyph there would mean the charge logic is wrong.
    r256_c = payload["cols"].index(next(c for c in payload["cols"] if c["label"] == "R256"))
    r256_sb = [x for x in payload["cells"] if x["c"] == r256_c and x["sb"]]
    print(f"  R256 salt-bridge cells          : {len(r256_sb)}   (expected 0; "
          f"dSASA {payload['cols'][r256_c]['dsasa']} A^2)")
    assert not r256_sb, "R256 should make no salt bridge"

    want_sb = {("R121", "E322"), ("K117", "E322"), ("R64", "E106")}
    got_sb = {(payload["rows"][x["r"]]["label"], payload["cols"][x["c"]]["label"])
              for x in payload["cells"] if x["sb"] and x["d"] <= 5.0}
    print(f"\n  salt bridges <=5.0 A            : {sorted(got_sb)}")
    assert got_sb == want_sb, f"salt bridges {got_sb}, expected {want_sb}"

    nlike = sum(1 for x in payload["cells"] if x["like"] and not x["sb"])
    print(f"  like-charge cells               : {nlike}")

    colnums = {c["resid"] for c in payload["cols"]}
    print(f"\n  D3 patch in columns             : {sorted(set(D3_PATCH) & colnums)}")
    print(f"  D1 patch in columns             : {sorted(set(D1_PATCH) & colnums)}")
    assert set(D3_PATCH) <= colnums and set(D1_PATCH) <= colnums, "patch coverage lost"

    no_dsasa = [r["label"] for r in payload["rows"] + payload["cols"] if r["dsasa"] is None]
    print(f"\n  in grid without dSASA ({len(no_dsasa):2d})       : {no_dsasa}")
    print(f"  BSA (from JSON)                 : {analysis['buried_surface_area_A2']:.0f} A^2")
    print(f"  H-bond-like / salt bridges      : {analysis['n_hbond_like']} / {analysis['n_saltbridge']}"
          f"   (atom pairs, at 5.0 A)")
    print("\n  all checks passed")


def write_html(payload, surface, analysis, config, output_dir):
    from html import escape

    ngl_dir = PROJECT_ROOT / "interface_app" / "static"
    ngl_js = (ngl_dir / "ngl.js").read_text()
    ngl_license = (ngl_dir / "License_for_NGL.txt").read_text()

    with Path(__file__).with_name("contact_map_template.html").open() as fh:
        tpl = fh.read()
    page = {"title": config["title"], "pdb": config["pdb"],
            "resolution": config["resolution"], "analysis_file": config["analysis"]}

    def inline(value):
        return json.dumps(value, separators=(",", ":"), ensure_ascii=False,
                          allow_nan=False).replace("<", "\\u003c")

    data = inline(payload)
    surface_data = inline(surface)
    analysis_data = inline(analysis)
    page_data = inline(page)
    # </script> inside a JSON string would end the block early; < is safe in JSON.
    pdb_data = inline((EXAMPLES / config["pdb"]).read_text())
    with (output_dir / config["html"]).open("w") as fh:
        fh.write(tpl.replace("__PAYLOAD__", data)
                 .replace("__SURFACE_DATA__", surface_data)
                 .replace("__ANALYSIS__", analysis_data)
                 .replace("__PAGE_DATA__", page_data)
                 .replace("__PAGE_TITLE__", escape(config["title"]))
                 .replace("__PDB_TEXT__", pdb_data)
                 .replace("__NGL_LICENSE__", ngl_license.replace("--", "—"))
                 .replace("__NGL_JS__", ngl_js))


def main():
    from zernike_convexity import compute

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", choices=COMPLEXES, default="il4",
                        help="ternary structure to render (default: il4)")
    parser.add_argument("--analysis-dir", type=Path, default=REFERENCE_DATA,
                        help="directory containing the interface-analysis JSON inputs")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs" / "standalone",
                        help="directory for generated HTML and surface reports")
    args = parser.parse_args()
    config = COMPLEXES[args.complex]
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    model = clean(EXAMPLES / config["pdb"])
    print(f"{config['pdb']} chains: "
          f"{[(c.get_id(), config['chain_names'][c.get_id()]) for c in model]}")
    with (args.analysis_dir / config["analysis"]).open() as fh:
        analysis = normalize_analysis(json.load(fh))
    payload = {}
    for pair in PAIRS:
        row_chain, col_chain = pair
        rows = res_list(model, row_chain)
        cols = res_list(model, col_chain)
        D = min_dist_matrix(rows, cols)
        dsasa, other_epitope, third = load_dsasa(analysis, pair)
        payload[pair] = build_payload(rows, cols, D, pair, dsasa, other_epitope, third,
                                      config)
        verify(rows, cols, D, payload[pair], analysis[pair], config)

    surface = compute(model, analysis, chain_names=config["chain_names"],
                      source=config["pdb"], source_resolution_A=config["resolution"])
    with (output_dir / config["surface_json"]).open("w") as fh:
        json.dump(surface["report"], fh, indent=2, ensure_ascii=False, allow_nan=False)
    write_html(payload, surface, analysis, config, output_dir)
    print(f"\nWrote {output_dir / config['html']}")


if __name__ == "__main__":
    main()
