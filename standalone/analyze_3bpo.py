#!/usr/bin/env python3
"""
Interface + target analysis for the IL-13 / IL-4Ra / IL-13Ra1 ternary complex
(PDB 3BPO, 3.00 A).  Chains: A = IL-13, B = IL-4Ra, C = IL-13Ra1.

The IL-13 system assembles in the REVERSE order to IL-4:
  IL-13 -> IL-13Ra1   is the high-affinity "driver"  (30 nM)
  IL-13/IL-13Ra1 -> IL-4Ra is the "trigger"          (20 nM)
  IL-13 -> IL-4Ra alone: no measurable affinity
so the weak-link logic used for IL-4 does NOT transfer, and has to be redone.

Outputs, per pairwise interface:
  BSA, interface residues, dSASA hotspots, H-bonds, salt bridges,
  epitope decomposition onto IL-13 secondary structure, and a check of which
  IL-13 element could be replaced by a peptide.
"""

import copy
import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, ShrakeRupley

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PDB = PROJECT_ROOT / "examples" / "3BPO.pdb"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "standalone"
CUT, TIGHT, HB, SALT = 5.0, 4.0, 3.5, 4.0

NAMES = {"A": "IL-13", "B": "IL-4Ra", "C": "IL-13Ra1"}

AA3TO1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}

ANION = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
CATION = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "HIS": ["ND1", "NE2"]}

# IL-13 secondary structure from 3BPO HELIX records
IL13_ELEMENTS = [
    ("helix A", 6, 21),
    ("loop AB", 22, 27),
    ("mini-helix A'", 28, 32),
    ("loop", 33, 43),
    ("helix B", 44, 50),
    ("loop BC", 51, 56),
    ("helix C", 57, 68),
    ("loop CD", 69, 90),
    ("helix D", 91, 110),
    ("C-term", 111, 112),
]


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


def sasa_of(model, chains):
    sub = copy.deepcopy(model)
    for ch in list(sub):
        if ch.get_id() not in chains:
            sub.detach_child(ch.get_id())
    ShrakeRupley(probe_radius=1.40, n_points=960).compute(sub, level="R")
    return {(ch.get_id(), r.get_id()[1]): r.sasa for ch in sub for r in ch}


def pair_contacts(model, c1, c2):
    a1 = [(r, a) for r in model[c1] for a in heavy(r)]
    a2 = [(r, a) for r in model[c2] for a in heavy(r)]
    x1 = np.array([a.coord for _, a in a1]); x2 = np.array([a.coord for _, a in a2])
    d = np.linalg.norm(x1[:, None] - x2[None], axis=-1)
    pairs = defaultdict(lambda: {"n": 0, "tight": 0, "min": 99.0, "hb": [], "sb": []})
    for i, j in np.argwhere(d <= CUT):
        r1, at1 = a1[i]; r2, at2 = a2[j]
        dist = float(d[i, j])
        k = (r1.get_id()[1], r1.get_resname(), r2.get_id()[1], r2.get_resname())
        p = pairs[k]
        p["n"] += 1
        p["tight"] += dist <= TIGHT
        p["min"] = min(p["min"], dist)
        if at1.element in "NO" and at2.element in "NO" and dist <= HB:
            p["hb"].append((at1.get_id(), at2.get_id(), round(dist, 2)))
        if dist <= SALT:
            n1, n2 = r1.get_resname(), r2.get_resname()
            if (at1.get_id() in ANION.get(n1, []) and at2.get_id() in CATION.get(n2, [])) or \
               (at1.get_id() in CATION.get(n1, []) and at2.get_id() in ANION.get(n2, [])):
                p["sb"].append((at1.get_id(), at2.get_id(), round(dist, 2)))
    return pairs


def main():
    m = clean(PDB)
    chains = [c.get_id() for c in m]
    print(f"3BPO chains: {[(c, NAMES[c]) for c in chains]}\n")

    iso = {c: sasa_of(m, [c]) for c in chains}
    report = {}

    for c1, c2 in itertools.combinations(chains, 2):
        nm = f"{NAMES[c1]}({c1}) - {NAMES[c2]}({c2})"
        pairs = pair_contacts(m, c1, c2)
        bound = sasa_of(m, [c1, c2])

        per = defaultdict(lambda: {"n": 0, "tight": 0, "min": 99.0, "hb": 0,
                                   "sb": 0, "p": []})
        for (r1, n1, r2, n2), p in pairs.items():
            for side, rid, rn, prid, prn in ((c1, r1, n1, r2, n2), (c2, r2, n2, r1, n1)):
                e = per[(side, rid, rn)]
                e["n"] += p["n"]; e["tight"] += p["tight"]
                e["min"] = min(e["min"], p["min"])
                e["hb"] += len(p["hb"]); e["sb"] += len(p["sb"])
                e["p"].append(f"{AA3TO1[prn]}{prid}")

        rows = []
        for (side, rid, rn), e in per.items():
            free = iso[side].get((side, rid), 0.0)
            ds = free - bound.get((side, rid), 0.0)
            rows.append({"chain": side, "resid": rid, "label": f"{AA3TO1[rn]}{rid}",
                         "dSASA": round(ds, 1), "n": e["n"], "tight": e["tight"],
                         "min": round(e["min"], 2), "hb": e["hb"], "sb": e["sb"],
                         "partners": sorted(set(e["p"]))})
        rows.sort(key=lambda r: -r["dSASA"])
        tot = sum(r["dSASA"] for r in rows if r["dSASA"] > 0)

        nhb = sum(len(p["hb"]) for p in pairs.values())
        nsb = sum(len(p["sb"]) for p in pairs.values())

        print("=" * 78)
        print(f"INTERFACE: {nm}")
        print("=" * 78)
        print(f"  BSA {tot/2:8.0f} A^2 | interface res {len(rows):3d} | "
              f"pairs {len(pairs):3d} | H-bonds {nhb:2d} | salt bridges {nsb:2d}")
        for side in (c1, c2):
            sub = [r for r in rows if r["chain"] == side and r["dSASA"] > 5]
            print(f"\n  --- {NAMES[side]} (chain {side}) epitope, {len(sub)} res dSASA>5 ---")
            print(f"  {'res':>6} {'dSASA':>7} {'ncont':>6} {'dmin':>6} {'HB':>3} {'SB':>3}  partners")
            for r in sub:
                print(f"  {r['label']:>6} {r['dSASA']:7.1f} {r['n']:6d} {r['min']:6.2f} "
                      f"{r['hb']:3d} {r['sb']:3d}  {','.join(r['partners'])}")
        print(f"\n  --- polar network ---")
        for (r1, n1, r2, n2), p in sorted(pairs.items(), key=lambda kv: kv[1]["min"]):
            if p["hb"] or p["sb"]:
                tag = "SALT" if p["sb"] else "HB  "
                det = p["sb"] or p["hb"]
                print(f"  {tag} {c1}:{AA3TO1[n1]}{r1:<5} .. {c2}:{AA3TO1[n2]}{r2:<5} "
                      + "; ".join(f"{a}-{b} {d}A" for a, b, d in det))
        print()
        report[f"{c1}{c2}"] = {"name": nm, "bsa": tot / 2, "rows": rows,
                              "nhb": nhb, "nsb": nsb}

    # ---------- IL-13 epitope decomposition ----------
    print("=" * 78)
    print("IL-13 EPITOPE DECOMPOSITION BY SECONDARY-STRUCTURE ELEMENT")
    print("=" * 78)
    dAB = {(r["chain"], r["resid"]): r["dSASA"] for r in report["AB"]["rows"]}
    dAC = {(r["chain"], r["resid"]): r["dSASA"] for r in report["AC"]["rows"]}
    t1 = sum(v for (c, i), v in dAB.items() if c == "A" and v > 0)
    t2 = sum(v for (c, i), v in dAC.items() if c == "A" and v > 0)
    print(f"  {'element':<16} {'range':>9} | {'vs IL-4Ra':>10} {'%':>5} | "
          f"{'vs IL-13Ra1':>12} {'%':>5}")
    for enm, lo, hi in IL13_ELEMENTS:
        s1 = sum(max(0, dAB.get(("A", i), 0)) for i in range(lo, hi + 1))
        s2 = sum(max(0, dAC.get(("A", i), 0)) for i in range(lo, hi + 1))
        print(f"  {enm:<16} {f'{lo}-{hi}':>9} | {s1:10.1f} {100*s1/t1:5.1f} | "
              f"{s2:12.1f} {100*s2/t2:5.1f}")
    print(f"  {'TOTAL':<16} {'':>9} | {t1:10.1f}       | {t2:12.1f}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / "interface_analysis_3BPO.json"
    with output.open("w") as fh:
        json.dump({k: {"name": v["name"], "bsa": v["bsa"], "nhb": v["nhb"],
                       "nsb": v["nsb"], "rows": v["rows"]} for k, v in report.items()},
                  fh, indent=2)
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
