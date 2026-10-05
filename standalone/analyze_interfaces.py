#!/usr/bin/env python3
"""
Interface analysis of the IL-4 / IL-4Ralpha / IL-13Ralpha1 ternary complex (PDB 3BPN).

For each pairwise interface (A-B, A-C, B-C) we compute:
  * heavy-atom contacts within CONTACT_CUTOFF
  * minimum interatomic distance per residue pair
  * polar contacts (H-bond-like N/O-N/O <= HBOND_CUTOFF) and salt bridges
  * per-residue buried SASA (dSASA = SASA_free - SASA_complex), Shrake-Rupley
  * total interface area and hotspot ranking

Chains: A = IL-4, B = IL-4Ralpha (ECD), C = IL-13Ralpha1 (ECD)
"""

import itertools
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Selection, ShrakeRupley

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PDB = PROJECT_ROOT / "examples" / "3BPN.pdb"
OUTPUT_DIR = PROJECT_ROOT / "outputs" / "standalone"
CONTACT_CUTOFF = 5.0   # heavy-atom, generous "in contact" shell
TIGHT_CUTOFF = 4.0     # packing contacts
HBOND_CUTOFF = 3.5     # N/O ... N/O
SALT_CUTOFF = 4.0      # charged group centroid atoms

CHAIN_NAMES = {
    "A": "IL-4",
    "B": "IL-4Ra",
    "C": "IL-13Ra1",
}

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

ANION_ATOMS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
CATION_ATOMS = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "HIS": ["ND1", "NE2"]}


def is_aa(res):
    return res.get_resname() in AA3TO1


def rlabel(res):
    """Chain-free label, e.g. GLU9 -> E9."""
    return f"{AA3TO1[res.get_resname()]}{res.get_id()[1]}"


def heavy_atoms(res):
    return [a for a in res if a.element != "H" and a.get_altloc() in (" ", "A")]


def load():
    st = PDBParser(QUIET=True).get_structure("3bpn", PDB)
    model = st[0]
    # strip waters / sugars / hetero for the protein-protein analysis
    for ch in model:
        for res in list(ch):
            if res.get_id()[0] != " " or not is_aa(res):
                ch.detach_child(res.get_id())
    return st, model


def residue_sasa(model, chains):
    """Shrake-Rupley SASA per residue for the sub-structure made of `chains`."""
    import copy

    sub = copy.deepcopy(model)
    for ch in list(sub):
        if ch.get_id() not in chains:
            sub.detach_child(ch.get_id())
    ShrakeRupley(probe_radius=1.40, n_points=960).compute(sub, level="R")
    return {
        (ch.get_id(), res.get_id()[1]): res.sasa
        for ch in sub
        for res in ch
    }


def pairwise_contacts(model, c1, c2):
    """Return dict[(res1,res2)] -> stats for every residue pair in contact."""
    a1 = [(res, a) for res in model[c1] for a in heavy_atoms(res)]
    a2 = [(res, a) for res in model[c2] for a in heavy_atoms(res)]
    x1 = np.array([a.coord for _, a in a1])
    x2 = np.array([a.coord for _, a in a2])

    d = np.linalg.norm(x1[:, None, :] - x2[None, :, :], axis=-1)
    idx = np.argwhere(d <= CONTACT_CUTOFF)

    pairs = defaultdict(lambda: {
        "n_contacts": 0, "n_tight": 0, "min_dist": 99.0,
        "hbonds": [], "saltbridges": [],
    })

    for i, j in idx:
        res1, at1 = a1[i]
        res2, at2 = a2[j]
        dist = float(d[i, j])
        key = (res1.get_id()[1], res1.get_resname(), res2.get_id()[1], res2.get_resname())
        p = pairs[key]
        p["n_contacts"] += 1
        if dist <= TIGHT_CUTOFF:
            p["n_tight"] += 1
        p["min_dist"] = min(p["min_dist"], dist)

        # polar / H-bond-like
        if at1.element in ("N", "O") and at2.element in ("N", "O") and dist <= HBOND_CUTOFF:
            p["hbonds"].append((at1.get_id(), at2.get_id(), round(dist, 2)))

        # salt bridge
        n1, n2 = res1.get_resname(), res2.get_resname()
        if dist <= SALT_CUTOFF:
            if (at1.get_id() in ANION_ATOMS.get(n1, []) and at2.get_id() in CATION_ATOMS.get(n2, [])) or \
               (at1.get_id() in CATION_ATOMS.get(n1, []) and at2.get_id() in ANION_ATOMS.get(n2, [])):
                p["saltbridges"].append((at1.get_id(), at2.get_id(), round(dist, 2)))

    return pairs


def analyze():
    st, model = load()
    chains = [c.get_id() for c in model]
    print(f"Chains present: {[(c, CHAIN_NAMES[c]) for c in chains]}\n")

    # SASA references
    sasa_iso = {c: residue_sasa(model, [c]) for c in chains}          # chain alone
    sasa_all = residue_sasa(model, chains)                            # full ternary complex
    sasa_pair = {}
    for c1, c2 in itertools.combinations(chains, 2):
        sasa_pair[(c1, c2)] = residue_sasa(model, [c1, c2])

    report = {}

    for c1, c2 in itertools.combinations(chains, 2):
        name = f"{CHAIN_NAMES[c1]}({c1}) - {CHAIN_NAMES[c2]}({c2})"
        pairs = pairwise_contacts(model, c1, c2)

        # per-residue aggregation
        per_res = defaultdict(lambda: {
            "partners": [], "n_contacts": 0, "n_tight": 0,
            "min_dist": 99.0, "n_hb": 0, "n_sb": 0,
        })
        for (r1, n1, r2, n2), p in pairs.items():
            for (side, rid, rn, prid, prn) in (
                (c1, r1, n1, r2, n2),
                (c2, r2, n2, r1, n1),
            ):
                e = per_res[(side, rid, rn)]
                e["partners"].append(f"{AA3TO1[prn]}{prid}")
                e["n_contacts"] += p["n_contacts"]
                e["n_tight"] += p["n_tight"]
                e["min_dist"] = min(e["min_dist"], p["min_dist"])
                e["n_hb"] += len(p["hbonds"])
                e["n_sb"] += len(p["saltbridges"])

        # dSASA using the isolated-chain vs the two-chain complex (pair-specific burial)
        rows = []
        for (side, rid, rn), e in per_res.items():
            free = sasa_iso[side].get((side, rid), 0.0)
            bound = sasa_pair[(c1, c2)].get((side, rid), 0.0)
            dsasa = free - bound
            rows.append({
                "chain": side,
                "protein": CHAIN_NAMES[side],
                "resid": rid,
                "resname": rn,
                "label": f"{AA3TO1[rn]}{rid}",
                "dSASA": round(dsasa, 1),
                "sasa_free": round(free, 1),
                "n_contacts": e["n_contacts"],
                "n_tight": e["n_tight"],
                "min_dist": round(e["min_dist"], 2),
                "n_hbond": e["n_hb"],
                "n_saltbridge": e["n_sb"],
                "partners": sorted(set(e["partners"])),
            })
        rows.sort(key=lambda r: (-r["dSASA"], -r["n_contacts"]))

        total_dsasa = sum(r["dSASA"] for r in rows if r["dSASA"] > 0)
        bsa = total_dsasa / 2.0

        n_hb = sum(len(p["hbonds"]) for p in pairs.values())
        n_sb = sum(len(p["saltbridges"]) for p in pairs.values())

        report[f"{c1}{c2}"] = {
            "name": name,
            "buried_surface_area_A2": round(bsa, 0),
            "total_dSASA_A2": round(total_dsasa, 0),
            "n_residue_pairs": len(pairs),
            "n_interface_residues": len(rows),
            "n_hbond_like": n_hb,
            "n_saltbridge": n_sb,
            "residues": rows,
            "pairs": [
                {
                    "r1": f"{c1}:{AA3TO1[n1]}{r1}",
                    "r2": f"{c2}:{AA3TO1[n2]}{r2}",
                    "min_dist": round(p["min_dist"], 2),
                    "n_contacts": p["n_contacts"],
                    "hbonds": p["hbonds"],
                    "saltbridges": p["saltbridges"],
                }
                for (r1, n1, r2, n2), p in sorted(
                    pairs.items(), key=lambda kv: kv[1]["min_dist"]
                )
            ],
        }

        # ---- print ----
        print("=" * 78)
        print(f"INTERFACE: {name}")
        print("=" * 78)
        print(f"  Buried surface area (BSA)      : {bsa:8.0f} A^2")
        print(f"  Total dSASA (both sides)       : {total_dsasa:8.0f} A^2")
        print(f"  Interface residues             : {len(rows)}")
        print(f"  Residue-residue contact pairs  : {len(pairs)}")
        print(f"  H-bond-like (N/O..N/O <=3.5A)  : {n_hb}")
        print(f"  Salt bridges (<=4.0A)          : {n_sb}")

        for side in (c1, c2):
            sub = [r for r in rows if r["chain"] == side and r["dSASA"] > 5]
            print(f"\n  --- {CHAIN_NAMES[side]} (chain {side}) epitope, {len(sub)} residues "
                  f"with dSASA>5 A^2 ---")
            print(f"  {'res':>7} {'dSASA':>7} {'ncont':>6} {'tight':>6} {'dmin':>6} "
                  f"{'HB':>3} {'SB':>3}  partners")
            for r in sub:
                print(f"  {r['label']:>7} {r['dSASA']:7.1f} {r['n_contacts']:6d} "
                      f"{r['n_tight']:6d} {r['min_dist']:6.2f} {r['n_hbond']:3d} "
                      f"{r['n_saltbridge']:3d}  {','.join(r['partners'])}")

        # polar detail
        print(f"\n  --- Polar network ({name}) ---")
        for pr in report[f'{c1}{c2}']["pairs"]:
            if pr["hbonds"] or pr["saltbridges"]:
                tag = "SALT" if pr["saltbridges"] else "HB  "
                det = pr["saltbridges"] if pr["saltbridges"] else pr["hbonds"]
                items = "; ".join(f"{a}-{b} {d}A" for a, b, d in det)
                print(f"  {tag} {pr['r1']:>12} .. {pr['r2']:<12} {items}")
        print()

    # ---- cross-interface: which IL-4 residues are shared / unique ----
    ab = {r["label"] for r in report["AB"]["residues"] if r["chain"] == "A" and r["dSASA"] > 5}
    ac = {r["label"] for r in report["AC"]["residues"] if r["chain"] == "A" and r["dSASA"] > 5}
    print("=" * 78)
    print("IL-4 EPITOPE PARTITIONING")
    print("=" * 78)
    print(f"  IL-4 residues engaging IL-4Ra only  ({len(ab - ac)}): {sorted(ab - ac, key=lambda s: int(s[1:]))}")
    print(f"  IL-4 residues engaging IL-13Ra1 only({len(ac - ab)}): {sorted(ac - ab, key=lambda s: int(s[1:]))}")
    print(f"  Shared by both receptors            ({len(ab & ac)}): {sorted(ab & ac, key=lambda s: int(s[1:]))}")

    report["summary"] = {
        "il4_epitope_IL4Ra": sorted(ab, key=lambda s: int(s[1:])),
        "il4_epitope_IL13Ra1": sorted(ac, key=lambda s: int(s[1:])),
        "shared": sorted(ab & ac, key=lambda s: int(s[1:])),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / "interface_analysis.json"
    with output.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nWrote {output}")


if __name__ == "__main__":
    analyze()
