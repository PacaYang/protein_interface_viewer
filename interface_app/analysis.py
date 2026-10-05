"""Generic pairwise protein-interface calculations.

The original project scripts were written around the A/B/C IL-4 complexes.
This module keeps their thresholds and interaction definitions while accepting
arbitrary chain identifiers, insertion codes, and chain counts.
"""

from __future__ import annotations

import base64
import copy
import itertools

import numpy as np
from Bio.PDB import ShrakeRupley
from scipy.spatial import cKDTree

from .constants import (
    ANION_ATOMS,
    CATION_ATOMS,
    FORMAL_CHARGE,
    atom_element,
    is_amino_acid,
    is_heavy_atom,
    residue_identity,
)
from .structure import LoadedStructure

CONTACT_CUTOFF_A = 5.0
GRID_CONTACT_CUTOFF_A = 6.0
TIGHT_CUTOFF_A = 4.0
HBOND_CUTOFF_A = 3.5
SALT_CUTOFF_A = 4.0
LIKE_CHARGE_CUTOFF_A = GRID_CONTACT_CUTOFF_A


def encode_pair(chain_a: str, chain_b: str) -> str:
    raw = f"{chain_a}\0{chain_b}".encode("utf-8")
    return base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def decode_pair(pair_id: str) -> tuple[str, str]:
    padded = pair_id + "=" * (-len(pair_id) % 4)
    raw = base64.urlsafe_b64decode(padded).decode("utf-8")
    return tuple(raw.split("\0", 1))  # type: ignore[return-value]


def protein_chains(model) -> list[str]:
    return [str(chain.id) for chain in model if any(is_amino_acid(r) for r in chain)]


def residues_for_chain(model, chain_id: str):
    if chain_id not in model.child_dict:
        raise ValueError(f"Protein chain {chain_id!r} was not found")
    return [r for r in model[chain_id] if is_amino_acid(r)]


def _atom_records(model, chain_id: str) -> list[dict]:
    records = []
    for residue in residues_for_chain(model, chain_id):
        ident = residue_identity(chain_id, residue)
        for atom in residue:
            if not is_heavy_atom(atom):
                continue
            records.append({
                "coord": np.asarray(atom.coord, dtype=np.float64),
                "residue_key": ident["key"],
                "residue": ident,
                "atom_name": atom.get_name().strip(),
                "element": atom_element(atom),
                "resname": residue.get_resname().upper(),
            })
    if not records:
        raise ValueError(f"Chain {chain_id!r} contains no supported heavy atoms")
    return records


def _filtered_copy(model, chains: list[str]):
    sub = copy.deepcopy(model)
    wanted = set(chains)
    for chain in list(sub):
        if str(chain.id) not in wanted:
            sub.detach_child(chain.id)
            continue
        for residue in list(chain):
            if not is_amino_acid(residue):
                chain.detach_child(residue.id)
    return sub


def residue_sasa(model, chains: list[str]) -> dict[tuple[str, str], float]:
    """Return residue SASA keyed by ``(chain_id, residue_key)``."""

    sub = _filtered_copy(model, chains)
    ShrakeRupley(probe_radius=1.40, n_points=960).compute(sub, level="R")
    values = {}
    for chain in sub:
        for residue in chain:
            ident = residue_identity(str(chain.id), residue)
            values[(str(chain.id), ident["key"])] = float(getattr(residue, "sasa", 0.0))
    return values


def _interaction_kind(a: dict, b: dict, distance: float) -> tuple[str | None, str | None, str | None]:
    hbond = None
    salt = None
    like = None
    if a["element"] in {"N", "O"} and b["element"] in {"N", "O"} and distance <= HBOND_CUTOFF_A:
        hbond = f"{a['atom_name']}-{b['atom_name']}"
    anion_a = a["atom_name"] in ANION_ATOMS.get(a["resname"], ())
    cation_a = a["atom_name"] in CATION_ATOMS.get(a["resname"], ())
    anion_b = b["atom_name"] in ANION_ATOMS.get(b["resname"], ())
    cation_b = b["atom_name"] in CATION_ATOMS.get(b["resname"], ())
    if distance <= SALT_CUTOFF_A and ((anion_a and cation_b) or (cation_a and anion_b)):
        salt = f"{a['atom_name']}-{b['atom_name']}"
    if distance <= LIKE_CHARGE_CUTOFF_A and ((anion_a and anion_b) or (cation_a and cation_b)):
        like = f"{a['atom_name']}-{b['atom_name']}"
    return hbond, salt, like


def pairwise_contacts(model, chain_a: str, chain_b: str, cutoff: float = GRID_CONTACT_CUTOFF_A) -> dict:
    """Compute atom contacts without constructing a dense atom-distance matrix."""

    first = _atom_records(model, chain_a)
    second = _atom_records(model, chain_b)
    coordinates = np.asarray([x["coord"] for x in second])
    tree = cKDTree(coordinates)
    pairs: dict[tuple[str, str], dict] = {}
    for atom_a in first:
        for index in tree.query_ball_point(atom_a["coord"], cutoff):
            atom_b = second[index]
            distance = float(np.linalg.norm(atom_a["coord"] - atom_b["coord"]))
            key = (atom_a["residue_key"], atom_b["residue_key"])
            stat = pairs.setdefault(key, {
                "residue_a": atom_a["residue"],
                "residue_b": atom_b["residue"],
                "n_contacts": 0,
                "n_tight": 0,
                "min_dist_A": 999.0,
                "hbond_like": [],
                "salt_bridges": [],
                "like_charges": [],
            })
            stat["n_contacts"] += 1
            stat["n_tight"] += int(distance <= TIGHT_CUTOFF_A)
            stat["min_dist_A"] = min(stat["min_dist_A"], distance)
            hbond, salt, like = _interaction_kind(atom_a, atom_b, distance)
            if hbond:
                stat["hbond_like"].append({"atoms": hbond, "distance_A": round(distance, 2)})
            if salt:
                stat["salt_bridges"].append({"atoms": salt, "distance_A": round(distance, 2)})
            if like:
                stat["like_charges"].append({"atoms": like, "distance_A": round(distance, 2)})
    return pairs


def _unique_interactions(items: list[dict]) -> list[dict]:
    seen = set()
    result = []
    for item in items:
        key = (item["atoms"], item["distance_A"])
        if key not in seen:
            seen.add(key)
            result.append(item)
    return result


def analyze_pair(model, chain_a: str, chain_b: str, *, chain_names: dict[str, str] | None = None) -> dict:
    if chain_a == chain_b:
        raise ValueError("An interface requires two different chains")
    chain_names = chain_names or {}
    pairs = pairwise_contacts(model, chain_a, chain_b)
    reported_pairs = [stat for stat in pairs.values() if stat["min_dist_A"] <= CONTACT_CUTOFF_A]
    sasa_a = residue_sasa(model, [chain_a])
    sasa_b = residue_sasa(model, [chain_b])
    sasa_bound = residue_sasa(model, [chain_a, chain_b])

    aggregates: dict[tuple[str, str], dict] = {}
    for stat in reported_pairs:
        for side, residue_key, ident, partner in (
            (chain_a, stat["residue_a"]["key"], stat["residue_a"], stat["residue_b"]),
            (chain_b, stat["residue_b"]["key"], stat["residue_b"], stat["residue_a"]),
        ):
            item = aggregates.setdefault((side, residue_key), {
                "identity": ident,
                "partners": set(),
                "n_contacts": 0,
                "n_tight": 0,
                "min_dist_A": 999.0,
                "n_hbond_like": 0,
                "n_salt_bridges": 0,
                "n_like_charges": 0,
            })
            item["partners"].add(partner["label"])
            item["n_contacts"] += stat["n_contacts"]
            item["n_tight"] += stat["n_tight"]
            item["min_dist_A"] = min(item["min_dist_A"], stat["min_dist_A"])
            item["n_hbond_like"] += len(stat["hbond_like"])
            item["n_salt_bridges"] += len(stat["salt_bridges"])
            item["n_like_charges"] += len(stat["like_charges"])

    free = {**sasa_a, **sasa_b}
    rows = []
    for key, item in aggregates.items():
        side, residue_key = key
        bound_value = sasa_bound.get((side, residue_key), 0.0)
        d_sasa = free.get((side, residue_key), 0.0) - bound_value
        ident = item["identity"]
        charge = FORMAL_CHARGE.get(ident["resname"], 0)
        rows.append({
            **ident,
            "chain_id": side,
            "dSASA_A2": round(float(d_sasa), 2),
            "sasa_free_A2": round(float(free.get((side, residue_key), 0.0)), 2),
            "sasa_bound_A2": round(float(bound_value), 2),
            "charge": charge,
            "partners": sorted(item["partners"]),
            "n_contacts": item["n_contacts"],
            "n_tight": item["n_tight"],
            "min_dist_A": round(float(item["min_dist_A"]), 2),
            "n_hbond_like": item["n_hbond_like"],
            "n_salt_bridges": item["n_salt_bridges"],
            "n_like_charges": item["n_like_charges"],
        })
    rows.sort(key=lambda x: (-x["dSASA_A2"], -x["n_contacts"], x["key"]))

    pair_rows = []
    for stat in sorted(pairs.values(), key=lambda x: x["min_dist_A"]):
        pair_rows.append({
            "residue_a": stat["residue_a"],
            "residue_b": stat["residue_b"],
            "min_dist_A": round(float(stat["min_dist_A"]), 2),
            "n_contacts": stat["n_contacts"],
            "n_tight": stat["n_tight"],
            "hbond_like": _unique_interactions(stat["hbond_like"]),
            "salt_bridges": _unique_interactions(stat["salt_bridges"]),
            "like_charges": _unique_interactions(stat["like_charges"]),
            "within_reported_cutoff": stat["min_dist_A"] <= CONTACT_CUTOFF_A,
        })

    # Keep the residue grid deterministic and in sequence order. The generated
    # contact-map pages use every <=6 Å residue pair, while dSASA is only
    # available for the historical <=5 Å interface rows.
    d_sasa_by_key = {
        (row["chain_id"], row["key"]): row["dSASA_A2"] for row in rows
    }
    contact_identities = {}
    row_keys, col_keys = set(), set()
    for stat in pairs.values():
        first, second = stat["residue_a"], stat["residue_b"]
        contact_identities[first["key"]] = first
        contact_identities[second["key"]] = second
        row_keys.add(first["key"])
        col_keys.add(second["key"])

    def identity_order(chain_id: str, keys: set[str]) -> list[str]:
        residues = residues_for_chain(model, chain_id)
        order = {
            residue_identity(chain_id, residue)["key"]: index
            for index, residue in enumerate(residues)
        }
        return sorted(keys, key=lambda key: order.get(key, len(order)))

    row_axis_keys = identity_order(chain_a, row_keys)
    col_axis_keys = identity_order(chain_b, col_keys)

    def axis_entry(key: str, chain_id: str) -> dict:
        ident = contact_identities[key]
        charge = FORMAL_CHARGE.get(ident["resname"], 0)
        residue_class = "pos" if charge > 0 else "neg" if charge < 0 else "neutral"
        return {
            **ident,
            "dSASA_A2": d_sasa_by_key.get((chain_id, key)),
            "charge": charge,
            "class": residue_class,
        }

    def sequence_breaks(axis: list[dict]) -> list[int]:
        return [
            index for index in range(1, len(axis))
            if axis[index]["number"] - axis[index - 1]["number"] > 1
        ]

    map_rows = [axis_entry(key, chain_a) for key in row_axis_keys]
    map_cols = [axis_entry(key, chain_b) for key in col_axis_keys]
    row_index = {key: index for index, key in enumerate(row_axis_keys)}
    col_index = {key: index for index, key in enumerate(col_axis_keys)}
    map_cells = []
    for stat in sorted(
        pairs.values(),
        key=lambda x: (x["min_dist_A"], x["residue_a"]["key"], x["residue_b"]["key"]),
    ):
        map_cells.append({
            "row": row_index[stat["residue_a"]["key"]],
            "col": col_index[stat["residue_b"]["key"]],
            "min_dist_A": round(float(stat["min_dist_A"]), 2),
            "n_contacts": stat["n_contacts"],
            "n_tight": stat["n_tight"],
            "hbond_like": _unique_interactions(stat["hbond_like"]),
            "salt_bridges": _unique_interactions(stat["salt_bridges"]),
            "like_charges": _unique_interactions(stat["like_charges"]),
            "within_reported_cutoff": stat["min_dist_A"] <= CONTACT_CUTOFF_A,
        })
    contact_map = {
        "row_chain": chain_a,
        "col_chain": chain_b,
        "row_name": chain_names.get(chain_a, f"Chain {chain_a}"),
        "col_name": chain_names.get(chain_b, f"Chain {chain_b}"),
        "cutoff_A": GRID_CONTACT_CUTOFF_A,
        "reported_cutoff_A": CONTACT_CUTOFF_A,
        "tight_cutoff_A": TIGHT_CUTOFF_A,
        "rows": map_rows,
        "cols": map_cols,
        "row_breaks": sequence_breaks(map_rows),
        "col_breaks": sequence_breaks(map_cols),
        "cells": map_cells,
    }

    total_dsasa = sum(max(0.0, free.get(key, 0.0) - sasa_bound.get(key, 0.0)) for key in free)
    bsa = total_dsasa / 2.0
    return {
        "id": encode_pair(chain_a, chain_b),
        "chain_a": chain_a,
        "chain_b": chain_b,
        "name_a": chain_names.get(chain_a, f"Chain {chain_a}"),
        "name_b": chain_names.get(chain_b, f"Chain {chain_b}"),
        "parameters": {
            "contact_cutoff_A": CONTACT_CUTOFF_A,
            "grid_contact_cutoff_A": GRID_CONTACT_CUTOFF_A,
            "tight_cutoff_A": TIGHT_CUTOFF_A,
            "hbond_cutoff_A": HBOND_CUTOFF_A,
            "salt_cutoff_A": SALT_CUTOFF_A,
            "sasa_probe_A": 1.40,
            "sasa_points": 960,
        },
        "buried_surface_area_A2": round(float(bsa), 2),
        "total_positive_dSASA_A2": round(float(total_dsasa), 2),
        "n_contact_residue_pairs": len(reported_pairs),
        "n_grid_residue_pairs": len(pairs),
        "n_interface_residues": len(rows),
        "n_hbond_like": sum(len(x["hbond_like"]) for x in pair_rows),
        "n_salt_bridges": sum(len(x["salt_bridges"]) for x in pair_rows),
        "n_like_charges": sum(len(x["like_charges"]) for x in pair_rows),
        "residues": rows,
        "contacts": pair_rows,
        "contact_map": contact_map,
    }


def analyze_loaded(loaded: LoadedStructure) -> dict:
    model = loaded.model
    chain_ids = protein_chains(model)
    chain_metadata = {item["id"]: item for item in loaded.metadata["chains"]}
    chain_names = {chain_id: chain_metadata[chain_id].get("name", f"Chain {chain_id}") for chain_id in chain_ids}
    pairs = [
        analyze_pair(model, a, b, chain_names=chain_names)
        for a, b in itertools.combinations(chain_ids, 2)
    ]
    return {
        "source": {
            "name": loaded.source_name,
            "format": loaded.source_format,
        },
        "parameters": {
            "interface_definition": "heavy-atom contacts <= 6.0 A; reported contact pairs <= 5.0 A",
            "bsa_definition": "[SASA(chain A) + SASA(chain B) - SASA(pair)] / 2 using all residues",
        },
        "metadata": loaded.metadata,
        "chains": [
            {k: v for k, v in item.items() if k != "residues"}
            for item in loaded.metadata["chains"]
        ],
        "pairs": pairs,
        "summary": {
            "protein_chain_count": len(chain_ids),
            "pair_count": len(pairs),
            "total_contact_pairs": sum(x["n_contact_residue_pairs"] for x in pairs),
        },
    }
