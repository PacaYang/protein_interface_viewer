"""Shared molecular constants and small serialization helpers."""

from __future__ import annotations

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common crystallographic substitutions. Their original residue name is
    # retained in reports while calculations use the corresponding amino acid.
    "MSE": "M", "SEC": "C", "PYL": "K",
}

ANION_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
CATION_ATOMS = {
    "LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2"),
    "HIS": ("ND1", "NE2"),
}
FORMAL_CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}


def is_amino_acid(residue) -> bool:
    """Return whether a Biopython residue is a supported protein residue."""

    return residue.get_id()[0] == " " and residue.get_resname().upper() in AA3TO1


def insertion_code(residue) -> str:
    code = residue.get_id()[2]
    return "" if code in (" ", "?", "") else str(code).strip()


def residue_identity(chain_id: str, residue) -> dict:
    """Stable, JSON-friendly residue identity preserving insertion codes."""

    name = residue.get_resname().upper()
    number = int(residue.get_id()[1])
    icode = insertion_code(residue)
    label = f"{AA3TO1.get(name, 'X')}{number}{icode}"
    return {
        "chain_id": str(chain_id),
        "number": number,
        "insertion_code": icode,
        "resname": name,
        "aa1": AA3TO1.get(name, "X"),
        "label": label,
        "key": f"{chain_id}:{number}{icode}",
    }


def residue_label(residue) -> str:
    name = residue.get_resname().upper()
    return f"{AA3TO1.get(name, 'X')}{int(residue.get_id()[1])}{insertion_code(residue)}"


def atom_element(atom) -> str:
    element = (getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element
    name = atom.get_name().strip()
    return name[0].upper() if name else ""


def is_heavy_atom(atom) -> bool:
    return atom_element(atom) != "H" and atom.get_altloc() in (" ", "A", "")
