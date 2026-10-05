"""Input loading, RCSB retrieval, and structure metadata normalization."""

from __future__ import annotations

import gzip
import io
import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser

from .constants import is_amino_acid, residue_identity

PDB_ID_RE = re.compile(r"^[A-Za-z0-9]{4}$")
# RCSB uses numeric assembly identifiers today.  Keep the validator a little
# broader for future identifiers while excluding path separators and query
# characters before interpolating the value into a download URL.
ASSEMBLY_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,31}$")


@dataclass
class LoadedStructure:
    """A parsed first-model structure and the text served to the viewer."""

    structure: object
    model: object
    source_name: str
    source_format: str
    source_text: str
    metadata: dict


def parse_structure(path: str | Path, *, source_name: str | None = None) -> LoadedStructure:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
        source_format = "mmcif"
    elif suffix == ".gz":
        text = gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
        is_mmcif = path.name.lower().endswith((".cif.gz", ".mmcif.gz"))
        parser = (MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
                  if is_mmcif else PDBParser(QUIET=True))
        structure = parser.get_structure("input", io.StringIO(text))
        return make_loaded(structure, text, source_name or path.name,
                           "mmcif" if is_mmcif else "pdb")
    else:
        parser = PDBParser(QUIET=True)
        source_format = "pdb"
    structure = parser.get_structure("input", str(path))
    return make_loaded(structure, path.read_text(errors="replace"), source_name or path.name, source_format)


def parse_text(text: str, *, source_name: str = "uploaded.pdb", source_format: str | None = None) -> LoadedStructure:
    fmt = source_format or ("mmcif" if source_name.lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")) else "pdb")
    if fmt == "mmcif":
        parser = MMCIFParser(QUIET=True, auth_chains=True, auth_residues=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", io.StringIO(text))
    return make_loaded(structure, text, source_name, fmt)


def make_loaded(structure, source_text: str, source_name: str, source_format: str) -> LoadedStructure:
    models = list(structure.get_models())
    if not models:
        raise ValueError("The coordinate file contains no models")
    model = models[0]
    chains = []
    for chain in model:
        residues = [r for r in chain if is_amino_acid(r)]
        if not residues:
            continue
        chains.append({
            "id": str(chain.id),
            "name": f"Chain {chain.id}",
            "residue_count": len(residues),
            "residues": [residue_identity(str(chain.id), r) for r in residues],
        })
    if not chains:
        raise ValueError("No supported protein chains were found in the coordinate file")
    metadata = {
        "model_count": len(models),
        "model_id": str(model.id),
        "chains": chains,
        "source_name": source_name,
        "source_format": source_format,
    }
    return LoadedStructure(structure, model, source_name, source_format, source_text, metadata)


def fetch_rcsb(pdb_id: str, *, assembly_id: str | None = None, timeout: float = 30.0) -> tuple[bytes, str, dict]:
    """Fetch an RCSB mmCIF entry or biological assembly.

    RCSB's assembly endpoint uses ``<id>-assembly<assembly>.cif``. The plain
    entry endpoint remains useful for deposited asymmetric-unit coordinates.
    """

    pdb_id = pdb_id.strip().upper()
    if not PDB_ID_RE.fullmatch(pdb_id):
        raise ValueError("PDB IDs must contain exactly four letters or digits")
    if assembly_id is not None:
        assembly_id = str(assembly_id).strip()
        if assembly_id and not ASSEMBLY_ID_RE.fullmatch(assembly_id):
            raise ValueError("Assembly IDs may contain only letters, digits, '.', '_' and '-'")
        if not assembly_id:
            assembly_id = None
    suffix = f"-assembly{assembly_id}" if assembly_id else ""
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}{suffix}.cif"
    request = urllib.request.Request(url, headers={"User-Agent": "interface-analysis-app/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = response.read()
    except urllib.error.HTTPError as exc:
        raise ValueError(f"RCSB could not retrieve {pdb_id}{suffix}: HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise ValueError(f"RCSB retrieval failed for {pdb_id}: {exc.reason}") from exc
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    name = f"{pdb_id}{suffix}.cif"
    return data, name, {"pdb_id": pdb_id, "assembly_id": assembly_id, "url": url}


def rcsb_entry_metadata(pdb_id: str, *, timeout: float = 15.0) -> dict:
    pdb_id = pdb_id.strip().upper()
    if not PDB_ID_RE.fullmatch(pdb_id):
        raise ValueError("PDB IDs must contain exactly four letters or digits")
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    request = urllib.request.Request(url, headers={"User-Agent": "interface-analysis-app/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode())
    except Exception as exc:  # metadata is optional; coordinate retrieval remains usable
        return {"metadata_warning": str(exc), "rcsb_id": pdb_id}
