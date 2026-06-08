"""Input/output helpers for oligo design workflows."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path


def read_designs(path: str | Path) -> tuple[dict[str, str], dict[str, str | None]]:
    """Read name/DNA or name/protein/DNA whitespace-delimited input."""
    designs: dict[str, str] = {}
    protein_seqs: dict[str, str | None] = {}
    with open(path, "r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                name, dna = parts
                protein = None
            elif len(parts) == 3:
                name, protein, dna = parts
            else:
                raise ValueError(
                    f"{path}:{line_number}: expected 2 or 3 columns, got {len(parts)}"
                )
            designs[name] = dna.upper()
            protein_seqs[name] = protein
    return designs, protein_seqs


def write_name_sequence_table(records: dict[str, str], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for name, seq in records.items():
            handle.write(f"{name}\t{seq}\n")


def read_subpool_barcodes(path: str | Path, frag_num: int) -> list[list[str]]:
    with open(path, "r", encoding="utf-8-sig") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"{path}: expected a header and at least one barcode row")
    rows = [line.split() for line in lines[1:]]
    barcode_columns = frag_num * 2
    for row in rows:
        if len(row) < barcode_columns + 1:
            raise ValueError(
                f"{path}: expected at least {barcode_columns + 1} columns per row"
            )
    return [[row[i + 1] for row in rows] for i in range(barcode_columns)]


def read_sequence_barcodes(path: str | Path, seq_num: int) -> list[list[str]]:
    with open(path, "r", encoding="utf-8-sig") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    rows = [re.split(r"[,\s]+", line)[1:] for line in lines[1:]]
    if len(rows) < seq_num:
        raise ValueError(f"{path}: expected at least {seq_num} sequence barcode rows")
    for row in rows[:seq_num]:
        if len(row) < 2:
            raise ValueError(f"{path}: each sequence barcode row needs two sequences")
    return rows[:seq_num]


def read_codon_table(path: str | Path) -> dict[str, list[str]]:
    codons: dict[str, list[str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            aa, codon, *_rest = line.strip().split()
            codons.setdefault(aa, []).append(codon.upper())
    return codons


def write_oligo_table(records: dict[str, dict[str, str]], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for frag_seq in records.values():
            for name, seq in frag_seq.items():
                handle.write(f"{name}, {seq}\n")


def write_fragment_table(records: dict[str, list[str]], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for name, fragments in records.items():
            handle.write(f"{name}, {json.dumps(fragments)}\n")


def read_fragment_table(path: str | Path) -> dict[str, list[str]]:
    """Read naked fragments from JSON list or historical Python-list format."""
    fragments: dict[str, list[str]] = {}
    with open(path, "r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            name, raw = line.strip().split(", ", 1)
            try:
                value = json.loads(raw)
            except json.JSONDecodeError:
                value = ast.literal_eval(raw)
            if not isinstance(value, list):
                raise ValueError(f"{path}:{line_number}: fragment value must be a list")
            fragments[name] = [str(item) for item in value]
    return fragments

