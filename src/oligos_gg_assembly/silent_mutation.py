"""Silent mutation of restriction enzyme sites."""

from __future__ import annotations

import argparse
from pathlib import Path

from tqdm import tqdm

from .io import read_designs, write_name_sequence_table
from .sequence_utils import find_new_codon, translate

RESTRICTION_SITES = {
    "NdeI": "CATATG",
    "XhoI": "CTCGAG",
    "HindIII": "AAGCTT",
    "EcoRI": "GAATTC",
    "NcoI": "CCATGG",
    "BamHI": "GGATCC",
    "BsaI5": "GGTCTC",
    "BsmBI5": "CGTCTC",
    "BsaI3": "GAGACC",
    "BsmBI3": "GAGACG",
    "NheI": "GCTAGC",
    "XbaI": "TCTAGA",
    "BbsI5": "GAAGAC",
    "BbsI3": "GTCTTC",
}

ENZYME_COMBINATIONS = {
    "BsaI": ("BsaI5", "BsaI3"),
    "BsmBI": ("BsmBI5", "BsmBI3"),
    "TypeII": ("BsaI5", "BsaI3", "BsmBI5", "BsmBI3"),
    "TypeI": ("NdeI", "XhoI", "HindIII", "EcoRI", "NcoI", "BamHI", "NheI", "XbaI"),
    "All": (
        "BsaI5", "BsaI3", "BsmBI5", "BsmBI3", "NdeI", "XhoI", "HindIII",
        "EcoRI", "NcoI", "BamHI", "NheI", "XbaI", "BbsI5", "BbsI3",
    ),
}


def mutate_restriction_sites(seq: str, enzyme_types: tuple[str, ...]) -> str:
    seq = seq.upper()
    selected = {
        enzyme: RESTRICTION_SITES[enzyme]
        for enzyme in enzyme_types
        if enzyme in RESTRICTION_SITES
    }

    enzyme_site_free = False
    while not enzyme_site_free:
        for enzyme, pattern in selected.items():
            while pattern in seq:
                site = seq.find(pattern)
                remainder = site % 3
                if remainder == 0:
                    old_codon = seq[site:site + 3]
                elif remainder == 1:
                    old_codon = seq[site - 1:site + 2] if enzyme == "NcoI" else seq[site + 2:site + 5]
                else:
                    old_codon = seq[site + 1:site + 4]

                aa = translate(old_codon)
                new_codon = find_new_codon(aa, old_codon)
                if new_codon is None:
                    raise ValueError(f"Cannot find synonymous codon for {old_codon}")

                if remainder == 0:
                    seq = seq[:site] + new_codon + seq[site + 3:]
                elif remainder == 1:
                    if enzyme == "NcoI":
                        seq = seq[:site - 1] + new_codon + seq[site + 2:]
                    else:
                        seq = seq[:site + 2] + new_codon + seq[site + 5:]
                else:
                    seq = seq[:site + 1] + new_codon + seq[site + 4:]

        if all(site not in seq for site in RESTRICTION_SITES.values()):
            enzyme_site_free = True
        else:
            selected = {
                enzyme: pattern
                for enzyme, pattern in RESTRICTION_SITES.items()
                if pattern in seq
            }
    return seq


def mutate_designs(input_path: str | Path, enzyme_type: str) -> dict[str, str]:
    designs, _protein_seqs = read_designs(input_path)
    try:
        enzyme_types = ENZYME_COMBINATIONS[enzyme_type]
    except KeyError as exc:
        raise ValueError(f"Unknown enzyme type: {enzyme_type}") from exc

    changed: dict[str, str] = {}
    for name, seq in tqdm(designs.items()):
        before = translate(seq)
        new_seq = mutate_restriction_sites(seq, enzyme_types)
        after = translate(new_seq)
        if before != after:
            raise ValueError(f"Translation changed for {name}")
        changed[name] = new_seq
    return changed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Silently mutate restriction enzyme sites.")
    parser.add_argument("--input", dest="input", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument(
        "--enzyme-type",
        default="All",
        choices=sorted(ENZYME_COMBINATIONS),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    changed = mutate_designs(args.input, args.enzyme_type)
    write_name_sequence_table(changed, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

