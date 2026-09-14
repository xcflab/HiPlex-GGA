"""Add Golden Gate adapters, subpool barcodes, and sequence barcodes."""

from __future__ import annotations

import argparse

from .io import (
    read_fragment_table,
    read_sequence_barcodes,
    read_subpool_barcodes,
    write_oligo_table,
)

HindIII = "AAGCTT"
NcoI = "CCATGG"
NdeI = "CATATG"
XhoI = "CTCGAG"
BamHI = "GGATCC"
EcoRI = "GAATTC"
NheI = "GCTAGC"
XbaI = "TCTAGA"
BsaI = "GGTCTCA"
BsmBI = "CGTCTCA"
BbsI = "GAAGACTG"
bsai = "AGAGACC"
bsmbi = "AGAGACG"
bbsi = "CAGTCTTC"

FILLER_50BP = "atgagccatattcaacgggaaacgtcttgctgtaat"
FILLER_100BP = (
    "atgagccatattcaacgggaaacgtcttgctgcgattaaattccaacatggatgctgatttatatgggtatataat"
)
ALPHABET = list("abcdefghijklmnopqrstuvwxyz")


def _inner_enzymes(index: int, scheme: str) -> tuple[str, str]:
    if scheme == "standard":
        if index == 3:
            return BsaI, bbsi
        if index == 4:
            return BbsI, bsai
    if scheme in {"standard", "alternating-bsmbi-bsai"}:
        return (BsmBI, bsai) if index % 2 == 0 else (BsaI, bsmbi)
    if scheme == "bsai":
        return BsaI, bsai
    raise ValueError(f"Unknown enzyme scheme: {scheme}")


def add_adapter_list(
    seq_name: str,
    frag_seq: list[str],
    subpool_barcodes: list[list[str]],
    subpool_index: int,
    frag_num: int,
    adapter_f: str,
    adapter_r: str,
    seq_barcode: list[str] | str,
    enzyme_scheme: str = "standard",
) -> dict[str, str]:
    """Add adapters and barcodes to a list of naked fragment sequences."""
    frag_ad_list: dict[str, str] = {}

    gene_5_en = NheI
    gene_3_en = BamHI
    seq_ba_5_en = NcoI
    seq_ba_3_en = XhoI

    for i in range(frag_num):
        name = f"{ALPHABET[i]}_{seq_name}"
        if frag_seq[i] == "":
            frag_ad_list[name] = "not_split"
            continue

        if seq_barcode:
            seq_barcode_mod_5prime = seq_ba_5_en + seq_barcode[0]
            seq_barcode_mod_3prime = seq_barcode[1] + seq_ba_3_en
        else:
            seq_barcode_mod_5prime = ""
            seq_barcode_mod_3prime = ""

        if i == 0:
            five_prime = (
                adapter_f
                + subpool_barcodes[2 * i][subpool_index]
                + seq_barcode_mod_5prime
                + gene_5_en
            )
            three_prime = bsai + subpool_barcodes[2 * i + 1][subpool_index] + adapter_r
            frag_ad = five_prime + frag_seq[i] + three_prime
            if len(frag_ad) < 211:
                frag_ad = five_prime + frag_seq[i] + bsai + FILLER_100BP + three_prime
            elif len(frag_ad) < 251:
                frag_ad = five_prime + frag_seq[i] + bsai + FILLER_50BP + three_prime
        elif i == frag_num - 1:
            enzyme = BsmBI if i % 2 == 0 else BsaI
            if enzyme_scheme == "bsai":
                enzyme = BsaI
            five_prime = adapter_f + subpool_barcodes[2 * i][subpool_index] + enzyme
            three_prime = (
                gene_3_en
                + seq_barcode_mod_3prime
                + subpool_barcodes[2 * i + 1][subpool_index]
                + adapter_r
            )
            frag_ad = five_prime + frag_seq[i] + three_prime
            if len(frag_ad) < 211:
                frag_ad = five_prime + FILLER_100BP + enzyme + frag_seq[i] + three_prime
            elif len(frag_ad) < 251:
                frag_ad = five_prime + FILLER_50BP + enzyme + frag_seq[i] + three_prime
        else:
            enzyme_5, enzyme_3 = _inner_enzymes(i, enzyme_scheme)
            five_prime = adapter_f + subpool_barcodes[2 * i][subpool_index] + enzyme_5
            three_prime = enzyme_3 + subpool_barcodes[2 * i + 1][subpool_index] + adapter_r
            frag_ad = five_prime + frag_seq[i] + three_prime
            if len(frag_ad) < 211:
                if i % 2 == 1:
                    frag_ad = five_prime + FILLER_100BP + enzyme_5 + frag_seq[i] + three_prime
                else:
                    frag_ad = five_prime + frag_seq[i] + enzyme_3 + FILLER_100BP + three_prime
            elif len(frag_ad) < 251:
                if i % 2 == 1:
                    frag_ad = five_prime + FILLER_50BP + enzyme_5 + frag_seq[i] + three_prime
                else:
                    frag_ad = five_prime + frag_seq[i] + enzyme_3 + FILLER_50BP + three_prime

        frag_ad_list[name] = frag_ad
    return frag_ad_list


def add_adapters_to_fragments(
    fragments: dict[str, list[str]],
    subpool_barcodes: list[list[str]],
    subpool_index: int,
    frag_num: int,
    adapter_f: str,
    adapter_r: str,
    sequence_barcodes: list[list[str]] | None = None,
    enzyme_scheme: str = "standard",
) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for seq_num, (name, frag_seq) in enumerate(fragments.items()):
        seq_barcode = "" if sequence_barcodes is None else sequence_barcodes[seq_num]
        result[name] = add_adapter_list(
            name,
            frag_seq,
            subpool_barcodes,
            subpool_index,
            frag_num,
            adapter_f,
            adapter_r,
            seq_barcode,
            enzyme_scheme=enzyme_scheme,
        )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add Golden Gate adapters to naked fragments.")
    parser.add_argument("--input", dest="input", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--subpool-barcodes", dest="subpool_barcodes", required=True)
    parser.add_argument("--sequence-barcodes", dest="sequence_barcodes")
    parser.add_argument("--subpool-index", dest="subpool_index", type=int, required=True)
    parser.add_argument("--frag-num", dest="frag_num", type=int, default=3)
    parser.add_argument("--adapter-f", dest="adapter_f", default="F" * 20)
    parser.add_argument("--adapter-r", dest="adapter_r", default="R" * 20)
    parser.add_argument(
        "--enzyme-scheme",
        choices=["standard", "alternating-bsmbi-bsai", "bsai"],
        default="standard",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    fragments = read_fragment_table(args.input)
    subpool_barcodes = read_subpool_barcodes(args.subpool_barcodes, args.frag_num)
    sequence_barcodes = (
        None
        if args.sequence_barcodes is None
        else read_sequence_barcodes(args.sequence_barcodes, len(fragments))
    )
    result = add_adapters_to_fragments(
        fragments,
        subpool_barcodes,
        args.subpool_index - 1,
        args.frag_num,
        args.adapter_f,
        args.adapter_r,
        sequence_barcodes,
        args.enzyme_scheme,
    )
    write_oligo_table(result, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

