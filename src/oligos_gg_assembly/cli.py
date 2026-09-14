"""Top-level command line interface for HiPlex-GGA."""

from __future__ import annotations

import argparse

from . import __version__
from . import add_adapters, silent_mutation, split_oligos


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hiplex-gga",
        description="HiPlex-GGA hierarchical multiplexed Golden Gate assembly tools.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    mutate = subparsers.add_parser("mutate-sites", help="Silently mutate restriction sites")
    mutate.add_argument("--input", "--input-list", "-input_list", dest="input", required=True)
    mutate.add_argument("--output", "-o", required=True)
    mutate.add_argument(
        "--enzyme-type",
        "-enzyme_type",
        default="All",
        choices=sorted(silent_mutation.ENZYME_COMBINATIONS),
    )
    mutate.set_defaults(func=lambda args: silent_mutation.main([
        "--input", args.input,
        "--output", args.output,
        "--enzyme-type", args.enzyme_type,
    ]))

    split = subparsers.add_parser("split", help="Split DNA sequences into HiPlex-GGA oligos")
    split_oligos.add_arguments(split)
    split.set_defaults(func=lambda args: split_oligos.run_from_args(args))

    adapters = subparsers.add_parser("add-adapters", help="Add adapters to HiPlex-GGA naked fragments")
    adapters.add_argument("--input", "--input-list", "-input_list", dest="input", required=True)
    adapters.add_argument("--output", "-o", required=True)
    adapters.add_argument("--subpool-barcodes", "--subpool_barcode_fname", dest="subpool_barcodes", required=True)
    adapters.add_argument("--sequence-barcodes", "--sequence_barcode_fname", dest="sequence_barcodes")
    adapters.add_argument("--subpool-index", "--subp_barc_index", dest="subpool_index", type=int, required=True)
    adapters.add_argument("--frag-num", "--frag_num", dest="frag_num", type=int, default=3)
    adapters.add_argument("--adapter-f", "--adapter_f", dest="adapter_f", default="F" * 20)
    adapters.add_argument("--adapter-r", "--adapter_r", dest="adapter_r", default="R" * 20)
    adapters.add_argument(
        "--enzyme-scheme",
        choices=["standard", "alternating-bsmbi-bsai", "bsai"],
        default="standard",
    )
    adapters.set_defaults(func=lambda args: add_adapters.main([
        "--input", args.input,
        "--output", args.output,
        "--subpool-barcodes", args.subpool_barcodes,
        "--subpool-index", str(args.subpool_index),
        "--frag-num", str(args.frag_num),
        "--adapter-f", args.adapter_f,
        "--adapter-r", args.adapter_r,
        "--enzyme-scheme", args.enzyme_scheme,
    ] + ([] if args.sequence_barcodes is None else ["--sequence-barcodes", args.sequence_barcodes])))
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

